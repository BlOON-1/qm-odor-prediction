"""Organize structure-disjoint split files and audit leakage risk."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from qmodor.descriptors.common import backup_existing, load_benchmark_table, write_csv, write_text


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_with_backup(source: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    backup_existing(dest)
    shutil.copy2(source, dest)


def _inspect_legacy_script(script_text: str) -> dict[str, Any]:
    lowered = script_text.lower()
    return {
        "expected_input_structure": "base_dir with fold_* directories, each containing fold_<xx>.csv and test_<xx>.csv, plus ECFP4_data.csv",
        "treats_fold_01_as_train": "train_src = fold_dir / f\"{fold_name}.csv\"" in script_text,
        "treats_test_01_as_test": "test_src = fold_dir / f\"test_{suffix}.csv\"" in script_text,
        "checks_cas_overlap": "overlap = train_cas & test_cas" in script_text,
        "removes_overlapping_cas_from_training_set": "train_cas = train_cas - overlap" in script_text,
        "expects_ecfp4_data_csv": "ecfp4_data.csv" in lowered,
    }


def _summarize_csv(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "file_name": path.name,
            "exists": False,
            "row_count": 0,
            "column_count": 0,
            "columns": "",
            "unique_cas_count": 0,
            "missing_cas_count": 0,
            "duplicate_cas_count": 0,
        }
    df = pd.read_csv(path)
    columns = [str(col) for col in df.columns.tolist()]
    cas_series = df["cas"].fillna("").astype(str).str.strip() if "cas" in df.columns else pd.Series([], dtype="string")
    nonempty_cas = cas_series[cas_series.ne("")] if "cas" in df.columns else pd.Series([], dtype="string")
    duplicate_cas_count = int(nonempty_cas.duplicated().sum()) if "cas" in df.columns else 0
    return {
        "file_name": path.name,
        "exists": True,
        "row_count": int(len(df)),
        "column_count": int(len(columns)),
        "columns": ";".join(columns),
        "unique_cas_count": int(nonempty_cas.nunique()) if "cas" in df.columns else 0,
        "missing_cas_count": int((cas_series.eq("")).sum()) if "cas" in df.columns else len(df),
        "duplicate_cas_count": duplicate_cas_count,
    }


def _read_override(root_candidates: list[Path]) -> dict[str, str] | None:
    for root in root_candidates:
        override_path = root / "split_role_override.yaml"
        if override_path.exists():
            payload = yaml.safe_load(override_path.read_text(encoding="utf-8")) or {}
            return {
                "path": str(override_path),
                "fold_01": str(payload.get("fold_01", "")),
                "test_01": str(payload.get("test_01", "")),
            }
    return None


def _upsert_yaml_section(config_path: Path, section_key: str, section_value: dict[str, Any]) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) if config_path.exists() else {}
    config[section_key] = section_value
    backup_existing(config_path)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8")


def organize_structure_disjoint_split(project_root: Path) -> dict[str, Any]:
    """Inspect, canonicalize, and audit the structure-disjoint split files."""

    workspace_root = project_root.parent
    split_script = workspace_root / "split_ecfp4_by_folds.py"
    fold_path = workspace_root / "fold_01.csv"
    test_path = workspace_root / "test_01.csv"

    splits_dir = project_root / "data" / "splits" / "structure_disjoint_splits"
    reports_dir = project_root / "reports"
    src_split_dir = project_root / "src" / "qmodor" / "splits"
    scripts_dir = project_root / "scripts"
    splits_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    src_split_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir.mkdir(parents=True, exist_ok=True)

    generated_files: list[Path] = []
    warnings: list[str] = []

    if split_script.exists():
        src_copy = src_split_dir / "legacy_split_ecfp4_by_folds.py"
        script_copy = scripts_dir / "legacy_split_ecfp4_by_folds.py"
        _copy_with_backup(split_script, src_copy)
        _copy_with_backup(split_script, script_copy)
        generated_files.extend([src_copy, script_copy])
        script_text = split_script.read_text(encoding="utf-8", errors="replace")
    else:
        script_text = ""
        warnings.append("CRITICAL: split_ecfp4_by_folds.py is missing.")

    script_inspection = _inspect_legacy_script(script_text) if script_text else {
        "expected_input_structure": "unavailable",
        "treats_fold_01_as_train": False,
        "treats_test_01_as_test": False,
        "checks_cas_overlap": False,
        "removes_overlapping_cas_from_training_set": False,
        "expects_ecfp4_data_csv": False,
    }

    fold_summary = _summarize_csv(fold_path)
    test_summary = _summarize_csv(test_path)

    script_convention_role = {"fold_01.csv": "train", "test_01.csv": "test"}
    user_statement_role = {"fold_01.csv": "test", "test_01.csv": "train"}
    inferred_by_row_count = {
        "fold_01.csv": "train" if fold_summary["row_count"] >= test_summary["row_count"] else "test",
        "test_01.csv": "test" if fold_summary["row_count"] >= test_summary["row_count"] else "train",
    }

    override = _read_override([workspace_root, project_root])
    final_role = dict(script_convention_role)
    final_role_source = "script_convention"
    if override:
        fold_role = override.get("fold_01", "").strip().lower()
        test_role = override.get("test_01", "").strip().lower()
        if fold_role in {"train", "test"} and test_role in {"train", "test"}:
            final_role = {"fold_01.csv": fold_role, "test_01.csv": test_role}
            final_role_source = "split_role_override"
    if user_statement_role != script_convention_role:
        warnings.append(
            "CRITICAL: user-stated split roles conflict with the legacy script convention; canonical outputs follow script convention unless split_role_override.yaml explicitly overrides it."
        )
    if final_role != script_convention_role:
        warnings.append("CRITICAL: final split role overrides the legacy script convention.")
    if final_role != inferred_by_row_count:
        warnings.append("CRITICAL: final split role conflicts with row-count evidence.")

    benchmark_info = None
    benchmark_cas_set: set[str] = set()
    try:
        benchmark_info = load_benchmark_table(project_root)
        benchmark_df = benchmark_info.df.copy()
        if "cas" in benchmark_df.columns:
            benchmark_cas_set = set(benchmark_df["cas"].fillna("").astype(str).str.strip())
            benchmark_cas_set.discard("")
    except Exception:
        benchmark_info = None

    ecfp4_index_path = project_root / "data" / "processed" / "ecfp4_fingerprint_index.csv"
    ecfp4_cas_set: set[str] = set()
    if ecfp4_index_path.exists():
        ecfp4_index = pd.read_csv(ecfp4_index_path)
        if "cas" in ecfp4_index.columns:
            ecfp4_cas_set = set(ecfp4_index["cas"].fillna("").astype(str).str.strip())
            ecfp4_cas_set.discard("")

    train_df = None
    test_df = None
    train_output = splits_dir / "structure_disjoint_train.csv"
    test_output = splits_dir / "structure_disjoint_test.csv"
    fold_train_output = splits_dir / "fold_01_train.csv"
    fold_test_output = splits_dir / "fold_01_test.csv"

    overlap_cas: set[str] = set()
    missing_in_benchmark_train: set[str] = set()
    missing_in_benchmark_test: set[str] = set()
    missing_in_ecfp4_train: set[str] = set()
    missing_in_ecfp4_test: set[str] = set()
    duplicate_train = 0
    duplicate_test = 0

    if fold_path.exists() and test_path.exists():
        fold_df = pd.read_csv(fold_path)
        test_input_df = pd.read_csv(test_path)
        if "cas" not in fold_df.columns or "cas" not in test_input_df.columns:
            warnings.append("CRITICAL: one or both split files are missing the `cas` column.")
        else:
            if final_role["fold_01.csv"] == "train":
                train_df = fold_df.copy()
                test_df = test_input_df.copy()
            else:
                train_df = test_input_df.copy()
                test_df = fold_df.copy()

            train_cas = train_df["cas"].fillna("").astype(str).str.strip()
            test_cas = test_df["cas"].fillna("").astype(str).str.strip()
            train_unique = set(train_cas[train_cas.ne("")])
            test_unique = set(test_cas[test_cas.ne("")])
            overlap_cas = train_unique & test_unique
            duplicate_train = int(train_cas[train_cas.ne("")].duplicated().sum())
            duplicate_test = int(test_cas[test_cas.ne("")].duplicated().sum())

            if overlap_cas:
                warnings.append(
                    f"CRITICAL: train/test CAS overlap detected ({len(overlap_cas)} CAS)."
                )
            if benchmark_cas_set:
                missing_in_benchmark_train = train_unique - benchmark_cas_set
                missing_in_benchmark_test = test_unique - benchmark_cas_set
            if ecfp4_cas_set:
                missing_in_ecfp4_train = train_unique - ecfp4_cas_set
                missing_in_ecfp4_test = test_unique - ecfp4_cas_set

            write_csv(train_df, train_output)
            write_csv(test_df, test_output)
            write_csv(train_df, fold_train_output)
            write_csv(test_df, fold_test_output)
            generated_files.extend([train_output, test_output, fold_train_output, fold_test_output])
    else:
        warnings.append("CRITICAL: fold_01.csv and/or test_01.csv is missing; final canonical split files were not created.")

    audit_rows = [
        {"section": "script_inspection", "metric": key, "value": json.dumps(value) if isinstance(value, (dict, list)) else value}
        for key, value in script_inspection.items()
    ]
    for summary in [fold_summary, test_summary]:
        for metric_key in ["row_count", "column_count", "columns", "unique_cas_count", "missing_cas_count", "duplicate_cas_count", "exists"]:
            audit_rows.append(
                {
                    "section": summary["file_name"],
                    "metric": metric_key,
                    "value": summary[metric_key],
                }
            )
    audit_rows.extend(
        [
            {"section": "role_detection", "metric": "script_convention_role", "value": json.dumps(script_convention_role)},
            {"section": "role_detection", "metric": "user_statement_role", "value": json.dumps(user_statement_role)},
            {"section": "role_detection", "metric": "inferred_by_row_count", "value": json.dumps(inferred_by_row_count)},
            {"section": "role_detection", "metric": "final_role", "value": json.dumps(final_role)},
            {"section": "role_detection", "metric": "final_role_source", "value": final_role_source},
            {"section": "split_quality", "metric": "train_test_cas_overlap_count", "value": len(overlap_cas)},
            {"section": "split_quality", "metric": "duplicate_cas_train", "value": duplicate_train},
            {"section": "split_quality", "metric": "duplicate_cas_test", "value": duplicate_test},
            {"section": "split_quality", "metric": "missing_in_benchmark_train", "value": len(missing_in_benchmark_train)},
            {"section": "split_quality", "metric": "missing_in_benchmark_test", "value": len(missing_in_benchmark_test)},
            {"section": "split_quality", "metric": "missing_in_ecfp4_index_train", "value": len(missing_in_ecfp4_train)},
            {"section": "split_quality", "metric": "missing_in_ecfp4_index_test", "value": len(missing_in_ecfp4_test)},
        ]
    )
    if train_df is not None and test_df is not None:
        audit_rows.extend(
            [
                {"section": "final_split", "metric": "train_rows", "value": int(len(train_df))},
                {"section": "final_split", "metric": "test_rows", "value": int(len(test_df))},
                {"section": "final_split", "metric": "train_unique_cas", "value": int(train_df["cas"].fillna("").astype(str).str.strip().replace("", pd.NA).dropna().nunique())},
                {"section": "final_split", "metric": "test_unique_cas", "value": int(test_df["cas"].fillna("").astype(str).str.strip().replace("", pd.NA).dropna().nunique())},
            ]
        )
    if benchmark_cas_set and train_df is not None and test_df is not None:
        final_union = set(train_df["cas"].fillna("").astype(str).str.strip()) | set(test_df["cas"].fillna("").astype(str).str.strip())
        final_union.discard("")
        audit_rows.append(
            {
                "section": "coverage",
                "metric": "benchmark_coverage_fraction",
                "value": len(final_union & benchmark_cas_set) / len(benchmark_cas_set) if benchmark_cas_set else "",
            }
        )

    audit_df = pd.DataFrame(audit_rows)
    audit_output = splits_dir / "structure_disjoint_split_audit.csv"
    write_csv(audit_df, audit_output)
    generated_files.append(audit_output)

    report_lines = [
        "# Structure-Disjoint Split Report",
        "",
        "## 1. Input files inspected",
        f"- `{split_script.name}` exists: `{split_script.exists()}`",
        f"- `{fold_path.name}` exists: `{fold_path.exists()}`",
        f"- `{test_path.name}` exists: `{test_path.exists()}`",
        "",
        "## 2. Role convention detected from split_ecfp4_by_folds.py",
        f"- expected input structure: {script_inspection['expected_input_structure']}",
        f"- treats fold_01.csv as train: `{script_inspection['treats_fold_01_as_train']}`",
        f"- treats test_01.csv as test: `{script_inspection['treats_test_01_as_test']}`",
        f"- checks CAS overlap: `{script_inspection['checks_cas_overlap']}`",
        f"- removes overlapping CAS from training set: `{script_inspection['removes_overlapping_cas_from_training_set']}`",
        f"- expects ECFP4_data.csv: `{script_inspection['expects_ecfp4_data_csv']}`",
        "",
        "## 3. Row counts and unique CAS counts",
        f"- fold_01.csv rows: `{fold_summary['row_count']}`, unique CAS: `{fold_summary['unique_cas_count']}`, missing CAS: `{fold_summary['missing_cas_count']}`, duplicate CAS: `{fold_summary['duplicate_cas_count']}`",
        f"- test_01.csv rows: `{test_summary['row_count']}`, unique CAS: `{test_summary['unique_cas_count']}`, missing CAS: `{test_summary['missing_cas_count']}`, duplicate CAS: `{test_summary['duplicate_cas_count']}`",
        "",
        "## 4. Train/test CAS overlap",
        f"- overlap count: `{len(overlap_cas)}`",
        "",
        "## 5. Missing CAS in benchmark or ECFP4 index if checked",
        f"- missing in benchmark train: `{len(missing_in_benchmark_train)}`",
        f"- missing in benchmark test: `{len(missing_in_benchmark_test)}`",
        f"- missing in ECFP4 index train: `{len(missing_in_ecfp4_train)}`",
        f"- missing in ECFP4 index test: `{len(missing_in_ecfp4_test)}`",
        "",
        "## 6. Final role convention used",
        f"- script_convention_role: `{script_convention_role}`",
        f"- user_statement_role: `{user_statement_role}`",
        f"- inferred_by_row_count: `{inferred_by_row_count}`",
        f"- final_role: `{final_role}`",
        f"- final_role_source: `{final_role_source}`",
        "",
        "## 7. Critical warnings",
    ]
    if warnings:
        report_lines.extend([f"- {warning}" for warning in warnings])
    else:
        report_lines.append("- None.")
    report_lines.extend(["", "## 8. Files generated"])
    report_lines.extend([f"- `{path.relative_to(project_root).as_posix()}`" for path in generated_files])

    report_output = reports_dir / "structure_disjoint_split_report.md"
    write_text("\n".join(report_lines) + "\n", report_output)
    generated_files.append(report_output)

    split_settings_path = project_root / "configs" / "split_settings.yaml"
    _upsert_yaml_section(
        split_settings_path,
        "structure_disjoint_split",
        {
            "source_files": {
                "original_script": "split_ecfp4_by_folds.py",
                "train_file": "data/splits/structure_disjoint_splits/structure_disjoint_train.csv",
                "test_file": "data/splits/structure_disjoint_splits/structure_disjoint_test.csv",
            },
            "similarity_metric": "tanimoto",
            "fingerprint": "ecfp4",
            "high_similarity_threshold": 0.8,
            "final_role_convention": "script_convention",
            "audit_report": "reports/structure_disjoint_split_report.md",
        },
    )

    return {
        "fold_summary": fold_summary,
        "test_summary": test_summary,
        "script_inspection": script_inspection,
        "script_convention_role": script_convention_role,
        "user_statement_role": user_statement_role,
        "inferred_by_row_count": inferred_by_row_count,
        "final_role": final_role,
        "warnings": warnings,
        "train_rows": int(len(train_df)) if train_df is not None else 0,
        "test_rows": int(len(test_df)) if test_df is not None else 0,
        "train_test_overlap": len(overlap_cas),
        "generated_files": generated_files,
    }
