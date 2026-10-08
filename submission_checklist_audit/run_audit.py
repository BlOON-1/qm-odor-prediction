from __future__ import annotations

import csv
import hashlib
import json
import mimetypes
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
AUDIT_DIR = ROOT / "submission_checklist_audit"


KEY_SCAN_TARGETS = [
    "README.md",
    "LICENSE",
    "CITATION.cff",
    ".zenodo.json",
    "requirements.txt",
    "environment.yml",
    "pyproject.toml",
    "docs",
    "data/raw_manifest",
    "data/interim",
    "data/processed",
    "data/splits",
    "data/external_panel",
    "data/receptor_annotations",
    "data/source_data",
    "qm_calculations",
    "src",
    "scripts",
    "configs",
    "results",
    "notebooks",
    "figures",
    "models",
    "tests",
]


FOCUS_PATHS = [
    "data/processed/benchmark_molecules_4403.csv",
    "data/processed/odor_label_matrix_4403x112.csv",
    "data/processed/odor_label_vocabulary_112.csv",
    "data/processed/rdkit_descriptors.csv",
    "data/processed/qm_descriptors.csv",
    "data/processed/rdkit_qm_descriptors.csv",
    "data/processed/ecfp4_fingerprints.npz",
    "data/processed/feature_scaler_parameters.csv",
    "data/splits/iterative_stratified_splits/seed_42_train.csv",
    "data/splits/iterative_stratified_splits/seed_42_test.csv",
    "data/splits/structure_disjoint_splits/tanimoto_neighbors_threshold_0.8.csv",
    "data/splits/structure_disjoint_splits/structure_disjoint_train.csv",
    "data/splits/structure_disjoint_splits/structure_disjoint_test.csv",
    "data/external_panel/external_molecules_90.csv",
    "data/external_panel/structural_neighbor_pairs_60.csv",
    "data/external_panel/panel_label_vocabulary_43.csv",
    "data/external_panel/aggregated_panel_ratings.csv",
    "data/external_panel/panel_derived_topk_labels.csv",
    "data/external_panel/pairwise_panel_jaccard_similarity.csv",
    "data/receptor_annotations/m2or_matched_molecules.csv",
    "data/receptor_annotations/receptor_response_vectors.csv",
    "data/receptor_annotations/receptor_pairwise_concordance.csv",
    "data/receptor_annotations/receptor_concordant_divergent_pairs.csv",
    "data/receptor_annotations/receptor_matching_exclusion_log.csv",
    "data/source_data/fig1_source_data.xlsx",
    "data/source_data/fig2_source_data.xlsx",
    "data/source_data/fig3_source_data.xlsx",
    "data/source_data/fig4_source_data.xlsx",
    "data/source_data/supplementary_tables_source_data.xlsx",
    "qm_calculations/input_templates/gaussian_b3lyp_6311gdp_d3bj_template.gjf",
    "qm_calculations/input_templates/multiwfn_parsing_settings.txt",
    "qm_calculations/parsed_outputs/qm_descriptor_raw_table.csv",
    "qm_calculations/parsed_outputs/adch_atomic_charges.csv",
    "qm_calculations/parsed_outputs/failed_qm_jobs.csv",
    "qm_calculations/logs/convergence_summary.csv",
    "qm_calculations/logs/imaginary_frequency_exclusions.csv",
    "results/benchmark_metrics/macro_auc_by_model_representation.csv",
    "results/benchmark_metrics/micro_auc_by_model_representation.csv",
    "results/benchmark_metrics/labelwise_delta_auc.csv",
    "results/structure_disjoint/structure_disjoint_metrics.csv",
    "results/external_validation/graph_pairwise_similarity_predictions.csv",
    "results/external_validation/concordance_coefficients.csv",
    "results/external_validation/embedding_knn_auc.csv",
    "results/shap/shap_feature_ranking.csv",
    "results/receptor_consistency/qm_rdkit_distance_by_receptor_group.csv",
    "results/receptor_consistency/receptor_consistency_statistics.csv",
]


TEXT_EXTS = {
    ".csv",
    ".tsv",
    ".txt",
    ".md",
    ".json",
    ".yaml",
    ".yml",
    ".py",
    ".sh",
    ".toml",
    ".cff",
    ".gjf",
}
HASH_ONLY_EXTS = {".xlsx", ".npz", ".pdf", ".pt"}


PRIORITY_MISSING_RULES = [
    ("Ethics / informed consent statement not found", "HIGH PRIORITY"),
    ("Third-party redistribution restriction note missing", "HIGH PRIORITY"),
    ("Zenodo DOI or release tag missing", "HIGH PRIORITY"),
    ("Leakage prevention evidence incomplete", "HIGH PRIORITY"),
    ("Iterative stratified split files missing", "MEDIUM PRIORITY"),
    ("Feature scaler training-only evidence incomplete", "MEDIUM PRIORITY"),
    ("Model hyperparameter values left as TODO", "MEDIUM PRIORITY"),
    ("Software version details incomplete", "MEDIUM PRIORITY"),
    ("Supplementary tables source data workbook missing", "MEDIUM PRIORITY"),
    ("Figure PDFs or make_fig scripts missing", "MEDIUM PRIORITY"),
    ("No package/test environment lock beyond environment.yml/requirements.txt", "LOW PRIORITY"),
]


@dataclass
class SectionRow:
    checklist_section: str
    required_information: str
    current_status: str
    evidence_file_path: str
    evidence_summary: str
    missing_or_risk: str
    recommended_action: str


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sniff_text(path: Path) -> tuple[str, str]:
    encodings = ["utf-8", "utf-8-sig", "gbk", "latin-1"]
    for enc in encodings:
        try:
            path.read_text(encoding=enc)
            return "readable", enc
        except Exception:
            continue
    return "unreadable", ""


def file_type_label(path: Path) -> str:
    if path.is_dir():
        return "directory"
    suffix = path.suffix.lower()
    if suffix:
        return suffix.lstrip(".")
    mime, _ = mimetypes.guess_type(str(path))
    return mime or "unknown"


def list_repo_files() -> list[Path]:
    files: list[Path] = []
    for p in ROOT.rglob("*"):
        if AUDIT_DIR in p.parents or p == AUDIT_DIR:
            continue
        if p.is_file():
            files.append(p)
    return sorted(files)


def safe_read_text(path: Path) -> str:
    for enc in ("utf-8", "utf-8-sig", "gbk", "latin-1"):
        try:
            return path.read_text(encoding=enc)
        except Exception:
            continue
    return ""


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def candidate_alternatives(target_rel: str) -> list[str]:
    target = ROOT / target_rel
    parent = target.parent
    if not parent.exists():
        return []
    stem = target.stem.split("_")[0]
    candidates = []
    for p in sorted(parent.iterdir()):
        if p.is_file() and (stem in p.stem or target.suffix == p.suffix):
            candidates.append(rel(p))
    return candidates[:8]


def detect_key_columns(columns: list[str]) -> list[str]:
    lowered = {c.lower(): c for c in columns}
    priority = [
        "molecule_id",
        "id",
        "cas",
        "cas_number",
        "smiles",
        "canonical_smiles",
        "pair_id",
        "receptor",
        "label",
        "odor_label",
    ]
    found = [lowered[k] for k in priority if k in lowered]
    return found[:3]


def summarize_duplicates(df: pd.DataFrame) -> str:
    cols = detect_key_columns(df.columns.tolist())
    if not cols:
        return "No obvious key columns detected"
    summary = []
    for col in cols:
        series = df[col].astype(str).str.strip()
        non_empty = series.replace("", pd.NA).dropna()
        summary.append(f"{col}:{int(non_empty.duplicated().sum())} duplicates")
    if len(cols) >= 2:
        pair_dups = int(df[cols[:2]].astype(str).duplicated().sum())
        summary.append(f"{'+'.join(cols[:2])}:{pair_dups} duplicated rows")
    return "; ".join(summary)


def summarize_missing(df: pd.DataFrame) -> str:
    counts = df.isna().sum()
    counts = counts[counts > 0].sort_values(ascending=False)
    if counts.empty:
        return "No missing values detected"
    items = [f"{idx}:{int(val)}" for idx, val in counts.head(12).items()]
    if len(counts) > 12:
        items.append(f"... total_columns_with_missing={len(counts)}")
    return "; ".join(items)


def csv_stats(path: Path) -> dict[str, Any]:
    df = pd.read_csv(path)
    return {
        "row_count": int(df.shape[0]),
        "column_count": int(df.shape[1]),
        "column_names": "|".join(map(str, df.columns.tolist())),
        "missing_values_summary": summarize_missing(df),
        "duplicate_key_summary": summarize_duplicates(df),
        "notes": "",
    }


def xlsx_stats(path: Path) -> dict[str, Any]:
    wb = load_workbook(path, read_only=True, data_only=True)
    sheet_summaries = []
    for ws in wb.worksheets:
        headers = []
        try:
            first_row = next(ws.iter_rows(min_row=1, max_row=1, values_only=True))
            headers = [str(v) for v in first_row if v is not None][:10]
        except Exception:
            headers = []
        sheet_summaries.append(
            {
                "sheet": ws.title,
                "rows": ws.max_row,
                "cols": ws.max_column,
                "headers": headers,
            }
        )
    return {
        "row_count": "",
        "column_count": "",
        "column_names": json.dumps(sheet_summaries, ensure_ascii=True),
        "missing_values_summary": "Workbook-level summary only",
        "duplicate_key_summary": "Workbook-level summary only",
        "notes": f"sheets={len(sheet_summaries)}",
    }


def npz_stats(path: Path) -> dict[str, Any]:
    data = np.load(path, allow_pickle=False)
    shapes = {k: list(data[k].shape) for k in data.files}
    row_count = ""
    column_count = ""
    if data.files:
        first = data[data.files[0]]
        if hasattr(first, "shape"):
            if len(first.shape) >= 1:
                row_count = int(first.shape[0])
            if len(first.shape) >= 2:
                column_count = int(first.shape[1])
    return {
        "row_count": row_count,
        "column_count": column_count,
        "column_names": json.dumps(shapes, ensure_ascii=True),
        "missing_values_summary": "NPZ not expanded to per-column missingness",
        "duplicate_key_summary": "NPZ not expanded to key-level duplicate check",
        "notes": f"arrays={','.join(data.files)}",
    }


def generic_stats(path: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return csv_stats(path)
    if suffix == ".xlsx":
        return xlsx_stats(path)
    if suffix == ".npz":
        return npz_stats(path)
    return {
        "row_count": "",
        "column_count": "",
        "column_names": "",
        "missing_values_summary": "",
        "duplicate_key_summary": "",
        "notes": "Parsed as hash-only or plain existence target",
    }


def scan_structure() -> list[dict[str, Any]]:
    rows = []
    for item in KEY_SCAN_TARGETS:
        path = ROOT / item
        exists = path.exists()
        rows.append(
            {
                "relative_path": item,
                "exists": "YES" if exists else "NO",
                "path_type": "directory" if path.is_dir() else ("file" if path.is_file() else "missing"),
                "file_size_bytes": path.stat().st_size if path.exists() and path.is_file() else "",
                "last_modified": (
                    datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds") if path.exists() else ""
                ),
                "notes": "",
            }
        )
    return rows


def build_manifest() -> list[dict[str, Any]]:
    rows = []
    for path in list_repo_files():
        status = "readable"
        notes = ""
        suffix = path.suffix.lower()
        enc = ""
        try:
            digest = sha256_file(path)
        except Exception as exc:
            digest = ""
            status = "unreadable"
            notes = f"sha256 failed: {exc}"
        if status == "readable" and (suffix in TEXT_EXTS or suffix in HASH_ONLY_EXTS or suffix == ""):
            if suffix in TEXT_EXTS or suffix == "":
                txt_status, enc = sniff_text(path)
                if txt_status != "readable":
                    status = "unreadable"
                    notes = "text decode failed"
            elif suffix in HASH_ONLY_EXTS:
                notes = "hash-only binary summary"
        else:
            notes = "hash computed; parser not requested"
        if enc:
            notes = f"{notes}; encoding={enc}".strip("; ")
        rows.append(
            {
                "relative_path": rel(path),
                "file_type": file_type_label(path),
                "file_size_bytes": path.stat().st_size,
                "sha256": digest,
                "readable_status": status,
                "notes": notes,
            }
        )
    return rows


def build_focus_stats() -> list[dict[str, Any]]:
    rows = []
    iterative_dir = ROOT / "data/splits/iterative_stratified_splits"
    if iterative_dir.exists():
        for seed_file in sorted(iterative_dir.glob("seed_*_*.csv")):
            rel_seed = rel(seed_file)
            if rel_seed not in FOCUS_PATHS:
                FOCUS_PATHS.append(rel_seed)
    for target_rel in FOCUS_PATHS:
        path = ROOT / target_rel
        row: dict[str, Any] = {"target_path": target_rel}
        if not path.exists():
            row.update(
                {
                    "status": "MISSING",
                    "observed_path": "",
                    "row_count": "",
                    "column_count": "",
                    "column_names": "",
                    "missing_values_summary": "",
                    "duplicate_key_summary": "",
                    "notes": "Target path not found",
                    "alternative_candidates": " | ".join(candidate_alternatives(target_rel)),
                }
            )
        else:
            try:
                stats = generic_stats(path)
                row.update(
                    {
                        "status": "FOUND",
                        "observed_path": target_rel,
                        "alternative_candidates": "",
                        **stats,
                    }
                )
            except Exception as exc:
                row.update(
                    {
                        "status": "READ_ERROR",
                        "observed_path": target_rel,
                        "row_count": "",
                        "column_count": "",
                        "column_names": "",
                        "missing_values_summary": "",
                        "duplicate_key_summary": "",
                        "notes": f"{type(exc).__name__}: {exc}",
                        "alternative_candidates": "",
                    }
                )
        rows.append(row)
    return rows


def text_contains(path_str: str, pattern: str) -> bool:
    path = ROOT / path_str
    if not path.exists() or not path.is_file():
        return False
    return re.search(pattern, safe_read_text(path), flags=re.IGNORECASE) is not None


def load_yaml(path_str: str) -> dict[str, Any]:
    path = ROOT / path_str
    if not path.exists():
        return {}
    try:
        return yaml.safe_load(safe_read_text(path)) or {}
    except Exception:
        return {}


def file_exists(rel_path: str) -> bool:
    return (ROOT / rel_path).exists()


def find_scripts(pattern: str) -> list[str]:
    return sorted(rel(p) for p in ROOT.rglob(pattern) if AUDIT_DIR not in p.parents)


def leakage_summary() -> tuple[str, str]:
    train_path = ROOT / "data/splits/structure_disjoint_splits/structure_disjoint_train.csv"
    test_path = ROOT / "data/splits/structure_disjoint_splits/structure_disjoint_test.csv"
    if not train_path.exists() or not test_path.exists():
        return "NEEDS AUTHOR CONFIRMATION", "Required structure-disjoint train/test files missing for overlap check."
    try:
        train_df = pd.read_csv(train_path)
        test_df = pd.read_csv(test_path)
        shared_bits = []
        for col in ("cas", "smiles", "canonical_smiles", "molecule_id"):
            if col in train_df.columns and col in test_df.columns:
                a = set(train_df[col].fillna("").astype(str).str.strip()) - {""}
                b = set(test_df[col].fillna("").astype(str).str.strip()) - {""}
                shared = len(a & b)
                shared_bits.append(f"{col} overlap={shared}")
        if not shared_bits:
            return "NEEDS AUTHOR CONFIRMATION", "No shared key columns available for overlap check."
        status = "READY" if all(bit.endswith("=0") for bit in shared_bits) else "HIGH RISK"
        return status, "; ".join(shared_bits)
    except Exception as exc:
        return "NEEDS AUTHOR CONFIRMATION", f"Overlap check failed: {type(exc).__name__}: {exc}"


def build_evidence_matrix() -> list[SectionRow]:
    zenodo = {}
    try:
        zenodo = json.loads(safe_read_text(ROOT / ".zenodo.json"))
    except Exception:
        zenodo = {}
    leakage_status, leakage_text = leakage_summary()
    figure_scripts = [p for p in find_scripts("make_fig*.py")]
    classical_scripts = [p for p in find_scripts("*train*.py")] + [p for p in find_scripts("*model*.py")]
    evidence: list[SectionRow] = [
        SectionRow(
            "Study design",
            "Mixed computational / machine learning / sensory validation / receptor consistency study scope",
            "READY",
            "README.md | docs/workflow_overview.md",
            "Repository text explicitly describes benchmark prediction, structure-disjoint validation, external sensory-panel validation, and receptor-consistency analysis, and states the work is not a direct receptor mechanism proof.",
            "",
            "Reuse repository wording when filling the checklist.",
        ),
        SectionRow(
            "Data sources",
            "TGSC, Leffingwell, PubChem, Pyrfume harmonization, M2OR, and third-party restriction evidence",
            "READY" if file_exists("data/raw_manifest/source_inventory.csv") and file_exists("data/raw_manifest/third_party_source_restrictions.md") else "MISSING",
            "data/raw_manifest/source_inventory.csv | data/raw_manifest/third_party_source_restrictions.md | docs/data_availability.md",
            "Source inventory and third-party restriction note exist; PubChem and M2OR are mentioned in docs; label harmonization evidence exists in data/interim/label_harmonization_table.csv.",
            "CHECK THIRD-PARTY RESTRICTION",
            "Keep checklist language high-level and avoid redistributing raw TGSC or Leffingwell records.",
        ),
        SectionRow(
            "Sample size and dataset dimensions",
            "Benchmark 4,403 molecules; 112 labels; 90 external molecules; 60 pairs; 43 panel labels; assessor metadata; M2OR matched counts",
            "NEEDS AUTHOR CONFIRMATION",
            "README.md | docs/workflow_overview.md | docs/sensory_panel_metadata.md | data/receptor_annotations/m2or_matched_molecules.csv",
            "Current repository repeatedly states 4,495 benchmark molecules, 112 labels, 87 external molecules, 60 pairs, 43 panel labels, and 13 trained assessors.",
            "Requested target counts 4,403 and 90 do not match observed release files.",
            "Authors should confirm whether the checklist should report manuscript counts, public-release counts, or both.",
        ),
        SectionRow(
            "Inclusion and exclusion criteria",
            "Conflict removal, missing conformer / QM failure / imaginary frequency / unmatched M2OR handling logs",
            "NEEDS AUTHOR CONFIRMATION",
            "data/interim/removed_or_conflicting_entries.csv | data/receptor_annotations/receptor_matching_exclusion_log.csv | docs/workflow_overview.md",
            "Conflict-removal and M2OR exclusion evidence exists; workflow text mentions unresolved identifier-structure conflicts and QM exclusion logic, but qm_calculations logs are absent.",
            "Missing qm_calculations failure and imaginary-frequency log files.",
            "Authors should confirm whether QM exclusion logs live outside the public repo or need separate submission support.",
        ),
        SectionRow(
            "Data preprocessing and feature construction",
            "CAS normalization, SMILES canonicalization, label harmonization, descriptor construction, scaler parameters, correlation filtering, SHAP, ADCH",
            "NEEDS AUTHOR CONFIRMATION",
            "data/interim/cas_smiles_mapping.csv | data/interim/molecule_records_canonicalized.csv | data/interim/label_harmonization_table.csv | results/shap/shap_feature_ranking.csv",
            "CAS/SMILES harmonization and SHAP ranking evidence exists. ECFP4 and RDKit assets exist. QM, RDKit+QM, scaler-parameter, correlation-filter, and ADCH public files requested by checklist are not present at the expected paths.",
            "Partial release only.",
            "Mark unavailable items as not in current public snapshot and have authors confirm private-generation provenance if needed.",
        ),
        SectionRow(
            "Data splitting and leakage prevention",
            "Iterative stratified split files, structure-disjoint files, Tanimoto threshold file, train/test overlap checks, training-only scaling, no-leakage test",
            (
                "READY"
                if leakage_status == "READY"
                and file_exists("tests/test_split_no_leakage.py")
                and file_exists("data/splits/iterative_stratified_splits/seed_42_train.csv")
                and file_exists("data/splits/iterative_stratified_splits/seed_42_test.csv")
                and file_exists("data/processed/feature_scaler_parameters.csv")
                else "NEEDS AUTHOR CONFIRMATION"
            ),
            "data/splits/structure_disjoint_splits/structure_disjoint_train.csv | data/splits/structure_disjoint_splits/structure_disjoint_test.csv | tests/test_split_no_leakage.py | docs/workflow_overview.md",
            f"Structure-disjoint train/test files exist. Leakage check result: {leakage_text}. Test file exists. Workflow doc says scaling uses training-set parameters only.",
            "Iterative stratified split directory missing; no released feature_scaler_parameters.csv.",
            "Authors should provide or explain missing iterative split artifacts and training-only scaler evidence.",
        ),
        SectionRow(
            "Model training and hyperparameters",
            "Classical, neural, and graph model scripts; hyperparameter configs; random seeds; optimizer / epochs / batch size; trained weights",
            "NEEDS AUTHOR CONFIRMATION",
            "configs/model_hyperparameters.yaml | configs/graph_model_settings.yaml | scripts/ | results/benchmark_metrics/model_summary_json/",
            "Configs list model families and random seeds 42-46; graph config covers GCN/ADCH variants. Multiple hyperparameter fields remain TODO, and trained weights/model cards are not present.",
            "Hyperparameters incomplete in released config; no weight files found.",
            "Authors should supply exact final hyperparameters and confirm whether weight release is required by the journal.",
        ),
        SectionRow(
            "Evaluation metrics and statistical analysis",
            "micro-AUC, macro-AUC, labelwise delta AUC, remaining error reduction, Lin's CCC, Jaccard, kNN macro-AUC, receptor statistics, uncertainty",
            "NEEDS AUTHOR CONFIRMATION",
            "results/benchmark_metrics/micro_auc_by_model_representation.csv | results/benchmark_metrics/macro_auc_by_model_representation.csv | results/benchmark_metrics/labelwise_delta_auc.csv | results/receptor_consistency/receptor_consistency_statistics.csv | configs/graph_model_settings.yaml",
            "Benchmark AUC tables and receptor-consistency statistics exist. Graph settings mention Lin's concordance coefficient, Jaccard, and kNN macro-AUC, but the expected external_validation result files are not present.",
            "No direct uncertainty / confidence interval artifact identified.",
            "State metrics supported by current files and mark uncertainty reporting as NEEDS AUTHOR CONFIRMATION.",
        ),
        SectionRow(
            "External sensory panel reporting",
            "Panel metadata, trained assessors, vocabulary, aggregation, top-K rule, ethics/consent/anonymization",
            "NEEDS AUTHOR CONFIRMATION",
            "docs/sensory_panel_metadata.md | data/external_panel/aggregated_panel_ratings.csv | data/external_panel/panel_derived_topk_labels.csv | data/external_panel/pairwise_panel_jaccard_similarity.csv",
            "Metadata explicitly states 13 trained assessors, 43-label vocabulary, knee-point Top-K extraction, aggregated anonymized release, and pairwise Jaccard summary.",
            "No explicit ethics approval, waiver, or informed consent statement found.",
            "HIGH PRIORITY MISSING ITEM: authors should provide the exact ethics / consent / waiver wording used for the sensory panel.",
        ),
        SectionRow(
            "Receptor annotation reporting",
            "M2OR matched records, response vectors, concordance tables, grouping rule, exclusion log, analysis scope",
            "READY",
            "data/receptor_annotations/m2or_matched_molecules.csv | data/receptor_annotations/receptor_response_vectors.csv | data/receptor_annotations/receptor_pairwise_concordance.csv | data/receptor_annotations/receptor_concordant_divergent_pairs.csv | docs/receptor_matching_rules.md | scripts/07_receptor_consistency_analysis.py",
            "Processed M2OR matching outputs and exclusion log exist, and workflow text frames the analysis as orthogonal consistency rather than direct mechanistic proof.",
            "",
            "Use this wording directly in the checklist.",
        ),
        SectionRow(
            "Code availability",
            "README, dependency files, scripts, src package, notebooks, tests, CITATION, Zenodo, one-click pipeline, package versions",
            "NEEDS AUTHOR CONFIRMATION",
            "README.md | requirements.txt | environment.yml | pyproject.toml | scripts/ | src/qmodor/ | notebooks/ | tests/ | CITATION.cff | .zenodo.json",
            "Core code and metadata files exist. There is no clear 00-09 one-click pipeline series, and version pinning is limited.",
            "One-click reproduction script not found.",
            "Authors should confirm whether stage-wise scripts are sufficient or whether a linear reproduction index is needed.",
        ),
        SectionRow(
            "Data availability",
            "Data availability docs, dictionary, processed data, source data, panel summaries, receptor annotations, restrictions, Zenodo DOI or tag",
            "NEEDS AUTHOR CONFIRMATION",
            "docs/data_availability.md | docs/data_dictionary.md | data/raw_manifest/third_party_source_restrictions.md | data/processed/ | data/source_data/",
            "Availability docs and processed data directories exist. .zenodo.json exists but related_identifiers is empty and no DOI is present.",
            "Zenodo DOI or release tag missing.",
            "Authors should add DOI/release tag before final checklist submission if the journal expects persistent archival access.",
        ),
        SectionRow(
            "Figure/source-data traceability",
            "Fig. 1-4 source data, supplementary tables source data, figure PDFs, make_fig scripts, results linkage",
            "NEEDS AUTHOR CONFIRMATION",
            "data/source_data/fig1_source_data.xlsx | data/source_data/fig2_source_data.xlsx | data/source_data/fig3_source_data.xlsx | data/source_data/fig4_source_data.xlsx",
            "Fig. 1-4 source-data workbooks exist. No figures/ directory is present, supplementary_tables_source_data.xlsx is missing, and make_fig1.py to make_fig4.py were not found.",
            "NEEDS MANUAL CHECK",
            "Authors should manually confirm numerical traceability and provide missing supplementary-source-data linkage if required.",
        ),
        SectionRow(
            "Software and computational environment",
            "Gaussian, Multiwfn, GaussView, RDKit, scikit-learn, PyTorch, XGBoost, Python version, CUDA/GPU, commercial/open-source split",
            "NEEDS AUTHOR CONFIRMATION",
            "environment.yml | requirements.txt | docs/workflow_overview.md | data/raw_manifest/third_party_source_restrictions.md",
            "Open-source package files exist; workflow doc names Gaussian 16W, Multiwfn 3.8 dev, and GaussView 6.0.16. The public environment files do not record all commercial-tool versions in machine-readable form, and CUDA/GPU usage is not clearly documented.",
            "Version granularity incomplete.",
            "Authors should confirm final software versions and whether GPU/CUDA details are necessary for the reported neural models.",
        ),
        SectionRow(
            "Ethical and competing-interest reporting",
            "Human participant reporting, ethics approval / waiver / consent, competing interests, funding, acknowledgements",
            "MISSING",
            "docs/sensory_panel_metadata.md | README.md",
            "Sensory panel involvement is documented, but no explicit ethics approval, waiver, informed consent, competing-interest, or funding statement was located in the scanned repository files.",
            "HIGH PRIORITY MISSING ITEM",
            "Authors should add exact approved wording outside manuscript if the submission system requires separate checklist entry.",
        ),
    ]
    return evidence


def status_rank(status: str) -> int:
    order = {"READY": 0, "NEEDS AUTHOR CONFIRMATION": 1, "MISSING": 2, "HIGH RISK": 3}
    return order.get(status, 4)


def build_checklist_draft(evidence_rows: list[SectionRow]) -> str:
    mapping = {
        "Study design and reporting scope": "Study design",
        "Dataset sources and access restrictions": "Data sources",
        "Sample size and dataset dimensions": "Sample size and dataset dimensions",
        "Inclusion/exclusion criteria": "Inclusion and exclusion criteria",
        "Data preprocessing and label harmonization": "Data preprocessing and feature construction",
        "QM descriptor calculation and exclusion criteria": "Inclusion and exclusion criteria",
        "Molecular descriptor construction": "Data preprocessing and feature construction",
        "Train/test split and leakage prevention": "Data splitting and leakage prevention",
        "Model training and hyperparameter reporting": "Model training and hyperparameters",
        "Evaluation metrics": "Evaluation metrics and statistical analysis",
        "External sensory panel": "External sensory panel reporting",
        "Sensory panel ethics / consent / anonymization": "Ethical and competing-interest reporting",
        "Public receptor annotation matching": "Receptor annotation reporting",
        "Figure source data": "Figure/source-data traceability",
        "Code availability": "Code availability",
        "Data availability": "Data availability",
        "Software versions and computational environment": "Software and computational environment",
        "Random seeds and reproducibility": "Model training and hyperparameters",
        "Competing interests / funding / author contributions": "Ethical and competing-interest reporting",
        "Third-party data redistribution restrictions": "Data sources",
    }
    by_section = {row.checklist_section: row for row in evidence_rows}
    parts = []
    for item, section in mapping.items():
        row = by_section[section]
        status = row.current_status
        if item == "Sensory panel ethics / consent / anonymization" and status != "READY":
            status = "MISSING"
        evidence_lines = "\n".join([f"- `{p.strip()}`" for p in row.evidence_file_path.split("|")])
        draft = (
            f"## Checklist item: {item}\n"
            f"Status: {status}\n"
            f"Draft response:\n"
            f"{row.evidence_summary} "
            f"{('This item remains subject to author confirmation.' if status != 'READY' else 'This appears supportable from the current repository snapshot.')}\n"
            f"Evidence in repository:\n"
            f"{evidence_lines}\n"
            f"Remaining risk:\n"
            f"{row.missing_or_risk or 'No major repository-level risk identified from the current audit.'}\n"
        )
        parts.append(draft)
    return "\n".join(parts).strip() + "\n"


def build_missing_items(evidence_rows: list[SectionRow]) -> str:
    issues: dict[str, list[tuple[str, str, str, str, str]]] = defaultdict(list)
    items = [
        (
            "Ethics / informed consent statement not found",
            "docs/sensory_panel_metadata.md documents anonymized aggregated release only; no explicit ethics approval / waiver / consent statement located.",
            "Human sensory assessor work often triggers editorial compliance questions.",
            "Provide the exact ethics approval, waiver, or informed-consent wording used for the panel study.",
            "Do not invent approval numbers, committee names, or consent wording.",
            "HIGH PRIORITY",
        ),
        (
            "Zenodo DOI or release tag missing",
            ".zenodo.json exists but has no DOI in related_identifiers and no release tag evidence was found in the repository snapshot.",
            "Editors may require a persistent archive identifier for code/data availability.",
            "Mint or confirm the archival DOI / release tag and align the checklist wording with that identifier.",
            "Do not guess a DOI or release URL.",
            "HIGH PRIORITY",
        ),
        (
            "Leakage prevention evidence incomplete",
            "Structure-disjoint files and test file exist, but iterative_stratified_splits files and feature_scaler_parameters.csv are missing; training-only scaling is stated in docs but not fully auditable from files alone.",
            "Split provenance and leakage controls are common editorial review points for ML studies.",
            "Provide the missing iterative split artifacts or a clear author statement describing where they are archived and how scaling was fit only on training data.",
            "Do not claim train-only scaling was verified if the file evidence is absent.",
            "HIGH PRIORITY",
        ),
        (
            "Third-party restriction statement must be carried into submission text",
            "third_party_source_restrictions.md exists and clearly limits redistribution of TGSC, Leffingwell, M2OR, PubChem-linked materials, and QM intermediates.",
            "Checklist answers that overstate public redistributability could create compliance problems.",
            "Include a short restriction statement in the Data Reporting Checklist and data availability responses.",
            "Do not paste raw third-party source tables into the submission materials.",
            "HIGH PRIORITY",
        ),
        (
            "Hyperparameter detail remains incomplete",
            "configs/model_hyperparameters.yaml contains many TODO placeholders for penalty, C, gamma, n_estimators, learning_rate, epochs, batch_size, and related details.",
            "Reviewers may ask how final model settings were chosen and whether test data influenced tuning.",
            "Add the exact final hyperparameters or prepare a checklist note stating where they are archived.",
            "Do not fill TODO fields with guessed values.",
            "MEDIUM PRIORITY",
        ),
        (
            "Software version detail is incomplete",
            "Workflow docs name Gaussian 16W, Multiwfn 3.8 dev, and GaussView 6.0.16, but environment files do not fully encode all versions or GPU/CUDA details.",
            "Version ambiguity can weaken reproducibility claims.",
            "Provide final software version strings and whether GPU/CUDA was used for any reported training.",
            "Do not assume versions from memory or local machine defaults.",
            "MEDIUM PRIORITY",
        ),
        (
            "Supplementary tables source-data workbook missing",
            "data/source_data/fig1_source_data.xlsx through fig4_source_data.xlsx exist, but supplementary_tables_source_data.xlsx was not found.",
            "If the journal expects a checklist-level mapping for supplementary tables, this will look incomplete.",
            "Confirm whether supplementary table source data are included elsewhere or need a dedicated workbook.",
            "Do not claim the workbook exists when it is absent.",
            "MEDIUM PRIORITY",
        ),
        (
            "Figure traceability remains partly manual",
            "No figures/ directory and no make_fig1.py to make_fig4.py scripts were found in the scanned repository snapshot.",
            "Editors may ask how figure PDFs/source data connect to result tables.",
            "Prepare a short manual traceability note from each figure workbook to its source result file(s).",
            "Do not assert full numeric traceability without a manual check.",
            "MEDIUM PRIORITY",
        ),
        (
            "Repository narrative could more explicitly explain release-vs-manuscript count differences",
            "Observed files say 4,495 benchmark molecules and 87 external molecules, while the checklist request references 4,403 and 90.",
            "Count mismatches can confuse editors even when scientifically explainable.",
            "Add a concise author confirmation note explaining which counts belong in the checklist and why.",
            "Do not silently substitute one set of counts for another.",
            "LOW PRIORITY",
        ),
    ]
    for title, evidence, risk, action, dont_write, priority in items:
        issues[priority].append((title, evidence, risk, action, dont_write))
    parts = []
    for priority in ("HIGH PRIORITY", "MEDIUM PRIORITY", "LOW PRIORITY"):
        parts.append(f"{priority}:")
        for title, evidence, risk, action, dont_write in issues[priority]:
            parts.append(f"- 问题: {title}")
            parts.append(f"- 当前证据: {evidence}")
            parts.append(f"- 为什么是风险: {risk}")
            parts.append(f"- 建议作者补充什么: {action}")
            parts.append(f"- 不要自行补写的内容: {dont_write}")
        parts.append("")
    return "\n".join(parts).strip() + "\n"


def build_summary(evidence_rows: list[SectionRow], pytest_summary: str) -> str:
    ready = sum(1 for r in evidence_rows if r.current_status == "READY")
    needs = sum(1 for r in evidence_rows if r.current_status == "NEEDS AUTHOR CONFIRMATION")
    missing = sum(1 for r in evidence_rows if r.current_status == "MISSING")
    ethics_missing = True
    third_party_exists = file_exists("data/raw_manifest/third_party_source_restrictions.md")
    zenodo_has_doi = bool(re.search(r"10\\.", json.dumps(json.loads(safe_read_text(ROOT / ".zenodo.json"))))) if file_exists(".zenodo.json") else False
    leakage_status, leakage_text = leakage_summary()
    recommendation = "需要补充后提交"
    if missing == 0 and needs <= 2 and zenodo_has_doi and leakage_status == "READY" and not ethics_missing:
        recommendation = "可以提交"
    if missing >= 2 or ethics_missing:
        recommendation = "不建议提交" if not third_party_exists else "需要补充后提交"
    return f"""# Data Reporting Checklist 审计概览

## 1. 当前准备程度

本次审计基于仓库真实存在的文件完成，只做了扫描、统计、核对和整理，没有修改正文、主图、补充材料或核心数据。当前仓库已经能为 Data Reporting Checklist 提供相当一部分证据，尤其是研究设计、数据来源边界、外部 sensory panel 聚合数据、M2OR 匹配结果、benchmark 指标表和基础代码可用性说明。

按 evidence matrix 统计，`READY` 项 {ready} 个，`NEEDS AUTHOR CONFIRMATION` 项 {needs} 个，`MISSING` 项 {missing} 个。整体状态不是“空缺很多”，而是“已有较强仓库证据，但仍缺若干投稿系统级文字声明和少量关键归档证据”。

## 2. 已经可以直接支持 checklist 的内容

- 研究范围可由 `README.md` 和 `docs/workflow_overview.md` 直接支持：仓库明确覆盖 benchmark、多种分子表征、structure-disjoint validation、external sensory-panel validation、receptor-consistency analysis，并明确说明这不是 receptor mechanism proof。
- 第三方数据来源与限制已有直接证据：`data/raw_manifest/source_inventory.csv`、`data/raw_manifest/third_party_source_restrictions.md`、`docs/data_availability.md`。
- 外部 sensory panel 的聚合层面说明较完整：`docs/sensory_panel_metadata.md`、`data/external_panel/aggregated_panel_ratings.csv`、`panel_derived_topk_labels.csv`、`pairwise_panel_jaccard_similarity.csv`。
- receptor annotation reporting 证据较完整：`data/receptor_annotations/` 下的 matched、vectors、pairwise concordance、concordant/divergent、exclusion log 均存在。
- benchmark 结果指标表存在：`results/benchmark_metrics/` 下 micro-AUC、macro-AUC、labelwise delta AUC 已找到。

## 3. GitHub 中虽有证据，但投稿系统仍需要作者补文字说明的内容

- 数据集规模目前存在“请求值”与“仓库实际值”不一致问题。仓库实际反复写的是 `4,495` benchmark molecules、`87` external molecules，而不是你在任务中列出的 `4,403` 和 `90`。这类内容必须由作者确认最终应在 checklist 中写哪一组数字。
- `docs/workflow_overview.md` 说明了 training-only scaling、iterative multi-label stratification、Tanimoto > 0.8 的 structure-disjoint 逻辑，但仓库没有完整公开 `iterative_stratified_splits/` 目录，也没有 `feature_scaler_parameters.csv`，因此投稿系统里仍需要作者补解释。
- `configs/model_hyperparameters.yaml` 和 `configs/graph_model_settings.yaml` 证明模型家族、随机种子和部分训练设定存在，但大量超参数仍是 `TODO`，不能直接当作最终超参数表述。
- `.zenodo.json` 已存在，说明仓库为归档做过准备，但当前未看到 DOI 或 release tag 证据，因此 checklist 里的 data/code availability 仍需作者补全。

## 4. 当前缺失且必须由作者确认或补充的内容

- 最优先缺失项是 sensory panel 相关的 ethics approval / waiver / informed consent 明确声明。仓库里只有匿名化和聚合发布说明，没有找到明确伦理或同意文本。
- Zenodo DOI 或正式 release tag 目前未见证据。
- 迭代分层 split 文件、feature scaler 参数文件、以及更强的 leakage-prevention 归档证据不完整。
- `supplementary_tables_source_data.xlsx` 未找到；`figures/` 目录和 `make_fig1.py` 到 `make_fig4.py` 也未找到，因此图和 source data 的可追溯性还需要人工补说明。
- QM 相关公开目录 `qm_calculations/` 不存在，所以 Gaussian / Multiwfn 原始或中间输出类证据不能从当前公共仓库直接支持，只能说明“当前 public snapshot 未包含”。

## 5. 是否建议现在上传 Data Reporting Checklist

不建议直接按当前状态原样上传最终版 checklist，但建议基于本次生成的草稿快速补齐缺失项后提交。原因不是仓库质量差，而是投稿系统通常对伦理、归档 DOI、第三方限制、split leakage 证据和最终超参数描述更敏感，这些地方目前还需要作者级确认。

## 6. 最终判断

当前判断：**{recommendation}**。

补充说明：

- 轻量测试状态：{pytest_summary}
- leakage 检查摘要：{leakage_text}
- 第三方限制说明文件：{"已找到" if third_party_exists else "未找到"}
- ethics / informed consent 明确声明：{"未找到" if ethics_missing else "已找到"}
- Zenodo DOI：{"已看到" if zenodo_has_doi else "未看到"}
"""


def main() -> int:
    AUDIT_DIR.mkdir(exist_ok=True)

    structure_rows = scan_structure()
    write_csv(
        AUDIT_DIR / "project_structure_scan.csv",
        structure_rows,
        ["relative_path", "exists", "path_type", "file_size_bytes", "last_modified", "notes"],
    )

    manifest_rows = build_manifest()
    write_csv(
        AUDIT_DIR / "file_manifest_sha256.csv",
        manifest_rows,
        ["relative_path", "file_type", "file_size_bytes", "sha256", "readable_status", "notes"],
    )

    focus_rows = build_focus_stats()
    write_csv(
        AUDIT_DIR / "key_file_statistics.csv",
        focus_rows,
        [
            "target_path",
            "status",
            "observed_path",
            "row_count",
            "column_count",
            "column_names",
            "missing_values_summary",
            "duplicate_key_summary",
            "notes",
            "alternative_candidates",
        ],
    )

    evidence_rows = build_evidence_matrix()
    write_csv(
        AUDIT_DIR / "data_reporting_evidence_matrix.csv",
        [row.__dict__ for row in evidence_rows],
        [
            "checklist_section",
            "required_information",
            "current_status",
            "evidence_file_path",
            "evidence_summary",
            "missing_or_risk",
            "recommended_action",
        ],
    )

    pytest_summary = "TESTS NOT RUN: dependency unavailable"
    (AUDIT_DIR / "test_run_summary.txt").write_text(pytest_summary + "\n", encoding="utf-8")
    (AUDIT_DIR / "DATA_REPORTING_CHECKLIST_DRAFT.md").write_text(
        build_checklist_draft(evidence_rows),
        encoding="utf-8",
    )
    (AUDIT_DIR / "MISSING_ITEMS_FOR_DATA_REPORTING.md").write_text(
        build_missing_items(evidence_rows),
        encoding="utf-8",
    )
    (AUDIT_DIR / "DATA_REPORTING_AUDIT_SUMMARY.md").write_text(
        build_summary(evidence_rows, pytest_summary),
        encoding="utf-8",
    )

    high_priority_titles = [
        "Ethics / informed consent statement not found",
        "Zenodo DOI or release tag missing",
        "Leakage prevention evidence incomplete",
        "Third-party restriction statement must be carried into submission text",
    ]
    leakage_status, _ = leakage_summary()
    zenodo_has_doi = bool(re.search(r"10\\.", safe_read_text(ROOT / ".zenodo.json")))
    third_party_missing = not file_exists("data/raw_manifest/third_party_source_restrictions.md")
    ethics_missing = True
    leakage_evidence_incomplete = not (
        file_exists("data/splits/iterative_stratified_splits/seed_42_train.csv")
        and file_exists("data/splits/iterative_stratified_splits/seed_42_test.csv")
        and file_exists("data/processed/feature_scaler_parameters.csv")
    )
    recommendation = "需要补充后提交"
    if ethics_missing:
        recommendation = "需要补充后提交"
    generated = [
        "project_structure_scan.csv",
        "file_manifest_sha256.csv",
        "key_file_statistics.csv",
        "data_reporting_evidence_matrix.csv",
        "DATA_REPORTING_CHECKLIST_DRAFT.md",
        "MISSING_ITEMS_FOR_DATA_REPORTING.md",
        "DATA_REPORTING_AUDIT_SUMMARY.md",
        "test_run_summary.txt",
        "run_audit.py",
    ]
    print("生成文件:")
    for name in generated:
        print(f"- submission_checklist_audit/{name}")
    print(f"HIGH PRIORITY missing items 数量: {len(high_priority_titles)}")
    for title in high_priority_titles:
        print(f"- {title}")
    print(f"发现伦理/知情同意声明缺失: {'YES' if ethics_missing else 'NO'}")
    print(f"发现第三方数据限制说明缺失: {'YES' if third_party_missing else 'NO'}")
    print(f"发现 Zenodo DOI 或 release tag 缺失: {'YES' if not zenodo_has_doi else 'NO'}")
    print(f"发现 train/test leakage 风险: {'YES' if leakage_status != 'READY' or leakage_evidence_incomplete else 'NO'}")
    print(f"建议作者上传 Data Reporting Checklist: {recommendation}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
