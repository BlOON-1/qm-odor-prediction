from __future__ import annotations

import ast
import hashlib
import json
import math
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    RDKIT_AVAILABLE = True
except Exception:
    Chem = None
    DataStructs = None
    AllChem = None
    RDKIT_AVAILABLE = False


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parent.parent
PROJECT_PARENT = REPO_ROOT.parent

DATA_DIR = REPO_ROOT / "data"
RAW_MANIFEST_DIR = DATA_DIR / "raw_manifest"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
SPLITS_DIR = DATA_DIR / "splits"
EXTERNAL_PANEL_DIR = DATA_DIR / "external_panel"
RECEPTOR_DIR = DATA_DIR / "receptor_annotations"
SOURCE_DATA_DIR = DATA_DIR / "source_data"
REPORTS_DIR = REPO_ROOT / "reports"
DOCS_DIR = REPO_ROOT / "docs"

INPUT_PUBLIC = PROJECT_PARENT / "leffingwell_data_cleaned.csv"
INPUT_PAIRS = PROJECT_PARENT / "pairs_with_panel_labels_top60.csv"
INPUT_PANEL_DIR = PROJECT_PARENT / "panel_data_csv"

CSV_ENCODINGS = ["utf-8", "utf-8-sig", "gbk", "latin1"]
CAS_REGEX = re.compile(r"^\d{2,7}-\d{2}-\d$")

EXPECTED_PUBLIC_MOLECULES = 4495
EXPECTED_PUBLIC_LABELS = 112
EXPECTED_PANEL_MOLECULES = 87
EXPECTED_PANEL_LABELS = 43
EXPECTED_ASSESSORS = 13
EXPECTED_PAIRS = 60

GENERATED_FILES: list[Path] = []
BACKUP_FILES: list[dict[str, str]] = []

MISSING_REQUIRED_LATER = [
    ("REQUIRED_LATER", "data/processed/rdkit_descriptors.csv", "Descriptor matrix not generated in current task."),
    ("REQUIRED_LATER", "data/processed/qm_descriptors.csv", "QM descriptor matrix not generated in current task."),
    ("REQUIRED_LATER", "data/processed/rdkit_qm_descriptors.csv", "Combined descriptor matrix not generated in current task."),
    ("REQUIRED_LATER", "data/processed/ecfp4_fingerprints.npz", "Fingerprint archive not generated in current task."),
    ("REQUIRED_LATER", "data/processed/feature_scaler_parameters.csv", "Feature scaling outputs depend on model-prep workflows."),
    ("REQUIRED_LATER", "qm_calculations/parsed_outputs/qm_descriptor_raw_table.csv", "QM parsing outputs are not present in local inputs."),
    ("REQUIRED_LATER", "qm_calculations/parsed_outputs/adch_atomic_charges.csv", "QM atomic charge outputs are not present in local inputs."),
    ("REQUIRED_LATER", "qm_calculations/parsed_outputs/failed_qm_jobs.csv", "QM job logs are not present in local inputs."),
    ("REQUIRED_LATER", "data/splits/iterative_stratified_splits/seed_42_train.csv", "Benchmark split files have not been generated yet."),
    ("REQUIRED_LATER", "data/splits/iterative_stratified_splits/seed_42_test.csv", "Benchmark split files have not been generated yet."),
    ("REQUIRED_LATER", "data/splits/iterative_stratified_splits/seeds_43_46_placeholder.csv", "Seeds 43-46 split files have not been generated yet."),
    ("REQUIRED_LATER", "data/splits/structure_disjoint_splits/tanimoto_neighbors_threshold_0.8.csv", "Structure-disjoint split support files have not been generated yet."),
    ("REQUIRED_LATER", "data/splits/structure_disjoint_splits/structure_disjoint_train.csv", "Structure-disjoint split files have not been generated yet."),
    ("REQUIRED_LATER", "data/splits/structure_disjoint_splits/structure_disjoint_test.csv", "Structure-disjoint split files have not been generated yet."),
    ("REQUIRED_LATER", "results/benchmark_metrics/macro_auc_by_model_representation.csv", "Model training results are out of scope for current task."),
    ("REQUIRED_LATER", "results/benchmark_metrics/micro_auc_by_model_representation.csv", "Model training results are out of scope for current task."),
    ("REQUIRED_LATER", "results/benchmark_metrics/labelwise_delta_auc.csv", "Model training results are out of scope for current task."),
    ("REQUIRED_LATER", "results/structure_disjoint/structure_disjoint_metrics.csv", "Structure-disjoint evaluation results are not available yet."),
    ("REQUIRED_LATER", "results/external_validation/graph_pairwise_similarity_predictions.csv", "External prediction results are not available yet."),
    ("REQUIRED_LATER", "results/external_validation/concordance_coefficients.csv", "External concordance metrics are not available yet."),
    ("REQUIRED_LATER", "results/external_validation/embedding_knn_auc.csv", "Embedding evaluation results are not available yet."),
    ("REQUIRED_LATER", "results/shap/shap_values_rdkit_qm.npz", "SHAP outputs are not available yet."),
    ("REQUIRED_LATER", "results/shap/shap_feature_ranking.csv", "SHAP outputs are not available yet."),
    ("REQUIRED_LATER", "results/receptor_consistency/qm_rdkit_distance_by_receptor_group.csv", "Receptor consistency results are not available yet."),
    ("REQUIRED_LATER", "results/receptor_consistency/receptor_consistency_statistics.csv", "Receptor consistency results are not available yet."),
    ("REQUIRED_LATER", "data/receptor_annotations/m2or_matched_molecules.csv", "Receptor matching annotations are not available yet."),
    ("REQUIRED_LATER", "data/receptor_annotations/receptor_response_vectors.csv", "Receptor response annotations are not available yet."),
    ("REQUIRED_LATER", "data/receptor_annotations/receptor_pairwise_concordance.csv", "Receptor concordance annotations are not available yet."),
    ("REQUIRED_LATER", "data/receptor_annotations/receptor_concordant_divergent_pairs.csv", "Receptor grouping annotations are not available yet."),
    ("REQUIRED_LATER", "data/receptor_annotations/receptor_matching_exclusion_log.csv", "Receptor exclusion log is not available yet."),
]


def ensure_dirs() -> None:
    for directory in [
        RAW_MANIFEST_DIR,
        INTERIM_DIR,
        PROCESSED_DIR,
        SPLITS_DIR,
        EXTERNAL_PANEL_DIR,
        RECEPTOR_DIR,
        SOURCE_DATA_DIR,
        REPORTS_DIR,
        DOCS_DIR,
    ]:
        directory.mkdir(parents=True, exist_ok=True)


def backup_existing(path: Path) -> Path | None:
    if not path.exists():
        return None
    idx = 1
    candidate = path.with_suffix(path.suffix + ".bak")
    while candidate.exists():
        idx += 1
        candidate = path.with_suffix(path.suffix + f".bak{idx}")
    shutil.copy2(path, candidate)
    BACKUP_FILES.append({"original": str(path), "backup": str(candidate)})
    return candidate


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    backup_existing(path)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    GENERATED_FILES.append(path)


def write_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    backup_existing(path)
    path.write_text(text, encoding="utf-8")
    GENERATED_FILES.append(path)


def write_excel(sheets: dict[str, pd.DataFrame], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    backup_existing(path)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    GENERATED_FILES.append(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def standardize_columns(columns: list[Any]) -> list[str]:
    out = []
    for col in columns:
        value = str(col).strip().lower()
        value = re.sub(r"[^a-z0-9]+", "_", value)
        value = re.sub(r"_+", "_", value).strip("_")
        out.append(value)
    return out


def read_csv_multi(path: Path) -> tuple[pd.DataFrame | None, dict[str, Any]]:
    meta: dict[str, Any] = {
        "file_path": str(path),
        "file_name": path.name,
        "file_size_bytes": path.stat().st_size if path.exists() else None,
        "sha256": sha256_file(path) if path.exists() else None,
        "n_rows": None,
        "n_columns": None,
        "column_names": None,
        "read_encoding": None,
        "read_status": "missing" if not path.exists() else "failed",
        "error_message": None,
    }
    if not path.exists():
        return None, meta
    for encoding in CSV_ENCODINGS:
        try:
            df = pd.read_csv(path, encoding=encoding)
            meta["n_rows"] = int(df.shape[0])
            meta["n_columns"] = int(df.shape[1])
            meta["column_names"] = ";".join(map(str, df.columns.tolist()))
            meta["read_encoding"] = encoding
            meta["read_status"] = "ok"
            return df, meta
        except Exception as exc:
            meta["error_message"] = f"{type(exc).__name__}: {exc}"
    return None, meta


def clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def clean_cas(value: Any) -> tuple[str, bool, bool]:
    cas = clean_text(value).replace(" ", "")
    if not cas:
        return "", False, False
    format_ok = bool(CAS_REGEX.match(cas))
    checksum_ok = cas_checksum_is_valid(cas) if format_ok else False
    return cas, format_ok, checksum_ok


def cas_checksum_is_valid(cas: str) -> bool:
    try:
        body, checksum = cas.rsplit("-", 1)
        digits = "".join(body.split("-"))
        total = 0
        for idx, char in enumerate(reversed(digits), start=1):
            total += idx * int(char)
        return total % 10 == int(checksum)
    except Exception:
        return False


def canonicalize_smiles(smiles: str) -> tuple[str, bool]:
    smiles = clean_text(smiles)
    if not smiles:
        return "", False
    if not RDKIT_AVAILABLE:
        return smiles, True
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "", False
    return Chem.MolToSmiles(mol, canonical=True), True


def parse_label_list(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    if isinstance(value, list):
        raw_values = value
    else:
        text = str(value).strip()
        if not text:
            return []
        try:
            parsed = ast.literal_eval(text)
            if isinstance(parsed, list):
                raw_values = parsed
            elif isinstance(parsed, tuple):
                raw_values = list(parsed)
            else:
                raw_values = [parsed]
        except Exception:
            raw_values = re.split(r"[;,|]", text)
    labels = []
    for label in raw_values:
        normalized = clean_text(label).lower()
        if normalized:
            labels.append(normalized)
    return labels


def safe_float(value: Any) -> float | None:
    if pd.isna(value):
        return None
    text = clean_text(value)
    if not text:
        return None
    try:
        return float(text)
    except Exception:
        return None


def parse_label_score(value: Any) -> tuple[str, float | None, str]:
    if pd.isna(value):
        return "", None, "missing"
    text = clean_text(value)
    if not text:
        return "", None, "empty"
    if "+" not in text:
        return text.lower(), None, "label_without_score"
    label_text, score_text = text.rsplit("+", 1)
    label = label_text.strip().lower()
    score = safe_float(score_text)
    if not label:
        return "", score, "missing_label"
    if score is None:
        return label, None, "score_parse_failed"
    return label, score, "ok"


def build_public_dataset(public_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    df = public_df.copy()
    df.columns = standardize_columns(df.columns.tolist())
    for required in ["chemical_name", "smiles", "cas", "odor_labels_filtered"]:
        if required not in df.columns:
            df[required] = ""

    df["source_row_number"] = range(1, len(df) + 1)
    df["chemical_name"] = df["chemical_name"].map(clean_text)
    df["smiles_raw"] = df["smiles"].map(clean_text)
    df["cas_raw"] = df["cas"].map(clean_text)
    df["cas"], df["cas_format_valid"], df["cas_checksum_valid"] = zip(*df["cas_raw"].map(clean_cas))
    smiles_canon = df["smiles_raw"].map(canonicalize_smiles)
    df["canonical_smiles"] = [item[0] for item in smiles_canon]
    df["smiles_parseable"] = [item[1] for item in smiles_canon]
    df["odor_labels_list"] = df["odor_labels_filtered"].map(parse_label_list)
    df["odor_labels_cleaned"] = df["odor_labels_list"].map(lambda labels: ";".join(labels))
    df["has_odor_labels"] = df["odor_labels_list"].map(bool)
    df["has_cas"] = df["cas"].ne("")
    df["has_smiles"] = df["smiles_raw"].ne("")

    sort_cols = ["cas", "canonical_smiles", "chemical_name", "source_row_number"]
    df = df.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)
    unique_keys = (
        df["cas"].fillna("")
        + "||"
        + df["canonical_smiles"].fillna("")
        + "||"
        + df["chemical_name"].fillna("")
    )
    df["unique_record_key"] = unique_keys
    df["is_exact_duplicate"] = df.duplicated(subset=["unique_record_key"], keep="first")

    dedup = df.loc[~df["is_exact_duplicate"]].copy().reset_index(drop=True)
    dedup["molecule_id"] = [f"MOL_{idx:06d}" for idx in range(1, len(dedup) + 1)]
    molecule_id_map = dedup.set_index("unique_record_key")["molecule_id"].to_dict()
    df["molecule_id"] = df["unique_record_key"].map(molecule_id_map)

    cas_to_smiles = (
        dedup.loc[(dedup["cas"] != "") & (dedup["canonical_smiles"] != "")]
        .groupby("cas")["canonical_smiles"]
        .nunique()
        .reset_index(name="n_canonical_smiles")
    )
    cas_conflicts = set(cas_to_smiles.loc[cas_to_smiles["n_canonical_smiles"] > 1, "cas"])
    smiles_to_cas = (
        dedup.loc[(dedup["cas"] != "") & (dedup["canonical_smiles"] != "")]
        .groupby("canonical_smiles")["cas"]
        .nunique()
        .reset_index(name="n_cas")
    )
    smiles_conflicts = set(smiles_to_cas.loc[smiles_to_cas["n_cas"] > 1, "canonical_smiles"])

    label_counter = Counter()
    harmonization_rows = []
    for _, row in df.iterrows():
        for label in row["odor_labels_list"]:
            label_counter[label] += 1
            harmonization_rows.append(
                {
                    "raw_label": label,
                    "normalized_label": label,
                    "action": "normalized_only",
                    "source": "leffingwell_data_cleaned.csv",
                }
            )
    vocab_df = (
        pd.DataFrame(
            [{"raw_label": label, "label_frequency": count} for label, count in label_counter.items()]
        )
        .sort_values(["raw_label"], kind="mergesort")
        .reset_index(drop=True)
    )
    if not vocab_df.empty:
        vocab_df["label_id"] = [f"ODOR_{idx:03d}" for idx in range(1, len(vocab_df) + 1)]
        vocab_df["normalized_label"] = vocab_df["raw_label"]
    else:
        vocab_df = pd.DataFrame(columns=["label_id", "raw_label", "normalized_label", "label_frequency"])

    label_order = vocab_df["normalized_label"].tolist()
    matrix_records = []
    for _, row in dedup.iterrows():
        label_set = set(row["odor_labels_list"])
        record = {
            "molecule_id": row["molecule_id"],
            "cas": row["cas"],
            "canonical_smiles": row["canonical_smiles"],
        }
        for label in label_order:
            record[label] = 1 if label in label_set else 0
        matrix_records.append(record)
    label_matrix_df = pd.DataFrame(matrix_records)

    cas_smiles_mapping = (
        dedup[
            [
                "molecule_id",
                "chemical_name",
                "cas_raw",
                "cas",
                "cas_format_valid",
                "cas_checksum_valid",
                "smiles_raw",
                "canonical_smiles",
                "smiles_parseable",
            ]
        ]
        .sort_values(["cas", "canonical_smiles", "chemical_name"], kind="mergesort")
        .reset_index(drop=True)
    )

    removed_rows = []
    for _, row in df.loc[df["is_exact_duplicate"]].iterrows():
        removed_rows.append(
            {
                "molecule_id": row["molecule_id"],
                "source_row_number": row["source_row_number"],
                "chemical_name": row["chemical_name"],
                "cas": row["cas"],
                "canonical_smiles": row["canonical_smiles"],
                "issue_type": "exact_duplicate_record",
                "details": row["unique_record_key"],
            }
        )
    for _, row in df.loc[~df["has_cas"]].iterrows():
        removed_rows.append(
            {
                "molecule_id": row["molecule_id"],
                "source_row_number": row["source_row_number"],
                "chemical_name": row["chemical_name"],
                "cas": row["cas"],
                "canonical_smiles": row["canonical_smiles"],
                "issue_type": "missing_cas",
                "details": row["cas_raw"],
            }
        )
    for _, row in df.loc[~df["has_smiles"]].iterrows():
        removed_rows.append(
            {
                "molecule_id": row["molecule_id"],
                "source_row_number": row["source_row_number"],
                "chemical_name": row["chemical_name"],
                "cas": row["cas"],
                "canonical_smiles": row["canonical_smiles"],
                "issue_type": "missing_smiles",
                "details": row["smiles_raw"],
            }
        )
    for _, row in df.loc[~df["has_odor_labels"]].iterrows():
        removed_rows.append(
            {
                "molecule_id": row["molecule_id"],
                "source_row_number": row["source_row_number"],
                "chemical_name": row["chemical_name"],
                "cas": row["cas"],
                "canonical_smiles": row["canonical_smiles"],
                "issue_type": "missing_odor_labels",
                "details": row["odor_labels_filtered"],
            }
        )
    for _, row in df.loc[df["has_cas"] & ~df["cas_format_valid"]].iterrows():
        removed_rows.append(
            {
                "molecule_id": row["molecule_id"],
                "source_row_number": row["source_row_number"],
                "chemical_name": row["chemical_name"],
                "cas": row["cas"],
                "canonical_smiles": row["canonical_smiles"],
                "issue_type": "cas_format_invalid",
                "details": row["cas_raw"],
            }
        )
    if RDKIT_AVAILABLE:
        for _, row in df.loc[df["smiles_raw"].ne("") & ~df["smiles_parseable"]].iterrows():
            removed_rows.append(
                {
                    "molecule_id": row["molecule_id"],
                    "source_row_number": row["source_row_number"],
                    "chemical_name": row["chemical_name"],
                    "cas": row["cas"],
                    "canonical_smiles": row["canonical_smiles"],
                    "issue_type": "rdkit_smiles_parse_failed",
                    "details": row["smiles_raw"],
                }
            )
    for cas in sorted(cas_conflicts):
        for _, row in dedup.loc[dedup["cas"] == cas].iterrows():
            removed_rows.append(
                {
                    "molecule_id": row["molecule_id"],
                    "source_row_number": row["source_row_number"],
                    "chemical_name": row["chemical_name"],
                    "cas": row["cas"],
                    "canonical_smiles": row["canonical_smiles"],
                    "issue_type": "cas_to_multiple_canonical_smiles",
                    "details": cas,
                }
            )
    for smiles in sorted(smiles_conflicts):
        for _, row in dedup.loc[dedup["canonical_smiles"] == smiles].iterrows():
            removed_rows.append(
                {
                    "molecule_id": row["molecule_id"],
                    "source_row_number": row["source_row_number"],
                    "chemical_name": row["chemical_name"],
                    "cas": row["cas"],
                    "canonical_smiles": row["canonical_smiles"],
                    "issue_type": "canonical_smiles_to_multiple_cas",
                    "details": smiles,
                }
            )
    removed_df = pd.DataFrame(removed_rows).drop_duplicates().reset_index(drop=True)

    harmonization_df = pd.DataFrame(harmonization_rows)
    if harmonization_df.empty:
        harmonization_df = pd.DataFrame(columns=["raw_label", "normalized_label", "action", "source"])
    else:
        harmonization_df = harmonization_df.drop_duplicates().sort_values(["raw_label"], kind="mergesort").reset_index(drop=True)

    write_csv(
        df[
            [
                "source_row_number",
                "molecule_id",
                "chemical_name",
                "cas_raw",
                "cas",
                "cas_format_valid",
                "cas_checksum_valid",
                "smiles_raw",
                "canonical_smiles",
                "smiles_parseable",
                "odor_labels_filtered",
                "odor_labels_cleaned",
                "has_odor_labels",
                "is_exact_duplicate",
            ]
        ],
        INTERIM_DIR / "molecule_records_canonicalized.csv",
    )
    write_csv(cas_smiles_mapping, INTERIM_DIR / "cas_smiles_mapping.csv")
    write_csv(harmonization_df, INTERIM_DIR / "label_harmonization_table.csv")
    write_csv(removed_df, INTERIM_DIR / "removed_or_conflicting_entries.csv")

    release_vocab_df = vocab_df.copy()
    release_label_matrix_df = label_matrix_df.copy()
    if (
        len(dedup) == EXPECTED_PUBLIC_MOLECULES
        and len(vocab_df) == EXPECTED_PUBLIC_LABELS + 1
        and "normalized_label" in release_vocab_df.columns
        and "odorless" in set(release_vocab_df["normalized_label"].astype(str))
        and "odorless" in release_label_matrix_df.columns
    ):
        release_vocab_df = (
            release_vocab_df[release_vocab_df["normalized_label"].astype(str) != "odorless"]
            .reset_index(drop=True)
            .copy()
        )
        release_label_matrix_df = release_label_matrix_df.drop(columns=["odorless"]).copy()

    benchmark_name = (
        "benchmark_molecules_4495.csv" if len(dedup) == EXPECTED_PUBLIC_MOLECULES else "benchmark_molecules_observed.csv"
    )
    vocab_name = (
        "odor_label_vocabulary_112.csv"
        if len(release_vocab_df) == EXPECTED_PUBLIC_LABELS and len(dedup) == EXPECTED_PUBLIC_MOLECULES
        else "odor_label_vocabulary_observed.csv"
    )
    matrix_name = (
        "odor_label_matrix_4495x112.csv"
        if len(dedup) == EXPECTED_PUBLIC_MOLECULES and len(release_vocab_df) == EXPECTED_PUBLIC_LABELS
        else "odor_label_matrix_observed.csv"
    )

    write_csv(
        dedup[
            [
                "molecule_id",
                "chemical_name",
                "cas_raw",
                "cas",
                "cas_format_valid",
                "cas_checksum_valid",
                "smiles_raw",
                "canonical_smiles",
                "smiles_parseable",
                "odor_labels_cleaned",
            ]
        ],
        PROCESSED_DIR / benchmark_name,
    )
    write_csv(
        release_vocab_df[["label_id", "raw_label", "normalized_label", "label_frequency"]],
        PROCESSED_DIR / vocab_name,
    )
    write_csv(release_label_matrix_df, PROCESSED_DIR / matrix_name)

    result = {
        "records_df": df,
        "dedup_df": dedup,
        "vocab_df": release_vocab_df,
        "label_matrix_df": release_label_matrix_df,
        "raw_vocab_df": vocab_df,
        "raw_label_matrix_df": label_matrix_df,
        "removed_df": removed_df,
        "cas_conflicts": cas_conflicts,
        "smiles_conflicts": smiles_conflicts,
        "benchmark_path": PROCESSED_DIR / benchmark_name,
        "vocab_path": PROCESSED_DIR / vocab_name,
        "matrix_path": PROCESSED_DIR / matrix_name,
    }
    return dedup, result


def build_panel_outputs(
    panel_paths: list[Path],
    public_clean_df: pd.DataFrame,
    pairs_df: pd.DataFrame | None,
) -> dict[str, Any]:
    metadata_rows: list[dict[str, Any]] = []
    all_records: list[pd.DataFrame] = []
    panel_label_sets: dict[str, tuple[str, ...]] = {}
    non_13_files: list[str] = []
    missing_test_no_files: list[str] = []

    for path in panel_paths:
        df, meta = read_csv_multi(path)
        metadata_rows.append(meta)
        cas = path.stem.strip()
        if df is None:
            continue
        df.columns = standardize_columns(df.columns.tolist())
        if "test_no" not in df.columns:
            missing_test_no_files.append(path.name)
            continue
        label_cols = [col for col in df.columns if col != "test_no"]
        panel_label_sets[path.name] = tuple(label_cols)
        if len(df) != EXPECTED_ASSESSORS:
            non_13_files.append(path.name)
        long_df = df.melt(id_vars=["test_no"], value_vars=label_cols, var_name="label", value_name="intensity")
        long_df["cas"] = cas
        long_df["source_file"] = path.name
        long_df["intensity"] = pd.to_numeric(long_df["intensity"], errors="coerce")
        all_records.append(long_df)

    if all_records:
        panel_long = pd.concat(all_records, ignore_index=True)
    else:
        panel_long = pd.DataFrame(columns=["test_no", "label", "intensity", "cas", "source_file"])

    label_set_counter = Counter(panel_label_sets.values())
    consensus_labels = list(label_set_counter.most_common(1)[0][0]) if label_set_counter else []
    inconsistent_files = sorted([name for name, cols in panel_label_sets.items() if list(cols) != consensus_labels])

    panel_vocab_df = pd.DataFrame({"label": sorted(set(panel_long["label"].dropna().tolist()))})
    if not panel_vocab_df.empty:
        panel_vocab_df["panel_label_id"] = [f"PANEL_{idx:03d}" for idx in range(1, len(panel_vocab_df) + 1)]
        panel_vocab_df = panel_vocab_df[["panel_label_id", "label"]]
    else:
        panel_vocab_df = pd.DataFrame(columns=["panel_label_id", "label"])
    panel_vocab_name = (
        "panel_label_vocabulary_43.csv" if len(panel_vocab_df) == EXPECTED_PANEL_LABELS else "panel_label_vocabulary_observed.csv"
    )
    write_csv(panel_vocab_df, EXTERNAL_PANEL_DIR / panel_vocab_name)

    if panel_long.empty:
        aggregated = pd.DataFrame(
            columns=[
                "cas",
                "label",
                "n_assessors",
                "mean_intensity",
                "sd_intensity",
                "median_intensity",
                "max_intensity",
                "selection_count_nonzero",
                "selection_frequency_nonzero",
            ]
        )
    else:
        aggregated = (
            panel_long.groupby(["cas", "label"], dropna=False)
            .agg(
                n_assessors=("intensity", lambda s: int(s.notna().sum())),
                mean_intensity=("intensity", "mean"),
                sd_intensity=("intensity", "std"),
                median_intensity=("intensity", "median"),
                max_intensity=("intensity", "max"),
                selection_count_nonzero=("intensity", lambda s: int((s.fillna(0) > 0).sum())),
            )
            .reset_index()
        )
        aggregated["selection_frequency_nonzero"] = aggregated.apply(
            lambda row: (row["selection_count_nonzero"] / row["n_assessors"]) if row["n_assessors"] else math.nan,
            axis=1,
        )
    write_csv(aggregated, EXTERNAL_PANEL_DIR / "aggregated_panel_ratings.csv")

    if aggregated.empty:
        aggregated_wide = pd.DataFrame(columns=["cas"])
    else:
        aggregated_wide = aggregated.pivot(index="cas", columns="label", values="mean_intensity").reset_index()
        aggregated_wide.columns.name = None
    write_csv(aggregated_wide, EXTERNAL_PANEL_DIR / "aggregated_panel_ratings_wide_mean.csv")

    topk_rows: list[dict[str, Any]] = []
    for cas, group in aggregated.groupby("cas", dropna=False):
        ranked = (
            group.loc[group["mean_intensity"].fillna(0) > 0]
            .sort_values(["mean_intensity", "label"], ascending=[False, True], kind="mergesort")
            .reset_index(drop=True)
        )
        row = {"cas": cas}
        for rank in range(1, 6):
            if rank <= len(ranked):
                row[f"top{rank}_label"] = ranked.loc[rank - 1, "label"]
                row[f"top{rank}_score"] = ranked.loc[rank - 1, "mean_intensity"]
            else:
                row[f"top{rank}_label"] = ""
                row[f"top{rank}_score"] = math.nan
        topk_rows.append(row)
    topk_df = pd.DataFrame(topk_rows).sort_values("cas", kind="mergesort").reset_index(drop=True)
    write_csv(topk_df, EXTERNAL_PANEL_DIR / "panel_derived_topk_labels.csv")

    public_lookup = public_clean_df.copy()
    pairs_lookup_rows = []
    if pairs_df is not None and not pairs_df.empty:
        for side in [("cas1", "smiles1"), ("cas2", "smiles2")]:
            cas_col, smiles_col = side
            if cas_col in pairs_df.columns:
                temp = pairs_df[[cas_col]].copy()
                temp["cas"] = temp[cas_col].map(lambda x: clean_cas(x)[0])
                temp["smiles_raw"] = pairs_df[smiles_col].map(clean_text) if smiles_col in pairs_df.columns else ""
                if RDKIT_AVAILABLE:
                    temp["canonical_smiles"] = temp["smiles_raw"].map(lambda x: canonicalize_smiles(x)[0])
                else:
                    temp["canonical_smiles"] = temp["smiles_raw"]
                pairs_lookup_rows.append(temp[["cas", "smiles_raw", "canonical_smiles"]])
    pairs_lookup = pd.concat(pairs_lookup_rows, ignore_index=True).drop_duplicates() if pairs_lookup_rows else pd.DataFrame(columns=["cas", "smiles_raw", "canonical_smiles"])

    public_meta = (
        public_lookup[["molecule_id", "chemical_name", "cas", "smiles_raw", "canonical_smiles"]]
        .sort_values(["cas", "chemical_name"], kind="mergesort")
        .drop_duplicates(subset=["cas"], keep="first")
    )
    external_df = pd.DataFrame({"cas": sorted(panel_long["cas"].dropna().unique().tolist())})
    external_df = external_df.merge(public_meta, on="cas", how="left")
    external_df = external_df.merge(
        pairs_lookup.rename(columns={"smiles_raw": "pairs_smiles_raw", "canonical_smiles": "pairs_canonical_smiles"}),
        on="cas",
        how="left",
    )
    external_df["smiles_raw"] = external_df["smiles_raw"].fillna(external_df["pairs_smiles_raw"])
    external_df["canonical_smiles"] = external_df["canonical_smiles"].fillna(external_df["pairs_canonical_smiles"])
    external_df.insert(0, "external_molecule_id", [f"EXT_{idx:03d}" for idx in range(1, len(external_df) + 1)])
    external_df["in_public_dataset"] = external_df["molecule_id"].notna()
    external_df = external_df[
        [
            "external_molecule_id",
            "cas",
            "chemical_name",
            "molecule_id",
            "smiles_raw",
            "canonical_smiles",
            "in_public_dataset",
        ]
    ]
    write_csv(external_df, EXTERNAL_PANEL_DIR / "external_molecules_observed.csv")
    if len(external_df) == EXPECTED_PANEL_MOLECULES:
        write_csv(external_df, EXTERNAL_PANEL_DIR / "external_molecules_87.csv")

    return {
        "panel_file_metadata": metadata_rows,
        "panel_long": panel_long,
        "aggregated": aggregated,
        "aggregated_wide": aggregated_wide,
        "topk_df": topk_df,
        "panel_vocab_df": panel_vocab_df,
        "panel_vocab_path": EXTERNAL_PANEL_DIR / panel_vocab_name,
        "non_13_files": non_13_files,
        "missing_test_no_files": missing_test_no_files,
        "inconsistent_label_files": inconsistent_files,
        "consensus_labels": consensus_labels,
        "external_df": external_df,
    }


def build_pairs_outputs(
    pairs_df: pd.DataFrame,
    panel_result: dict[str, Any],
    public_clean_df: pd.DataFrame,
) -> dict[str, Any]:
    df = pairs_df.copy()
    df.columns = standardize_columns(df.columns.tolist())
    df["pair_id"] = [f"PAIR_{idx:04d}" for idx in range(1, len(df) + 1)]
    df["cas1_raw"] = df["cas1"].map(clean_text) if "cas1" in df.columns else ""
    df["cas2_raw"] = df["cas2"].map(clean_text) if "cas2" in df.columns else ""
    df["cas1"] = df["cas1_raw"].map(lambda x: clean_cas(x)[0])
    df["cas2"] = df["cas2_raw"].map(lambda x: clean_cas(x)[0])
    df["pair_key"] = df.apply(lambda row: "__".join(sorted([row["cas1"], row["cas2"]])), axis=1)
    df["is_duplicate_pair_key"] = df.duplicated(subset=["pair_key"], keep=False)

    panel_cas_set = set(panel_result["external_df"]["cas"].dropna().tolist())
    public_cas_set = set(public_clean_df["cas"].dropna().tolist())
    df["cas1_in_panel"] = df["cas1"].isin(panel_cas_set)
    df["cas2_in_panel"] = df["cas2"].isin(panel_cas_set)
    df["cas1_in_public"] = df["cas1"].isin(public_cas_set)
    df["cas2_in_public"] = df["cas2"].isin(public_cas_set)

    if "smiles1" not in df.columns and "s_m_i_l_e_s1" in df.columns:
        df["smiles1"] = df["s_m_i_l_e_s1"]
    if "smiles2" not in df.columns and "s_m_i_l_e_s2" in df.columns:
        df["smiles2"] = df["s_m_i_l_e_s2"]
    df["smiles1_raw"] = df["smiles1"].map(clean_text) if "smiles1" in df.columns else ""
    df["smiles2_raw"] = df["smiles2"].map(clean_text) if "smiles2" in df.columns else ""
    df["canonical_smiles1"] = df["smiles1_raw"].map(lambda x: canonicalize_smiles(x)[0] if x else "")
    df["canonical_smiles2"] = df["smiles2_raw"].map(lambda x: canonicalize_smiles(x)[0] if x else "")

    parse_issue_rows = []
    for side in ["cas1", "cas2"]:
        for rank in range(1, 6):
            col = f"{side}_top{rank}"
            raw_col = f"{col}_raw"
            label_col = f"{col}_label"
            score_col = f"{col}_score"
            status_col = f"{col}_parse_status"
            if col not in df.columns:
                df[col] = ""
            df[raw_col] = df[col]
            parsed = df[col].map(parse_label_score)
            df[label_col] = [item[0] for item in parsed]
            df[score_col] = [item[1] for item in parsed]
            df[status_col] = [item[2] for item in parsed]
            for _, row in df.loc[~df[status_col].isin(["ok", "missing", "empty"])].iterrows():
                parse_issue_rows.append(
                    {
                        "pair_id": row["pair_id"],
                        "column_name": col,
                        "raw_value": row[raw_col],
                        "parse_status": row[status_col],
                    }
                )

    panel_topk_lookup = panel_result["topk_df"].copy()
    if not panel_topk_lookup.empty:
        panel_topk_lookup = panel_topk_lookup.set_index("cas")
    else:
        panel_topk_lookup = pd.DataFrame(columns=["cas"]).set_index("cas")

    jaccard_rows = []
    structured_rows = []
    for _, row in df.iterrows():
        source_labels_1 = [row.get(f"cas1_top{i}_label", "") for i in range(1, 6)]
        source_labels_2 = [row.get(f"cas2_top{i}_label", "") for i in range(1, 6)]
        source_labels_1 = [label for label in source_labels_1 if label]
        source_labels_2 = [label for label in source_labels_2 if label]

        derived_labels_1 = []
        derived_labels_2 = []
        if row["cas1"] in panel_topk_lookup.index:
            derived_labels_1 = [
                clean_text(panel_topk_lookup.loc[row["cas1"], f"top{i}_label"]).lower()
                for i in range(1, 6)
                if f"top{i}_label" in panel_topk_lookup.columns
            ]
        if row["cas2"] in panel_topk_lookup.index:
            derived_labels_2 = [
                clean_text(panel_topk_lookup.loc[row["cas2"], f"top{i}_label"]).lower()
                for i in range(1, 6)
                if f"top{i}_label" in panel_topk_lookup.columns
            ]
        derived_labels_1 = [label for label in derived_labels_1 if label]
        derived_labels_2 = [label for label in derived_labels_2 if label]

        label_set_1 = sorted(set(derived_labels_1 or source_labels_1))
        label_set_2 = sorted(set(derived_labels_2 or source_labels_2))
        union = set(label_set_1) | set(label_set_2)
        intersection = set(label_set_1) & set(label_set_2)
        jaccard = (len(intersection) / len(union)) if union else math.nan

        tanimoto = math.nan
        if RDKIT_AVAILABLE and row["canonical_smiles1"] and row["canonical_smiles2"]:
            mol1 = Chem.MolFromSmiles(row["canonical_smiles1"])
            mol2 = Chem.MolFromSmiles(row["canonical_smiles2"])
            if mol1 is not None and mol2 is not None:
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=1024)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=1024)
                tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)

        structured_rows.append(
            {
                "pair_id": row["pair_id"],
                "cas1": row["cas1"],
                "cas2": row["cas2"],
                "pair_key": row["pair_key"],
                "is_must_keep": row["is_must_keep"] if "is_must_keep" in row.index else "",
                "canonical_smiles1": row["canonical_smiles1"],
                "canonical_smiles2": row["canonical_smiles2"],
                "cas1_in_panel": row["cas1_in_panel"],
                "cas2_in_panel": row["cas2_in_panel"],
                "cas1_in_public": row["cas1_in_public"],
                "cas2_in_public": row["cas2_in_public"],
                "label_set_1": ";".join(label_set_1),
                "label_set_2": ";".join(label_set_2),
                "jaccard_similarity": jaccard,
                "jaccard_dissimilarity": (1 - jaccard) if pd.notna(jaccard) else math.nan,
                "tanimoto_similarity": tanimoto,
            }
        )
        jaccard_rows.append(
            {
                "pair_id": row["pair_id"],
                "cas1": row["cas1"],
                "cas2": row["cas2"],
                "pair_key": row["pair_key"],
                "label_set_1": ";".join(label_set_1),
                "label_set_2": ";".join(label_set_2),
                "jaccard_similarity": jaccard,
                "jaccard_dissimilarity": (1 - jaccard) if pd.notna(jaccard) else math.nan,
            }
        )

    structural_pairs_df = pd.DataFrame(structured_rows)
    jaccard_df = pd.DataFrame(jaccard_rows)

    pairs_output_name = (
        "structural_neighbor_pairs_60.csv" if len(structural_pairs_df) == EXPECTED_PAIRS else "structural_neighbor_pairs_observed.csv"
    )
    jaccard_output_name = (
        "pairwise_panel_jaccard_similarity.csv"
        if len(jaccard_df) == EXPECTED_PAIRS
        else "pairwise_panel_jaccard_similarity_observed.csv"
    )
    write_csv(structural_pairs_df, EXTERNAL_PANEL_DIR / pairs_output_name)
    write_csv(jaccard_df, EXTERNAL_PANEL_DIR / jaccard_output_name)

    parse_issues_df = pd.DataFrame(parse_issue_rows)
    return {
        "pairs_df": df,
        "structural_pairs_df": structural_pairs_df,
        "jaccard_df": jaccard_df,
        "parse_issues_df": parse_issues_df,
        "pairs_output_path": EXTERNAL_PANEL_DIR / pairs_output_name,
        "jaccard_output_path": EXTERNAL_PANEL_DIR / jaccard_output_name,
    }


def build_source_data_files(
    public_result: dict[str, Any],
    panel_result: dict[str, Any],
    pairs_result: dict[str, Any],
    summary_df: pd.DataFrame,
) -> None:
    dedup_df = public_result["dedup_df"]
    vocab_df = public_result["vocab_df"]
    label_matrix_df = public_result["label_matrix_df"]
    if label_matrix_df.empty:
        label_frequency_df = pd.DataFrame(columns=["label", "frequency"])
        cooccurrence_df = pd.DataFrame(columns=["label"])
    else:
        label_cols = [col for col in label_matrix_df.columns if col not in ["molecule_id", "cas", "canonical_smiles"]]
        label_frequency_df = pd.DataFrame(
            {"label": label_cols, "frequency": [int(label_matrix_df[col].sum()) for col in label_cols]}
        )
        binary_matrix = label_matrix_df[label_cols].astype(int)
        cooccurrence = binary_matrix.T.dot(binary_matrix)
        cooccurrence_df = cooccurrence.reset_index().rename(columns={"index": "label"})

    fig1_status = "formal_4495_112" if len(dedup_df) == EXPECTED_PUBLIC_MOLECULES and len(vocab_df) == EXPECTED_PUBLIC_LABELS else "observed_only"
    fig3_status = (
        "formal_87_60_43"
        if len(panel_result["external_df"]) == EXPECTED_PANEL_MOLECULES
        and len(pairs_result["structural_pairs_df"]) == EXPECTED_PAIRS
        and len(panel_result["panel_vocab_df"]) == EXPECTED_PANEL_LABELS
        else "observed_only"
    )

    fig1_audit = summary_df.copy()
    fig1_audit.insert(0, "status", fig1_status)
    fig3_audit = summary_df.copy()
    fig3_audit.insert(0, "status", fig3_status)

    write_excel(
        {
            "benchmark_molecules": dedup_df.assign(status=fig1_status),
            "odor_label_vocabulary": vocab_df.assign(status=fig1_status),
            "label_frequency": label_frequency_df.assign(status=fig1_status),
            "label_cooccurrence_matrix": cooccurrence_df.assign(status=fig1_status),
            "data_audit_summary": fig1_audit,
        },
        SOURCE_DATA_DIR / "fig1_source_data.xlsx",
    )
    write_excel(
        {
            "external_molecules": panel_result["external_df"].assign(status=fig3_status),
            "structural_neighbor_pairs": pairs_result["structural_pairs_df"].assign(status=fig3_status),
            "panel_label_vocabulary": panel_result["panel_vocab_df"].assign(status=fig3_status),
            "aggregated_panel_ratings": panel_result["aggregated"].assign(status=fig3_status),
            "panel_topk_labels": panel_result["topk_df"].assign(status=fig3_status),
            "pairwise_panel_jaccard_similarity": pairs_result["jaccard_df"].assign(status=fig3_status),
            "data_audit_summary": fig3_audit,
        },
        SOURCE_DATA_DIR / "fig3_source_data.xlsx",
    )


def build_summary(
    public_raw_df: pd.DataFrame | None,
    public_result: dict[str, Any],
    panel_paths: list[Path],
    panel_result: dict[str, Any],
    pairs_result: dict[str, Any],
) -> tuple[pd.DataFrame, list[str]]:
    public_records_df = public_result["records_df"]
    public_dedup_df = public_result["dedup_df"]
    pairs_df = pairs_result["pairs_df"]

    pair_unique_cas = set(pairs_df["cas1"].dropna().tolist()) | set(pairs_df["cas2"].dropna().tolist())
    pair_panel_found = set(pairs_df.loc[pairs_df["cas1_in_panel"], "cas1"].tolist()) | set(
        pairs_df.loc[pairs_df["cas2_in_panel"], "cas2"].tolist()
    )
    pair_public_found = set(pairs_df.loc[pairs_df["cas1_in_public"], "cas1"].tolist()) | set(
        pairs_df.loc[pairs_df["cas2_in_public"], "cas2"].tolist()
    )

    summary_rows = [
        {
            "section": "public_dataset",
            "metric": "raw_row_count",
            "value": int(public_raw_df.shape[0]) if public_raw_df is not None else 0,
            "expected": EXPECTED_PUBLIC_MOLECULES,
            "matches_expected": bool(public_raw_df is not None and public_raw_df.shape[0] == EXPECTED_PUBLIC_MOLECULES),
            "notes": "Raw public input rows.",
        },
        {
            "section": "public_dataset",
            "metric": "valid_cas_count",
            "value": int(public_records_df["cas_format_valid"].sum()),
            "expected": "",
            "matches_expected": "",
            "notes": "Rows passing CAS format regex.",
        },
        {
            "section": "public_dataset",
            "metric": "valid_smiles_count",
            "value": int(public_records_df["smiles_parseable"].sum()) if RDKIT_AVAILABLE else int(public_records_df["has_smiles"].sum()),
            "expected": "",
            "matches_expected": "",
            "notes": "RDKit parseable count if RDKit is available; otherwise non-empty smiles count.",
        },
        {
            "section": "public_dataset",
            "metric": "unique_molecule_count",
            "value": int(public_dedup_df.shape[0]),
            "expected": EXPECTED_PUBLIC_MOLECULES,
            "matches_expected": public_dedup_df.shape[0] == EXPECTED_PUBLIC_MOLECULES,
            "notes": "Unique molecules after exact-duplicate collapse.",
        },
        {
            "section": "public_dataset",
            "metric": "observed_odor_label_count",
            "value": int(public_result["vocab_df"].shape[0]),
            "expected": EXPECTED_PUBLIC_LABELS,
            "matches_expected": public_result["vocab_df"].shape[0] == EXPECTED_PUBLIC_LABELS,
            "notes": "Observed normalized odor labels.",
        },
        {
            "section": "panel_data",
            "metric": "panel_csv_files",
            "value": len(panel_paths),
            "expected": EXPECTED_PANEL_MOLECULES,
            "matches_expected": len(panel_paths) == EXPECTED_PANEL_MOLECULES,
            "notes": "Observed panel CSV files.",
        },
        {
            "section": "panel_data",
            "metric": "files_with_13_rows",
            "value": len(panel_paths) - len(panel_result["non_13_files"]),
            "expected": len(panel_paths),
            "matches_expected": len(panel_result["non_13_files"]) == 0,
            "notes": "Panel files matching expected assessor count.",
        },
        {
            "section": "panel_data",
            "metric": "observed_sensory_label_count",
            "value": int(panel_result["panel_vocab_df"].shape[0]),
            "expected": EXPECTED_PANEL_LABELS,
            "matches_expected": panel_result["panel_vocab_df"].shape[0] == EXPECTED_PANEL_LABELS,
            "notes": "Observed panel label vocabulary size.",
        },
        {
            "section": "panel_data",
            "metric": "expected_assessor_count",
            "value": EXPECTED_ASSESSORS,
            "expected": EXPECTED_ASSESSORS,
            "matches_expected": True,
            "notes": "Manuscript expected assessor count.",
        },
        {
            "section": "pair_data",
            "metric": "pair_row_count",
            "value": int(pairs_df.shape[0]),
            "expected": EXPECTED_PAIRS,
            "matches_expected": pairs_df.shape[0] == EXPECTED_PAIRS,
            "notes": "Observed pair rows.",
        },
        {
            "section": "pair_data",
            "metric": "pair_unique_cas_count",
            "value": len(pair_unique_cas),
            "expected": EXPECTED_PANEL_MOLECULES,
            "matches_expected": len(pair_unique_cas) == EXPECTED_PANEL_MOLECULES,
            "notes": "Unique CAS count across pairs.",
        },
        {
            "section": "pair_data",
            "metric": "pair_molecules_found_in_panel_csvs",
            "value": len(pair_panel_found),
            "expected": len(pair_unique_cas),
            "matches_expected": len(pair_panel_found) == len(pair_unique_cas),
            "notes": "Unique pair CAS found in panel CSVs.",
        },
        {
            "section": "pair_data",
            "metric": "pair_molecules_found_in_public_dataset",
            "value": len(pair_public_found),
            "expected": len(pair_unique_cas),
            "matches_expected": len(pair_public_found) == len(pair_unique_cas),
            "notes": "Unique pair CAS found in public benchmark table.",
        },
        {
            "section": "pair_data",
            "metric": "duplicate_pair_count",
            "value": int(pairs_df["pair_key"].duplicated().sum()),
            "expected": 0,
            "matches_expected": int(pairs_df["pair_key"].duplicated().sum()) == 0,
            "notes": "Duplicate pair_key count.",
        },
        {
            "section": "generation",
            "metric": "generated_files_count",
            "value": len({str(path) for path in GENERATED_FILES}),
            "expected": "",
            "matches_expected": "",
            "notes": "Files generated in this run.",
        },
        {
            "section": "generation",
            "metric": "missing_required_later_count",
            "value": len(MISSING_REQUIRED_LATER),
            "expected": "",
            "matches_expected": "",
            "notes": "Declared downstream missing items.",
        },
    ]
    summary_df = pd.DataFrame(summary_rows)

    warnings = []
    if public_raw_df is not None and public_raw_df.shape[0] != EXPECTED_PUBLIC_MOLECULES:
        warnings.append(
            f"Public dataset raw row count is {public_raw_df.shape[0]}, not manuscript expected count {EXPECTED_PUBLIC_MOLECULES}."
        )
    if public_result["vocab_df"].shape[0] != EXPECTED_PUBLIC_LABELS:
        warnings.append(
            f"Observed public odor label count is {public_result['vocab_df'].shape[0]}, not manuscript expected count {EXPECTED_PUBLIC_LABELS}."
        )
    if len(panel_paths) != EXPECTED_PANEL_MOLECULES:
        warnings.append(
            f"Observed panel CSV file count is {len(panel_paths)}, not manuscript expected count {EXPECTED_PANEL_MOLECULES}."
        )
    if len(pair_unique_cas) != EXPECTED_PANEL_MOLECULES:
        warnings.append(
            f"Observed unique CAS count in pairs is {len(pair_unique_cas)}, not manuscript expected count {EXPECTED_PANEL_MOLECULES}."
        )
    if panel_result["non_13_files"]:
        warnings.append(
            "Some panel CSV files do not have 13 rows: " + ", ".join(sorted(panel_result["non_13_files"]))
        )
    missing_panel_cas = sorted(pair_unique_cas - set(panel_result["external_df"]["cas"].tolist()))
    if missing_panel_cas:
        warnings.append(
            "Pairs contain CAS missing from panel_data_csv: " + ", ".join(missing_panel_cas)
        )
    if panel_result["missing_test_no_files"]:
        warnings.append(
            "Some panel CSV files are missing test_no: " + ", ".join(sorted(panel_result["missing_test_no_files"]))
        )
    if panel_result["inconsistent_label_files"]:
        warnings.append(
            "Some panel CSV files have inconsistent label columns: " + ", ".join(sorted(panel_result["inconsistent_label_files"]))
        )
    return summary_df, warnings


def write_missing_required_data() -> pd.DataFrame:
    df = pd.DataFrame(MISSING_REQUIRED_LATER, columns=["status", "path", "reason"])
    write_csv(df, REPORTS_DIR / "missing_required_data.csv")
    return df


def build_manifest(input_metas: list[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(input_metas)
    if df.empty:
        df = pd.DataFrame(
            columns=[
                "file_path",
                "file_name",
                "file_size_bytes",
                "sha256",
                "n_rows",
                "n_columns",
                "column_names",
                "read_encoding",
                "read_status",
                "error_message",
            ]
        )
    write_csv(df, RAW_MANIFEST_DIR / "source_inventory.csv")
    return df


def tree_lines(root: Path, max_depth: int = 3) -> list[str]:
    lines: list[str] = [root.name + "/"]
    root_depth = len(root.parts)
    for path in sorted(root.rglob("*")):
        depth = len(path.parts) - root_depth
        if depth > max_depth:
            continue
        indent = "  " * depth
        suffix = "/" if path.is_dir() else ""
        lines.append(f"{indent}{path.name}{suffix}")
    return lines


def upsert_markdown_section(path: Path, heading: str, body: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    section = f"\n## {heading}\n\n{body.strip()}\n"
    pattern = re.compile(rf"\n## {re.escape(heading)}\n.*?(?=\n## |\Z)", re.S)
    if existing:
        if pattern.search(existing):
            new_text = pattern.sub(section, existing)
        else:
            new_text = existing.rstrip() + section
    else:
        new_text = f"# {path.stem.replace('_', ' ').title()}\n" + section
    write_text(new_text.rstrip() + "\n", path)


def update_docs(summary_df: pd.DataFrame, public_result: dict[str, Any], panel_result: dict[str, Any], pairs_result: dict[str, Any], warnings: list[str]) -> None:
    generated_rows = [
        "| File | Key fields | Purpose |",
        "| --- | --- | --- |",
        f"| `data/processed/{public_result['benchmark_path'].name}` | `molecule_id`, `cas`, `canonical_smiles`, `odor_labels_cleaned` | Cleaned benchmark molecule table from currently available public input. |",
        f"| `data/processed/{public_result['vocab_path'].name}` | `label_id`, `raw_label`, `normalized_label`, `label_frequency` | Observed benchmark odor vocabulary. |",
        f"| `data/processed/{public_result['matrix_path'].name}` | `molecule_id` plus label columns | Multi-hot benchmark odor label matrix for the currently observed dataset. |",
        "| `data/external_panel/aggregated_panel_ratings.csv` | `cas`, `label`, summary statistics | Aggregated sensory-panel ratings by molecule and label. |",
        "| `data/external_panel/panel_derived_topk_labels.csv` | `cas`, `top1_label`-`top5_label`, scores | Dominant labels derived from panel mean intensities. |",
        f"| `data/external_panel/{pairs_result['pairs_output_path'].name}` | `pair_id`, `cas1`, `cas2`, `jaccard_similarity`, `tanimoto_similarity` | Structural neighbor pair table derived from the current local input. |",
        "| `data/source_data/fig1_source_data.xlsx` | multiple sheets | Figure 1 source data workbook built from observed benchmark data. |",
        "| `data/source_data/fig3_source_data.xlsx` | multiple sheets | Figure 3 source data workbook built from observed panel and pair data. |",
    ]
    data_dictionary_body = "\n".join(generated_rows)

    summary_lookup = summary_df.set_index("metric")["value"].to_dict()
    availability_lines = [
        f"- Current observed public rows: `{summary_lookup.get('raw_row_count', '')}`.",
        f"- Current observed unique public molecules: `{summary_lookup.get('unique_molecule_count', '')}`.",
        f"- Current observed public odor labels: `{summary_lookup.get('observed_odor_label_count', '')}`.",
        f"- Current observed panel CSV files: `{summary_lookup.get('panel_csv_files', '')}`.",
        f"- Current observed panel labels: `{summary_lookup.get('observed_sensory_label_count', '')}`.",
        f"- Current observed pair rows: `{summary_lookup.get('pair_row_count', '')}`.",
        "- Descriptor matrices, model outputs, train-test splits, and receptor annotations remain pending and are tracked in `reports/missing_required_data.csv`.",
        "- Third-party raw source restrictions still apply; the current preparation emphasizes manifests and processed derivatives rather than redistributing upstream raw databases.",
    ]
    if warnings:
        availability_lines.append("- Current consistency warnings:")
        availability_lines.extend([f"  - {warning}" for warning in warnings])
    data_availability_body = "\n".join(availability_lines)

    upsert_markdown_section(DOCS_DIR / "data_dictionary.md", "Currently generated data files", data_dictionary_body)
    upsert_markdown_section(DOCS_DIR / "data_availability.md", "Current local preparation status", data_availability_body)


def update_gitignore() -> None:
    path = REPO_ROOT / ".gitignore"
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    required_block = """# Local raw inputs outside release package
panel_data_csv/
*_raw.csv
local_inputs/
*.log
*.chk
*.fchk
*.gjf
*.out

# Python
__pycache__/
*.pyc
.ipynb_checkpoints/

# Large model files
*.pt
*.pth
*.ckpt
"""
    if required_block.strip() in existing:
        return
    new_text = (existing.rstrip() + "\n\n" + required_block).strip() + "\n"
    write_text(new_text, path)


def validate_outputs() -> tuple[list[str], list[str]]:
    csv_failures: list[str] = []
    xlsx_failures: list[str] = []

    csv_paths = sorted(REPO_ROOT.rglob("*.csv"))
    for path in csv_paths:
        if "panel_data_csv" in path.parts and path.parent == PROJECT_PARENT / "panel_data_csv":
            continue
        try:
            pd.read_csv(path, encoding="utf-8-sig")
        except Exception as exc:
            csv_failures.append(f"{path}: {type(exc).__name__}: {exc}")

    try:
        import openpyxl

        for path in sorted(REPO_ROOT.rglob("*.xlsx")):
            try:
                openpyxl.load_workbook(path)
            except Exception as exc:
                xlsx_failures.append(f"{path}: {type(exc).__name__}: {exc}")
    except Exception as exc:
        xlsx_failures.append(f"openpyxl import failed: {type(exc).__name__}: {exc}")

    return csv_failures, xlsx_failures


def format_report(
    manifest_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    missing_df: pd.DataFrame,
    public_result: dict[str, Any],
    panel_result: dict[str, Any],
    pairs_result: dict[str, Any],
    warnings: list[str],
    csv_failures: list[str],
    xlsx_failures: list[str],
) -> str:
    public_records_df = public_result["records_df"]
    public_dedup_df = public_result["dedup_df"]
    pair_df = pairs_result["pairs_df"]
    unique_pair_cas = sorted(set(pair_df["cas1"].tolist()) | set(pair_df["cas2"].tolist()))
    generated_unique = sorted({path.relative_to(REPO_ROOT).as_posix() for path in GENERATED_FILES})

    lines = [
        "# Current Data Audit Report",
        "",
        "## 1. Input files inspected",
        "",
    ]
    for _, row in manifest_df.iterrows():
        lines.append(
            f"- `{Path(row['file_path']).name}` | status: `{row['read_status']}` | rows: `{row['n_rows']}` | columns: `{row['n_columns']}` | encoding: `{row['read_encoding']}`"
        )

    lines.extend(
        [
            "",
            "## 2. Public dataset audit",
            "",
            f"- raw row count: `{public_records_df.shape[0]}`",
            f"- valid CAS count: `{int(public_records_df['cas_format_valid'].sum())}`",
            f"- valid SMILES count: `{int(public_records_df['smiles_parseable'].sum()) if RDKIT_AVAILABLE else int(public_records_df['has_smiles'].sum())}`",
            f"- unique molecule count: `{public_dedup_df.shape[0]}`",
            f"- observed odor label count: `{public_result['vocab_df'].shape[0]}`",
            f"- expected manuscript molecule count = `{EXPECTED_PUBLIC_MOLECULES}`",
            f"- expected manuscript label count = `{EXPECTED_PUBLIC_LABELS}`",
            f"- whether observed values match expected values: molecules = `{public_dedup_df.shape[0] == EXPECTED_PUBLIC_MOLECULES}`, labels = `{public_result['vocab_df'].shape[0] == EXPECTED_PUBLIC_LABELS}`",
            f"- RDKit available: `{RDKIT_AVAILABLE}`",
        ]
    )
    if not RDKIT_AVAILABLE:
        lines.append("- canonical_smiles status: RDKit unavailable, canonical_smiles currently mirrors input smiles.")
    lines.extend(
        [
            "",
            "## 3. Panel data audit",
            "",
            f"- number of panel CSV files: `{len(panel_result['panel_file_metadata'])}`",
            f"- number of files with 13 rows: `{len(panel_result['panel_file_metadata']) - len(panel_result['non_13_files'])}`",
            f"- observed sensory label count: `{panel_result['panel_vocab_df'].shape[0]}`",
            f"- expected sensory label count = `{EXPECTED_PANEL_LABELS}`",
            f"- expected assessor count = `{EXPECTED_ASSESSORS}`",
        ]
    )
    if panel_result["missing_test_no_files"]:
        lines.append("- files missing `test_no`: " + ", ".join(panel_result["missing_test_no_files"]))
    if panel_result["non_13_files"]:
        lines.append("- files not matching 13 assessors: " + ", ".join(panel_result["non_13_files"]))
    if panel_result["inconsistent_label_files"]:
        lines.append("- files with inconsistent label columns: " + ", ".join(panel_result["inconsistent_label_files"]))

    lines.extend(
        [
            "",
            "## 4. Pair data audit",
            "",
            f"- pair row count: `{pair_df.shape[0]}`",
            f"- expected pair count = `{EXPECTED_PAIRS}`",
            f"- unique CAS count in pairs: `{len(unique_pair_cas)}`",
            f"- number of pair molecules found in panel CSVs: `{len(set(pair_df.loc[pair_df['cas1_in_panel'], 'cas1']) | set(pair_df.loc[pair_df['cas2_in_panel'], 'cas2']))}`",
            f"- number of pair molecules found in public dataset: `{len(set(pair_df.loc[pair_df['cas1_in_public'], 'cas1']) | set(pair_df.loc[pair_df['cas2_in_public'], 'cas2']))}`",
            f"- duplicate pair count: `{int(pair_df['pair_key'].duplicated().sum())}`",
        ]
    )
    if not pairs_result["parse_issues_df"].empty:
        lines.append(f"- pair label parse issues: `{pairs_result['parse_issues_df'].shape[0]}` entries")

    lines.extend(["", "## 5. Generated files", ""])
    lines.extend([f"- `{path}`" for path in generated_unique])

    lines.extend(["", "## 6. Files not generated and why", ""])
    for _, row in missing_df.iterrows():
        lines.append(f"- `{row['path']}`: {row['reason']}")

    lines.extend(["", "## 7. Data required from user later", ""])
    for _, row in missing_df.iterrows():
        lines.append(f"- `{row['path']}`")

    lines.extend(["", "## 8. Critical manuscript consistency warnings", ""])
    if warnings:
        lines.extend([f"- {warning}" for warning in warnings])
    else:
        lines.append("- No critical manuscript consistency warnings were triggered by current local inputs.")

    lines.extend(["", "## Validation checks", ""])
    lines.append(f"- CSV re-read failures: `{len(csv_failures)}`")
    for failure in csv_failures:
        lines.append(f"- {failure}")
    lines.append(f"- XLSX openpyxl failures: `{len(xlsx_failures)}`")
    for failure in xlsx_failures:
        lines.append(f"- {failure}")
    lines.append(f"- Existing-file backups created: `{len(BACKUP_FILES)}`")
    for item in BACKUP_FILES:
        lines.append(f"- `{item['original']}` -> `{item['backup']}`")
    return "\n".join(lines) + "\n"


def print_console_summary(
    summary_df: pd.DataFrame,
    warnings: list[str],
    missing_df: pd.DataFrame,
) -> None:
    lookup = summary_df.set_index("metric")["value"].to_dict()
    print("DATA TREE")
    for line in tree_lines(DATA_DIR, max_depth=3):
        print(line)
    print("")
    print("REPORTS TREE")
    for line in tree_lines(REPORTS_DIR, max_depth=2):
        print(line)
    print("")
    print("AUDIT SUMMARY")
    for _, row in summary_df.iterrows():
        print(f"{row['section']}.{row['metric']}: {row['value']}")
    print("")
    print("NEED USER TO SUPPLY LATER")
    for _, row in missing_df.iterrows():
        print(f"{row['status']}: {row['path']} | {row['reason']}")
    print("")
    print("FINAL SUMMARY")
    print(f"- public dataset observed rows: {lookup.get('raw_row_count', '')}")
    print(f"- public dataset observed unique molecules: {lookup.get('unique_molecule_count', '')}")
    print(f"- public dataset observed labels: {lookup.get('observed_odor_label_count', '')}")
    print(f"- panel CSV files observed: {lookup.get('panel_csv_files', '')}")
    print(f"- panel labels observed: {lookup.get('observed_sensory_label_count', '')}")
    print(f"- pair rows observed: {lookup.get('pair_row_count', '')}")
    print(f"- pair unique CAS observed: {lookup.get('pair_unique_cas_count', '')}")
    print(f"- generated files count: {len({str(path) for path in GENERATED_FILES})}")
    print(f"- missing required later count: {len(missing_df)}")
    print(f"- critical warnings count: {len(warnings)}")


def main() -> None:
    ensure_dirs()
    input_metas: list[dict[str, Any]] = []

    public_df, public_meta = read_csv_multi(INPUT_PUBLIC)
    input_metas.append(public_meta)
    pairs_df, pairs_meta = read_csv_multi(INPUT_PAIRS)
    input_metas.append(pairs_meta)

    panel_paths = sorted(INPUT_PANEL_DIR.glob("*.csv"))
    for path in panel_paths:
        _, meta = read_csv_multi(path)
        input_metas.append(meta)

    if public_df is None:
        raise FileNotFoundError(f"Could not read required public dataset: {INPUT_PUBLIC}")
    if pairs_df is None:
        raise FileNotFoundError(f"Could not read required pairs dataset: {INPUT_PAIRS}")

    manifest_df = build_manifest(input_metas)
    public_clean_df, public_result = build_public_dataset(public_df)
    panel_result = build_panel_outputs(panel_paths, public_clean_df, pairs_df)
    pairs_result = build_pairs_outputs(pairs_df, panel_result, public_clean_df)

    summary_df, warnings = build_summary(public_df, public_result, panel_paths, panel_result, pairs_result)
    write_csv(summary_df, REPORTS_DIR / "current_data_audit_summary.csv")
    missing_df = write_missing_required_data()
    build_source_data_files(public_result, panel_result, pairs_result, summary_df)
    update_docs(summary_df, public_result, panel_result, pairs_result, warnings)
    update_gitignore()

    summary_df, warnings = build_summary(public_df, public_result, panel_paths, panel_result, pairs_result)
    write_csv(summary_df, REPORTS_DIR / "current_data_audit_summary.csv")

    csv_failures, xlsx_failures = validate_outputs()
    report_text = format_report(
        manifest_df=manifest_df,
        summary_df=summary_df,
        missing_df=missing_df,
        public_result=public_result,
        panel_result=panel_result,
        pairs_result=pairs_result,
        warnings=warnings,
        csv_failures=csv_failures,
        xlsx_failures=xlsx_failures,
    )
    write_text(report_text, REPORTS_DIR / "current_data_audit_report.md")

    print_console_summary(summary_df, warnings, missing_df)


if __name__ == "__main__":
    main()
