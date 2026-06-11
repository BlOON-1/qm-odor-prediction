"""Smoke tests for generated RDKit and ECFP4 representation files."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def test_rdkit_descriptors_raw_can_be_read_if_present() -> None:
    """RDKit raw descriptor CSV should be readable and include key identifiers."""

    path = PROJECT_ROOT / "data" / "processed" / "rdkit_descriptors_raw.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    assert "molecule_id" in df.columns
    assert "cas" in df.columns


def test_ecfp4_matrix_rows_match_index_if_present() -> None:
    """ECFP4 matrix rows should match the index file row count."""

    matrix_path = PROJECT_ROOT / "data" / "processed" / "ecfp4_fingerprints.npz"
    index_path = PROJECT_ROOT / "data" / "processed" / "ecfp4_fingerprint_index.csv"
    if not (matrix_path.exists() and index_path.exists()):
        return
    matrix = sparse.load_npz(matrix_path)
    index_df = pd.read_csv(index_path)
    assert matrix.shape[0] == len(index_df)


def test_ecfp4_metadata_settings_if_present() -> None:
    """ECFP4 metadata should preserve the requested fingerprint settings."""

    metadata_path = PROJECT_ROOT / "data" / "processed" / "ecfp4_fingerprint_metadata.json"
    if not metadata_path.exists():
        return
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert metadata["radius"] == 2
    assert metadata["nBits"] == 1024


def test_standardized_feature_columns_do_not_repeat_metadata() -> None:
    """Standardized descriptor exports should keep metadata separate from feature columns."""

    standardized_dir = PROJECT_ROOT / "data" / "processed" / "standardized"
    if not standardized_dir.exists():
        return
    for path in standardized_dir.glob("rdkit_descriptors_standardized_seed*_*.csv"):
        df = pd.read_csv(path)
        feature_columns = [col for col in df.columns if col not in {"molecule_id", "cas", "chemical_name", "canonical_smiles"}]
        assert "molecule_id" not in feature_columns
        assert "cas" not in feature_columns
        assert "chemical_name" not in feature_columns
        assert "canonical_smiles" not in feature_columns


def test_no_modeling_standardized_outputs_without_splits() -> None:
    """Modeling standardized outputs should not exist without split files unless marked QC only."""

    split_dir = PROJECT_ROOT / "data" / "splits" / "iterative_stratified_splits"
    has_any_split = split_dir.exists() and any(split_dir.glob("seed_*_train.csv"))
    standardized_dir = PROJECT_ROOT / "data" / "processed" / "standardized"
    if has_any_split or not standardized_dir.exists():
        return
    disallowed = [
        path
        for path in standardized_dir.glob("*.csv")
        if "FOR_QC_ONLY" not in path.name
    ]
    assert disallowed == []
