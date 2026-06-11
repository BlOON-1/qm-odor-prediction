"""Smoke tests for structure-disjoint split leakage checks."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def test_structure_disjoint_split_no_cas_overlap_if_present() -> None:
    train_path = PROJECT_ROOT / "data" / "splits" / "structure_disjoint_splits" / "structure_disjoint_train.csv"
    test_path = PROJECT_ROOT / "data" / "splits" / "structure_disjoint_splits" / "structure_disjoint_test.csv"
    if not (train_path.exists() and test_path.exists()):
        return
    train_df = pd.read_csv(train_path)
    test_df = pd.read_csv(test_path)
    assert "cas" in train_df.columns
    assert "cas" in test_df.columns
    train_cas = set(train_df["cas"].fillna("").astype(str).str.strip()) - {""}
    test_cas = set(test_df["cas"].fillna("").astype(str).str.strip()) - {""}
    assert len(train_cas & test_cas) == 0
    assert int(train_df["cas"].fillna("").astype(str).str.strip().replace("", pd.NA).dropna().duplicated().sum()) == 0
    assert int(test_df["cas"].fillna("").astype(str).str.strip().replace("", pd.NA).dropna().duplicated().sum()) == 0

