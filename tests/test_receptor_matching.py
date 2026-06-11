"""Smoke tests for receptor annotation outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def test_receptor_pairwise_concordance_schema_if_present() -> None:
    path = PROJECT_ROOT / "data" / "receptor_annotations" / "receptor_pairwise_concordance.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    assert "pair_id" in df.columns
    assert "cas_a" in df.columns
    assert "cas_b" in df.columns
    assert "n_jointly_tested_receptors" in df.columns
    valid = df["differential_response_proportion"].dropna()
    assert ((valid >= 0) & (valid <= 1)).all()
    allowed = {
        "concordant",
        "divergent",
        "indeterminate",
        "insufficient_joint_receptor_data",
        "no_matched_or_data",
    }
    assert set(df["receptor_response_group"].dropna().unique()).issubset(allowed)


def test_grouped_pairs_only_concordant_or_divergent_if_present() -> None:
    path = PROJECT_ROOT / "data" / "receptor_annotations" / "receptor_concordant_divergent_pairs.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    assert set(df["receptor_response_group"].dropna().unique()).issubset({"concordant", "divergent"})

