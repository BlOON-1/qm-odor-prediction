"""Optional receptor-group descriptor distance comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from qmodor.descriptors.common import write_csv, write_text


def _find_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def compute_receptor_descriptor_distances(project_root: Path) -> dict[str, Any]:
    """Compute descriptor distances for receptor-classified pairs when standardized descriptors exist."""

    results_dir = project_root / "results" / "receptor_consistency"
    reports_dir = project_root / "reports"
    results_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    pair_path = project_root / "data" / "receptor_annotations" / "receptor_concordant_divergent_pairs.csv"
    rdkit_candidates = sorted((project_root / "data" / "processed" / "standardized").glob("rdkit_descriptors_standardized_seed*_train.csv"))
    qm_candidates = sorted((project_root / "data" / "processed").glob("qm_descriptors*.csv")) + sorted((project_root / "data" / "processed" / "standardized").glob("qm_descriptors_standardized*.csv"))

    rdkit_path = rdkit_candidates[0] if rdkit_candidates else None
    qm_path = qm_candidates[0] if qm_candidates else None

    missing = []
    if not pair_path.exists():
        missing.append(str(pair_path.relative_to(project_root)).replace("\\", "/"))
    if rdkit_path is None:
        missing.append("data/processed/standardized/rdkit_descriptors_standardized_seed*_train.csv")
    if qm_path is None:
        missing.append("data/processed/qm_descriptors.csv or data/processed/standardized/qm_descriptors_standardized*.csv")

    if missing:
        blocked_report = "\n".join(
            [
                "# Receptor Distance Comparison Blocked",
                "",
                "Descriptor-distance comparison was not generated because one or more required descriptor files were missing.",
                "",
                "Missing requirements:",
                *[f"- {item}" for item in missing],
            ]
        )
        blocked_path = reports_dir / "receptor_distance_comparison_blocked_missing_descriptors.md"
        write_text(blocked_report + "\n", blocked_path)
        return {
            "generated": False,
            "missing": missing,
            "blocked_report": blocked_path,
            "output_files": [blocked_path],
        }

    pair_df = pd.read_csv(pair_path)
    rdkit_df = pd.read_csv(rdkit_path)
    qm_df = pd.read_csv(qm_path)
    if "cas" not in rdkit_df.columns or "cas" not in qm_df.columns:
        blocked_path = reports_dir / "receptor_distance_comparison_blocked_missing_descriptors.md"
        write_text(
            "# Receptor Distance Comparison Blocked\n\nRequired descriptor files exist, but one or more files do not contain a `cas` column.\n",
            blocked_path,
        )
        return {
            "generated": False,
            "missing": ["cas column in descriptor files"],
            "blocked_report": blocked_path,
            "output_files": [blocked_path],
        }

    meta_cols = {"molecule_id", "cas", "chemical_name", "smiles", "canonical_smiles"}
    rdkit_features = [col for col in rdkit_df.columns if col not in meta_cols]
    qm_features = [col for col in qm_df.columns if col not in meta_cols]
    rdkit_lookup = rdkit_df.drop_duplicates(subset=["cas"], keep="first").set_index("cas")
    qm_lookup = qm_df.drop_duplicates(subset=["cas"], keep="first").set_index("cas")

    rows = []
    for _, row in pair_df.iterrows():
        cas_a = row["cas_a"]
        cas_b = row["cas_b"]
        if cas_a not in rdkit_lookup.index or cas_b not in rdkit_lookup.index or cas_a not in qm_lookup.index or cas_b not in qm_lookup.index:
            continue
        rdkit_dist = float(np.linalg.norm(rdkit_lookup.loc[cas_a, rdkit_features].to_numpy(dtype=float) - rdkit_lookup.loc[cas_b, rdkit_features].to_numpy(dtype=float)))
        qm_dist = float(np.linalg.norm(qm_lookup.loc[cas_a, qm_features].to_numpy(dtype=float) - qm_lookup.loc[cas_b, qm_features].to_numpy(dtype=float)))
        rows.append(
            {
                "pair_id": row["pair_id"],
                "cas_a": cas_a,
                "cas_b": cas_b,
                "receptor_response_group": row["receptor_response_group"],
                "n_jointly_tested_receptors": row["n_jointly_tested_receptors"],
                "qm_distance": qm_dist,
                "rdkit_distance": rdkit_dist,
                "descriptor_source_note": f"RDKit from {rdkit_path.name}; QM from {qm_path.name}",
            }
        )

    distance_df = pd.DataFrame(rows)
    summary_df = (
        distance_df.groupby("receptor_response_group", dropna=False)
        .agg(
            n_pairs=("pair_id", "count"),
            mean_qm_distance=("qm_distance", "mean"),
            mean_rdkit_distance=("rdkit_distance", "mean"),
            median_qm_distance=("qm_distance", "median"),
            median_rdkit_distance=("rdkit_distance", "median"),
        )
        .reset_index()
        if not distance_df.empty
        else pd.DataFrame(columns=["receptor_response_group", "n_pairs", "mean_qm_distance", "mean_rdkit_distance", "median_qm_distance", "median_rdkit_distance"])
    )

    distance_output = results_dir / "qm_rdkit_distance_by_receptor_group.csv"
    summary_output = results_dir / "qm_rdkit_distance_group_summary.csv"
    write_csv(distance_df, distance_output)
    write_csv(summary_df, summary_output)
    return {
        "generated": True,
        "distance_df": distance_df,
        "summary_df": summary_df,
        "output_files": [distance_output, summary_output],
    }

