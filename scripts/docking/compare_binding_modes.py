#!/usr/bin/env python3
"""Compare eligible cross-CAS ligand binding modes in a fixed receptor frame."""
from __future__ import annotations

import argparse
import math
import os
from collections import defaultdict
from pathlib import Path

import numpy as np

from pipeline_common import analysis_run_dir, dock_run_dir, parse_pdbqt_atoms, read_tsv, update_stage_status, write_tsv
from build_receptor_ligand_index import SKIP_FIELDS


COMPARISON_FIELDS = [
    "comparison_id", "receptor", "canonical_receptor_key", "cas_A", "cas_B", "pocket_index",
    "pocket_A", "pocket_B", "comparison_status", "centroid_distance_A",
    "principal_axis_angle_deg", "residue_jaccard", "interaction_fingerprint_similarity",
    "hbond_residue_jaccard", "hydrophobic_residue_jaccard", "pi_residue_jaccard",
    "salt_bridge_residue_jaccard", "median_energy_A", "median_energy_B",
    "median_energy_difference", "mean_energy_A", "mean_energy_B", "mean_energy_difference",
    "dominant_cluster_fraction_A", "dominant_cluster_fraction_B",
    "dominant_cluster_fraction_difference", "shared_residue_count", "union_residue_count",
    "shared_interaction_count", "union_interaction_count", "same_spatial_region",
    "representative_pose_A", "representative_pose_B", "message",
]


def jaccard(a, b):
    union = set(a) | set(b)
    return 1.0 if not union else len(set(a) & set(b)) / len(union)


def geometry(path):
    atoms = [atom for atom in parse_pdbqt_atoms(path) if atom["heavy"]]
    if len(atoms) < 2:
        raise ValueError("representative pose has fewer than two heavy atoms")
    coords = np.asarray([[atom["x"], atom["y"], atom["z"]] for atom in atoms], dtype=float)
    center = coords.mean(axis=0)
    centered = coords - center
    if float(np.linalg.norm(centered)) == 0.0:
        raise ValueError("representative pose has no defined principal axis")
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    return center, vh[0]


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _fmt(value, digits=6):
    return "NA" if value is None else ("{:." + str(digits) + "f}").format(value)


def _feature_parts(feature):
    parts = str(feature).split(":", 3)
    if len(parts) != 4:
        return None
    interaction_type, chain, residue_name, residue_number = parts
    return interaction_type, "{}:{}:{}".format(chain, residue_name, residue_number)


def _merge_skips(path: Path, new_rows) -> None:
    rows = read_tsv(path) + list(new_rows)
    unique = {}
    for row in rows:
        key = tuple(row.get(field, "") for field in SKIP_FIELDS)
        unique[key] = row
    write_tsv(path, [unique[key] for key in sorted(unique)], SKIP_FIELDS)


def _spatial_classification(centroid, residue_similarity, fp_similarity):
    centroid_limit = float(os.environ.get("COMPARISON_CENTROID_THRESHOLD_A", "5.0"))
    residue_limit = float(os.environ.get("COMPARISON_RESIDUE_JACCARD_THRESHOLD", "0.3"))
    fp_limit = float(os.environ.get("COMPARISON_FINGERPRINT_THRESHOLD", "0.3"))
    if residue_similarity is None or fp_similarity is None:
        return "NOT_ASSESSED_PARTIAL_INTERACTION_DATA"
    if centroid <= centroid_limit and residue_similarity >= residue_limit and fp_similarity >= fp_limit:
        return "WORKFLOW_SIMILAR_REGION"
    return "WORKFLOW_DISTINCT_OR_INDETERMINATE_REGION"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-run", type=Path, default=None)
    parser.add_argument("--dock-run", type=Path, default=None)
    args = parser.parse_args()
    analysis = args.analysis_run or analysis_run_dir()
    dock = args.dock_run or dock_run_dir()
    comparison_root = analysis / "comparisons"

    index_rows = read_tsv(comparison_root / "receptor_ligand_index.tsv")
    receptor_groups = read_tsv(comparison_root / "receptor_ligand_groups.tsv")
    manifests = read_tsv(dock / "job_manifest.tsv")
    if not manifests or not index_rows or not receptor_groups:
        raise RuntimeError("receptor index/group tables or docking job_manifest.tsv are missing or empty")
    eligibility = [row for row in read_tsv(comparison_root / "comparison_eligibility.tsv")
                   if row.get("eligible") == "1"]
    reps = [row for row in read_tsv(analysis / "clustering" / "representatives.tsv")
            if row.get("role") == "dominant_cluster_medoid"]
    clusters = read_tsv(analysis / "clustering" / "cluster_summary.tsv")
    assignments = read_tsv(analysis / "clustering" / "pose_cluster_assignments.tsv")
    energy_rows = read_tsv(comparison_root / "energy_summary.tsv")
    fingerprints = read_tsv(analysis / "interactions" / "residue_interaction_fingerprint.tsv")
    # Loaded intentionally as part of the comparison evidence bundle.  The
    # medoid-specific metrics come from fingerprints, while this table remains
    # available for traceability and final reporting.
    read_tsv(analysis / "interactions" / "residue_interaction_frequency.tsv")
    plip_status = read_tsv(analysis / "interactions" / "plip_execution_status.tsv")

    rep_by_key = {(row.get("cas"), row.get("receptor"), row.get("pocket")): row for row in reps}
    cluster_by_key = {(row.get("cas"), row.get("receptor"), row.get("pocket")): row
                      for row in clusters if row.get("cluster_id") == "cluster_001"}
    assignment_keys = {(row.get("cas"), row.get("receptor"), row.get("pocket"), row.get("pose_id"))
                       for row in assignments}
    energy_by_key = {(row.get("cas"), row.get("receptor"), row.get("pocket")): row for row in energy_rows}

    features = defaultdict(set)
    for row in fingerprints:
        if row.get("present", "1") in {"0", "0.0", "false", "False"}:
            continue
        features[(row.get("cas"), row.get("receptor"), row.get("pocket"), row.get("pose_id"))].add(row.get("feature", ""))
    successful_plip = {
        (row.get("cas"), row.get("receptor"), row.get("pocket"), row.get("pose_id"))
        for row in plip_status
        if row.get("status") in {"SUCCESS", "SKIPPED_COMPLETE"}
    }

    comparisons, runtime_skips, residue_matrix, fp_matrix = [], [], [], []
    exceptions = 0
    for eligible in eligibility:
        canonical = eligible["canonical_receptor_key"]
        display = eligible.get("receptor_display_name", "")
        cas_a, cas_b = eligible["cas_A"], eligible["cas_B"]
        pocket_index = eligible["pocket_index"]
        receptor_a, receptor_b = eligible["receptor_name_A"], eligible["receptor_name_B"]
        pocket_a, pocket_b = eligible["pocket_A"], eligible["pocket_B"]
        key_a = (cas_a, receptor_a, pocket_a)
        key_b = (cas_b, receptor_b, pocket_b)
        rep_a, rep_b = rep_by_key.get(key_a), rep_by_key.get(key_b)
        cluster_a, cluster_b = cluster_by_key.get(key_a), cluster_by_key.get(key_b)
        missing = []
        if not rep_a or not cluster_a:
            missing.append("{} dominant cluster/medoid".format(cas_a))
        if not rep_b or not cluster_b:
            missing.append("{} dominant cluster/medoid".format(cas_b))
        if rep_a and assignments and (cas_a, receptor_a, pocket_a, rep_a.get("pose_id")) not in assignment_keys:
            missing.append("{} medoid assignment".format(cas_a))
        if rep_b and assignments and (cas_b, receptor_b, pocket_b, rep_b.get("pose_id")) not in assignment_keys:
            missing.append("{} medoid assignment".format(cas_b))
        if missing:
            runtime_skips.append({
                "canonical_receptor_key": canonical, "receptor_display_name": display,
                "cas_A": cas_a, "cas_B": cas_b, "pocket_index": pocket_index,
                "status": "SKIPPED_CLUSTER_RESULT_MISSING", "reason": "Required clustering result is missing",
                "details": "; ".join(missing),
            })
            continue
        try:
            path_a, path_b = Path(rep_a["representative_pdbqt"]), Path(rep_b["representative_pdbqt"])
            center_a, axis_a = geometry(path_a)
            center_b, axis_b = geometry(path_b)
            centroid = float(np.linalg.norm(center_a - center_b))
            cosine = min(1.0, max(0.0, abs(float(np.dot(axis_a, axis_b)))))
            angle = math.degrees(math.acos(cosine))

            pose_key_a = key_a + (rep_a.get("pose_id"),)
            pose_key_b = key_b + (rep_b.get("pose_id"),)
            plip_a = pose_key_a in successful_plip or pose_key_a in features
            plip_b = pose_key_b in successful_plip or pose_key_b in features
            interaction_a = features[pose_key_a] if plip_a else set()
            interaction_b = features[pose_key_b] if plip_b else set()
            type_residues_a, type_residues_b = defaultdict(set), defaultdict(set)
            for feature in interaction_a:
                parsed = _feature_parts(feature)
                if parsed:
                    type_residues_a[parsed[0]].add(parsed[1])
            for feature in interaction_b:
                parsed = _feature_parts(feature)
                if parsed:
                    type_residues_b[parsed[0]].add(parsed[1])
            residues_a = set().union(*type_residues_a.values()) if type_residues_a else set()
            residues_b = set().union(*type_residues_b.values()) if type_residues_b else set()
            interaction_complete = plip_a and plip_b
            residue_similarity = jaccard(residues_a, residues_b) if interaction_complete else None
            fp_similarity = jaccard(interaction_a, interaction_b) if interaction_complete else None

            energy_a, energy_b = energy_by_key.get(key_a, {}), energy_by_key.get(key_b, {})
            median_a = _float(energy_a.get("median_energy_kcal_mol"))
            median_b = _float(energy_b.get("median_energy_kcal_mol"))
            mean_a = _float(energy_a.get("mean_energy_kcal_mol"))
            mean_b = _float(energy_b.get("mean_energy_kcal_mol"))
            fraction_a = _float(cluster_a.get("cluster_fraction"))
            fraction_b = _float(cluster_b.get("cluster_fraction"))
            energy_complete = None not in (median_a, median_b, mean_a, mean_b)
            cluster_fraction_complete = None not in (fraction_a, fraction_b)
            if not interaction_complete:
                comparison_status = "PARTIAL_INTERACTION_DATA"
                message = "Geometry and available energy/cluster metrics completed; medoid PLIP data missing for {}".format(
                    ",".join(cas for cas, ok in ((cas_a, plip_a), (cas_b, plip_b)) if not ok))
            elif not energy_complete or not cluster_fraction_complete:
                comparison_status = "PARTIAL_ENERGY_OR_CLUSTER_DATA"
                message = "Geometry and interaction metrics completed; energy or cluster fraction is incomplete"
            else:
                comparison_status = "SUCCESS"
                message = "Eligible cross-ligand comparison completed"

            def typed_similarity(types):
                if not interaction_complete:
                    return None
                left = set().union(*(type_residues_a[t] for t in types))
                right = set().union(*(type_residues_b[t] for t in types))
                return jaccard(left, right)

            comparisons.append({
                "comparison_id": eligible["comparison_id"], "receptor": display,
                "canonical_receptor_key": canonical, "cas_A": cas_a, "cas_B": cas_b,
                "pocket_index": pocket_index, "pocket_A": pocket_a, "pocket_B": pocket_b,
                "comparison_status": comparison_status, "centroid_distance_A": _fmt(centroid, 3),
                "principal_axis_angle_deg": _fmt(angle, 3), "residue_jaccard": _fmt(residue_similarity),
                "interaction_fingerprint_similarity": _fmt(fp_similarity),
                "hbond_residue_jaccard": _fmt(typed_similarity(("hydrogen_bond",))),
                "hydrophobic_residue_jaccard": _fmt(typed_similarity(("hydrophobic_contact",))),
                "pi_residue_jaccard": _fmt(typed_similarity(("pi_stacking", "pi_cation"))),
                "salt_bridge_residue_jaccard": _fmt(typed_similarity(("salt_bridge",))),
                "median_energy_A": _fmt(median_a, 3), "median_energy_B": _fmt(median_b, 3),
                "median_energy_difference": _fmt(abs(median_a - median_b) if energy_complete else None, 3),
                "mean_energy_A": _fmt(mean_a, 3), "mean_energy_B": _fmt(mean_b, 3),
                "mean_energy_difference": _fmt(abs(mean_a - mean_b) if energy_complete else None, 3),
                "dominant_cluster_fraction_A": _fmt(fraction_a),
                "dominant_cluster_fraction_B": _fmt(fraction_b),
                "dominant_cluster_fraction_difference": _fmt(abs(fraction_a - fraction_b) if cluster_fraction_complete else None),
                "shared_residue_count": len(residues_a & residues_b) if interaction_complete else "NA",
                "union_residue_count": len(residues_a | residues_b) if interaction_complete else "NA",
                "shared_interaction_count": len(interaction_a & interaction_b) if interaction_complete else "NA",
                "union_interaction_count": len(interaction_a | interaction_b) if interaction_complete else "NA",
                "same_spatial_region": _spatial_classification(centroid, residue_similarity, fp_similarity),
                "representative_pose_A": str(path_a), "representative_pose_B": str(path_b),
                "message": message,
            })
            residue_matrix.append({"comparison_id": eligible["comparison_id"], "canonical_receptor_key": canonical,
                                   "pocket_index": pocket_index, "cas_A": cas_a, "cas_B": cas_b,
                                   "residue_jaccard": _fmt(residue_similarity)})
            fp_matrix.append({"comparison_id": eligible["comparison_id"], "canonical_receptor_key": canonical,
                              "pocket_index": pocket_index, "cas_A": cas_a, "cas_B": cas_b,
                              "fingerprint_similarity": _fmt(fp_similarity)})
        except Exception as exc:
            exceptions += 1
            runtime_skips.append({
                "canonical_receptor_key": canonical, "receptor_display_name": display,
                "cas_A": cas_a, "cas_B": cas_b, "pocket_index": pocket_index,
                "status": "FAILED_COMPARISON_EXCEPTION", "reason": "Comparison program exception",
                "details": repr(exc),
            })

    write_tsv(comparison_root / "binding_mode_comparison.tsv", comparisons, COMPARISON_FIELDS)
    write_tsv(comparison_root / "residue_similarity_matrix.tsv", residue_matrix,
              ["comparison_id", "canonical_receptor_key", "pocket_index", "cas_A", "cas_B", "residue_jaccard"])
    write_tsv(comparison_root / "interaction_fingerprint_similarity_matrix.tsv", fp_matrix,
              ["comparison_id", "canonical_receptor_key", "pocket_index", "cas_A", "cas_B", "fingerprint_similarity"])
    _merge_skips(comparison_root / "comparison_skips.tsv", runtime_skips)

    project = Path(__file__).resolve().parent.parent
    qm_status = "AVAILABLE_NOT_INFERRED" if (project / "qm_descriptors.tsv").is_file() and (project / "comparison_design.tsv").is_file() else "NOT_PROVIDED"
    (comparison_root / "qm_analysis_status.txt").write_text("QM_ANALYSIS_STATUS={}\n".format(qm_status), encoding="utf-8")
    if exceptions:
        status = "WARNING"
    elif not eligibility:
        status = "SUCCESS_NO_ELIGIBLE_PAIRS"
    elif any(row["comparison_status"].startswith("PARTIAL") for row in comparisons) or runtime_skips:
        status = "WARNING"
    else:
        status = "SUCCESS"
    update_stage_status(analysis, "comparison", status,
                        "{} manifest rows; {} receptor groups; {} eligible pocket comparisons; {} results; {} runtime skips; {} exceptions; QM_ANALYSIS_STATUS={}".format(
                            len(manifests), len(receptor_groups), len(eligibility), len(comparisons), len(runtime_skips), exceptions, qm_status))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
