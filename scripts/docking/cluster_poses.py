#!/usr/bin/env python3
"""Global fixed-receptor-frame RMSD clustering for each CAS/receptor/pocket."""
from __future__ import annotations

import argparse
import os
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from pipeline_common import (
    analysis_run_dir,
    dock_run_dir,
    map_sdf_heavy_atoms_to_pdbqt,
    parse_pdbqt_atoms,
    read_tsv,
    update_stage_status,
    write_tsv,
)


def fixed_frame_rmsd(a, b, permutations):
    return min(float(np.sqrt(np.mean(np.sum((a - b[list(perm), :]) ** 2, axis=1)))) for perm in permutations)


def stable_labels(raw_labels, energies):
    groups = defaultdict(list)
    for i, label in enumerate(raw_labels): groups[int(label)].append(i)
    order = sorted(groups, key=lambda label: (-len(groups[label]), min(energies[i] for i in groups[label]), label))
    remap = {old: rank + 1 for rank, old in enumerate(order)}
    return np.asarray([remap[int(label)] for label in raw_labels], dtype=int)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-run", type=Path, default=None)
    parser.add_argument("--primary-cutoff", type=float, default=float(os.environ.get("RMSD_PRIMARY_CUTOFF", "2.0")))
    parser.add_argument("--cutoffs", nargs="*", type=float,
                        default=[float(v) for v in os.environ.get("RMSD_SENSITIVITY_CUTOFFS", "1.5 2.0 2.5").split()])
    args = parser.parse_args()
    top_clusters = max(1, int(os.environ.get("TOP_CLUSTERS_FOR_DETAILED_ANALYSIS", "3")))
    analysis = args.analysis_run or analysis_run_dir()
    root = analysis / "clustering"
    matrix_root = root / "pairwise_rmsd"
    rep_root = root / "representatives"
    matrix_root.mkdir(parents=True, exist_ok=True)
    rows = [r for r in read_tsv(analysis / "docking_post_analysis" / "pose_inventory.tsv") if r["extraction_status"] == "SUCCESS"]
    manifests = {(r["cas"], r["receptor"], r["pocket"]): r for r in read_tsv(dock_run_dir() / "job_manifest.tsv")}
    groups = defaultdict(list)
    for row in rows: groups[(row["cas"], row["receptor"], row["pocket"])].append(row)
    assignments, summaries, sensitivity, representatives, failures, mapping_audit = [], [], [], [], [], []

    for key, poses in sorted(groups.items()):
        cas, receptor, pocket = key
        poses.sort(key=lambda r: r["pose_id"])
        try:
            n = len(poses)
            matrix_path = matrix_root / ("{}__{}__{}.npz".format(cas, receptor, pocket))
            expected_pose_ids = np.asarray([p["pose_id"] for p in poses])
            matrix = None
            if os.environ.get("SKIP_EXISTING", "1") == "1" and matrix_path.is_file():
                try:
                    saved = np.load(str(matrix_path), allow_pickle=False)
                    candidate = saved["rmsd"]
                    saved_ids = saved["pose_ids"].astype(str)
                    if (candidate.shape == (n, n) and np.array_equal(saved_ids, expected_pose_ids) and
                            np.all(np.isfinite(candidate)) and np.allclose(candidate, candidate.T) and
                            np.allclose(np.diag(candidate), 0.0)):
                        matrix = candidate
                except Exception:
                    matrix = None
            if matrix is None:
                atom_sets = [parse_pdbqt_atoms(r["pose_pdbqt"]) for r in poses]
                coords = [np.asarray([[a["x"], a["y"], a["z"]] for a in atoms if a["heavy"]], dtype=float) for atoms in atom_sets]
                elements = [[a["element"] for a in atoms if a["heavy"]] for atoms in atom_sets]
                if any(e != elements[0] for e in elements[1:]):
                    raise ValueError("heavy-atom element/order mismatch")
                if key not in manifests:
                    raise ValueError("missing job manifest row for {}".format(key))
                sdf = Path(manifests[key]["ligand_sdf"])
                mapping = map_sdf_heavy_atoms_to_pdbqt(sdf, atom_sets[0])
                perms = mapping["pdbqt_permutations"]
                mapping_audit.append({
                    "cas": cas,
                    "receptor": receptor,
                    "pocket": pocket,
                    "sdf": str(sdf),
                    "heavy_atom_count": len(elements[0]),
                    "mapping_method": mapping["mapping_algorithm"],
                    "topology_mapping_count": mapping["topology_mapping_count"],
                    "atom_type_valid_mapping_count": mapping["atom_type_valid_mapping_count"],
                    "stereo_valid_mapping_count": mapping["stereo_valid_mapping_count"],
                    "stereochemistry_constraint_count": mapping["stereochemistry_constraint_count"],
                    "symmetry_permutation_count": mapping["symmetry_permutation_count"],
                    "source_to_pdbqt": ",".join(
                        str(index) for index in mapping["source_to_pdbqt"]
                    ),
                    "connectivity_threshold": mapping["connectivity_threshold"],
                    "bond_length_rms_A": "{:.6f}".format(mapping["bond_length_rms_A"]),
                    "aromatic_mismatches": mapping["aromatic_mismatches"],
                    "coordinate_fitting": "NONE",
                    "status": "SUCCESS",
                })
                matrix = np.zeros((n, n), dtype=float)
                for i in range(n):
                    for j in range(i + 1, n):
                        matrix[i, j] = matrix[j, i] = fixed_frame_rmsd(coords[i], coords[j], perms)
                np.savez_compressed(str(matrix_path), pose_ids=expected_pose_ids, rmsd=matrix)
            condensed = squareform(matrix, checks=True)
            tree = linkage(condensed, method="average") if n > 1 else None
            energies = np.asarray([float(p["binding_energy_kcal_mol"]) for p in poses])
            labels_by_cutoff = {}
            for cutoff in sorted(set(args.cutoffs + [args.primary_cutoff])):
                raw = fcluster(tree, t=cutoff, criterion="distance") if tree is not None else np.ones(1, dtype=int)
                labels = stable_labels(raw, energies)
                labels_by_cutoff[cutoff] = labels
                counts = [int(np.sum(labels == label)) for label in sorted(set(labels))]
                sensitivity.append({"cas": cas, "receptor": receptor, "pocket": pocket, "cutoff_A": cutoff,
                                    "cluster_count": len(counts), "dominant_cluster_size": max(counts),
                                    "dominant_cluster_fraction": "{:.6f}".format(max(counts) / n)})
            labels = labels_by_cutoff[args.primary_cutoff]
            chosen = {}
            for label in sorted(set(labels)):
                idxs = np.where(labels == label)[0]
                intra = matrix[np.ix_(idxs, idxs)]
                mean_dist = intra.mean(axis=1)
                medoid = int(idxs[int(np.argmin(mean_dist))])
                low = int(idxs[int(np.argmin(energies[idxs]))])
                avg_pair = float(intra[np.triu_indices(len(idxs), 1)].mean()) if len(idxs) > 1 else 0.0
                summaries.append({"cas": cas, "receptor": receptor, "pocket": pocket, "cluster_id": "cluster_{:03d}".format(label),
                                  "pose_count": len(idxs), "cluster_fraction": "{:.6f}".format(len(idxs) / n),
                                  "medoid_pose_id": poses[medoid]["pose_id"], "lowest_energy_pose_id": poses[low]["pose_id"],
                                  "lowest_energy_kcal_mol": "{:.3f}".format(energies[low]),
                                  "mean_energy_kcal_mol": "{:.3f}".format(float(np.mean(energies[idxs]))),
                                  "median_energy_kcal_mol": "{:.3f}".format(float(np.median(energies[idxs]))),
                                  "sd_energy_kcal_mol": "{:.3f}".format(float(np.std(energies[idxs], ddof=1)) if len(idxs) > 1 else 0.0),
                                  "mean_intra_cluster_rmsd_A": "{:.3f}".format(avg_pair)})
                chosen[label] = medoid
            for index, pose in enumerate(poses):
                assignments.append({"pose_id": pose["pose_id"], "cas": cas, "receptor": receptor, "pocket": pocket,
                                    "cluster_id": "cluster_{:03d}".format(labels[index]),
                                    "binding_energy_kcal_mol": pose["binding_energy_kcal_mol"],
                                    "is_cluster_medoid": int(chosen[labels[index]] == index),
                                    "is_global_lowest_energy": int(index == int(np.argmin(energies))),
                                    "pose_pdbqt": pose["pose_pdbqt"]})
            out_dir = rep_root / cas / receptor / pocket
            out_dir.mkdir(parents=True, exist_ok=True)
            selections = [("global_lowest_energy", int(np.argmin(energies))), ("dominant_cluster_medoid", chosen[1])]
            for label in range(2, min(top_clusters, max(chosen)) + 1):
                if label in chosen: selections.append(("cluster_{:03d}_medoid".format(label), chosen[label]))
            for role, index in selections:
                target = out_dir / (role + ".pdbqt")
                source = Path(poses[index]["pose_pdbqt"])
                if not target.is_file() or target.stat().st_size != source.stat().st_size:
                    shutil.copy2(str(source), str(target))
                representatives.append({"cas": cas, "receptor": receptor, "pocket": pocket, "role": role,
                                        "pose_id": poses[index]["pose_id"], "binding_energy_kcal_mol": energies[index],
                                        "representative_pdbqt": str(target), "source_pose_pdbqt": str(source)})
        except Exception as exc:
            failures.append({"cas": cas, "receptor": receptor, "pocket": pocket, "message": repr(exc)})

    write_tsv(root / "pose_cluster_assignments.tsv", assignments,
              ["pose_id", "cas", "receptor", "pocket", "cluster_id", "binding_energy_kcal_mol", "is_cluster_medoid", "is_global_lowest_energy", "pose_pdbqt"])
    write_tsv(root / "cluster_summary.tsv", summaries,
              ["cas", "receptor", "pocket", "cluster_id", "pose_count", "cluster_fraction", "medoid_pose_id", "lowest_energy_pose_id", "lowest_energy_kcal_mol", "mean_energy_kcal_mol", "median_energy_kcal_mol", "sd_energy_kcal_mol", "mean_intra_cluster_rmsd_A"])
    write_tsv(root / "cutoff_sensitivity.tsv", sensitivity,
              ["cas", "receptor", "pocket", "cutoff_A", "cluster_count", "dominant_cluster_size", "dominant_cluster_fraction"])
    write_tsv(root / "representatives.tsv", representatives,
              ["cas", "receptor", "pocket", "role", "pose_id", "binding_energy_kcal_mol", "representative_pdbqt", "source_pose_pdbqt"])
    write_tsv(root / "clustering_failures.tsv", failures, ["cas", "receptor", "pocket", "message"])
    write_tsv(root / "atom_mapping_audit.tsv", mapping_audit,
              ["cas", "receptor", "pocket", "sdf", "heavy_atom_count", "mapping_method",
               "topology_mapping_count", "atom_type_valid_mapping_count", "stereo_valid_mapping_count",
               "stereochemistry_constraint_count", "symmetry_permutation_count", "source_to_pdbqt",
               "connectivity_threshold",
               "bond_length_rms_A", "aromatic_mismatches", "coordinate_fitting", "status"])
    status = "SUCCESS" if groups and not failures else ("WARNING" if representatives else "FAIL")
    update_stage_status(analysis, "clustering", status, "{} representatives; {} failed groups".format(len(representatives), len(failures)), 0 if representatives else 1)
    return 0 if representatives else 1


if __name__ == "__main__":
    raise SystemExit(main())
