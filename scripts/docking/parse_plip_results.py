#!/usr/bin/env python3
"""Normalize PLIP XML and calculate pose-based residue/interaction frequencies."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

from pipeline_common import analysis_run_dir, read_tsv, write_tsv

TYPE_MAP = {
    "hydrogen_bond": "hydrogen_bond", "hydrophobic_interaction": "hydrophobic_contact",
    "pi_stack": "pi_stacking", "pi_cation_interaction": "pi_cation", "salt_bridge": "salt_bridge",
    "halogen_bond": "halogen_bond", "metal_complex": "metal_interaction",
}
FILES = {
    "hydrogen_bond": "hydrogen_bonds.tsv", "hydrophobic_contact": "hydrophobic_contacts.tsv",
    "pi_stacking": "pi_stacking.tsv", "pi_cation": "pi_cation.tsv", "salt_bridge": "salt_bridges.tsv",
    "halogen_bond": "halogen_bonds.tsv", "metal_interaction": "metal_interactions.tsv",
}
INTERACTION_FIELDS = ["cas", "receptor", "pocket", "kind", "role", "pose_id", "interaction_type",
                      "residue_name", "residue_number", "residue_chain", "ligand_residue", "distance_A", "source_xml"]


def local(tag): return tag.rsplit("}", 1)[-1]


def children(node):
    result = {}
    for child in node.iter():
        if child is node: continue
        value = (child.text or "").strip()
        if value and local(child.tag) not in result: result[local(child.tag)] = value
    return result


def first(data, names, default=""):
    for name in names:
        if data.get(name, "") != "": return data[name]
    return default


def parse_one(row):
    root = ET.parse(row["plip_xml"]).getroot()
    records = []
    for node in root.iter():
        kind = TYPE_MAP.get(local(node.tag))
        if not kind: continue
        data = children(node)
        record = {key: row.get(key, "") for key in ["cas", "receptor", "pocket", "kind", "role", "pose_id"]}
        record.update({
            "interaction_type": kind,
            "residue_name": first(data, ["restype", "resname", "restype_l"]),
            "residue_number": first(data, ["resnr", "resnum", "resnr_l"]),
            "residue_chain": first(data, ["reschain", "chain", "reschain_l"]),
            "ligand_residue": first(data, ["restype_lig", "resnr_lig", "ligand"]),
            "distance_A": first(data, ["dist", "distance", "dist_h-a", "centdist"]),
            "source_xml": row["plip_xml"],
        })
        records.append(record)
    return records


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--analysis-run", type=Path, default=None); args = parser.parse_args()
    analysis = args.analysis_run or analysis_run_dir()
    root = analysis / "interactions"
    statuses = read_tsv(root / "plip_execution_status.tsv")
    successful = [r for r in statuses if r["status"] in {"SUCCESS", "SKIPPED_COMPLETE"} and r["plip_xml"]]
    records, parse_failures = [], []
    for row in successful:
        try: records.extend(parse_one(row))
        except Exception as exc: parse_failures.append(dict(row, parse_error=repr(exc)))
    for kind, filename in FILES.items():
        write_tsv(root / filename, [r for r in records if r["interaction_type"] == kind], INTERACTION_FIELDS)
    write_tsv(root / "plip_parse_failures.tsv", parse_failures,
              ["cas", "receptor", "pocket", "kind", "role", "pose_id", "plip_xml", "parse_error"])

    requested_by_group, successful_by_group = defaultdict(set), defaultdict(set)
    for row in read_tsv(analysis / "complexes" / "complex_inventory.tsv"):
        if row["kind"] == "all_pose":
            requested_by_group[(row["cas"], row["receptor"], row["pocket"])].add(row["pose_id"])
    for row in statuses:
        if row["kind"] != "all_pose": continue
        key = (row["cas"], row["receptor"], row["pocket"])
        if row["status"] in {"SUCCESS", "SKIPPED_COMPLETE"} and not any(f["plip_xml"] == row["plip_xml"] for f in parse_failures):
            successful_by_group[key].add(row["pose_id"])
    unique_interactions = set()
    for r in records:
        if r["kind"] != "all_pose": continue
        unique_interactions.add((r["cas"], r["receptor"], r["pocket"], r["pose_id"], r["interaction_type"], r["residue_chain"], r["residue_number"], r["residue_name"]))
    counts = Counter((a, b, c, itype, chain, number, name) for a, b, c, pose, itype, chain, number, name in unique_interactions)
    frequency = []
    for (cas, receptor, pocket, itype, chain, number, name), count in sorted(counts.items()):
        key = (cas, receptor, pocket); denominator = len(successful_by_group[key])
        frequency.append({"cas": cas, "receptor": receptor, "pocket": pocket, "interaction_type": itype,
                          "residue_chain": chain, "residue_number": number, "residue_name": name,
                          "pose_occurrence_count": count, "requested_pose_count": len(requested_by_group[key]),
                          "successful_plip_pose_count": denominator, "frequency_denominator": denominator,
                          "frequency": "{:.6f}".format(count / denominator) if denominator else ""})
    write_tsv(root / "residue_interaction_frequency.tsv", frequency,
              ["cas", "receptor", "pocket", "interaction_type", "residue_chain", "residue_number", "residue_name", "pose_occurrence_count", "requested_pose_count", "successful_plip_pose_count", "frequency_denominator", "frequency"])
    fingerprints = []
    for item in sorted(unique_interactions):
        cas, receptor, pocket, pose, itype, chain, number, name = item
        fingerprints.append({"cas": cas, "receptor": receptor, "pocket": pocket, "pose_id": pose,
                             "feature": "{}:{}:{}:{}".format(itype, chain, name, number), "present": 1})
    write_tsv(root / "residue_interaction_fingerprint.tsv", fingerprints,
              ["cas", "receptor", "pocket", "pose_id", "feature", "present"])
    summary = []
    interaction_counts = Counter((r["cas"], r["receptor"], r["pocket"], r["kind"], r["role"], r["pose_id"], r["interaction_type"]) for r in records)
    for key, count in sorted(interaction_counts.items()):
        summary.append(dict(zip(["cas", "receptor", "pocket", "kind", "role", "pose_id", "interaction_type"], key), interaction_count=count))
    write_tsv(root / "interaction_summary.tsv", summary,
              ["cas", "receptor", "pocket", "kind", "role", "pose_id", "interaction_type", "interaction_count"])
    return 0


if __name__ == "__main__": raise SystemExit(main())
