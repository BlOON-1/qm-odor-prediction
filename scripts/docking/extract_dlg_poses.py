#!/usr/bin/env python3
"""Extract each AutoDock4 DLG model with its own energy and shard provenance."""
from __future__ import annotations

import argparse
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from pipeline_common import analysis_run_dir, dock_run_dir, parse_pdbqt_atoms, read_tsv, update_stage_status, write_tsv

ENERGY_RE = re.compile(r"Estimated Free Energy of Binding\s*=\s*([-+]?\d+(?:\.\d+)?)\s*kcal/mol", re.I)
MODEL_RE = re.compile(r"^MODEL\s+(\d+)", re.I)
DOCKED_RE = re.compile(r"^DOCKED:\s?(.*)$")
FIELDS = ["pose_id", "cas", "receptor", "pocket", "shard_index", "model_index",
          "binding_energy_kcal_mol", "seed1", "seed2", "source_dlg", "pose_pdbqt",
          "atom_count", "heavy_atom_count", "extraction_status", "message"]


def parse_dlg(path: Path) -> Tuple[List[Dict[str, object]], List[str]]:
    """Parse only energy lines that occur inside their corresponding DOCKED model."""
    poses: List[Dict[str, object]] = []
    errors: List[str] = []
    current: List[str] = []
    model_index = None
    saw_docked = False

    def finish() -> None:
        nonlocal current, model_index
        if not current:
            return
        energies = []
        for line in current:
            match = ENERGY_RE.search(line)
            if match:
                energies.append(float(match.group(1)))
        if model_index is None:
            errors.append("model block without MODEL number")
        elif len(energies) != 1:
            errors.append("model {} has {} binding-energy records".format(model_index, len(energies)))
        else:
            coordinate_lines = [line for line in current if line.startswith(("ATOM  ", "HETATM"))]
            if not coordinate_lines:
                errors.append("model {} has no coordinates".format(model_index))
            else:
                poses.append({"model_index": model_index, "energy": energies[0], "lines": list(current)})
        current = []
        model_index = None

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            match = DOCKED_RE.match(raw.rstrip("\r\n"))
            if not match:
                continue
            saw_docked = True
            line = match.group(1)
            model_match = MODEL_RE.match(line)
            if model_match:
                finish()
                model_index = int(model_match.group(1))
                current = [line]
            elif current:
                current.append(line)
                if line.strip().upper() == "ENDMDL":
                    finish()
    finish()
    if not saw_docked:
        errors.append("no DOCKED: records found; unsupported or incomplete DLG")
    return poses, errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dock-run", type=Path, default=None)
    parser.add_argument("--analysis-run", type=Path, default=None)
    args = parser.parse_args()
    dock = args.dock_run or dock_run_dir()
    analysis = args.analysis_run or analysis_run_dir()
    out_root = analysis / "docking_post_analysis"
    inventory_path = out_root / "pose_inventory.tsv"
    failures_path = out_root / "pose_extraction_failures.tsv"
    out_root.mkdir(parents=True, exist_ok=True)

    queue = read_tsv(dock / "docking_shard_queue.tsv")
    results = {row["task_id"]: row for row in read_tsv(dock / "docking_shard_results.tsv")}
    manifests = {(r["cas"], r["receptor"], r["pocket"]): r for r in read_tsv(dock / "job_manifest.tsv")}
    if not queue or not manifests:
        raise RuntimeError("missing/empty docking shard queue or job manifest under {}".format(dock))

    grouped = defaultdict(list)
    failure_rows = []
    for shard in queue:
        key = (shard["cas"], shard["receptor"], shard["pocket"])
        result = results.get(shard["task_id"], {})
        if result.get("status") not in {"SUCCESS", "SKIPPED_COMPLETE"}:
            failure_rows.append({"cas": key[0], "receptor": key[1], "pocket": key[2],
                                 "shard_index": shard["shard_index"], "source_dlg": shard["dlg"],
                                 "message": "shard status is {}".format(result.get("status", "MISSING"))})
            continue
        dlg = Path(shard["dlg"])
        if not dlg.is_file() or dlg.stat().st_size == 0:
            failure_rows.append({"cas": key[0], "receptor": key[1], "pocket": key[2],
                                 "shard_index": shard["shard_index"], "source_dlg": str(dlg),
                                 "message": "successful shard has missing/empty DLG"})
            continue
        poses, errors = parse_dlg(dlg)
        for error in errors:
            failure_rows.append({"cas": key[0], "receptor": key[1], "pocket": key[2],
                                 "shard_index": shard["shard_index"], "source_dlg": str(dlg), "message": error})
        for pose in poses:
            pose.update({"shard": shard, "dlg": str(dlg)})
            grouped[key].append(pose)

    inventory = []
    for key, manifest in sorted(manifests.items()):
        cas, receptor, pocket = key
        poses = sorted(grouped.get(key, []), key=lambda p: (int(p["shard"]["shard_index"]), int(p["model_index"])))
        signature = None
        expected = sum(int(s["runs_requested"]) for s in queue if (s["cas"], s["receptor"], s["pocket"]) == key)
        pose_dir = out_root / "extracted_poses" / cas / receptor / pocket
        pose_dir.mkdir(parents=True, exist_ok=True)
        for ordinal, pose in enumerate(poses, 1):
            pose_id = "pose_{:04d}".format(ordinal)
            pose_path = pose_dir / (pose_id + ".pdbqt")
            content = "\n".join(pose["lines"]) + "\n"
            if not (pose_path.is_file() and pose_path.stat().st_size and pose_path.read_text(errors="replace") == content):
                # pathlib.Path.write_text() did not gain the ``newline``
                # argument until Python 3.10.  The production environment is
                # Python 3.9, while Path.open() supports newline explicitly.
                with pose_path.open("w", encoding="utf-8", newline="\n") as handle:
                    handle.write(content)
            status, message = "SUCCESS", ""
            try:
                atoms = parse_pdbqt_atoms(pose_path)
                this_signature = [(a["name"], a["element"]) for a in atoms]
                if signature is None:
                    signature = this_signature
                elif this_signature != signature:
                    raise ValueError("atom name/element/order differs from the first pose")
            except Exception as exc:
                atoms, status, message = [], "FAILED", str(exc)
                failure_rows.append({"cas": cas, "receptor": receptor, "pocket": pocket,
                                     "shard_index": pose["shard"]["shard_index"],
                                     "source_dlg": pose["dlg"], "message": message})
            inventory.append({
                "pose_id": pose_id, "cas": cas, "receptor": receptor, "pocket": pocket,
                "shard_index": pose["shard"]["shard_index"], "model_index": pose["model_index"],
                "binding_energy_kcal_mol": "{:.3f}".format(pose["energy"]),
                "seed1": pose["shard"]["seed1"], "seed2": pose["shard"]["seed2"],
                "source_dlg": pose["dlg"], "pose_pdbqt": str(pose_path),
                "atom_count": len(atoms), "heavy_atom_count": sum(bool(a["heavy"]) for a in atoms),
                "extraction_status": status, "message": message,
            })
        successful = sum(1 for r in inventory if (r["cas"], r["receptor"], r["pocket"]) == key and r["extraction_status"] == "SUCCESS")
        if successful != expected:
            failure_rows.append({"cas": cas, "receptor": receptor, "pocket": pocket, "shard_index": "ALL",
                                 "source_dlg": "", "message": "expected {} poses; validated {}".format(expected, successful)})

    write_tsv(inventory_path, inventory, FIELDS)
    write_tsv(failures_path, failure_rows, ["cas", "receptor", "pocket", "shard_index", "source_dlg", "message"])
    failed_poses = sum(r["extraction_status"] != "SUCCESS" for r in inventory)
    status = "SUCCESS" if not failure_rows and not failed_poses else "WARNING"
    update_stage_status(analysis, "extraction", status,
                        "{} validated poses; {} extraction/QC findings".format(len(inventory) - failed_poses, len(failure_rows)))
    return 0 if inventory else 1


if __name__ == "__main__":
    raise SystemExit(main())
