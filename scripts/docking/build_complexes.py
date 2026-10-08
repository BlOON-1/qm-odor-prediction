#!/usr/bin/env python3
"""Build standard receptor-ligand PDB complexes through strict topology/coordinate mapping."""
from __future__ import annotations

import argparse
import os
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

from pipeline_common import (
    analysis_run_dir,
    dock_run_dir,
    map_sdf_heavy_atoms_to_pdbqt,
    parse_pdbqt_atoms,
    read_tsv,
    update_stage_status,
    write_tsv,
)


MAPPING_CACHE = {}


def load_template_with_mapped_atoms(sdf: Path, pose_atoms):
    """Restore SDF topology while accepting a validated PDBQT atom reorder."""
    heavy_pose_atoms = [atom for atom in pose_atoms if bool(atom["heavy"])]
    cache_key = (
        str(sdf.resolve()),
        tuple((str(atom["element"]).upper(), str(atom.get("atom_type", "")).upper())
              for atom in heavy_pose_atoms),
    )
    mapping = MAPPING_CACHE.get(cache_key)
    if mapping is None:
        mapping = map_sdf_heavy_atoms_to_pdbqt(sdf, pose_atoms)
        MAPPING_CACHE[cache_key] = mapping

    # Start from the hydrogen-stripped SDF molecule so its original bond orders
    # are retained.  Coordinates are assigned by the validated source->PDBQT
    # mapping, not by coincidental file order.
    editable = Chem.RWMol(Chem.Mol(mapping["mol"]))
    coordinate_atoms = [
        heavy_pose_atoms[mapping["source_to_pdbqt"][source_index]]
        for source_index in range(editable.GetNumAtoms())
    ]
    pdbqt_to_source = {
        pdbqt_index: source_index
        for source_index, pdbqt_index in enumerate(mapping["source_to_pdbqt"])
    }

    # AutoDock PDBQT normally retains only polar hydrogens.  Attach each one to
    # a uniquely nearest heavy atom and preserve its docked coordinate.  An
    # ambiguous or non-covalent distance is rejected rather than guessed.
    for hydrogen in (atom for atom in pose_atoms if not bool(atom["heavy"])):
        distances = []
        for target_index, heavy in enumerate(heavy_pose_atoms):
            distance = float(np.linalg.norm(np.asarray([
                float(hydrogen["x"]) - float(heavy["x"]),
                float(hydrogen["y"]) - float(heavy["y"]),
                float(hydrogen["z"]) - float(heavy["z"]),
            ])))
            distances.append((distance, target_index))
        distances.sort()
        if not distances or not (0.60 <= distances[0][0] <= 1.35):
            raise ValueError(
                "ATOM_MAPPING_FAILED: PDBQT hydrogen has no heavy-atom neighbour at a covalent distance"
            )
        if len(distances) > 1 and distances[1][0] - distances[0][0] < 0.15:
            raise ValueError(
                "ATOM_MAPPING_FAILED: PDBQT hydrogen attachment is ambiguous"
            )
        parent_source = pdbqt_to_source[distances[0][1]]
        hydrogen_index = editable.AddAtom(Chem.Atom(1))
        editable.AddBond(parent_source, hydrogen_index, Chem.BondType.SINGLE)
        coordinate_atoms.append(hydrogen)

    mol = editable.GetMol()
    if mol.GetNumAtoms() != len(pose_atoms) or len(coordinate_atoms) != len(pose_atoms):
        raise ValueError("ATOM_MAPPING_FAILED: mapped ligand atom count differs from PDBQT")
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for index, atom in enumerate(coordinate_atoms):
        conformer.SetAtomPosition(index, Point3D(float(atom["x"]), float(atom["y"]), float(atom["z"])))
    mol.RemoveAllConformers()
    mol.AddConformer(conformer, assignId=True)
    return mol, mapping


def receptor_records(path: Path):
    records = []
    max_serial = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM  "):
                records.append(line.rstrip("\r\n"))
                try: max_serial = max(max_serial, int(line[6:11]))
                except ValueError: pass
            elif line.startswith("TER"):
                records.append(line.rstrip("\r\n"))
    if not any(line.startswith("ATOM  ") for line in records):
        raise ValueError("complete receptor has no ATOM records")
    return records, max_serial


def ligand_pdb_lines(mol, serial_offset):
    block = Chem.MolToPDBBlock(mol, flavor=4)
    output = []
    for line in block.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            old_serial = int(line[6:11])
            serial = old_serial + serial_offset
            line = "HETATM{:5d}{}{}{}{:4d}{}".format(serial, line[11:17], "LIG", " Z", 1, line[26:])
            output.append(line)
        elif line.startswith("CONECT"):
            numbers = [int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip()]
            if numbers:
                output.append("CONECT" + "".join("{:5d}".format(n + serial_offset) for n in numbers))
    return output


def validate_existing(path: Path, expected_elements) -> bool:
    if not path.is_file() or path.stat().st_size == 0: return False
    elements = []
    with path.open("r", errors="replace") as handle:
        for line in handle:
            if line.startswith("HETATM") and line[17:20].strip() == "LIG" and line[21:22] == "Z":
                element = line[76:78].strip().upper()
                if not element:
                    element = "".join(character for character in line[12:16] if character.isalpha()).upper()[:1]
                elements.append(element)
    return Counter(elements) == Counter(str(element).upper() for element in expected_elements)


def build_one(source_pdbqt: Path, sdf: Path, receptor: Path, pocket_pdb: Path, output: Path, skip_existing: bool):
    pose_atoms = parse_pdbqt_atoms(source_pdbqt)
    expected_elements = [atom["element"] for atom in pose_atoms]
    if skip_existing and validate_existing(output, expected_elements):
        return "SKIPPED_COMPLETE", "validated existing complex"
    mol, mapping = load_template_with_mapped_atoms(sdf, pose_atoms)
    receptor_lines, max_serial = receptor_records(receptor)

    pocket_atoms = parse_pdbqt_atoms(pocket_pdb) if pocket_pdb.suffix.lower() == ".pdbqt" else None
    pocket_coords = []
    if pocket_atoms:
        pocket_coords = [[a["x"], a["y"], a["z"]] for a in pocket_atoms]
    else:
        with pocket_pdb.open("r", errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM  ", "HETATM")):
                    try: pocket_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except ValueError: pass
    if not pocket_coords: raise ValueError("pocket PDB has no coordinates")
    center = np.mean(np.asarray(pocket_coords), axis=0)
    ligand_center = np.mean(np.asarray([[a["x"], a["y"], a["z"]] for a in pose_atoms if a["heavy"]]), axis=0)
    if float(np.linalg.norm(center - ligand_center)) > 40.0:
        raise ValueError("POSE_OUTSIDE_POCKET_REGION: ligand centroid is over 40 A from pocket center")

    ligand_lines = ligand_pdb_lines(mol, max_serial)
    if sum(line.startswith("HETATM") for line in ligand_lines) != len(pose_atoms):
        raise ValueError("ATOM_MAPPING_FAILED: output ligand atom count changed")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    content = "\n".join(receptor_lines + ligand_lines + ["END"]) + "\n"
    # Path.write_text(newline=...) is unavailable on Python 3.9.  Keep the
    # explicit Unix newline policy by opening the temporary file directly.
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(content)
    temporary.replace(output)
    if not validate_existing(output, expected_elements):
        raise ValueError("ATOM_MAPPING_FAILED: written complex failed ligand element/count validation")
    return "SUCCESS", (
        "topology-mapped SDF-to-PDBQT atoms; coordinates validated; "
        "{} symmetry permutations; no coordinate fitting".format(
            mapping["symmetry_permutation_count"]
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-run", type=Path, default=None)
    parser.add_argument("--dock-run", type=Path, default=None)
    args = parser.parse_args()
    analysis = args.analysis_run or analysis_run_dir()
    dock = args.dock_run or dock_run_dir()
    skip = os.environ.get("SKIP_EXISTING", "1") == "1"
    manifests = {(r["cas"], r["receptor"], r["pocket"]): r for r in read_tsv(dock / "job_manifest.tsv")}
    reps = read_tsv(analysis / "clustering" / "representatives.tsv")
    poses = [r for r in read_tsv(analysis / "docking_post_analysis" / "pose_inventory.tsv") if r["extraction_status"] == "SUCCESS"]
    requested = []
    for row in reps:
        requested.append((row["cas"], row["receptor"], row["pocket"], "representative", row["role"], row["pose_id"], Path(row["representative_pdbqt"])))
    if os.environ.get("RUN_PLIP_ALL_POSES", "1") == "1":
        for row in poses:
            requested.append((row["cas"], row["receptor"], row["pocket"], "all_pose", row["pose_id"], row["pose_id"], Path(row["pose_pdbqt"])))
    inventory = []
    for cas, receptor, pocket, kind, role, pose_id, source in requested:
        meta = manifests.get((cas, receptor, pocket))
        if not meta:
            status, message, output = "FAILED", "missing docking manifest row", Path("")
        else:
            subdir = "representative" if kind == "representative" else "all_poses"
            filename = (role + "_complex.pdb") if kind == "representative" else (pose_id + "_complex.pdb")
            output = analysis / "complexes" / cas / receptor / pocket / subdir / filename
            try:
                status, message = build_one(source, Path(meta["ligand_sdf"]), Path(meta["full_receptor_pdb"]), Path(meta["pocket_pdb"]), output, skip)
            except Exception as exc:
                status, message = "FAILED", str(exc)
        inventory.append({"cas": cas, "receptor": receptor, "pocket": pocket, "kind": kind, "role": role,
                          "pose_id": pose_id, "source_pdbqt": str(source), "complex_pdb": str(output),
                          "status": status, "message": message})
    write_tsv(analysis / "complexes" / "complex_inventory.tsv", inventory,
              ["cas", "receptor", "pocket", "kind", "role", "pose_id", "source_pdbqt", "complex_pdb", "status", "message"])
    failures = sum(row["status"] == "FAILED" for row in inventory)
    rep_failures = sum(row["status"] == "FAILED" and row["kind"] == "representative" for row in inventory)
    status = "SUCCESS" if inventory and failures == 0 else ("WARNING" if inventory and rep_failures == 0 else "FAIL")
    update_stage_status(analysis, "complexes", status, "{} complexes; {} failures".format(len(inventory) - failures, failures), 0 if rep_failures == 0 and inventory else 1)
    return 0 if inventory and rep_failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
