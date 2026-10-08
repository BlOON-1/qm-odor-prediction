#!/usr/bin/env python3
"""Shared, Python-3.9-compatible helpers for the docking post-analysis modules."""
from __future__ import annotations

import csv
import hashlib
import math
import os
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

CAS_RE = re.compile(r"^[0-9]+-[0-9]+-[0-9]+$")
PDB_RECORDS = ("ATOM  ", "HETATM")


def read_tsv(path: os.PathLike) -> List[Dict[str, str]]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: os.PathLike, rows: Iterable[Dict[str, object]], fields: Sequence[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_field(row.get(key, "")) for key in fields})


def clean_field(value: object) -> str:
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def md5_file(path: os.PathLike) -> str:
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_pdbqt_atoms(path: os.PathLike) -> List[Dict[str, object]]:
    atoms = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(PDB_RECORDS):
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError as exc:
                raise ValueError("invalid coordinate in {}: {!r}".format(path, line.rstrip())) from exc
            fields = line.split()
            atom_type = fields[-1] if fields else ""
            # AutoDock types are not chemical symbols: A is aromatic carbon,
            # NA is acceptor nitrogen, OA is acceptor oxygen, and HD is H.
            ad_element = {
                "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
                "OA": "O", "OS": "O", "O": "O", "SA": "S", "S": "S",
                "HD": "H", "HS": "H", "H": "H", "F": "F", "CL": "CL",
                "BR": "BR", "I": "I", "P": "P", "SI": "SI", "MG": "MG",
                "ZN": "ZN", "FE": "FE", "MN": "MN", "CU": "CU", "CA": "CA",
            }
            normalized_type = re.sub(r"[^A-Za-z]", "", atom_type).upper()
            element = ad_element.get(normalized_type)
            if not element:
                raw_element = re.sub(r"[^A-Za-z]", "", line[76:78]).upper()
                element = raw_element or re.sub(r"[^A-Za-z]", "", line[12:16]).upper()[:1]
            try:
                partial_charge = float(fields[-2])
            except (IndexError, ValueError):
                partial_charge = None
            atoms.append({
                "line": line.rstrip("\r\n"), "name": line[12:16].strip(), "element": element,
                "atom_type": normalized_type, "partial_charge": partial_charge,
                "x": x, "y": y, "z": z, "heavy": element != "H",
            })
    if not atoms:
        raise ValueError("no PDBQT atoms found in {}".format(path))
    return atoms


def _enumerate_graph_mappings(source_adj, target_adj, source_elements, target_elements, max_mappings=4096):
    """Enumerate element-labelled, induced graph isomorphisms.

    Each returned tuple maps an SDF-heavy-atom index to a PDBQT-heavy-atom
    index.  Degree and both edges and non-edges are checked, so a mere element
    permutation cannot pass as an atom mapping.
    """
    count = len(source_elements)
    domains = []
    for source_index in range(count):
        domains.append([
            target_index for target_index in range(count)
            if (source_elements[source_index] == target_elements[target_index]
                and len(source_adj[source_index]) == len(target_adj[target_index]))
        ])
    if any(not domain for domain in domains):
        return [], False
    order = sorted(range(count), key=lambda i: (len(domains[i]), -len(source_adj[i]), i))
    assigned = {}
    used = set()
    results = []
    truncated = False

    def visit(depth):
        nonlocal truncated
        if len(results) >= max_mappings:
            truncated = True
            return
        if depth == count:
            results.append(tuple(assigned[i] for i in range(count)))
            return
        source_index = order[depth]
        for target_index in domains[source_index]:
            if target_index in used:
                continue
            compatible = True
            for other_source, other_target in assigned.items():
                if ((other_source in source_adj[source_index]) !=
                        (other_target in target_adj[target_index])):
                    compatible = False
                    break
            if not compatible:
                continue
            assigned[source_index] = target_index
            used.add(target_index)
            visit(depth + 1)
            used.remove(target_index)
            del assigned[source_index]
            if truncated:
                return

    visit(0)
    return results, truncated


def _chemical_automorphisms(mol, max_mappings=4096):
    """Return deterministic full-chemistry automorphisms of an RDKit molecule."""
    count = mol.GetNumAtoms()
    matches = mol.GetSubstructMatches(
        mol,
        uniquify=False,
        useChirality=True,
        maxMatches=max_mappings + 1,
    )
    if len(matches) > max_mappings:
        raise ValueError(
            "ATOM_MAPPING_FAILED: more than {} chemical symmetry mappings; "
            "mapping is not safely enumerable".format(max_mappings)
        )
    automorphisms = {tuple(match) for match in matches}
    automorphisms.add(tuple(range(count)))
    return tuple(sorted(automorphisms))


def _pdbqt_atom_type_compatible(source_atom, target_atom):
    """Apply only AutoDock atom-type constraints that are chemically exact.

    In particular, AutoDock type ``A`` is aromatic carbon whereas ``C`` is
    non-aromatic carbon.  Other AutoDock types can merge multiple valence or
    charge states, so they are retained as audit evidence rather than being
    over-interpreted as formal bond orders or formal charges.
    """
    atom_type = str(target_atom.get("atom_type", "")).upper()
    if source_atom.GetAtomicNum() == 6 and atom_type in {"A", "C"}:
        return bool(source_atom.GetIsAromatic()) == (atom_type == "A")
    return True


def _mapped_stereochemistry_is_preserved(mol, source_coords, target_coords, mapping):
    """Validate tetrahedral and double-bond stereo in mapped PDBQT geometry.

    PDBQT omits stereo labels, but docking retains their three-dimensional
    geometry.  This test uses that local handedness/configuration as chemical
    evidence; it never aligns conformers or compares their global positions.
    """
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    specified_centres = [
        atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    ]
    if specified_centres:
        mapped = Chem.Mol(mol)
        for atom in mapped.GetAtoms():
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        conformer = Chem.Conformer(mapped.GetNumAtoms())
        for source_index in range(mapped.GetNumAtoms()):
            point = target_coords[mapping[source_index]]
            conformer.SetAtomPosition(
                source_index, Point3D(float(point[0]), float(point[1]), float(point[2]))
            )
        mapped.RemoveAllConformers()
        mapped.AddConformer(conformer, assignId=True)
        Chem.AssignAtomChiralTagsFromStructure(mapped, confId=0, replaceExistingTags=True)
        for source_index in specified_centres:
            if (mapped.GetAtomWithIdx(source_index).GetChiralTag() !=
                    mol.GetAtomWithIdx(source_index).GetChiralTag()):
                return False

    # RDKit stores the two reference substituents for a specified E/Z bond.
    # Compare their local same-side/opposite-side geometry before and after
    # mapping.  A nearly collinear/undefined target arrangement is not guessed.
    for bond in mol.GetBonds():
        if bond.GetStereo() in {Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY}:
            continue
        stereo_atoms = tuple(bond.GetStereoAtoms())
        if len(stereo_atoms) != 2:
            continue
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()

        def stereo_metric(coords):
            axis = [coords[end][i] - coords[begin][i] for i in range(3)]
            axis_sq = sum(value * value for value in axis)
            if axis_sq <= 1.0e-12:
                return None
            left = [coords[stereo_atoms[0]][i] - coords[begin][i] for i in range(3)]
            right = [coords[stereo_atoms[1]][i] - coords[end][i] for i in range(3)]
            left_dot = sum(left[i] * axis[i] for i in range(3)) / axis_sq
            right_dot = sum(right[i] * axis[i] for i in range(3)) / axis_sq
            left = [left[i] - left_dot * axis[i] for i in range(3)]
            right = [right[i] - right_dot * axis[i] for i in range(3)]
            scale = math.sqrt(
                sum(value * value for value in left) *
                sum(value * value for value in right)
            )
            if scale <= 1.0e-12:
                return None
            return sum(left[i] * right[i] for i in range(3)) / scale

        mapped_coords = [target_coords[mapping[i]] for i in range(mol.GetNumAtoms())]
        source_metric = stereo_metric(source_coords)
        target_metric = stereo_metric(mapped_coords)
        if (source_metric is None or target_metric is None or
                abs(source_metric) < 1.0e-4 or abs(target_metric) < 1.0e-4 or
                source_metric * target_metric < 0.0):
            return False
    return True


def map_sdf_heavy_atoms_to_pdbqt(sdf: os.PathLike, pdbqt_atoms) -> Dict[str, object]:
    """Map reordered PDBQT heavy atoms onto the original SDF topology.

    PDBQT does not carry explicit bond orders or stereo labels.  Connectivity
    is therefore inferred from covalent distances and required to be exactly
    graph-isomorphic to the SDF heavy-atom graph.  Candidates are then filtered
    by exact AutoDock aromatic-carbon typing and by preservation of specified
    tetrahedral/double-bond stereochemistry in the local PDBQT geometry.

    Remaining candidates are assessed by SDF-versus-PDBQT covalent-bond length
    consistency, which supplies bond-order geometry evidence after all discrete
    chemistry constraints.  A deterministic choice is made only when all
    near-best candidates differ by full-chemistry graph automorphisms; otherwise
    ambiguity is reported instead of selecting the first or nearest mapping.

    No rotation, translation, or pose-to-pose fitting is performed.
    """
    from rdkit import Chem

    sdf = Path(sdf)
    supplier = Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True)
    original = next((mol for mol in supplier if mol is not None), None)
    if original is None:
        raise ValueError("ATOM_MAPPING_FAILED: RDKit cannot read SDF {}".format(sdf))
    mol = Chem.RemoveHs(original)
    heavy_atoms = [atom for atom in pdbqt_atoms if bool(atom["heavy"])]
    source_elements = [atom.GetSymbol().upper() for atom in mol.GetAtoms()]
    target_elements = [str(atom["element"]).upper() for atom in heavy_atoms]
    if len(source_elements) != len(target_elements):
        raise ValueError(
            "ATOM_MAPPING_FAILED: SDF has {} heavy atoms but PDBQT has {}".format(
                len(source_elements), len(target_elements)
            )
        )
    if Counter(source_elements) != Counter(target_elements):
        raise ValueError(
            "ATOM_MAPPING_FAILED: heavy-element counts differ: SDF={} PDBQT={}".format(
                dict(sorted(Counter(source_elements).items())),
                dict(sorted(Counter(target_elements).items())),
            )
        )

    count = len(source_elements)
    source_adj = [set() for _ in range(count)]
    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        source_adj[begin].add(end)
        source_adj[end].add(begin)

    periodic_table = Chem.GetPeriodicTable()
    atomic_numbers = [periodic_table.GetAtomicNumber(element.title()) for element in target_elements]
    if any(number <= 0 for number in atomic_numbers):
        raise ValueError("ATOM_MAPPING_FAILED: PDBQT contains an unsupported element")
    radii = [float(periodic_table.GetRcovalent(number)) for number in atomic_numbers]
    target_coords = [
        (float(atom["x"]), float(atom["y"]), float(atom["z"])) for atom in heavy_atoms
    ]

    attempts = []
    mappings = []
    selected_threshold = None
    # 1.20 is the normal covalent-radius tolerance.  Nearby values accommodate
    # rounded PDBQT coordinates while still requiring exact SDF graph topology.
    for threshold in (1.20, 1.15, 1.25, 1.10, 1.30):
        target_adj = [set() for _ in range(count)]
        for i in range(count):
            for j in range(i + 1, count):
                distance = math.sqrt(sum((target_coords[i][axis] - target_coords[j][axis]) ** 2 for axis in range(3)))
                if 0.40 < distance <= threshold * (radii[i] + radii[j]):
                    target_adj[i].add(j)
                    target_adj[j].add(i)
        edge_count = sum(len(neighbours) for neighbours in target_adj) // 2
        attempts.append("{:.2f}:{}edges".format(threshold, edge_count))
        candidate_mappings, truncated = _enumerate_graph_mappings(
            source_adj, target_adj, source_elements, target_elements
        )
        if truncated:
            raise ValueError(
                "ATOM_MAPPING_FAILED: more than 4096 topology mappings; mapping is not safely unique"
            )
        if candidate_mappings:
            mappings = candidate_mappings
            selected_threshold = threshold
            break
    if not mappings:
        raise ValueError(
            "ATOM_MAPPING_FAILED: inferred PDBQT connectivity is not isomorphic to the SDF "
            "(SDF bonds={}; attempts={})".format(mol.GetNumBonds(), ",".join(attempts))
        )

    conformer = mol.GetConformer() if mol.GetNumConformers() else None
    source_coords = []
    if conformer is not None:
        source_coords = [
            (
                float(conformer.GetAtomPosition(i).x),
                float(conformer.GetAtomPosition(i).y),
                float(conformer.GetAtomPosition(i).z),
            )
            for i in range(count)
        ]

    type_valid_mappings = [
        mapping for mapping in mappings
        if all(
            _pdbqt_atom_type_compatible(
                mol.GetAtomWithIdx(source_index), heavy_atoms[target_index]
            )
            for source_index, target_index in enumerate(mapping)
        )
    ]
    if not type_valid_mappings:
        raise ValueError(
            "ATOM_MAPPING_FAILED: no topology mapping preserves exact AutoDock aromatic atom typing"
        )

    stereo_constraint_count = sum(
        atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        for atom in mol.GetAtoms()
    ) + sum(
        bond.GetStereo() not in {Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY}
        for bond in mol.GetBonds()
    )
    if stereo_constraint_count and conformer is None:
        raise ValueError(
            "ATOM_MAPPING_FAILED: SDF specifies stereochemistry but has no coordinates for validation"
        )
    stereo_valid_mappings = [
        mapping for mapping in type_valid_mappings
        if (not stereo_constraint_count or _mapped_stereochemistry_is_preserved(
            mol, source_coords, target_coords, mapping
        ))
    ]
    if not stereo_valid_mappings:
        raise ValueError(
            "ATOM_MAPPING_FAILED: no topology mapping preserves specified stereochemistry"
        )

    def mapping_score(mapping):
        aromatic_mismatches = 0
        for source_index, target_index in enumerate(mapping):
            atom_type = str(heavy_atoms[target_index].get("atom_type", "")).upper()
            source_atom = mol.GetAtomWithIdx(source_index)
            if source_atom.GetSymbol().upper() == "C" and atom_type:
                if bool(source_atom.GetIsAromatic()) != (atom_type == "A"):
                    aromatic_mismatches += 1
        squared = []
        if conformer is not None:
            for bond in mol.GetBonds():
                i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                pi, pj = conformer.GetAtomPosition(i), conformer.GetAtomPosition(j)
                source_distance = math.sqrt((pi.x - pj.x) ** 2 + (pi.y - pj.y) ** 2 + (pi.z - pj.z) ** 2)
                ti, tj = target_coords[mapping[i]], target_coords[mapping[j]]
                target_distance = math.sqrt(sum((ti[axis] - tj[axis]) ** 2 for axis in range(3)))
                squared.append((source_distance - target_distance) ** 2)
        bond_rms = math.sqrt(sum(squared) / len(squared)) if squared else 0.0
        return aromatic_mismatches, bond_rms, tuple(mapping)

    ranked = sorted((mapping_score(mapping), mapping) for mapping in stereo_valid_mappings)
    best_score, best_mapping = ranked[0]
    chemical_automorphisms = _chemical_automorphisms(mol)
    chemical_automorphism_set = set(chemical_automorphisms)
    inverse_best = {target: source for source, target in enumerate(best_mapping)}
    for score, candidate in ranked[1:]:
        if score[0] != best_score[0] or score[1] - best_score[1] > 0.03:
            break
        relative = tuple(inverse_best[candidate[source]] for source in range(count))
        if relative not in chemical_automorphism_set:
            raise ValueError(
                "ATOM_MAPPING_FAILED: multiple chemically non-equivalent mappings remain ambiguous "
                "(best bond RMS={:.4f} A)".format(best_score[1])
            )

    # Convert SDF-index automorphisms into permutations of PDBQT coordinate
    # indices, which fixed-frame RMSD can apply without coordinate fitting.
    pdbqt_permutations = set()
    for automorphism in chemical_automorphisms:
        permutation = list(range(count))
        for source_index in range(count):
            permutation[best_mapping[source_index]] = best_mapping[automorphism[source_index]]
        pdbqt_permutations.add(tuple(permutation))

    return {
        "mol": mol,
        "source_to_pdbqt": tuple(best_mapping),
        "pdbqt_permutations": tuple(sorted(pdbqt_permutations)),
        "topology_mapping_count": len(mappings),
        "atom_type_valid_mapping_count": len(type_valid_mappings),
        "stereo_valid_mapping_count": len(stereo_valid_mappings),
        "stereochemistry_constraint_count": stereo_constraint_count,
        "symmetry_permutation_count": len(pdbqt_permutations),
        "connectivity_threshold": selected_threshold,
        "bond_length_rms_A": best_score[1],
        "aromatic_mismatches": best_score[0],
        "mapping_algorithm": "chemical_graph_v2_stereo_symmetry",
    }


def discover_cas_inputs(project_root: os.PathLike) -> Tuple[List[Dict[str, str]], List[str]]:
    """Discover inputs case-insensitively while rejecting all ambiguous basenames."""
    root = Path(project_root)
    records: List[Dict[str, str]] = []
    errors: List[str] = []

    def unique_file(directory: Path, expected: str, label: str) -> Path:
        matches = sorted(
            (p for p in directory.iterdir() if p.is_file() and p.name.lower() == expected.lower() and p.stat().st_size),
            key=lambda p: p.name,
        )
        if len(matches) != 1:
            raise ValueError("{} expected {!s}: {} non-empty case-insensitive matches ({})".format(
                label, directory / expected, len(matches), ", ".join(p.name for p in matches)
            ))
        return matches[0]

    for cas_dir in sorted(p for p in root.iterdir() if p.is_dir() and CAS_RE.match(p.name)):
        cas = cas_dir.name
        try:
            ligand = unique_file(cas_dir, cas + ".sdf", "ligand")
        except (OSError, ValueError) as exc:
            errors.append("{}: {}".format(cas, exc))
            continue
        receptor_dirs = sorted(p for p in cas_dir.iterdir() if p.is_dir())
        receptor_counts: Dict[str, int] = {}
        for directory in receptor_dirs:
            receptor_counts[directory.name.lower()] = receptor_counts.get(directory.name.lower(), 0) + 1
        for receptor_dir in receptor_dirs:
            receptor = receptor_dir.name
            if receptor_counts[receptor.lower()] > 1:
                errors.append("{}: ambiguous receptor directories differing only by case: {}".format(
                    cas, ", ".join(p.name for p in receptor_dirs if p.name.lower() == receptor.lower())
                ))
                continue
            pockets = sorted(
                p for p in receptor_dir.iterdir()
                if p.is_dir() and p.name.lower().startswith((receptor + "_P_").lower())
            )
            if not pockets:
                errors.append("{}/{}: no pocket directory matching {}_P_*".format(cas, receptor, receptor))
                continue
            pocket_counts: Dict[str, int] = {}
            for directory in pockets:
                pocket_counts[directory.name.lower()] = pocket_counts.get(directory.name.lower(), 0) + 1
            for pocket_dir in pockets:
                pocket = pocket_dir.name
                if pocket_counts[pocket.lower()] > 1:
                    errors.append("{}/{}: ambiguous pocket directories differing only by case: {}".format(
                        cas, receptor, ", ".join(p.name for p in pockets if p.name.lower() == pocket.lower())
                    ))
                    continue
                try:
                    full_pdb = unique_file(pocket_dir, receptor + ".pdb", "complete receptor")
                    pocket_pdb = unique_file(pocket_dir, pocket + ".pdb", "pocket")
                except (OSError, ValueError) as exc:
                    errors.append("{}/{}/{}: {}".format(cas, receptor, pocket, exc))
                    continue
                records.append({
                    "cas": cas, "receptor": receptor, "pocket": pocket,
                    "ligand_sdf": str(ligand.resolve()),
                    "full_receptor_pdb": str(full_pdb.resolve()),
                    "pocket_pdb": str(pocket_pdb.resolve()),
                    "receptor_md5": md5_file(full_pdb),
                })
    return records, errors


def analysis_run_dir() -> Path:
    project = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parent.parent))
    pipeline_id = os.environ.get("PIPELINE_ID")
    if not pipeline_id:
        raise RuntimeError("PIPELINE_ID is required")
    return Path(os.environ.get("PIPELINE_RUN_DIR", project / "_full_docking_analysis" / "runs" / pipeline_id))


def dock_run_dir() -> Path:
    project = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parent.parent))
    tag = os.environ.get("RUN_TAG") or os.environ.get("PIPELINE_ID")
    if not tag:
        raise RuntimeError("RUN_TAG or PIPELINE_ID is required")
    root = Path(os.environ.get("DOCK_RESULT_ROOT", project / "_autodock4_batch"))
    return root / "runs" / tag


def update_stage_status(run_dir: Path, stage: str, status: str, message: str, exit_code: int = 0) -> None:
    from datetime import datetime
    import os
    target = run_dir / "stage_status" / (stage + ".tsv")
    temporary = target.with_name(".{}.tmp.{}".format(target.name, os.getpid()))
    write_tsv(temporary, [{
        "stage": stage, "status": status, "exit_code": exit_code,
        "message": message, "updated_at": datetime.now().astimezone().isoformat(),
    }], ["stage", "status", "exit_code", "message", "updated_at"])
    os.replace(str(temporary), str(target))
