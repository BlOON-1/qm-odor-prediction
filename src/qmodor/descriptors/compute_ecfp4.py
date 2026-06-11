"""Compute ECFP4 fingerprints from the benchmark molecule table."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .common import (
    InputTableInfo,
    load_benchmark_table,
    require_rdkit,
    resolve_smiles_column,
    write_csv,
    write_json,
)


LOGGER = logging.getLogger(__name__)


def compute_ecfp4(
    project_root: Path,
    input_path: str | None = None,
    radius: int = 2,
    n_bits: int = 1024,
    use_chirality: bool = False,
    use_bond_types: bool = True,
) -> dict[str, Any]:
    """Compute ECFP4 fingerprints and fingerprint-level audit files."""

    chem, _, rdkit_pkg = require_rdkit()
    from rdkit import DataStructs
    from rdkit.Chem import AllChem

    input_info: InputTableInfo = load_benchmark_table(project_root, input_path)
    benchmark_df = resolve_smiles_column(input_info.df)

    processed_dir = project_root / "data" / "processed"
    reports_dir = project_root / "reports"
    processed_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    invalid_rows: list[dict[str, Any]] = []
    index_rows: list[dict[str, Any]] = []
    density_rows: list[dict[str, Any]] = []
    fingerprint_rows = []

    for _, row in benchmark_df.iterrows():
        source_smiles = str(row["smiles_input"]).strip()
        mol = chem.MolFromSmiles(source_smiles) if source_smiles else None
        if mol is None:
            invalid_rows.append(
                {
                    "molecule_id": row["molecule_id"],
                    "cas": row["cas"],
                    "chemical_name": row["chemical_name"],
                    "smiles_input": source_smiles,
                    "failure_reason": "invalid_or_missing_smiles",
                }
            )
            continue

        canonical_smiles = chem.MolToSmiles(mol, canonical=True)
        bitvect = AllChem.GetMorganFingerprintAsBitVect(
            mol,
            radius=radius,
            nBits=n_bits,
            useChirality=use_chirality,
            useBondTypes=use_bond_types,
        )
        arr = np.zeros((n_bits,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(bitvect, arr)
        fingerprint_rows.append(arr)
        row_index = len(fingerprint_rows) - 1
        on_bits = int(bitvect.GetNumOnBits())
        index_rows.append(
            {
                "molecule_id": row["molecule_id"],
                "cas": row["cas"],
                "chemical_name": row["chemical_name"],
                "canonical_smiles": canonical_smiles,
                "row_index": row_index,
            }
        )
        density_rows.append(
            {
                "molecule_id": row["molecule_id"],
                "cas": row["cas"],
                "canonical_smiles": canonical_smiles,
                "on_bits_count": on_bits,
                "bit_density": on_bits / n_bits,
            }
        )

    if fingerprint_rows:
        matrix = sparse.csr_matrix(np.vstack(fingerprint_rows).astype("uint8"))
    else:
        matrix = sparse.csr_matrix((0, n_bits), dtype="uint8")

    npz_output = processed_dir / "ecfp4_fingerprints.npz"
    index_output = processed_dir / "ecfp4_fingerprint_index.csv"
    metadata_output = processed_dir / "ecfp4_fingerprint_metadata.json"
    invalid_output = reports_dir / "ecfp4_invalid_smiles.csv"
    density_output = reports_dir / "ecfp4_density_summary.csv"

    from .common import backup_existing

    backup_existing(npz_output)
    sparse.save_npz(npz_output, matrix, compressed=True)
    write_csv(
        pd.DataFrame(
            index_rows,
            columns=["molecule_id", "cas", "chemical_name", "canonical_smiles", "row_index"],
        ),
        index_output,
    )
    write_csv(
        pd.DataFrame(
            invalid_rows,
            columns=["molecule_id", "cas", "chemical_name", "smiles_input", "failure_reason"],
        ),
        invalid_output,
    )
    write_csv(
        pd.DataFrame(
            density_rows,
            columns=["molecule_id", "cas", "canonical_smiles", "on_bits_count", "bit_density"],
        ),
        density_output,
    )
    metadata = {
        "radius": radius,
        "nBits": n_bits,
        "useChirality": use_chirality,
        "useBondTypes": use_bond_types,
        "rdkit_version": rdkit_pkg.__version__,
        "input_file": input_info.path.relative_to(project_root).as_posix(),
        "n_molecules_total": int(len(benchmark_df)),
        "n_molecules_valid": int(matrix.shape[0]),
        "n_molecules_invalid": int(len(invalid_rows)),
    }
    write_json(metadata, metadata_output)

    LOGGER.info(
        "Computed ECFP4 fingerprints for %s valid molecules; matrix shape=%s.",
        matrix.shape[0],
        matrix.shape,
    )
    return {
        "input_info": input_info,
        "matrix_shape": tuple(int(x) for x in matrix.shape),
        "invalid_df": pd.DataFrame(invalid_rows),
        "index_df": pd.DataFrame(index_rows),
        "metadata": metadata,
        "output_files": [
            npz_output,
            index_output,
            metadata_output,
            invalid_output,
            density_output,
        ],
    }
