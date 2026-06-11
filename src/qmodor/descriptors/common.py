"""Shared utilities for RDKit-based representation generation."""

from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


class DependencyError(RuntimeError):
    """Raised when an optional scientific dependency is unavailable."""


@dataclass
class InputTableInfo:
    """Metadata describing the selected benchmark molecule table."""

    path: Path
    df: pd.DataFrame
    sha256: str


def ensure_parent(path: Path) -> None:
    """Create the parent directory of a path if it does not exist."""

    path.parent.mkdir(parents=True, exist_ok=True)


def backup_existing(path: Path) -> Path | None:
    """Back up an existing file before writing a replacement."""

    if not path.exists():
        return None
    idx = 1
    candidate = path.with_suffix(path.suffix + ".bak")
    while candidate.exists():
        idx += 1
        candidate = path.with_suffix(path.suffix + f".bak{idx}")
    shutil.copy2(path, candidate)
    return candidate


def write_csv(df: pd.DataFrame, path: Path) -> Path | None:
    """Write a dataframe as UTF-8-SIG CSV with automatic backup."""

    ensure_parent(path)
    backup_path = backup_existing(path)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return backup_path


def write_json(payload: dict[str, Any], path: Path) -> Path | None:
    """Write a JSON payload with backup support."""

    ensure_parent(path)
    backup_path = backup_existing(path)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return backup_path


def write_text(text: str, path: Path) -> Path | None:
    """Write a text file with backup support."""

    ensure_parent(path)
    backup_path = backup_existing(path)
    path.write_text(text, encoding="utf-8")
    return backup_path


def compute_sha256(path: Path) -> str:
    """Compute a SHA256 digest for a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_benchmark_table(project_root: Path, explicit_input: str | None = None) -> InputTableInfo:
    """Load the preferred benchmark molecule table according to project rules."""

    if explicit_input:
        input_path = (project_root / explicit_input).resolve() if not Path(explicit_input).is_absolute() else Path(explicit_input)
        candidates = [input_path]
    else:
        candidates = [
            project_root / "data" / "processed" / "benchmark_molecules_4403.csv",
            project_root / "data" / "processed" / "benchmark_molecules_observed.csv",
        ]
    for candidate in candidates:
        if candidate.exists():
            df = pd.read_csv(candidate)
            return InputTableInfo(path=candidate, df=df, sha256=compute_sha256(candidate))
    searched = "\n".join(str(path) for path in candidates)
    raise FileNotFoundError(
        "No benchmark molecule table was found. Checked:\n" + searched
    )


def resolve_smiles_column(df: pd.DataFrame) -> pd.DataFrame:
    """Resolve a canonical SMILES source column while preserving identifiers."""

    required_columns = ["molecule_id", "cas"]
    for column in required_columns:
        if column not in df.columns:
            raise ValueError(f"Required column `{column}` is missing from benchmark table.")

    working = df.copy()
    if "chemical_name" not in working.columns:
        working["chemical_name"] = ""
    if "canonical_smiles" not in working.columns:
        working["canonical_smiles"] = ""
    if "smiles" not in working.columns:
        working["smiles"] = ""
    if "smiles_raw" not in working.columns:
        working["smiles_raw"] = ""

    working["chemical_name"] = working["chemical_name"].fillna("").astype(str)
    working["canonical_smiles"] = working["canonical_smiles"].fillna("").astype(str)
    working["smiles"] = working["smiles"].fillna("").astype(str)
    working["smiles_raw"] = working["smiles_raw"].fillna("").astype(str)
    working["smiles_input"] = working["canonical_smiles"]
    fallback = working["smiles_input"].eq("")
    working.loc[fallback, "smiles_input"] = working.loc[fallback, "smiles"]
    fallback = working["smiles_input"].eq("")
    working.loc[fallback, "smiles_input"] = working.loc[fallback, "smiles_raw"]
    return working


def utc_now_iso() -> str:
    """Return a UTC timestamp for audit metadata."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_rdkit() -> tuple[Any, Any, Any]:
    """Import RDKit or raise a clear installation error."""

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        import rdkit

        return Chem, Descriptors, rdkit
    except Exception as exc:  # pragma: no cover - import guard
        raise DependencyError(
            "RDKit is required for descriptor generation but is not available. "
            "Please install the project environment from environment.yml."
        ) from exc


def require_sklearn() -> Any:
    """Import StandardScaler or raise a clear installation error."""

    try:
        from sklearn.preprocessing import StandardScaler

        return StandardScaler
    except Exception as exc:  # pragma: no cover - import guard
        raise DependencyError(
            "scikit-learn is required for descriptor standardization but is not available. "
            "Please install the project environment from environment.yml."
        ) from exc


def clean_numeric_frame(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Replace inf values with NaN in a numeric dataframe and count replacements."""

    cleaned = df.replace([np.inf, -np.inf], np.nan)
    inf_counts = pd.Series(0, index=df.columns, dtype="int64")
    for column in df.columns:
        values = df[column].to_numpy()
        inf_counts.loc[column] = int(np.isinf(values).sum()) if np.issubdtype(values.dtype, np.number) else 0
    return cleaned, inf_counts

