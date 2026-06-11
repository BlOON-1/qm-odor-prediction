"""Compatibility wrapper for receptor concordance preparation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .prepare_m2or_receptor_annotations import prepare_m2or_receptor_annotations


def compute_receptor_pair_concordance(project_root: Path) -> dict[str, Any]:
    """Compute pair-level receptor concordance tables from M2OR pairwise data."""

    return prepare_m2or_receptor_annotations(project_root)

