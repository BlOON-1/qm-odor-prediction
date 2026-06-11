"""Entry point for organizing structure-disjoint split files."""

from __future__ import annotations

import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qmodor.splits.organize_structure_disjoint_split import organize_structure_disjoint_split


def main() -> int:
    result = organize_structure_disjoint_split(PROJECT_ROOT)
    print(f"split train rows: {result['train_rows']}")
    print(f"split test rows: {result['test_rows']}")
    print(f"train/test CAS overlap: {result['train_test_overlap']}")
    print(f"final role: {result['final_role']}")
    print("warnings:")
    for warning in result["warnings"]:
        print(f"- {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

