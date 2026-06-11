# Structure-Disjoint Split Report

## 1. Input files inspected
- `split_ecfp4_by_folds.py` exists: `True`
- `fold_01.csv` exists: `True`
- `test_01.csv` exists: `True`

## 2. Role convention detected from split_ecfp4_by_folds.py
- expected input structure: base_dir with fold_* directories, each containing fold_<xx>.csv and test_<xx>.csv, plus ECFP4_data.csv
- treats fold_01.csv as train: `True`
- treats test_01.csv as test: `True`
- checks CAS overlap: `True`
- removes overlapping CAS from training set: `True`
- expects ECFP4_data.csv: `True`

## 3. Row counts and unique CAS counts
- fold_01.csv rows: `2142`, unique CAS: `2142`, missing CAS: `0`, duplicate CAS: `0`
- test_01.csv rows: `1047`, unique CAS: `1047`, missing CAS: `0`, duplicate CAS: `0`

## 4. Train/test CAS overlap
- overlap count: `0`

## 5. Missing CAS in benchmark or ECFP4 index if checked
- missing in benchmark train: `0`
- missing in benchmark test: `0`
- missing in ECFP4 index train: `0`
- missing in ECFP4 index test: `0`

## 6. Final role convention used
- script_convention_role: `{'fold_01.csv': 'train', 'test_01.csv': 'test'}`
- user_statement_role: `{'fold_01.csv': 'test', 'test_01.csv': 'train'}`
- inferred_by_row_count: `{'fold_01.csv': 'train', 'test_01.csv': 'test'}`
- final_role: `{'fold_01.csv': 'train', 'test_01.csv': 'test'}`
- final_role_source: `script_convention`

## 7. Critical warnings
- CRITICAL: user-stated split roles conflict with the legacy script convention; canonical outputs follow script convention unless split_role_override.yaml explicitly overrides it.

## 8. Files generated
- `src/qmodor/splits/legacy_split_ecfp4_by_folds.py`
- `scripts/legacy_split_ecfp4_by_folds.py`
- `data/splits/structure_disjoint_splits/structure_disjoint_train.csv`
- `data/splits/structure_disjoint_splits/structure_disjoint_test.csv`
- `data/splits/structure_disjoint_splits/fold_01_train.csv`
- `data/splits/structure_disjoint_splits/fold_01_test.csv`
- `data/splits/structure_disjoint_splits/structure_disjoint_split_audit.csv`
