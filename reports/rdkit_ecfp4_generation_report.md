# RDKit and ECFP4 Generation Report

## 1. Input table used
- `data/processed/benchmark_molecules_observed.csv`

## 2. Molecule count
- total molecules: `4495`

## 3. Valid molecule count
- valid molecules: `4494`

## 4. Invalid SMILES count
- invalid SMILES: `1`

## 5. RDKit descriptor count before filtering
- descriptor count before filtering: `217`

## 6. RDKit descriptor count after missing/constant filtering
- descriptor count after missing/constant filtering: `191`

## 7. Whether split files were found
- split files found: `False`

## 8. Whether standardized descriptors were generated
- standardized descriptors generated: `False`

## 9. ECFP4 settings: radius=2, nBits=1024
- radius: `2`
- nBits: `1024`
- matrix shape: `(4494, 1024)`

## 10. Files generated
- `data/processed/ecfp4_fingerprint_index.csv`
- `data/processed/ecfp4_fingerprint_metadata.json`
- `data/processed/ecfp4_fingerprints.npz`
- `data/processed/rdkit_descriptors_raw.csv`
- `data/processed/rdkit_descriptors_unstandardized_qc.csv`
- `reports/ecfp4_density_summary.csv`
- `reports/ecfp4_invalid_smiles.csv`
- `reports/rdkit_descriptor_audit.csv`
- `reports/rdkit_descriptor_constant_columns.csv`
- `reports/rdkit_descriptor_missingness.csv`
- `reports/rdkit_invalid_smiles.csv`
- `reports/rdkit_standardization_blocked_no_splits.md`
- `reports/rdkit_ecfp4_generation_report.md`

## 11. Files not generated and why
- data/processed/standardized/rdkit_descriptors_standardized_seed42_train.csv and related seed-specific outputs: split files were not available.
- data/processed/feature_scaler_parameters_seed*.csv: scaler fitting on all molecules was not performed to avoid data leakage.

## 12. Critical warnings
- Invalid or missing SMILES were encountered for 1 molecules.
- Train/test split files were not found, so finalized standardized RDKit descriptor matrices were not generated.
