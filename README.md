# Quantum-mechanical descriptors improve odor prediction for structurally similar molecules

## Paper

This repository supports the paper _Quantum-mechanical descriptors improve odor prediction for structurally similar molecules_.

## Repository overview

The repository packages the released data products, analysis scripts, configuration files, and supporting documentation used for benchmark odor prediction, structure-disjoint evaluation, external sensory-panel validation, receptor-consistency analysis, and figure source-data export.

QM descriptors are used here as receptor-related molecular proxies that preserve interaction-relevant physicochemical variation among structurally similar molecules. The repository does not claim that QM descriptors directly explain olfactory receptor mechanisms or receptor combinations.

## Main contents

- `data/processed/`: benchmark molecule table, 112-label odor matrix, label vocabulary, RDKit descriptor tables, and ECFP4 fingerprints.
- `data/external_panel/`: 87 external molecules, 60 structural-neighbor pairs, aggregated panel ratings, and the 43-label panel vocabulary.
- `data/receptor_annotations/` and `results/receptor_consistency/`: processed M2OR matching outputs and receptor-consistency summary tables.
- `data/source_data/`: source data files for Fig. 1-4.
- `results/`: released benchmark metrics, structure-disjoint summaries, SHAP ranking output, and model-summary JSON files.
- `scripts/`, `notebooks/`, `configs/`, and `docs/`: analysis-stage code, configuration, and reproducibility notes.

## Quick start

Use the environment file provided in this repository.

```bash
conda env create -f environment.yml
conda activate qmodor
pip install -e .
```

## Reproducing analyses

Analysis scripts and notebooks are organized by analysis stage rather than by a single one-click entry point.

- `scripts/prepare_available_data.py` prepares released benchmark, panel, and source-data tables from curated inputs.
- `scripts/01_compute_rdkit_features.py` computes RDKit descriptor features.
- `scripts/03_make_splits.py` and `src/qmodor/splits/` cover split-related utilities.
- `scripts/07_receptor_consistency_analysis.py` and `notebooks/05_receptor_consistency_check.ipynb` cover receptor-consistency outputs.
- `notebooks/01_dataset_overview.ipynb` through `notebooks/04_external_panel_validation.ipynb` summarize the released benchmark, descriptor, benchmark-metric, and external-panel assets.
- `docs/workflow_overview.md` describes the full repository workflow and reproducibility boundaries.

## Data notes

The released benchmark package contains `benchmark_molecules_4495.csv`, `odor_label_vocabulary_112.csv`, and `odor_label_matrix_4495x112.csv`. The external validation package contains `external_molecules_87.csv`, `structural_neighbor_pairs_60.csv`, and `panel_label_vocabulary_43.csv`.

Some upstream information used during dataset construction comes from The Good Scents Company and Leffingwell and may be subject to provider-specific redistribution terms. See `data/raw_manifest/` and `docs/data_availability.md` for release notes and restrictions.

Public olfactory receptor response annotations come from M2OR. Processed receptor matching and consistency tables are released under `data/receptor_annotations/` and `results/receptor_consistency/`.

## Citation

Please cite the associated paper and the repository metadata in `CITATION.cff` when using this package.

## License

This repository is distributed under the terms described in `LICENSE`.
