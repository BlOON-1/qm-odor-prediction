# Data Availability

## 1. Overview

This repository supports a manuscript on quantum-mechanical descriptor-enhanced molecular representations for odor prediction. Code, documentation, and selected processed outputs are organized for public release with explicit attention to reproducibility and third-party redistribution boundaries. The manuscript-aligned release contains 4,495 curated benchmark molecules, 112 standardized odor labels used for modeling, 87 structurally similar external molecules arranged into 60 comparison pairs, and a 43-label sensory vocabulary used for aggregated panel analysis.

## 2. Data categories targeted for release

The repository is intended to distribute the following processed or derived materials, subject to final curation and file-size considerations:

- processed benchmark molecule table
- odor label vocabulary and 112-dimensional label matrix
- label harmonization table
- RDKit descriptors
- ECFP4 fingerprints
- train-test split indices
- structure-disjoint split files
- external 87-molecule sensory-panel summary data
- 60 structural-neighbor pair list
- 43-label sensory vocabulary
- panel-derived pairwise Jaccard similarity
- M2OR matched receptor annotations
- receptor-response concordance and divergence tables
- source data for main and supplementary figures where available

Public release content therefore focuses on processed and harmonized outputs for which redistribution is believed to be appropriate, or otherwise provides instructions for reconstruction from legally obtained source access where needed.

## 3. Third-party source restrictions

Some upstream information used during dataset construction originates from third-party sources, including The Good Scents Company, Leffingwell, PubChem, and M2OR, and may be subject to separate provider-specific usage terms or redistribution restrictions. This repository does not imply that third-party raw source records can be freely redistributed.

Users are responsible for complying with the terms of the original data providers and any applicable software or database licenses.

## 4. What is not included in the current public package

The following items are not included in the current public repository snapshot:

- raw proprietary source records if redistribution is restricted
- commercial software binaries
- non-anonymized assessor-level sensory records
- large raw Gaussian output and QM intermediate files generated with commercial software

## 5. Archive strategy

- GitHub hosts code, environment files, usage notes, and lightweight processed tables intended for non-commercial research and educational use.
- Zenodo or an equivalent release archive should be used for versioned public release, archival citation, and larger processed outputs where appropriate.
- Release packages should distinguish repository-generated processed data from third-party source materials and should document any omissions required by redistribution boundaries.

## 6. License and citation note

Repository code, documentation, and processed outputs are released for non-commercial research and educational use only. Users should cite both the associated manuscript and the software archive when those records become available.

## Current local preparation status

- Manuscript-aligned benchmark size: `4495` molecules.
- Manuscript-aligned target odor vocabulary: `112` labels.
- Manuscript-aligned external validation set: `87` molecules arranged into `60` structural-neighbor pairs.
- Manuscript-aligned sensory vocabulary: `43` labels.
- Current checked raw benchmark snapshot contains `4495` rows in `data/processed/benchmark_molecules_observed.csv`.
- Current checked raw odor-vocabulary snapshot contains `113` rows in `data/processed/odor_label_vocabulary_observed.csv`.
- The current manuscript-aligned 112-label release excludes the `odorless` label from the raw 113-label snapshot so that the released benchmark label matrix matches the model-training vocabulary used in the archived summary JSON files.
- QM descriptor tables and QM intermediate outputs are not included in the current repository because they depend on commercial software and restricted local workflows.
