# Workflow Overview

## ASCII pipeline

```text
raw source inventory
  -> identifier standardization
  -> label harmonization
  -> descriptor generation
  -> feature preprocessing
  -> split construction
  -> model training
  -> benchmark evaluation
  -> external sensory-panel validation
  -> receptor-consistency analysis
  -> figure/source-data generation
```

## 1. Scientific question

The repository supports a study asking whether quantum-mechanical descriptors provide useful molecular information beyond standard cheminformatic representations for multi-label odor prediction and for discrimination among structurally similar molecules. A secondary goal is to test whether descriptor-space relationships are consistent with external sensory judgments and with public receptor-response annotations used as receptor-related molecular proxies. These analyses are interpretation-limited and are not presented as direct evidence of olfactory receptor mechanism or receptor combinatorial coding.

## 2. Repository-level workflow

The full workflow proceeds from source record curation through processed benchmark construction, descriptor generation, model development, external validation, orthogonal consistency analysis, and export of source data for figures and tables. The repository is designed to expose reproducible intermediate products where redistribution boundaries permit public release.

## 3. Dataset construction

Benchmark molecules are assembled from The Good Scents Company and Leffingwell, aligned by CAS number, and secondarily verified against molecular structure and canonical SMILES. Records with unresolved identifier-structure conflicts are excluded. The manuscript-aligned release contains a standardized 112-label multi-hot target matrix for 4,495 odorant molecules. The raw checked local snapshot still contains 113 labels; the public model-training release excludes `odorless` so that the released label matrix matches the archived benchmark-summary JSON files used for Fig. 2.

## 4. Molecular representation generation

The study compares ECFP4 fingerprints, RDKit-computed 2D descriptors, QM descriptors, and fused RDKit+QM descriptors. The QM workflow uses PubChem 3D starting structures, GaussView inspection, Gaussian 16W geometry optimization and frequency analysis at B3LYP/6-311+G(2d,p) with D3(BJ), and Multiwfn parsing. Descriptor categories cover electronic energetics, electrostatic landscape, polarizability, and aromaticity. All continuous features are intended to be standardized using training-set parameters only.

### Descriptor generation details

- RDKit continuous descriptors are generated from canonical SMILES.
- ECFP4 fingerprints are generated as 1024-bit Morgan fingerprints with radius 2.
- ECFP4 fingerprints are not z-score standardized.
- Continuous descriptor scaling is fit on training data only.

## 5. Model evaluation

Models include logistic regression, support vector machine, random forest, XGBoost, convolutional neural network, multilayer perceptron, graph convolutional network, an ADCH-only graph model, and a GCN+ADCH graph model. Benchmark evaluation uses iterative multi-label stratification with a 7:3 train-test split and repeated random seeds 42 to 46. A separate structure-disjoint setting is defined using ECFP4 Tanimoto similarity, with high-similarity structural neighbors identified by similarity greater than 0.8.

## 6. External validation

External validation uses 87 structurally similar fragrance molecules arranged into 60 near-neighbor pairs. A 43-label sensory vocabulary derived from the benchmark vocabulary is applied by 13 trained assessors under standardized test conditions. Responses are aggregated into molecule-level sensory profiles, dominant label sets, and pairwise Jaccard similarity summaries for perceptual comparison.

## 7. Receptor-consistency analysis

Public M2OR olfactory-receptor response annotations are matched by CAS number and curated for duplicate consistency. Only jointly tested receptors are retained when comparing molecules within a pair. Pair-level concordance or divergence is then summarized from matched receptor subsets and compared with descriptor-space distances. This analysis is an orthogonal consistency analysis that examines whether interaction-relevant physicochemical variation is consistent with public receptor-response annotations; it is not a receptor-resolved mechanism study and does not directly model receptor combinations.

## 8. Figure reproduction

Main and supplementary figure source data are exported from processed descriptor matrices, model outputs, split files, panel summaries, and receptor-consistency tables where those assets are available in the current public package.

## 9. Reproducibility boundaries

Reproducibility is bounded by third-party data access terms, external commercial software dependencies, and the practical size of QM intermediate outputs. The repository therefore emphasizes processed releases, documented reconstruction rules, explicit separation between redistributable derived materials and restricted upstream content, and non-commercial reuse of the released package.
