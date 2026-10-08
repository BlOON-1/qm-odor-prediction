# Final data-integrity audit

**Audit date:** 2026-10-08  
**Scope:** 340 release-candidate files under data/, results/, scripts/, src/, configs/, docs/, notebooks/, tests/ and root release metadata. Excluded internal audit/checklist folders and Git internals. No source data were modified.

## Overall decision: BLOCKED

The internally cross-checkable benchmark, triplet, and representative-structure subsets pass their stated checks, but the release is not ready as a complete paper data package: QM/ADCH/model-output assets are absent, panel raw ratings are absent, and portable docking metadata still contains absolute local paths. Git state could not be read (`git status` and `git rev-parse` returned “not a git repository”), despite a `.git` directory being present.

## File status counts

| PASS | MISSING | FAIL | UNVERIFIED | NOT_FOR_RELEASE | total |
|---:|---:|---:|---:|---:|---:|
| 157 | 0 (file-system manifest cannot enumerate absent files) | 0 | 183 | 0 | 340 |

`MISSING` planned assets are listed in the issue document rather than fabricated as manifest rows.

## Reproduction and consistency tests

| Test | status | evidence |
|---|---|---|
| Benchmark 4495×112 matrix | PASS | 4495 molecule IDs/order, 112 labels, and binary matrix are mutually consistent; odorless is absent from the vocabulary. |
| ECFP4 index/array | PASS | ECFP4 CSR matrix is 4494×1024 and index has 4494 rows; the sole omitted MOL_003034 has explicitly unparseable SMILES, matching metadata (4495 total/1 invalid). |
| External sensory-panel structure/Jaccard | PASS | 87 unique panel molecules, 60 valid pairs, 43 terms; all 60 published Jaccard values recompute exactly. Individual raw panel ratings are not present. |
| Triplet response counts | PASS | 30 triplets × 13 unique assessors = 390 valid codes; every summary count recomputes. Original Excel is not included for external source comparison. |
| Docking representative PDB index | PASS | 36 receptor-pocket jobs; 143 indexed structures have matching SHA-256 and coordinates. Required two core representatives/job=72 present; optional cluster medoids=71. |
| PLIP relationships/frequencies | PASS | 3743 PLIP tasks (3743 success); all interaction-summary pose keys map to successful tasks; all frequency rows recompute from stated denominators. |
| Model aggregate metrics | FAIL | [Errno 2] No such file or directory: 'D:\\Desktop\\project\\qm-odor-prediction\\results\\structure_disjoint\\fig2E_structure_disjoint_audit.csv' |

## Findings and unresolved boundaries

| status | item | evidence | effect |
|---|---|---|---|
| FAIL | `Model aggregate metrics` | [Errno 2] No such file or directory: 'D:\\Desktop\\project\\qm-odor-prediction\\results\\structure_disjoint\\fig2E_structure_disjoint_audit.csv' | May block the linked scientific result |
| UNVERIFIED | `Raw panel records` | Only aggregate 87×43 panel table is released; no 87×43×13 assessor-level panel matrix was found. | Figure 3 raw-rating provenance cannot be independently checked. |
| UNVERIFIED | `QM/ADCH descriptors and prediction probabilities` | QM-valued structure-disjoint tables (3,189 rows) and aggregate model/SHAP summaries exist, but no standalone full QM/ADCH/fused feature matrices, trained weights, or per-molecule probabilities were found. | Aggregate benchmark means can be checked, but electronic-feature models cannot be fully rerun. |
| NOT_FOR_RELEASE | `Docking provenance tables` | Several docking TSVs may retain absolute external paths. This is a portability/privacy disclosure risk. | Redact paths or publish portable replacements only. |

## Interpretation of collection versions

The formal model dataset is 4,495 molecules × 112 labels: `odorless` is not in the released 112-term vocabulary. “Observed” files are retained as a different collection state and are not silently forced to 4,495 rows. No 4,403-model package was found.

## Release conditions

Only manifest rows marked PASS and free of release-policy issues are candidates for direct publication. PASS proves the explicit checks above, not scientific validity of docking or a third-party redistribution license. Model/QM claims remain **BLOCKED** pending release or documented exclusion of the missing data products.
