# Minimum reproducible release

## Main conclusions that can be recalculated now

- The released 4,495 × 112 benchmark label matrix and ECFP4 indexing consistency.
- All 60 published panel-pair Jaccard values from their published label sets.
- All 30 triplet summary counts from the 390 anonymized responses.
- Representative PDB file integrity against the provided SHA-256 index, and PLIP frequency arithmetic (although docking metadata itself needs portable redaction before release).

## Conditions only, not currently reproducible

- QM/ADCH and fused-feature modeling claims, prediction probabilities, model metrics, SHAP analyses, and seed-level split/model reruns: required feature/output assets are missing.
- Figure 3 assessor-level analyses: only aggregate panel values are present.
- Docking validation claims for 8F76/8UXY/9WG4: validation records are absent.

## Extra assets for a complete rerun

Release ID-aligned QM/ADCH/fusion matrices, calculation-failure log, split seeds, model code/weights or deterministic training instructions, per-molecule predictions/labels, SHAP outputs, anonymized panel raw ratings, and validation inputs/outputs. Preserve third-party raw materials only where redistribution rights permit; otherwise provide rebuild instructions and immutable source/version hashes.
