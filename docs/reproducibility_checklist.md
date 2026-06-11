# Reproducibility Checklist

This checklist is intended to support strict, audit-oriented review of computational and release readiness.

## 1. Environment

- [ ] environment created from `environment.yml`
- [ ] package versions recorded
- [ ] external software dependencies documented separately where required

## 2. Data integrity

- [ ] molecule IDs are stable across tables
- [ ] identifier-structure conflicts are excluded or explicitly documented
- [ ] label vocabulary version is fixed for the release

## 3. Descriptor generation

- [ ] no test-set information is used during preprocessing
- [ ] scaler parameters estimated from training set only
- [ ] QM descriptor extraction settings documented

## 4. Split construction

- [ ] iterative stratification split files saved
- [ ] structure-disjoint split files saved
- [ ] random seeds 42–46 recorded

## 5. Model training

- [ ] representation type is recorded for every trained model
- [ ] hyperparameter search space or fixed settings are archived
- [ ] training logs or summary metrics are exportable from the released workflow

## 6. Benchmark evaluation

- [ ] benchmark metric definitions are documented
- [ ] repeated-seed aggregation procedure is documented
- [ ] source data for benchmark comparison figures are exported

## 7. External validation

- [ ] external panel data are anonymized and aggregated
- [ ] external molecule and pair identifiers are stable across summary tables
- [ ] panel-derived top-label extraction rule is documented

## 8. Receptor-consistency analysis

- [ ] M2OR matching exclusions documented
- [ ] jointly tested receptor filtering rule applied consistently
- [ ] concordant, divergent, and intermediate classification logic recorded

## 9. Figure generation

- [ ] all main-figure source data exported
- [ ] supplementary figure and table source data exported where applicable
- [ ] figure scripts reference only released or reconstructable inputs

## 10. Release readiness

- [ ] Zenodo DOI added after release
- [ ] third-party restrictions documented
- [ ] license notices checked against included content categories
