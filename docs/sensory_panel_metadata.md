# Sensory Panel Metadata

## 1. Purpose of external panel

The external sensory panel provides an out-of-benchmark perceptual comparison for structurally similar fragrance molecules. Its role is to assess whether molecular representations and model-derived comparisons remain consistent with aggregated human sensory judgments in a controlled external setting.

## 2. External molecule set

The manuscript-aligned external set contains 87 fragrance molecules selected for structural relevance to the study objectives. These molecules are handled as a distinct evaluation resource and are not treated as interchangeable with the benchmark training corpus.

## 3. Structural-neighbor pair selection

Sixty near-neighbor pairs are selected from the external molecule set to emphasize perceptually challenging comparisons among structurally similar molecules. Pair definitions and selection metadata are released as a processed table.

## 4. Sensory vocabulary

The panel uses a 43-label sensory vocabulary derived from the benchmark 112-label model vocabulary. This reduced vocabulary is intended to support consistent assessor use while preserving alignment with the benchmark semantic space where possible.

## 5. Assessor panel

The panel consists of 13 trained assessors. Testing is conducted under standardized test conditions with controlled presentation procedures. Public release materials document panel size, vocabulary, and aggregation strategy, but not assessor-identifying information.

## 6. Aggregation of responses

For each molecule, assessors select applicable labels and assign intensity scores under the standardized protocol. Responses are aggregated across assessors to produce anonymized molecule-level summaries such as mean label intensity and label selection frequency.

## 7. Top-K dominant-label extraction

Dominant labels are derived from ordered mean-intensity profiles. Labels are retained in descending order until a knee point is reached, yielding a panel-derived top-K label set for each molecule. This summary representation is used for pairwise perceptual similarity analyses.

## 8. Pairwise perceptual similarity

Pairwise perceptual similarity is defined as the Jaccard similarity between panel-derived label sets for the two molecules in a pair. Release tables provide pair identifiers, derived label sets or references to them, and the corresponding Jaccard similarity values for the 60 comparison pairs.

## 9. Privacy and ethics-related handling

Only anonymized and aggregated summaries are intended for public release. Non-anonymized assessor-level records are not planned for repository distribution. Release materials remain focused on reproducible molecule-level outcomes while minimizing disclosure of personally identifiable or sensitive participation information.
