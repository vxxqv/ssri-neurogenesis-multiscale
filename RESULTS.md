# Results in plain language

## The question

Counting new cells is useful, but a new neuron must also survive, integrate, and function. We tested how often those two ideas agree during simulated SSRI treatment, then asked whether public gene-expression data support the same distinction.

## How the study works

The lineage model follows six cell states for 56 days. Each treated system is paired with its untreated version. We sampled 2,048 plausible parameter sets and ran 16 paired stochastic replicates per set. A separate deterministic solver, exact attribution, global sensitivity analysis, and held-out prediction tested the model from different angles.

Six public transcriptomic studies provided independent biological context. Gene programs came from Gene Ontology before the main tests. Sample labels were permuted exactly where possible, and false-discovery rates were calculated across all tests and within each dataset.

## What we found

Cell number and the functional index were related, with Pearson r = 0.785, but not interchangeable. In 10-fold held-out prediction, cell number alone explained 47.0% of functional variation. Cell composition explained 59.2%, while process-aware features explained 86.4%.

The mismatch probability peaked at day 14. This suggests that timing matters: an early rise in cells may appear before their later integration and efficacy are clear.

The strongest independent result appeared in GSE309750 mossy cells. Fluoxetine produced a large integration-over-extent shift, Hedges g = 3.76, exact p = 0.00117, within-dataset q = 0.00668. The direction remained positive under three scoring methods, after removal of every single sample, and in both halves of all 200 random gene splits.

In GSE43261 dorsal dentate gyrus, responders also showed a positive integration-over-extent shift relative to resistant animals, Hedges g = 1.86. Its exact p was 0.00909 and global q was 0.0418, but within-dataset q was 0.183. We therefore treat this result as exploratory, not confirmed.

The broader datasets showed strong context dependence across cell type, brain region, experimental system, and SSRI. There was no consistent class-wide increase in the integration program. This argues against describing all SSRIs or all hippocampal samples with one transcriptomic effect.

## What this does not prove

The model is a controlled thought experiment, not a complete brain. The public analyses reuse animal data collected for other questions. They do not establish causality, human neurogenesis, or clinical response. The integration-over-extent score was developed after inspecting the first program results, so it is openly labeled exploratory even though its robustness checks were strong.

## Public studies

The analysis credits Rayan et al. ([GSE197622](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197622), [article](https://doi.org/10.1038/s41380-022-01725-1)), Samuels et al. ([GSE43261](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43261), [article](https://doi.org/10.1371/journal.pone.0085136)), Oh et al. ([GSE309750](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309750), [article](https://doi.org/10.1038/s41380-026-03461-2)), Rayan et al. ([GSE222756](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222756), [article](https://doi.org/10.1038/s41380-024-02619-0)), Demin et al. ([GSE205325](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205325), [article](https://doi.org/10.1038/s41598-022-22688-x)), and Yamamoto and Abe ([GSE292948](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292948), [article](https://doi.org/10.1016/j.isci.2025.113800)).
