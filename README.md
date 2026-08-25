# Cell number incompletely predicts functional neurogenesis during SSRI treatment

Reproducible code, derived results, and figures for a multiscale computational study of adult dentate-gyrus neurogenesis during SSRI treatment. The model is rodent-informed and does not estimate human neurogenesis or clinical antidepressant response.

## What is included

- A six-state C++17 Gillespie lineage model
- Paired treatment-control uncertainty and sensitivity analyses
- Replicate-aware reanalysis of GSE197622
- Model and transcriptomic cross-layer integration
- Locked result tables, figure data, and reproducibility tests
- Matplotlib, Seaborn, R, and ggplot2 figure code

## Main findings

- 1,500 parameter sets with four stochastic replicates per condition
- Extent versus functional neurogenesis index: Pearson r = 0.7749
- 51 of 1,500 systems showed the prespecified extent-function mismatch
- No transcriptomic test survived global BH FDR below 0.10 across 207 tests
- Seven-process cross-layer association: rho = -0.50, permutation p = 0.268

## Reproduce on Windows

Requirements: Python 3.12+, R 4.3+, a C++17 compiler, Git, and PowerShell 7.

```powershell
python -m pip install -r requirements.txt
Rscript .\scripts\install_r_packages.R
Set-ExecutionPolicy -Scope Process Bypass
.\scripts\run_analysis.ps1
python .\tests\test_invariants.py
```

The workflow downloads the audited public input, compiles the simulator, runs the analyses, regenerates the figures, and checks the locked numerical results.

## Public input data

The transcriptomic analysis uses replicate-level TP10K matrices from Rayan et al., *Molecular Psychiatry* 27 (2022), 4510-4525, [doi:10.1038/s41380-022-01725-1](https://doi.org/10.1038/s41380-022-01725-1), GEO [GSE197622](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197622).

The workflow retrieves the matrices from the [source repository](https://github.com/arulrayan/Integrative-multi-omics-landscape-of-fluoxetine-action-across-27-brain-regions) at audited commit `c3ae3f6d2a98a81d6b0208f19b70fcc56dc50c5f`. They are not redistributed here.

## Citation and license

Citation metadata are provided in `CITATION.cff`. Version 1.0.0 is archived at [Zenodo](https://doi.org/10.5281/zenodo.22096470). Original code, derived tables, and figures are available under the MIT License. The license does not cover separately downloaded source data.
