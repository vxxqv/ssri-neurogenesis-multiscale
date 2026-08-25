# Cell number incompletely predicts functional neurogenesis during SSRI treatment

Reproducible code and derived results for an original computational study of adult dentate-gyrus neurogenesis during SSRI treatment.

## Scientific question

Does a fluoxetine-like increase in late neurogenic-lineage cell number necessarily imply an increase in functional integration? The project combines:

1. a six-state C++ Gillespie model of adult dentate-gyrus neurogenesis;
2. paired treatment-control uncertainty experiments and sensitivity analysis;
3. an independent, replicate-aware reanalysis of GSE197622; and
4. a process-level model-versus-omics convergence test.

The model is rodent-informed and does not estimate human neurogenesis or clinical antidepressant response.

## Locked findings

- 1,500 parameter sets, four stochastic replicates per condition, 56 simulated days.
- Delta extent versus delta FNI: Pearson r = 0.7749; Spearman rho = 0.7100.
- Prespecified mismatch: 51/1,500 systems (3.4%).
- No module-by-cell-type test survived global BH FDR < 0.10 across 207 tests.
- Seven-process cross-layer association: rho = -0.50; permutation p = 0.268.

## Repository map

```text
config/       frozen parameter priors and gene modules
cpp/          ISO C++17 stochastic simulation engine
python/       design, analysis, integration, and figure code
results/      locked machine-readable outputs and figure data
figures/      publication-quality PNG files plus graphical abstract
scripts/      Windows PowerShell reproduction helpers
tests/        release invariants and numerical consistency checks
```

## Public input data

The transcriptomic reanalysis uses the source authors' replicate-level TP10K matrices from:

- Rayan et al., *Molecular Psychiatry* 27 (2022), 4510-4525. DOI: 10.1038/s41380-022-01725-1
- GEO accession: [GSE197622](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197622)
- Upstream data materials used to retrieve the GSE197622 matrices, not this study's repository: <https://github.com/arulrayan/Integrative-multi-omics-landscape-of-fluoxetine-action-across-27-brain-regions>
- Audited source commit: `c3ae3f6d2a98a81d6b0208f19b70fcc56dc50c5f`

The upstream normalized matrices are not redistributed here because no explicit reuse license was visible in that repository at analysis time. The reproduction script clones the public source and checks out the audited commit.

## Reproduce on Windows

Requirements: Python 3.12+, the packages in `requirements.txt`, R 4.3+ with ggplot2, patchwork, and scales, a C++17 compiler (`g++`), Git, and PowerShell 7.

```powershell
python -m pip install -r requirements.txt
Rscript .\scripts\install_r_packages.R
Set-ExecutionPolicy -Scope Process Bypass
.\scripts\run_analysis.ps1
python .\tests\test_invariants.py
```

The full run creates the design, compiles the simulator, executes all stochastic experiments, reanalyzes the public matrices, performs cross-layer testing, and regenerates every figure. It is deterministic under the recorded seeds. Runtime depends on CPU speed.

If Rscript is not on `PATH`, provide its full Windows path. To use an already cloned upstream repository:

```powershell
.\scripts\run_analysis.ps1 -Rscript "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -RayanRepoPath "D:\data\rayan2022" -SkipFetch
```

## Statistical decisions

- Biological replicate, not cell, is the transcriptomic inference unit.
- Treatment effects are Hedges' g with exact permutations when feasible.
- BH correction is applied once across all 207 region-by-cell-type module tests.
- Transcriptomic data are an independent convergence layer and do not tune model priors.
- Sobol total-order estimates are screening values; PRCC is the primary sensitivity analysis.

## Figures and source data

- Figure 1 and the graphical abstract use Matplotlib vector patches.
- Figure 2 uses Seaborn and Matplotlib with `phase_space.csv` and `model_summary.json`.
- Figure 3 uses R and ggplot2 with `prcc.csv`, `sobol_indices.csv`, and `structural_contrasts.csv`.
- Figure 4 uses Seaborn and Matplotlib with `module_effects.csv`, `cross_layer_processes.csv`, and `cross_layer_summary.json`.

## Release links

- GitHub: <https://github.com/vxxqv/ssri-neurogenesis-multiscale>
- Archived release: pending Zenodo DOI

The archived release DOI will be added after the Zenodo record is published.

## Citation

Use `CITATION.cff` when citing this software. Cite the source Rayan et al. dataset separately when using its expression matrices.

## License

The original code, derived tables, and original figures in this repository are released under the MIT License. This license does not apply to the separately downloaded upstream GSE197622 source files.
