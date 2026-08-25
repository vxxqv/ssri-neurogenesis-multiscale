param(
    [string]$Python = "python",
    [string]$Cxx = "g++",
    [string]$Git = "git",
    [string]$Rscript = "Rscript",
    [string]$RayanRepoPath = "external\rayan2022",
    [switch]$SkipFetch
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = Split-Path -Parent $PSScriptRoot
$BuildDir = Join-Path $ProjectRoot "build"
$ResultsDir = Join-Path $ProjectRoot "results"
$FigureDir = Join-Path $ProjectRoot "figures"
$SourceUrl = "https://github.com/arulrayan/Integrative-multi-omics-landscape-of-fluoxetine-action-across-27-brain-regions.git"
$SourceCommit = "c3ae3f6d2a98a81d6b0208f19b70fcc56dc50c5f"

if (-not [System.IO.Path]::IsPathRooted($RayanRepoPath)) {
    $RayanRepoPath = Join-Path $ProjectRoot $RayanRepoPath
}

New-Item -ItemType Directory -Force -Path $BuildDir | Out-Null
New-Item -ItemType Directory -Force -Path $ResultsDir | Out-Null
New-Item -ItemType Directory -Force -Path $FigureDir | Out-Null

if (-not $SkipFetch) {
    if (-not (Test-Path -LiteralPath $RayanRepoPath)) {
        & $Git clone $SourceUrl $RayanRepoPath
    }
    & $Git -C $RayanRepoPath checkout $SourceCommit
}

$Dorsal = Join-Path $RayanRepoPath "single-cell DEGs\Single-cell collated tp10k_dorDG.txt"
$Ventral = Join-Path $RayanRepoPath "single-cell DEGs\Single-cell collated tp10k_venDG.txt"
if (-not (Test-Path -LiteralPath $Dorsal) -or -not (Test-Path -LiteralPath $Ventral)) {
    throw "Rayan et al. TP10K files were not found under $RayanRepoPath"
}

$Design = Join-Path $ResultsDir "design.csv"
$Simulation = Join-Path $ResultsDir "simulation_results.csv"
$Simulator = Join-Path $BuildDir "ssa_batch.exe"

& $Cxx -O3 -std=c++17 (Join-Path $ProjectRoot "cpp\ssa_batch.cpp") -o $Simulator
& $Python (Join-Path $ProjectRoot "python\prepare_design.py") --priors (Join-Path $ProjectRoot "config\parameter_priors.csv") --out $Design --seed 20260825 --phase-n 1500 --structural-n 300 --sobol-n 128
& $Simulator $Design $Simulation
& $Python (Join-Path $ProjectRoot "python\analyze_model.py") --design $Design --simulation $Simulation --priors (Join-Path $ProjectRoot "config\parameter_priors.csv") --outdir (Join-Path $ResultsDir "model")
& $Python (Join-Path $ProjectRoot "python\analyze_transcriptomics.py") --dorsal $Dorsal --ventral $Ventral --modules (Join-Path $ProjectRoot "config\gene_modules.csv") --outdir (Join-Path $ResultsDir "omics") --seed 20260825
& $Python (Join-Path $ProjectRoot "python\cross_layer.py") --sobol (Join-Path $ResultsDir "model\sobol_indices.csv") --omics (Join-Path $ResultsDir "omics\process_omics_evidence.csv") --outdir (Join-Path $ResultsDir "integration") --seed 20260825
& $Python (Join-Path $ProjectRoot "python\make_figures.py") --results $ResultsDir --outdir $FigureDir --rscript $Rscript
& $Python (Join-Path $ProjectRoot "tests\test_invariants.py")

Write-Host "Analysis reproduced successfully in $ProjectRoot"
