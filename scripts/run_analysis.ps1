param(
    [string]$Python = "python",
    [string]$Cxx = "g++",
    [string]$Rscript = "Rscript",
    [string]$DataPath = "",
    [int]$Workers = 8,
    [switch]$SkipFetch
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Confirm-Success {
    if ($LASTEXITCODE -ne 0) { throw "A program failed with exit code $LASTEXITCODE" }
}

$ProjectRoot = Split-Path -Parent $PSScriptRoot
$PublicData = if ($DataPath) { $DataPath } else { Join-Path $ProjectRoot "external\public_data" }
$WorkDir = Join-Path $ProjectRoot "work"
$BuildDir = Join-Path $ProjectRoot "build"
$ResultsDir = Join-Path $ProjectRoot "results"
$ModelWork = Join-Path $WorkDir "model"
$TranscriptWork = Join-Path $WorkDir "transcriptomics"
$ModelResults = Join-Path $ResultsDir "model"
$TranscriptResults = Join-Path $ResultsDir "transcriptomics"
$FigureDir = Join-Path $ProjectRoot "figures"
$RLibrary = Join-Path $ProjectRoot ".r-lib"
$RCache = Join-Path $WorkDir "r-packages"

New-Item -ItemType Directory -Force -Path $PublicData,$WorkDir,$BuildDir,$ResultsDir,$ModelWork,$TranscriptWork,$ModelResults,$TranscriptResults,$FigureDir | Out-Null

if (-not $SkipFetch) {
    & $Python (Join-Path $ProjectRoot "python\fetch_public_data.py") --outdir $PublicData
    Confirm-Success
}

& $Python (Join-Path $ProjectRoot "python\build_gene_programs.py") --programs (Join-Path $ProjectRoot "config\go_programs.csv") --obo (Join-Path $PublicData "go-basic.obo") --mouse-gaf (Join-Path $PublicData "mgi.gaf.gz") --rat-gaf (Join-Path $PublicData "rgd.gaf.gz") --out (Join-Path $PublicData "gene_programs.tsv")
Confirm-Success
& $Python (Join-Path $ProjectRoot "python\install_r_packages.py") --rscript $Rscript --library $RLibrary --cache $RCache
Confirm-Success

$Design = Join-Path $WorkDir "design.csv"
$Aggregate = Join-Path $WorkDir "aggregate.csv"
$Replicates = Join-Path $WorkDir "replicates.csv"
$Deterministic = Join-Path $WorkDir "deterministic.csv"
$Simulator = Join-Path $BuildDir "ssa_batch.exe"

& $Cxx -O3 -std=c++17 (Join-Path $ProjectRoot "cpp\ssa_batch.cpp") -o $Simulator
Confirm-Success
& $Python (Join-Path $ProjectRoot "python\prepare_design.py") --priors (Join-Path $ProjectRoot "config\parameter_priors.csv") --out $Design --seed 20260825 --phase-n 2048 --phase-reps 16 --factorial-n 512 --factorial-reps 0 --sobol-n 2048 --sobol-reps 0 --temporal-n 512 --temporal-reps 4
Confirm-Success
& $Python (Join-Path $ProjectRoot "python\solve_deterministic.py") --design $Design --out $Deterministic
Confirm-Success
& $Python (Join-Path $ProjectRoot "python\run_simulation_parallel.py") --simulator $Simulator --design $Design --aggregate $Aggregate --replicates $Replicates --workdir (Join-Path $WorkDir "workers") --workers $Workers
Confirm-Success
& $Python (Join-Path $ProjectRoot "python\analyze_model.py") --design $Design --simulation $Aggregate --replicates $Replicates --deterministic $Deterministic --priors (Join-Path $ProjectRoot "config\parameter_priors.csv") --outdir $ModelWork --seed 20260825
Confirm-Success
& $Python (Join-Path $ProjectRoot "python\validate_model.py") --design $Design --shapley (Join-Path $ModelWork "channel_shapley.csv") --deterministic $Deterministic --out (Join-Path $ModelWork "numerical_validation.json")
Confirm-Success

& $Python (Join-Path $ProjectRoot "python\analyze_public_transcriptomes.py") --data $PublicData --rayan-dorsal (Join-Path $PublicData "GSE197622_tp10k_dorDG.txt") --rayan-ventral (Join-Path $PublicData "GSE197622_tp10k_venDG.txt") --programs (Join-Path $PublicData "gene_programs.tsv") --mouse-gtf (Join-Path $PublicData "Mus_musculus.GRCm39.116.gtf.gz") --rat-gtf (Join-Path $PublicData "Rattus_norvegicus.GRCr8.116.gtf.gz") --gpl1261 (Join-Path $PublicData "GPL1261.annot.gz") --outdir $TranscriptWork
Confirm-Success
& $Python (Join-Path $ProjectRoot "python\audit_relative_programs.py") --data $PublicData --outdir $TranscriptWork
Confirm-Success

$ModelFiles = @("model_summary.json","numerical_validation.json","phase_space.csv","predictive_information.csv","predictive_predictions.csv","channel_shapley_summary.csv","channel_pair_interaction_summary.csv","sobol_indices.csv","sobol_convergence.csv","replicate_convergence.csv","temporal_summary.csv","temporal_sets.csv","deterministic_validation.csv")
$TranscriptFiles = @("dataset_qc.csv","program_effects.csv","gene_effects_top.csv","meta_analysis.csv","cross_ssri_context.csv","regional_drug_concordance.csv","transcriptome_summary.json","relative_score_robustness.csv","relative_score_gene_splits.csv")
foreach ($Name in $ModelFiles) { Copy-Item -LiteralPath (Join-Path $ModelWork $Name) -Destination (Join-Path $ModelResults $Name) -Force }
foreach ($Name in $TranscriptFiles) { Copy-Item -LiteralPath (Join-Path $TranscriptWork $Name) -Destination (Join-Path $TranscriptResults $Name) -Force }

& $Python (Join-Path $ProjectRoot "python\make_figures.py") --results $ResultsDir --outdir $FigureDir --rscript $Rscript --r-library $RLibrary
Confirm-Success
& $Python (Join-Path $ProjectRoot "tests\test_invariants.py")
Confirm-Success

Write-Host "Analysis reproduced successfully."
