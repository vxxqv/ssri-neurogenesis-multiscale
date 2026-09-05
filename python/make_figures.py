from __future__ import annotations

from pathlib import Path
import argparse
import json
import subprocess
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, FancyArrowPatch


INK = "#171717"
GRAY = "#696969"
ORANGE = "#D9572B"


def style() -> None:
    sns.set_theme(style="whitegrid")
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "axes.titlecolor": INK,
        "text.color": INK,
        "xtick.color": GRAY,
        "ytick.color": GRAY,
        "grid.color": "#DDDDDD",
        "grid.linewidth": 0.55,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def panel(ax: plt.Axes, label: str, title: str) -> None:
    ax.text(0.0, 1.04, f"{label}  {title}", transform=ax.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold")


def save(fig: plt.Figure, outdir: Path, name: str) -> None:
    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(outdir / f"{name}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def figure1(outdir: Path) -> None:
    fig = plt.figure(figsize=(12, 7.2), layout="constrained")
    grid = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0])
    ax = fig.add_subplot(grid[0])
    ax.set_xlim(-0.6, 5.7)
    ax.set_ylim(-1.15, 1.05)
    ax.axis("off")
    panel(ax, "A", "Lineage model and treatment channels")
    nodes = ["Q", "A", "P", "N", "M", "G"]
    labels = ["Quiescent\nstem cell", "Activated\nprecursor", "Proliferating\nprogenitor", "Neuroblast", "Maturing\nneuron", "Integrated\ngranule cell"]
    channels = ["Activation", "Proliferation", "Maturation", "Integration", "Efficacy", "Survival"]
    for i in range(5):
        ax.add_patch(FancyArrowPatch((i + 0.24, 0), (i + 0.76, 0), arrowstyle="-|>", mutation_scale=10, linewidth=1.4, color=INK))
    for i, (node, label) in enumerate(zip(nodes, labels)):
        face = ORANGE if node == "G" else "white"
        color = "white" if node == "G" else INK
        ax.add_patch(Circle((i, 0), 0.23, facecolor=face, edgecolor=INK if node != "G" else ORANGE, linewidth=1.5))
        ax.text(i, 0, node, ha="center", va="center", fontsize=14, fontweight="bold", color=color)
        ax.text(i, -0.37, label, ha="center", va="top", fontsize=8.5, color=GRAY)
    channel_x = [0.15, 1.15, 2.15, 3.15, 5.0, 4.35]
    for x, label in zip(channel_x, channels):
        ax.text(x, 0.58, label, ha="center", va="center", fontsize=8.5, color=ORANGE, fontweight="bold")
        ax.plot([x, x], [0.47, 0.28], color=ORANGE, linewidth=1.1)
    ax.text(0, -0.98, "Numerical extent", fontsize=10, fontweight="bold")
    ax.text(1.13, -0.98, "M + G", fontsize=15, fontweight="bold")
    ax.text(3.1, -0.98, "Functional index", fontsize=10, fontweight="bold", color=ORANGE)
    ax.text(4.36, -0.98, "Σ integrated cells × efficacy", fontsize=13, fontweight="bold", color=ORANGE)

    ax = fig.add_subplot(grid[1])
    panel(ax, "B", "Independent evidence structure")
    rows = ["GSE197622", "GSE43261", "GSE309750", "GSE222756", "GSE205325", "GSE292948"]
    cols = ["Cell type", "Treatment\nresponse", "Hippocampal\ncircuit", "Brain\nregions", "Stress\ncontext", "Three\nSSRIs"]
    values = np.eye(6)
    cmap = LinearSegmentedColormap.from_list("evidence", ["white", ORANGE])
    sns.heatmap(values, ax=ax, cmap=cmap, vmin=0, vmax=1, cbar=False, linewidths=1.2, linecolor="white", xticklabels=cols, yticklabels=rows)
    ax.tick_params(axis="x", rotation=0, labelsize=8.5)
    ax.tick_params(axis="y", rotation=0, labelsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for i, text in enumerate(["Pseudobulk cell map", "Responder versus resistant", "Granule and mossy cells", "10 regions, 3 antidepressants", "Stress plus fluoxetine", "Fluoxetine, sertraline, citalopram"]):
        ax.text(6.18, i + 0.5, text, va="center", fontsize=8.5, color=GRAY, clip_on=False)
    save(fig, outdir, "Fig1_study_architecture")


def figure2(results: Path, outdir: Path) -> None:
    phase = pd.read_csv(results / "model" / "phase_space.csv")
    temporal = pd.read_csv(results / "model" / "temporal_summary.csv")
    prediction = pd.read_csv(results / "model" / "predictive_information.csv")
    validation = pd.read_csv(results / "model" / "deterministic_validation.csv")
    summary = json.loads((results / "model" / "model_summary.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    ax = axes[0, 0]
    mean_mismatch = (phase.delta_extent > 0) & (phase.delta_fni <= 0)
    ax.scatter(phase.loc[~mean_mismatch, "delta_extent"], phase.loc[~mean_mismatch, "delta_fni"], s=12, color="#BDBDBD", alpha=0.55, linewidth=0)
    ax.scatter(phase.loc[mean_mismatch, "delta_extent"], phase.loc[mean_mismatch, "delta_fni"], s=25, color=ORANGE, alpha=0.9, linewidth=0)
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.axvline(0, color=INK, linewidth=0.8)
    x = phase.delta_extent.to_numpy()
    y = phase.delta_fni.to_numpy()
    slope, intercept = np.polyfit(x, y, 1)
    line = np.linspace(np.quantile(x, 0.005), np.quantile(x, 0.995), 100)
    ax.plot(line, slope * line + intercept, color=INK, linewidth=1.2)
    ax.set_xlim(np.quantile(x, 0.005), np.quantile(x, 0.995))
    ax.set_ylim(np.quantile(y, 0.005), np.quantile(y, 0.995))
    ax.set_xlabel("Change in numerical extent")
    ax.set_ylabel("Change in functional index")
    ax.text(0.04, 0.92, f"r = {summary['pearson_extent_fni']:.3f}", transform=ax.transAxes, fontweight="bold")
    ax.text(0.04, 0.84, f"{mean_mismatch.sum()} mean mismatches", transform=ax.transAxes, color=ORANGE, fontweight="bold")
    panel(ax, "A", "Paired stochastic outcomes")

    ax = axes[0, 1]
    ax.plot(temporal.t_end, temporal.mean_mismatch_probability * 100, color=ORANGE, marker="o", linewidth=2.2, markersize=6)
    peak = temporal.loc[temporal.mean_mismatch_probability.idxmax()]
    ax.scatter([peak.t_end], [peak.mean_mismatch_probability * 100], s=90, facecolor="white", edgecolor=ORANGE, linewidth=2, zorder=4)
    ax.text(peak.t_end + 2, peak.mean_mismatch_probability * 100, f"Peak {peak.mean_mismatch_probability * 100:.1f}%", va="center", color=ORANGE, fontweight="bold")
    ax.set_xlabel("Treatment day")
    ax.set_ylabel("Replicate mismatch probability (%)")
    ax.set_xticks(temporal.t_end)
    ax.set_ylim(0, max(16, temporal.mean_mismatch_probability.max() * 120))
    panel(ax, "B", "Mismatch is most likely early")

    ax = axes[1, 0]
    labels = {"extent_only": "Extent only", "cell_composition": "Maturing + integrated", "process_aware": "Process aware"}
    prediction = prediction.assign(label=prediction.model.map(labels)).sort_values("r2")
    ypos = np.arange(len(prediction))
    ax.hlines(ypos, prediction.r2_low, prediction.r2_high, color=GRAY, linewidth=2)
    colors = [ORANGE if model == "process_aware" else GRAY for model in prediction.model]
    ax.scatter(prediction.r2, ypos, s=75, color=colors, zorder=3)
    for y0, value, model in zip(ypos, prediction.r2, prediction.model):
        ax.text(value + 0.025, y0, f"{value:.3f}", va="center", fontweight="bold", color=ORANGE if model == "process_aware" else INK)
    ax.set_yticks(ypos, prediction.label)
    ax.set_xlim(0.35, 1.0)
    ax.set_xlabel("Held-out R² for functional change")
    panel(ax, "C", "Information beyond cell count")

    ax = axes[1, 1]
    sample = validation.sample(900, random_state=20260825)
    ax.scatter(sample.delta_fni_ode, sample.delta_fni_ssa, s=12, color="#AFAFAF", alpha=0.55, linewidth=0)
    limits = np.quantile(np.r_[sample.delta_fni_ode, sample.delta_fni_ssa], [0.005, 0.995])
    ax.plot(limits, limits, color=ORANGE, linewidth=1.6)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel("Deterministic functional change")
    ax.set_ylabel("Mean stochastic functional change")
    ax.text(0.04, 0.92, f"r = {summary['ode_ssa_fni_r']:.3f}", transform=ax.transAxes, fontweight="bold")
    ax.text(0.04, 0.84, f"MAE = {summary['ode_ssa_fni_mae']:.2f}", transform=ax.transAxes, color=GRAY)
    panel(ax, "D", "Independent numerical agreement")
    save(fig, outdir, "Fig2_model_results")


def figure4(results: Path, outdir: Path) -> None:
    programs = pd.read_csv(results / "transcriptomics" / "program_effects.csv")
    regional = pd.read_csv(results / "transcriptomics" / "regional_drug_concordance.csv")
    cross = pd.read_csv(results / "transcriptomics" / "cross_ssri_context.csv")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), layout="constrained")
    key = [
        ("GSE309750", "mossy_cells_Fluoxetine_vs_Vehicle", "Mossy cells, fluoxetine"),
        ("GSE43261", "dorsal_Responder_vs_Resistant", "Dorsal DG, responder versus resistant"),
        ("GSE205325", "stress_Fluoxetine_vs_stress", "Stressed hippocampus, fluoxetine"),
        ("GSE292948", "vivo_Fluoxetine_vs_Control", "In vivo cortex, fluoxetine"),
        ("GSE222756", "dorDG_FT_vs_Control", "Dorsal DG atlas, fluoxetine"),
        ("GSE222756", "venDG_FT_vs_Control", "Ventral DG atlas, fluoxetine"),
    ]
    ax = axes[0, 0]
    ordinary, exclusive, labels = [], [], []
    for dataset, contrast, label in key:
        block = programs[(programs.dataset == dataset) & (programs.contrast == contrast)].set_index("program")
        ordinary.append(block.loc["integration_bias", "hedges_g"])
        exclusive.append(block.loc["exclusive_integration_bias", "hedges_g"])
        labels.append(label)
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color=INK, linewidth=0.8)
    for yi, a, b in zip(y, ordinary, exclusive):
        ax.plot([a, b], [yi, yi], color="#BEBEBE", linewidth=1.4)
    ax.scatter(ordinary, y, color=ORANGE, s=55, label="All ontology genes", zorder=3)
    ax.scatter(exclusive, y, facecolor="white", edgecolor=INK, linewidth=1.3, s=55, label="Nonoverlapping genes", zorder=3)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Integration bias, Hedges’ g")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    panel(ax, "A", "Integration beyond numerical programs")

    ax = axes[0, 1]
    heat = regional[regional.program == "integration_bias"].set_index("context")[["FT", "Bupropion", "Desipramine"]]
    order = ["dorDG", "venDG", "BLA", "CGC", "ILC", "LConly", "mPOA", "NACShell", "PLC", "Raphe"]
    heat = heat.loc[order]
    cmap = LinearSegmentedColormap.from_list("signed", [INK, "white", ORANGE])
    sns.heatmap(heat, ax=ax, cmap=cmap, center=0, annot=True, fmt=".1f", linewidths=0.8, linecolor="white", cbar_kws={"label": "Integration bias (g)", "shrink": 0.75})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticklabels(["Fluoxetine", "Bupropion", "Desipramine"], rotation=22, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    panel(ax, "B", "Region and antidepressant specificity")

    ax = axes[1, 0]
    subset = cross[cross.program == "integration_bias"].copy()
    for _, row in subset.iterrows():
        ax.plot([0, 1], [row.vitro, row.vivo], color="#BDBDBD", linewidth=1.7, zorder=1)
        ax.scatter(0, row.vitro, facecolor="white", edgecolor=INK, s=65, linewidth=1.3, zorder=2)
        ax.scatter(1, row.vivo, color=ORANGE, s=65, zorder=2)
        ax.text(1.05, row.vivo, row.treated, va="center", fontsize=8.5)
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.set_xticks([0, 1], ["Cultured neurons", "In vivo cortex"])
    ax.set_xlim(-0.2, 1.45)
    ax.set_ylabel("Integration bias, Hedges’ g")
    panel(ax, "C", "The cellular context changes SSRI effects")

    ax = axes[1, 1]
    wide = programs[programs.program.isin(["integration_bias", "exclusive_integration_bias"])].pivot_table(index=["dataset", "contrast"], columns="program", values="hedges_g").dropna()
    ax.scatter(wide.integration_bias, wide.exclusive_integration_bias, s=18, color="#AFAFAF", alpha=0.65, linewidth=0)
    highlight = wide.loc[[index for index in wide.index if "mossy_cells" in index[1] or "dorsal_Responder_vs_Resistant" in index[1]]]
    ax.scatter(highlight.integration_bias, highlight.exclusive_integration_bias, s=55, color=ORANGE, zorder=3)
    limits = np.quantile(np.r_[wide.integration_bias, wide.exclusive_integration_bias], [0.01, 0.99])
    padding = (limits[1] - limits[0]) * 0.06
    display_limits = [limits[0] - padding, limits[1] + padding]
    ax.plot(display_limits, display_limits, color=INK, linewidth=1)
    ax.set_xlim(display_limits)
    ax.set_ylim(display_limits)
    rho = wide.integration_bias.rank().corr(wide.exclusive_integration_bias.rank())
    ax.text(0.04, 0.92, f"Spearman ρ = {rho:.3f}", transform=ax.transAxes, fontweight="bold")
    ax.set_xlabel("All-gene integration bias")
    ax.set_ylabel("Nonoverlapping-gene integration bias")
    panel(ax, "D", "Ontology overlap does not drive the pattern")
    save(fig, outdir, "Fig4_transcriptomic_evidence")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--r-library", required=True)
    args = parser.parse_args()
    results = Path(args.results)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    style()
    figure1(outdir)
    figure2(results, outdir)
    figure4(results, outdir)
    script = Path(__file__).with_name("make_mechanism_figure.R")
    subprocess.run([args.rscript, str(script), str(results), str(outdir), str(Path(args.r_library).resolve())], check=True)


if __name__ == "__main__":
    main()
