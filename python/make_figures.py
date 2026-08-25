from __future__ import annotations

from pathlib import Path
import argparse
import json
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle


INK = "#161616"
TEXT = "#404040"
GRID = "#DEDAD3"
PAPER = "#FFFFFF"
TEAL = "#008B74"
ORANGE = "#D75B28"
PURPLE = "#7655B5"
GOLD = "#D99A18"
ROSE = "#B5486D"
PALE_TEAL = "#E6F3F0"
PALE_ORANGE = "#F9ECE5"
PALE_PURPLE = "#F0ECF7"


def configure_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "figure.facecolor": PAPER,
            "axes.facecolor": PAPER,
            "savefig.facecolor": PAPER,
            "font.family": "DejaVu Sans",
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK,
            "axes.titlecolor": INK,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "axes.linewidth": 0.9,
            "figure.dpi": 120,
        }
    )


def require_columns(frame: pd.DataFrame, name: str, columns: set[str]) -> None:
    missing = sorted(columns.difference(frame.columns))
    if missing:
        raise ValueError(f"{name} is missing required columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError(f"{name} is empty")


def save_figure(figure: plt.Figure, outdir: Path, stem: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    figure.savefig(outdir / f"{stem}.png", dpi=600, facecolor=PAPER)
    figure.savefig(outdir / f"{stem}.tif", dpi=600, facecolor=PAPER, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(figure)


def rounded_box(axis: plt.Axes, x: float, y: float, width: float, height: float, edge: str, fill: str = PAPER, radius: float = 0.02, linewidth: float = 1.6) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        transform=axis.transAxes,
        facecolor=fill,
        edgecolor=edge,
        linewidth=linewidth,
    )
    axis.add_patch(patch)
    return patch


def arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = TEXT, width: float = 1.5, connection: str = "arc3") -> None:
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=axis.transAxes,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=width,
            color=color,
            connectionstyle=connection,
        )
    )


def add_panel_label(axis: plt.Axes, label: str, x: float = -0.08, y: float = 1.05) -> None:
    axis.text(x, y, label, transform=axis.transAxes, fontsize=13, fontweight="bold", color=INK, va="top")


def figure1(outdir: Path) -> None:
    figure = plt.figure(figsize=(10, 6), constrained_layout=False)
    axis = figure.add_axes([0.055, 0.06, 0.89, 0.89])
    axis.set_axis_off()
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_autoscale_on(False)
    axis.text(0.00, 0.98, "A", fontsize=13, fontweight="bold", transform=axis.transAxes, va="top")
    axis.text(0.045, 0.98, "Lineage model", fontsize=12, fontweight="bold", transform=axis.transAxes, va="top")
    states = [
        ("Q", "Quiescent\nstem cell"),
        ("A", "Activated\nprecursor"),
        ("P", "Proliferating\nprogenitor"),
        ("N", "Neuroblast"),
        ("M", "Maturing\nneuron"),
        ("G", "Integrated\ngranule cell"),
    ]
    x_positions = np.linspace(0.08, 0.92, len(states))
    axis.plot([0.08, 0.92], [0.79, 0.79], transform=axis.transAxes, color=INK, linewidth=1.4, solid_capstyle="round")
    for index, (symbol, label) in enumerate(states):
        x = float(x_positions[index])
        color = ORANGE if index == len(states) - 1 else INK
        fill = PALE_ORANGE if index == len(states) - 1 else PAPER
        axis.scatter([x], [0.79], transform=axis.transAxes, s=760, facecolor=fill, edgecolor=color, linewidth=1.8, zorder=4)
        axis.text(x, 0.79, symbol, transform=axis.transAxes, ha="center", va="center", fontsize=13, fontweight="bold", color=INK, zorder=5)
        axis.text(x, 0.685, label, transform=axis.transAxes, ha="center", va="center", fontsize=8.8, color=TEXT, linespacing=1.0)
        if index > 0:
            axis.plot([x, x], [0.755, 0.61], transform=axis.transAxes, color="#A6A6A6", linestyle=(0, (2, 3)), linewidth=0.9, zorder=1)
            axis.scatter([x], [0.605], transform=axis.transAxes, s=15, color="#888888", marker="v", zorder=3)
    axis.text(0.50, 0.575, "Stage specific loss", transform=axis.transAxes, fontsize=8.8, color=TEXT, ha="center")
    axis.plot([0.00, 1.00], [0.535, 0.535], transform=axis.transAxes, color=GRID, linewidth=1.0)
    axis.text(0.00, 0.49, "B", fontsize=13, fontweight="bold", transform=axis.transAxes, va="top")
    axis.text(0.045, 0.49, "Study logic and molecular audit", fontsize=12, fontweight="bold", transform=axis.transAxes, va="top")
    axis.text(0.08, 0.385, "PAIRED DESIGN", transform=axis.transAxes, fontsize=8.1, fontweight="bold", color=TEXT, ha="center")
    axis.text(0.08, 0.335, "Control and SSRI", transform=axis.transAxes, fontsize=10.5, fontweight="bold", color=INK, ha="center")
    axis.text(0.08, 0.295, "1,500 parameter sets", transform=axis.transAxes, fontsize=8.6, color=TEXT, ha="center")
    arrow(axis, (0.18, 0.335), (0.30, 0.335), INK, 1.2)
    axis.add_patch(Rectangle((0.31, 0.235), 0.34, 0.20, transform=axis.transAxes, facecolor=PAPER, edgecolor=INK, linewidth=1.3))
    axis.plot([0.31, 0.65], [0.335, 0.335], transform=axis.transAxes, color=GRID, linewidth=1.0)
    axis.plot([0.31, 0.65], [0.235, 0.235], transform=axis.transAxes, color=ORANGE, linewidth=3.2, solid_capstyle="butt")
    axis.text(0.48, 0.395, "Numerical extent", transform=axis.transAxes, fontsize=10.2, fontweight="bold", color=INK, ha="center")
    axis.text(0.48, 0.355, "M + G", transform=axis.transAxes, fontsize=9.2, color=INK, ha="center")
    axis.text(0.48, 0.295, "Functional index", transform=axis.transAxes, fontsize=10.2, fontweight="bold", color=INK, ha="center")
    axis.text(0.48, 0.255, "integration and efficacy", transform=axis.transAxes, fontsize=9.2, color=ORANGE, ha="center")
    arrow(axis, (0.66, 0.335), (0.72, 0.335), INK, 1.2)
    axis.plot([0.735, 0.735], [0.245, 0.425], transform=axis.transAxes, color=ORANGE, linewidth=3.0, solid_capstyle="round")
    axis.text(0.765, 0.395, "JOINT INTERPRETATION", transform=axis.transAxes, fontsize=8.1, fontweight="bold", color=TEXT, ha="left")
    axis.text(0.765, 0.345, "Concordance or mismatch", transform=axis.transAxes, fontsize=10.5, fontweight="bold", color=INK, ha="left")
    axis.text(0.765, 0.285, "51 of 1,500 mismatched", transform=axis.transAxes, fontsize=9.5, color=ORANGE, fontweight="bold", ha="left")
    axis.plot([0.08, 0.92], [0.145, 0.145], transform=axis.transAxes, color=GRID, linewidth=1.0)
    axis.text(0.08, 0.095, "INDEPENDENT MOLECULAR AUDIT", transform=axis.transAxes, fontsize=8.1, fontweight="bold", color=TEXT, ha="left")
    axis.text(0.325, 0.095, "GSE197622", transform=axis.transAxes, fontsize=10.2, fontweight="bold", color=INK, ha="left")
    arrow(axis, (0.44, 0.095), (0.55, 0.095), TEXT, 1.0)
    axis.text(0.58, 0.095, "No module survived global FDR correction", transform=axis.transAxes, fontsize=9.5, color=TEXT, ha="left", va="center")
    save_figure(figure, outdir, "Fig1_architecture")


def figure2(outdir: Path, phase: pd.DataFrame, summary: dict) -> None:
    require_columns(phase, "phase_space.csv", {"delta_extent", "delta_fni", "decoupled"})
    phase = phase.copy()
    phase["decoupled"] = phase["decoupled"].astype(bool)
    x_low, x_high = phase["delta_extent"].quantile([0.005, 0.995])
    y_low, y_high = phase["delta_fni"].quantile([0.005, 0.995])
    clipped = phase.assign(
        display_extent=phase["delta_extent"].clip(x_low, x_high),
        display_fni=phase["delta_fni"].clip(y_low, y_high),
    )
    figure = plt.figure(figsize=(12, 7.5))
    marginal_x = figure.add_axes([0.085, 0.79, 0.625, 0.14])
    main = figure.add_axes([0.085, 0.19, 0.625, 0.58], sharex=marginal_x)
    marginal_y = figure.add_axes([0.72, 0.19, 0.055, 0.58], sharey=main)
    summary_axis = figure.add_axes([0.81, 0.16, 0.16, 0.77])
    other = clipped.loc[~clipped["decoupled"]]
    mismatch = clipped.loc[clipped["decoupled"]]
    sns.scatterplot(data=other, x="display_extent", y="display_fni", ax=main, s=18, color=PURPLE, alpha=0.25, edgecolor=None, rasterized=True)
    sns.regplot(data=clipped, x="display_extent", y="display_fni", ax=main, scatter=False, ci=95, color=TEAL, line_kws={"linewidth": 2.0}, truncate=False)
    sns.scatterplot(data=mismatch, x="display_extent", y="display_fni", ax=main, s=42, marker="X", color=ORANGE, edgecolor=PAPER, linewidth=0.5, zorder=5)
    sns.histplot(data=clipped, x="display_extent", bins=30, ax=marginal_x, color=PURPLE, alpha=0.78, edgecolor=PAPER, linewidth=0.4)
    sns.histplot(data=clipped, y="display_fni", bins=30, ax=marginal_y, color=TEAL, alpha=0.78, edgecolor=PAPER, linewidth=0.4)
    main.axvline(0, color=TEXT, linewidth=1.0, linestyle=(0, (4, 3)))
    main.axhline(0, color=TEXT, linewidth=1.0, linestyle=(0, (4, 3)))
    main.set_xlabel("Change in numerical extent (model cells)", fontsize=11, labelpad=8)
    main.set_ylabel("Change in functional neurogenesis index", fontsize=11, labelpad=8)
    main.grid(True, color=GRID, linewidth=0.7)
    main.spines[["top", "right"]].set_visible(False)
    marginal_x.set_axis_off()
    marginal_y.set_axis_off()
    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PURPLE, alpha=0.45, markersize=7, label="Other parameter sets"),
        Line2D([0], [0], marker="X", color="none", markerfacecolor=ORANGE, markeredgecolor=ORANGE, markersize=7, label="Mismatch"),
        Line2D([0], [0], color=TEAL, linewidth=2, label="Linear trend with 95% CI"),
    ]
    figure.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.08, 0.055), frameon=False, fontsize=8.5, ncol=3, borderaxespad=0, columnspacing=1.6, handlelength=2.8)
    summary_axis.set_axis_off()
    summary_axis.set_xlim(0, 1)
    summary_axis.set_ylim(0, 1)
    metrics = [
        ("Pearson r", f"{summary['pearson_delta_extent_fni']:.3f}", TEAL),
        ("Spearman rho", f"{summary['spearman_delta_extent_fni']:.3f}", PURPLE),
        ("Mismatch", f"{summary['decoupled_count']} / {summary['n_phase_parameter_sets']}", ORANGE),
        ("Mismatch rate", f"{100 * summary['decoupled_fraction']:.1f}%", ROSE),
    ]
    for index, (label, value, color) in enumerate(metrics):
        y = 0.87 - index * 0.18
        summary_axis.text(0.03, y + 0.035, label, transform=summary_axis.transAxes, fontsize=8.8, color=TEXT, va="center")
        summary_axis.text(0.03, y - 0.025, value, transform=summary_axis.transAxes, fontsize=15, fontweight="bold", color=color, va="center")
        if index < len(metrics) - 1:
            summary_axis.plot([0.03, 0.97], [y - 0.09, y - 0.09], transform=summary_axis.transAxes, color=GRID, linewidth=0.9)
    summary_axis.text(0.02, 0.02, "Axes show the 0.5th to 99.5th\npercentiles. All values remain in\nthe source table.", transform=summary_axis.transAxes, fontsize=7.6, color=TEXT, va="bottom", linespacing=1.35)
    figure.text(0.06, 0.955, "A   Paired simulation phase space", fontsize=12.5, fontweight="bold", color=INK)
    figure.text(0.79, 0.955, "B   Model summary", fontsize=12.5, fontweight="bold", color=INK)
    save_figure(figure, outdir, "Fig2_phase_space")


def figure4(outdir: Path, effects: pd.DataFrame, cross: pd.DataFrame, cross_summary: dict) -> None:
    require_columns(effects, "module_effects.csv", {"region", "cell_type", "process", "hedges_g"})
    require_columns(cross, "cross_layer_processes.csv", {"process", "model_ST", "omics_evidence"})
    effects = effects.copy()
    effects["region_short"] = effects["region"].map({"dorsal_DG": "D", "ventral_DG": "V"}).fillna(effects["region"])
    effects["profile"] = effects["region_short"] + " · " + effects["cell_type"].str.replace("_", " ", regex=False)
    process_order = ["activation", "activity", "differentiation", "energy", "integration", "maturation", "niche", "proliferation", "survival"]
    profile_order = effects["profile"].drop_duplicates().tolist()
    matrix = effects.pivot_table(index="process", columns="profile", values="hedges_g", aggfunc="first").reindex(index=process_order, columns=profile_order)
    figure = plt.figure(figsize=(12, 7.5))
    grid = figure.add_gridspec(1, 2, width_ratios=[3.4, 1.25], left=0.08, right=0.97, top=0.90, bottom=0.27, wspace=0.22)
    heat_axis = figure.add_subplot(grid[0, 0])
    rank_axis = figure.add_subplot(grid[0, 1])
    diverging = LinearSegmentedColormap.from_list("effect", [PURPLE, "#C8BCE0", PAPER, "#F2C6B3", ORANGE], N=256)
    limit = max(1.0, float(np.nanmax(np.abs(matrix.to_numpy()))))
    sns.heatmap(matrix, ax=heat_axis, cmap=diverging, center=0, vmin=-limit, vmax=limit, linewidths=0.5, linecolor=PAPER, cbar_kws={"label": "Hedges' g", "shrink": 0.72, "pad": 0.02})
    heat_axis.set_xlabel("")
    heat_axis.set_ylabel("")
    heat_axis.set_xticklabels(heat_axis.get_xticklabels(), rotation=63, ha="right", rotation_mode="anchor", fontsize=7.2)
    heat_axis.set_yticklabels([label.get_text().capitalize() for label in heat_axis.get_yticklabels()], rotation=0, fontsize=8.5)
    heat_axis.tick_params(length=0)
    heat_axis.set_title("A   Transcriptomic effects", loc="left", fontsize=11.5, fontweight="bold", pad=14, color=INK)
    cross = cross.copy()
    cross["model_rank"] = cross["model_ST"].rank(ascending=False, method="average")
    cross["omics_rank"] = cross["omics_evidence"].rank(ascending=False, method="average")
    cross = cross.sort_values("model_rank", ascending=False).reset_index(drop=True)
    y_positions = np.arange(len(cross))
    for y, row in zip(y_positions, cross.itertuples()):
        line_color = ORANGE if abs(row.model_rank - row.omics_rank) >= 3 else GRID
        rank_axis.plot([row.model_rank, row.omics_rank], [y, y], color=line_color, linewidth=2.2, solid_capstyle="round", zorder=1)
    rank_axis.scatter(cross["model_rank"], y_positions, s=58, color=TEAL, edgecolor=PAPER, linewidth=0.8, label="Model rank", zorder=3)
    rank_axis.scatter(cross["omics_rank"], y_positions, s=62, facecolor=PAPER, edgecolor=ORANGE, linewidth=1.8, label="Molecular rank", zorder=3)
    rank_axis.set_yticks(y_positions, [value.capitalize() for value in cross["process"]], fontsize=8.5)
    rank_axis.set_xticks(range(1, len(cross) + 1))
    rank_axis.set_xlim(0.5, len(cross) + 0.5)
    rank_axis.set_xlabel("Rank, 1 = strongest", fontsize=9.5)
    rank_axis.set_title("B   Process ranks", loc="left", fontsize=11.5, fontweight="bold", pad=14, color=INK)
    rank_axis.text(1.0, 1.015, f"rho = {cross_summary['spearman_rho']:.2f}   p = {cross_summary['permutation_p_two_sided']:.3f}", transform=rank_axis.transAxes, ha="right", va="bottom", fontsize=8.7, color=TEXT)
    rank_axis.grid(axis="x", color=GRID)
    rank_axis.grid(axis="y", visible=False)
    rank_axis.spines[["top", "right", "left"]].set_visible(False)
    rank_axis.tick_params(axis="y", length=0)
    rank_axis.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, fontsize=8.3)
    figure.text(0.08, 0.11, "No module test survived global BH FDR q < 0.10.", fontsize=9.0, fontweight="bold", color=INK)
    figure.text(0.08, 0.075, "Purple indicates lower and orange indicates higher module scores with fluoxetine.", fontsize=9.2, color=TEXT)
    save_figure(figure, outdir, "Fig4_omics_convergence")


def graphical_abstract(outdir: Path, summary: dict) -> None:
    figure = plt.figure(figsize=(12, 5.25))
    axis = figure.add_axes([0.03, 0.05, 0.94, 0.90])
    axis.set_axis_off()
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_autoscale_on(False)
    steps = [
        ("1", "Sample", "1,500 biological\nparameter sets", TEAL, PALE_TEAL),
        ("2", "Simulate", "Paired control and\nfluoxetine-like arms", PURPLE, PALE_PURPLE),
        ("3", "Separate", "Numerical extent\nfrom functional index", ORANGE, PALE_ORANGE),
        ("4", "Challenge", "Independent GSE197622\nprocess ranking", ROSE, "#F8EAF0"),
    ]
    x_positions = [0.0, 0.255, 0.51, 0.765]
    for index, (number, title, detail, color, fill) in enumerate(steps):
        x = x_positions[index]
        rounded_box(axis, x, 0.66, 0.205, 0.29, color, fill, 0.03, 1.6)
        axis.add_patch(Circle((x + 0.035, 0.88), 0.025, transform=axis.transAxes, facecolor=color, edgecolor=color))
        axis.text(x + 0.035, 0.88, number, transform=axis.transAxes, ha="center", va="center", color=PAPER, fontsize=9.5, fontweight="bold")
        axis.text(x + 0.075, 0.88, title, transform=axis.transAxes, fontsize=12, fontweight="bold", color=INK, va="center")
        axis.text(x + 0.025, 0.775, detail, transform=axis.transAxes, fontsize=10, color=TEXT, va="center", linespacing=1.35)
        if index < len(steps) - 1:
            arrow(axis, (x + 0.21, 0.805), (x_positions[index + 1] - 0.008, 0.805), color, 1.7)
    rounded_box(axis, 0.10, 0.08, 0.34, 0.22, TEAL, PAPER, 0.03, 1.8)
    rounded_box(axis, 0.56, 0.08, 0.34, 0.22, ORANGE, PAPER, 0.03, 1.8)
    axis.text(0.13, 0.24, "STRONG OVERALL RELATIONSHIP", transform=axis.transAxes, fontsize=9.2, fontweight="bold", color=TEAL)
    axis.text(0.13, 0.155, f"Pearson r = {summary['pearson_delta_extent_fni']:.3f}", transform=axis.transAxes, fontsize=18, fontweight="bold", color=INK)
    axis.text(0.59, 0.24, "PRESPECIFIED MISMATCH", transform=axis.transAxes, fontsize=9.2, fontweight="bold", color=ORANGE)
    axis.text(0.59, 0.155, f"{summary['decoupled_count']} / {summary['n_phase_parameter_sets']}   ({100 * summary['decoupled_fraction']:.1f}%)", transform=axis.transAxes, fontsize=18, fontweight="bold", color=INK)
    axis.text(0.59, 0.105, "Extent increased without positive FNI", transform=axis.transAxes, fontsize=9.5, color=TEXT)
    axis.text(0.50, 0.19, "BUT", transform=axis.transAxes, fontsize=10, fontweight="bold", color=PURPLE, ha="center")
    save_figure(figure, outdir, "Graphical_Abstract")


def run_ggplot_figure(results: Path, outdir: Path, rscript: str, r_library: str | None) -> None:
    script = Path(__file__).with_name("make_sensitivity_figure.R")
    environment = os.environ.copy()
    if r_library:
        environment["R_LIBS_USER"] = r_library
    subprocess.run([rscript, str(script), str(results), str(outdir)], check=True, env=environment)


def load_inputs(results: Path) -> dict[str, object]:
    return {
        "phase": pd.read_csv(results / "model" / "phase_space.csv"),
        "summary": json.loads((results / "model" / "model_summary.json").read_text(encoding="utf-8")),
        "effects": pd.read_csv(results / "omics" / "module_effects.csv"),
        "cross": pd.read_csv(results / "integration" / "cross_layer_processes.csv"),
        "cross_summary": json.loads((results / "integration" / "cross_layer_summary.json").read_text(encoding="utf-8")),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT_EXE", "Rscript"))
    parser.add_argument("--r-library", default=os.environ.get("R_LIBS_USER"))
    arguments = parser.parse_args()
    configure_style()
    results = Path(arguments.results).resolve()
    outdir = Path(arguments.outdir).resolve()
    inputs = load_inputs(results)
    figure1(outdir)
    figure2(outdir, inputs["phase"], inputs["summary"])
    run_ggplot_figure(results, outdir, arguments.rscript, arguments.r_library)
    figure4(outdir, inputs["effects"], inputs["cross"], inputs["cross_summary"])
    graphical_abstract(outdir, inputs["summary"])
    print(f"Wrote Matplotlib, Seaborn, and ggplot2 figures to {outdir}")


if __name__ == "__main__":
    main()
