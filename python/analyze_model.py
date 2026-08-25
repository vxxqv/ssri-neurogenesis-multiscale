from __future__ import annotations

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd


def ranks(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    out = np.empty(len(x), dtype=float)
    out[order] = np.arange(len(x), dtype=float)
    values, inverse, counts = np.unique(x, return_inverse=True, return_counts=True)
    if np.any(counts > 1):
        sums = np.bincount(inverse, weights=out)
        out = sums[inverse] / counts[inverse]
    return out


def corr(x: np.ndarray, y: np.ndarray) -> float:
    x = x - x.mean(); y = y - y.mean()
    den = np.sqrt(np.sum(x*x) * np.sum(y*y))
    return float(np.sum(x*y) / den) if den > 0 else float("nan")


def prcc(x: np.ndarray, y: np.ndarray, names: list[str]) -> pd.DataFrame:
    ranked_x = np.column_stack([ranks(x[:, j]) for j in range(x.shape[1])])
    ranked_y = ranks(y)
    z = np.column_stack([ranked_x, ranked_y])
    c = np.corrcoef(z, rowvar=False)
    precision = np.linalg.pinv(c)
    target = len(names)
    vals = []
    for j, name in enumerate(names):
        value = -precision[j, target] / np.sqrt(precision[j, j] * precision[target, target])
        vals.append((name, float(value)))
    return pd.DataFrame(vals, columns=["parameter", "prcc"])


def sobol_indices(merged: pd.DataFrame, parameters: list[str], outcome: str) -> pd.DataFrame:
    ya = merged.loc[merged.design_group == "sobol_A", outcome].to_numpy()
    yb = merged.loc[merged.design_group == "sobol_B", outcome].to_numpy()
    variance = np.var(np.concatenate([ya, yb]), ddof=1)
    records = []
    for parameter in parameters:
        yab = merged.loc[merged.design_group == f"sobol_AB_{parameter}", outcome].to_numpy()
        first = 1.0 - np.mean((yb - yab) ** 2) / (2.0 * variance)
        total = np.mean((ya - yab) ** 2) / (2.0 * variance)
        records.append({"parameter": parameter, "outcome": outcome, "S1": first, "ST": total})
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", required=True)
    parser.add_argument("--simulation", required=True)
    parser.add_argument("--priors", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    design = pd.read_csv(args.design)
    sim = pd.read_csv(args.simulation)
    merged = design.merge(sim, on=["set_id", "design_group", "model"], validate="one_to_one")
    if not np.array_equal(merged.reps.to_numpy(), merged.n_rep.to_numpy()):
        raise RuntimeError("replicate count mismatch between design and simulation output")
    priors = pd.read_csv(args.priors)
    parameters = priors.parameter.tolist()

    if int(merged.failed_reps.sum()) != 0:
        raise RuntimeError("simulation failures detected")

    phase = merged.loc[merged.design_group == "phase"].copy()
    phase["decoupled"] = (phase.delta_extent > 0) & (phase.delta_fni <= 0)
    phase.to_csv(outdir / "phase_space.csv", index=False)

    summary = {
        "n_phase_parameter_sets": int(len(phase)),
        "n_stochastic_replicates_per_set": int(phase.n_rep.iloc[0]),
        "pearson_delta_extent_fni": corr(phase.delta_extent.to_numpy(), phase.delta_fni.to_numpy()),
        "spearman_delta_extent_fni": corr(ranks(phase.delta_extent.to_numpy()), ranks(phase.delta_fni.to_numpy())),
        "decoupled_count": int(phase.decoupled.sum()),
        "decoupled_fraction": float(phase.decoupled.mean()),
        "median_delta_extent": float(phase.delta_extent.median()),
        "median_delta_fni": float(phase.delta_fni.median()),
        "delta_extent_q025": float(phase.delta_extent.quantile(0.025)),
        "delta_extent_q975": float(phase.delta_extent.quantile(0.975)),
        "delta_fni_q025": float(phase.delta_fni.quantile(0.025)),
        "delta_fni_q975": float(phase.delta_fni.quantile(0.975)),
    }
    (outdir / "model_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    structural = merged.loc[merged.design_group.str.startswith("structural_")].copy()
    structural["positive_extent"] = structural.delta_extent > 0
    structural["positive_fni"] = structural.delta_fni > 0
    structural["decoupled"] = structural.positive_extent & ~structural.positive_fni
    structural_summary = structural.groupby("model", as_index=False).agg(
        n=("set_id", "count"),
        median_delta_extent=("delta_extent", "median"),
        median_delta_fni=("delta_fni", "median"),
        decoupled_fraction=("decoupled", "mean"),
        positive_extent_fraction=("positive_extent", "mean"),
        positive_fni_fraction=("positive_fni", "mean"),
    )
    structural_summary.to_csv(outdir / "structural_contrasts.csv", index=False)

    sobol = pd.concat([
        sobol_indices(merged, parameters, "delta_extent"),
        sobol_indices(merged, parameters, "delta_fni"),
    ], ignore_index=True)
    sobol.to_csv(outdir / "sobol_indices.csv", index=False)

    x = phase[parameters].to_numpy()
    prcc_extent = prcc(x, phase.delta_extent.to_numpy(), parameters).assign(outcome="delta_extent")
    prcc_fni = prcc(x, phase.delta_fni.to_numpy(), parameters).assign(outcome="delta_fni")
    pd.concat([prcc_extent, prcc_fni], ignore_index=True).to_csv(outdir / "prcc.csv", index=False)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
