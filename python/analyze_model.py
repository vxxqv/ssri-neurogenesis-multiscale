from __future__ import annotations

from pathlib import Path
import argparse
import itertools
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import beta
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold


CHANNELS = ["activation", "proliferation", "maturation", "integration", "efficacy", "survival"]


def correlation(x: pd.Series, y: pd.Series, method: str) -> float:
    return float(x.corr(y, method=method))


def posterior_probability(successes: pd.Series, trials: pd.Series) -> np.ndarray:
    return beta.sf(0.5, successes.to_numpy() + 0.5, trials.to_numpy() - successes.to_numpy() + 0.5)


def phase_table(replicates: pd.DataFrame, design: pd.DataFrame, n_reps: int | None = None) -> pd.DataFrame:
    frame = replicates.loc[replicates.design_group == "phase"].copy()
    if n_reps is not None:
        frame = frame.loc[frame.replicate <= n_reps]
    frame["mismatch"] = (frame.delta_extent > 0) & (frame.delta_fni <= 0)
    frame["aligned_gain"] = (frame.delta_extent > 0) & (frame.delta_fni > 0)
    for state in "QAPNMG":
        frame[f"delta_{state}"] = frame[f"treat_{state}"] - frame[f"control_{state}"]
    named = {
        "n_rep": ("replicate", "count"),
        "delta_extent": ("delta_extent", "mean"),
        "sd_delta_extent": ("delta_extent", "std"),
        "delta_fni": ("delta_fni", "mean"),
        "sd_delta_fni": ("delta_fni", "std"),
        "mismatch_reps": ("mismatch", "sum"),
        "aligned_gain_reps": ("aligned_gain", "sum"),
    }
    for state in "QAPNMG":
        named[f"delta_{state}"] = (f"delta_{state}", "mean")
    summary = frame.groupby(["set_id", "base_id"], as_index=False).agg(**named)
    summary["mismatch_probability"] = (summary.mismatch_reps + 0.5) / (summary.n_rep + 1.0)
    summary["aligned_gain_probability"] = (summary.aligned_gain_reps + 0.5) / (summary.n_rep + 1.0)
    summary["posterior_mismatch_gt_half"] = posterior_probability(summary.mismatch_reps, summary.n_rep)
    summary["posterior_aligned_gt_half"] = posterior_probability(summary.aligned_gain_reps, summary.n_rep)
    summary["response_class"] = np.select(
        [summary.posterior_mismatch_gt_half >= 0.95, summary.posterior_aligned_gt_half >= 0.95],
        ["robust_mismatch", "robust_aligned_gain"],
        default="uncertain_or_other",
    )
    phase_design = design.loc[design.design_group == "phase"].drop(columns=["reps", "seed", "t_end", "design_group", "model"])
    return phase_design.merge(summary, on=["set_id", "base_id"], validate="one_to_one")


def variance_decomposition(replicates: pd.DataFrame, outcome: str) -> dict[str, float]:
    frame = replicates.loc[replicates.design_group == "phase", ["set_id", outcome]]
    grouped = frame.groupby("set_id")[outcome]
    within = float(grouped.var().mean())
    means = grouped.mean()
    n = int(grouped.size().iloc[0])
    between = max(0.0, float(means.var(ddof=1)) - within / n)
    total = within + between
    return {
        "within_stochastic_variance": within,
        "between_parameter_variance": between,
        "intraclass_correlation": between / total if total > 0 else float("nan"),
    }


def shapley_values(factorial: pd.DataFrame, outcome: str) -> pd.DataFrame:
    n = len(CHANNELS)
    weights = {size: math.factorial(size) * math.factorial(n - size - 1) / math.factorial(n) for size in range(n)}
    records = []
    for base_id, group in factorial.groupby("base_id", sort=False):
        values = {tuple(int(value) for value in row.model): getattr(row, outcome) for row in group.itertuples(index=False)}
        if len(values) != 2**n:
            raise RuntimeError(f"incomplete factorial for {base_id}")
        for index, channel in enumerate(CHANNELS):
            contribution = 0.0
            others = [j for j in range(n) if j != index]
            for size in range(n):
                for subset in itertools.combinations(others, size):
                    off = [0] * n
                    for j in subset:
                        off[j] = 1
                    on = off.copy()
                    on[index] = 1
                    contribution += weights[size] * (values[tuple(on)] - values[tuple(off)])
            records.append({"base_id": base_id, "channel": channel, "outcome": outcome, "shapley": contribution})
    return pd.DataFrame(records)


def sobol_tables(simulation: pd.DataFrame, parameters: list[str], outcome: str, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    groups = {name: group.set_index("base_id")[outcome].sort_index() for name, group in simulation.groupby("design_group") if name.startswith("sobol_")}
    ya = groups["sobol_A"].to_numpy()
    yb = groups["sobol_B"].to_numpy()

    def estimate(indices: np.ndarray, parameter: str) -> tuple[float, float]:
        a, b = ya[indices], yb[indices]
        ab = groups[f"sobol_AB_{parameter}"].to_numpy()[indices]
        variance = np.var(np.concatenate([a, b]), ddof=1)
        return 1.0 - np.mean((b - ab) ** 2) / (2.0 * variance), np.mean((a - ab) ** 2) / (2.0 * variance)

    final = []
    convergence = []
    all_indices = np.arange(len(ya))
    for parameter in parameters:
        s1, st = estimate(all_indices, parameter)
        boot = np.array([estimate(rng.integers(0, len(ya), len(ya)), parameter) for _ in range(500)])
        final.append({
            "parameter": parameter,
            "outcome": outcome,
            "S1": s1,
            "S1_low": np.quantile(boot[:, 0], 0.025),
            "S1_high": np.quantile(boot[:, 0], 0.975),
            "ST": st,
            "ST_low": np.quantile(boot[:, 1], 0.025),
            "ST_high": np.quantile(boot[:, 1], 0.975),
        })
        for size in (128, 256, 512, len(ya)):
            if size <= len(ya):
                prefix_s1, prefix_st = estimate(np.arange(size), parameter)
                convergence.append({"parameter": parameter, "outcome": outcome, "base_n": size, "S1": prefix_s1, "ST": prefix_st})
    return pd.DataFrame(final), pd.DataFrame(convergence)


def predictive_information(phase: pd.DataFrame, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    feature_sets = {
        "extent_only": ["delta_extent"],
        "cell_composition": ["delta_M", "delta_G"],
        "process_aware": ["delta_M", "delta_G", "tx_integration", "tx_eff", "eff_mean", "eff_cv", "tx_survival"],
    }
    folds = KFold(n_splits=min(10, len(phase)), shuffle=True, random_state=20260904)
    predictions = []
    metrics = []
    y = phase.delta_fni.to_numpy()
    for name, columns in feature_sets.items():
        predicted = np.full(len(phase), np.nan)
        for fold, (train, test) in enumerate(folds.split(phase), 1):
            model = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=20, l2_regularization=1.0, random_state=fold)
            model.fit(phase.iloc[train][columns], y[train])
            predicted[test] = model.predict(phase.iloc[test][columns])
        if np.isnan(predicted).any():
            raise RuntimeError("missing out of fold predictions")
        boot = []
        for _ in range(1000):
            index = rng.integers(0, len(y), len(y))
            boot.append(r2_score(y[index], predicted[index]))
        metrics.append({
            "model": name,
            "features": ";".join(columns),
            "r2": r2_score(y, predicted),
            "r2_low": np.quantile(boot, 0.025),
            "r2_high": np.quantile(boot, 0.975),
            "mae": mean_absolute_error(y, predicted),
        })
        predictions.append(pd.DataFrame({"set_id": phase.set_id, "model": name, "observed": y, "predicted": predicted}))
    return pd.DataFrame(metrics), pd.concat(predictions, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", required=True)
    parser.add_argument("--simulation", required=True)
    parser.add_argument("--replicates", required=True)
    parser.add_argument("--deterministic", required=True)
    parser.add_argument("--priors", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=20260904)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    design = pd.read_csv(args.design, dtype={"model": str}, low_memory=False)
    simulation = pd.read_csv(args.simulation, dtype={"model": str})
    replicates = pd.read_csv(args.replicates, dtype={"model": str})
    deterministic = pd.read_csv(args.deterministic, dtype={"model": str})
    parameters = pd.read_csv(args.priors).parameter.tolist()
    if simulation.failed_reps.sum() or replicates.failed.sum():
        raise RuntimeError("simulation failures detected")

    phase = phase_table(replicates, design)
    phase.to_csv(outdir / "phase_space.csv", index=False)
    convergence = []
    for n_reps in (4, 8, 12, 16, 24):
        if n_reps <= replicates.loc[replicates.design_group == "phase", "replicate"].max():
            table = phase_table(replicates, design, n_reps)
            convergence.append({
                "replicates": n_reps,
                "pearson": correlation(table.delta_extent, table.delta_fni, "pearson"),
                "spearman": correlation(table.delta_extent, table.delta_fni, "spearman"),
                "robust_mismatch_fraction": float((table.response_class == "robust_mismatch").mean()),
                "replicate_mismatch_fraction": float(table.mismatch_reps.sum() / table.n_rep.sum()),
                "mean_mismatch_probability": float(table.mismatch_probability.mean()),
            })
    pd.DataFrame(convergence).to_csv(outdir / "replicate_convergence.csv", index=False)

    phase_sim = simulation.loc[simulation.design_group == "phase"]
    validation = phase_sim.merge(deterministic.loc[deterministic.design_group == "phase"], on=["set_id", "base_id", "design_group", "t_end"], suffixes=("_ssa", "_ode"), validate="one_to_one")
    validation[["set_id", "base_id", "delta_extent_ssa", "delta_extent_ode", "delta_fni_ssa", "delta_fni_ode"]].to_csv(outdir / "deterministic_validation.csv", index=False)

    factorial = deterministic.loc[deterministic.design_group == "factorial"].copy()
    shapley = pd.concat([shapley_values(factorial, "delta_extent"), shapley_values(factorial, "delta_fni")], ignore_index=True)
    shapley.to_csv(outdir / "channel_shapley.csv", index=False)
    shapley_summary = shapley.groupby(["channel", "outcome"], as_index=False).agg(mean=("shapley", "mean"), median=("shapley", "median"), q025=("shapley", lambda x: x.quantile(0.025)), q975=("shapley", lambda x: x.quantile(0.975)), positive_fraction=("shapley", lambda x: (x > 0).mean()))
    shapley_summary.to_csv(outdir / "channel_shapley_summary.csv", index=False)
    interactions = []
    for base_id, group in factorial.groupby("base_id", sort=False):
        values = {row.model: row for row in group.itertuples(index=False)}
        for left, right in itertools.combinations(range(len(CHANNELS)), 2):
            empty = "0" * len(CHANNELS)
            a = list(empty); a[left] = "1"; a = "".join(a)
            b = list(empty); b[right] = "1"; b = "".join(b)
            ab = list(empty); ab[left] = "1"; ab[right] = "1"; ab = "".join(ab)
            for outcome in ("delta_extent", "delta_fni"):
                value = getattr(values[ab], outcome) - getattr(values[a], outcome) - getattr(values[b], outcome) + getattr(values[empty], outcome)
                interactions.append({"base_id": base_id, "channel_a": CHANNELS[left], "channel_b": CHANNELS[right], "outcome": outcome, "interaction": value})
    interaction_frame = pd.DataFrame(interactions)
    interaction_frame.to_csv(outdir / "channel_pair_interactions.csv", index=False)
    interaction_frame.groupby(["channel_a", "channel_b", "outcome"], as_index=False).agg(mean=("interaction", "mean"), median=("interaction", "median"), q025=("interaction", lambda x: x.quantile(0.025)), q975=("interaction", lambda x: x.quantile(0.975))).to_csv(outdir / "channel_pair_interaction_summary.csv", index=False)

    sobol, sobol_convergence = [], []
    for outcome in ("delta_extent", "delta_fni"):
        final, trace = sobol_tables(deterministic, parameters, outcome, rng)
        sobol.append(final)
        sobol_convergence.append(trace)
    pd.concat(sobol, ignore_index=True).to_csv(outdir / "sobol_indices.csv", index=False)
    pd.concat(sobol_convergence, ignore_index=True).to_csv(outdir / "sobol_convergence.csv", index=False)

    information, predictions = predictive_information(phase, rng)
    information.to_csv(outdir / "predictive_information.csv", index=False)
    predictions.to_csv(outdir / "predictive_predictions.csv", index=False)

    temporal_rep = replicates.loc[replicates.design_group.str.startswith("temporal_")].copy()
    temporal_rep["mismatch"] = (temporal_rep.delta_extent > 0) & (temporal_rep.delta_fni <= 0)
    temporal_set = temporal_rep.groupby(["set_id", "base_id", "t_end"], as_index=False).agg(delta_extent=("delta_extent", "mean"), delta_fni=("delta_fni", "mean"), mismatch_probability=("mismatch", "mean"))
    temporal_set.to_csv(outdir / "temporal_sets.csv", index=False)
    temporal_set.groupby("t_end", as_index=False).agg(median_delta_extent=("delta_extent", "median"), extent_q025=("delta_extent", lambda x: x.quantile(0.025)), extent_q975=("delta_extent", lambda x: x.quantile(0.975)), median_delta_fni=("delta_fni", "median"), fni_q025=("delta_fni", lambda x: x.quantile(0.025)), fni_q975=("delta_fni", lambda x: x.quantile(0.975)), mean_mismatch_probability=("mismatch_probability", "mean")).to_csv(outdir / "temporal_summary.csv", index=False)

    summary = {
        "phase_parameter_sets": int(len(phase)),
        "replicates_per_phase_set": int(phase.n_rep.iloc[0]),
        "pearson_extent_fni": correlation(phase.delta_extent, phase.delta_fni, "pearson"),
        "spearman_extent_fni": correlation(phase.delta_extent, phase.delta_fni, "spearman"),
        "robust_mismatch_count": int((phase.response_class == "robust_mismatch").sum()),
        "robust_mismatch_fraction": float((phase.response_class == "robust_mismatch").mean()),
        "robust_aligned_count": int((phase.response_class == "robust_aligned_gain").sum()),
        "replicate_mismatch_fraction": float(phase.mismatch_reps.sum() / phase.n_rep.sum()),
        "mean_mismatch_probability": float(phase.mismatch_probability.mean()),
        "median_delta_extent": float(phase.delta_extent.median()),
        "median_delta_fni": float(phase.delta_fni.median()),
        "extent_variance": variance_decomposition(replicates, "delta_extent"),
        "fni_variance": variance_decomposition(replicates, "delta_fni"),
        "ode_ssa_extent_r": correlation(validation.delta_extent_ssa, validation.delta_extent_ode, "pearson"),
        "ode_ssa_fni_r": correlation(validation.delta_fni_ssa, validation.delta_fni_ode, "pearson"),
        "ode_ssa_extent_mae": float(mean_absolute_error(validation.delta_extent_ode, validation.delta_extent_ssa)),
        "ode_ssa_fni_mae": float(mean_absolute_error(validation.delta_fni_ode, validation.delta_fni_ssa)),
    }
    (outdir / "model_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
