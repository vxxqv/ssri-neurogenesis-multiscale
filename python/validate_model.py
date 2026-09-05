from __future__ import annotations

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from solve_deterministic import solve


def relative_error(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.abs(left - right) / np.maximum(np.abs(right), 1e-9)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", required=True)
    parser.add_argument("--shapley", required=True)
    parser.add_argument("--deterministic", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    design = pd.read_csv(args.design, dtype={"model": str}, low_memory=False)
    sample = design.groupby("design_group", group_keys=False).head(32).reset_index(drop=True)
    errors = []
    for _, group in sample.groupby("t_end"):
        coarse = solve(group.reset_index(drop=True), True, 0.25)
        fine = solve(group.reset_index(drop=True), True, 0.125)
        errors.append(relative_error(coarse, fine))
    error = np.vstack(errors)
    deterministic = pd.read_csv(args.deterministic, dtype={"model": str})
    null = deterministic.loc[(deterministic.design_group == "factorial") & (deterministic.model == "000000")]
    shapley = pd.read_csv(args.shapley)
    full = deterministic.loc[(deterministic.design_group == "factorial") & (deterministic.model == "111111")].set_index("base_id")
    efficiency = []
    for outcome in ("delta_extent", "delta_fni"):
        summed = shapley.loc[shapley.outcome == outcome].groupby("base_id").shapley.sum()
        efficiency.extend((summed - full[outcome]).abs().tolist())
    report = {
        "rk4_median_relative_error": float(np.median(error)),
        "rk4_max_relative_error": float(np.max(error)),
        "null_max_absolute_extent": float(null.delta_extent.abs().max()),
        "null_max_absolute_fni": float(null.delta_fni.abs().max()),
        "shapley_max_efficiency_error": float(np.max(efficiency)),
    }
    if report["rk4_max_relative_error"] > 1e-4 or report["null_max_absolute_extent"] > 1e-9 or report["null_max_absolute_fni"] > 1e-9 or report["shapley_max_efficiency_error"] > 1e-8:
        raise RuntimeError(json.dumps(report))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
