from __future__ import annotations

from pathlib import Path
import argparse
import numpy as np
import pandas as pd


FIXED = {
    "t_end": 56.0,
    "Q0": 60,
    "A0": 4,
    "P0": 8,
    "N0": 10,
    "M0": 15,
    "G0": 12,
    "tau": 14.0,
}


def latin_hypercube(n: int, d: int, rng: np.random.Generator) -> np.ndarray:
    u = rng.random((n, d))
    out = np.empty_like(u)
    for j in range(d):
        out[:, j] = (rng.permutation(n) + u[:, j]) / n
    return out


def transform(unit: np.ndarray, priors: pd.DataFrame) -> np.ndarray:
    out = np.empty_like(unit)
    for j, row in priors.reset_index(drop=True).iterrows():
        lo, hi = float(row.low), float(row.high)
        if row.scale == "log":
            out[:, j] = np.exp(np.log(lo) + unit[:, j] * (np.log(hi) - np.log(lo)))
        else:
            out[:, j] = lo + unit[:, j] * (hi - lo)
    return out


def frame(values: np.ndarray, priors: pd.DataFrame, group: str, model: str, reps: int) -> pd.DataFrame:
    df = pd.DataFrame(values, columns=priors.parameter)
    df.insert(0, "model", model)
    df.insert(0, "design_group", group)
    df.insert(0, "set_id", "")
    df.insert(3, "reps", reps)
    for key, value in FIXED.items():
        df[key] = value
    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--priors", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--phase-n", type=int, default=1500)
    parser.add_argument("--structural-n", type=int, default=300)
    parser.add_argument("--sobol-n", type=int, default=128)
    args = parser.parse_args()

    priors = pd.read_csv(args.priors)
    rng = np.random.default_rng(args.seed)
    d = len(priors)
    pieces: list[pd.DataFrame] = []

    phase = transform(latin_hypercube(args.phase_n, d, rng), priors)
    pieces.append(frame(phase, priors, "phase", "full", 4))

    structural = transform(latin_hypercube(args.structural_n, d, rng), priors)
    for model in ("baseline", "proliferation", "maturation", "integration", "full"):
        pieces.append(frame(structural.copy(), priors, f"structural_{model}", model, 3))

    a_unit = rng.random((args.sobol_n, d))
    b_unit = rng.random((args.sobol_n, d))
    a = transform(a_unit, priors)
    b = transform(b_unit, priors)
    pieces.append(frame(a, priors, "sobol_A", "full", 2))
    pieces.append(frame(b, priors, "sobol_B", "full", 2))
    for j, parameter in enumerate(priors.parameter):
        ab = a.copy()
        ab[:, j] = b[:, j]
        pieces.append(frame(ab, priors, f"sobol_AB_{parameter}", "full", 2))

    design = pd.concat(pieces, ignore_index=True)
    design["set_id"] = [f"S{i:06d}" for i in range(len(design))]
    design["seed"] = [args.seed + (i + 1) * 1000003 for i in range(len(design))]

    ordered = [
        "set_id", "design_group", "model", "reps", "seed", "t_end", "Q0", "A0", "P0", "N0", "M0", "G0",
        *priors.parameter.tolist(), "tau"
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    design[ordered].to_csv(out, index=False)
    print(f"wrote {len(design)} design rows to {out}")


if __name__ == "__main__":
    main()
