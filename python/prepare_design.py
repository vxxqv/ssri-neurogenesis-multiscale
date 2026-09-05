from __future__ import annotations

from pathlib import Path
import argparse
import itertools
import numpy as np
import pandas as pd


FIXED = {"Q0": 60, "A0": 4, "P0": 8, "N0": 10, "M0": 15, "G0": 12, "tau": 14.0}
CHANNELS = ["activation", "proliferation", "maturation", "integration", "efficacy", "survival"]


def latin_hypercube(n: int, d: int, rng: np.random.Generator) -> np.ndarray:
    unit = rng.random((n, d))
    for j in range(d):
        unit[:, j] = (rng.permutation(n) + unit[:, j]) / n
    return unit


def transform(unit: np.ndarray, priors: pd.DataFrame) -> np.ndarray:
    values = np.empty_like(unit)
    for j, row in priors.reset_index(drop=True).iterrows():
        low, high = float(row.low), float(row.high)
        if row.scale == "log":
            values[:, j] = np.exp(np.log(low) + unit[:, j] * (np.log(high) - np.log(low)))
        else:
            values[:, j] = low + unit[:, j] * (high - low)
    return values


def base_frame(values: np.ndarray, priors: pd.DataFrame, prefix: str) -> pd.DataFrame:
    frame = pd.DataFrame(values, columns=priors.parameter)
    frame.insert(0, "base_id", [f"{prefix}{i:05d}" for i in range(len(frame))])
    for key, value in FIXED.items():
        frame[key] = value
    return frame


def scenario(frame: pd.DataFrame, group: str, model: str, reps: int, t_end: float, bits: tuple[int, ...]) -> pd.DataFrame:
    out = frame.copy()
    out.insert(1, "design_group", group)
    out.insert(2, "model", model)
    out.insert(3, "reps", reps)
    out.insert(4, "t_end", t_end)
    for channel, bit in zip(CHANNELS, bits):
        out[f"{channel}_on"] = bit
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--priors", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--phase-n", type=int, default=2048)
    parser.add_argument("--phase-reps", type=int, default=16)
    parser.add_argument("--factorial-n", type=int, default=512)
    parser.add_argument("--factorial-reps", type=int, default=0)
    parser.add_argument("--sobol-n", type=int, default=2048)
    parser.add_argument("--sobol-reps", type=int, default=0)
    parser.add_argument("--temporal-n", type=int, default=512)
    parser.add_argument("--temporal-reps", type=int, default=4)
    args = parser.parse_args()

    priors = pd.read_csv(args.priors)
    rng = np.random.default_rng(args.seed)
    full = tuple([1] * len(CHANNELS))
    pieces = []
    phase = base_frame(transform(latin_hypercube(args.phase_n, len(priors), rng), priors), priors, "P")
    pieces.append(scenario(phase, "phase", "full", args.phase_reps, 56.0, full))
    temporal = phase.iloc[: args.temporal_n].copy()
    for day in (7.0, 14.0, 28.0, 42.0, 56.0):
        pieces.append(scenario(temporal, f"temporal_{int(day)}", "full", args.temporal_reps, day, full))
    factorial = base_frame(transform(latin_hypercube(args.factorial_n, len(priors), rng), priors), priors, "F")
    for bits in itertools.product((0, 1), repeat=len(CHANNELS)):
        label = "".join(str(bit) for bit in bits)
        pieces.append(scenario(factorial, "factorial", label, args.factorial_reps, 56.0, bits))
    a_unit = rng.random((args.sobol_n, len(priors)))
    b_unit = rng.random((args.sobol_n, len(priors)))
    a = base_frame(transform(a_unit, priors), priors, "S")
    b = base_frame(transform(b_unit, priors), priors, "S")
    pieces.append(scenario(a, "sobol_A", "full", args.sobol_reps, 56.0, full))
    pieces.append(scenario(b, "sobol_B", "full", args.sobol_reps, 56.0, full))
    b_values = transform(b_unit, priors)
    for j, parameter in enumerate(priors.parameter):
        ab_values = transform(a_unit.copy(), priors)
        ab_values[:, j] = b_values[:, j]
        ab = base_frame(ab_values, priors, "S")
        pieces.append(scenario(ab, f"sobol_AB_{parameter}", "full", args.sobol_reps, 56.0, full))
    design = pd.concat(pieces, ignore_index=True)
    design.insert(0, "set_id", [f"S{i:07d}" for i in range(len(design))])
    seed_codes = pd.factorize(design.base_id, sort=True)[0]
    design.insert(5, "seed", args.seed + (seed_codes + 1) * 1000003)
    ordered = [
        "set_id", "base_id", "design_group", "model", "reps", "seed", "t_end",
        "Q0", "A0", "P0", "N0", "M0", "G0", *priors.parameter.tolist(), "tau",
        *[f"{channel}_on" for channel in CHANNELS],
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    design[ordered].to_csv(out, index=False)
    print(f"wrote {len(design)} design rows")


if __name__ == "__main__":
    main()
