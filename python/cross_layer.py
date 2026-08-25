from __future__ import annotations

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd


PARAMETER_TO_PROCESS = {
    "k_qa": "activation", "k_aq": "activation", "tx_activation": "activation",
    "b_a": "proliferation", "b_p": "proliferation", "tx_prolif": "proliferation",
    "k_pn": "differentiation",
    "k_nm": "maturation", "tx_maturation": "maturation",
    "k_mg": "integration", "tx_integration": "integration",
    "eff_mean": "activity", "eff_cv": "activity", "tx_eff": "activity",
    "d_a": "survival", "d_p": "survival", "d_n": "survival", "d_m": "survival", "d_g": "survival", "tx_survival": "survival",
}


def ranks(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    out = np.empty(len(x), dtype=float); out[order] = np.arange(len(x), dtype=float)
    return out


def rho(x: np.ndarray, y: np.ndarray) -> float:
    rx, ry = ranks(x), ranks(y)
    rx -= rx.mean(); ry -= ry.mean()
    return float(np.sum(rx*ry) / np.sqrt(np.sum(rx*rx)*np.sum(ry*ry)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sobol", required=True)
    parser.add_argument("--omics", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=20260825)
    args = parser.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    sobol = pd.read_csv(args.sobol)
    sobol = sobol.loc[sobol.outcome == "delta_fni"].copy()
    sobol["process"] = sobol.parameter.map(PARAMETER_TO_PROCESS)
    model = sobol.dropna(subset=["process"]).groupby("process", as_index=False).agg(
        model_ST=("ST", "max"), model_ST_mean=("ST", "mean"), n_parameters=("parameter", "count")
    )
    omics = pd.read_csv(args.omics)
    merged = model.merge(omics, on="process", how="inner")
    merged["omics_evidence"] = merged.median_abs_g
    merged.to_csv(outdir / "cross_layer_processes.csv", index=False)
    observed = rho(merged.model_ST.to_numpy(), merged.omics_evidence.to_numpy())
    rng = np.random.default_rng(args.seed)
    permutations = np.array([rho(merged.model_ST.to_numpy(), rng.permutation(merged.omics_evidence.to_numpy())) for _ in range(20000)])
    p = float((1 + np.sum(np.abs(permutations) >= abs(observed))) / (len(permutations) + 1))
    summary = {"n_mapped_processes": int(len(merged)), "spearman_rho": observed, "permutation_p_two_sided": p, "permutations": int(len(permutations))}
    (outdir / "cross_layer_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
