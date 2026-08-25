from __future__ import annotations

from pathlib import Path
from itertools import combinations
import argparse
import json
import math
import numpy as np
import pandas as pd


def parse_column(name: str) -> tuple[str, str, str]:
    parts = name.rsplit("_", 2)
    if len(parts) != 3 or parts[1] not in {"FT", "Sham"}:
        raise ValueError(f"unexpected pseudobulk column: {name}")
    return parts[0], parts[1], parts[2]


def standardized_effect(ft: np.ndarray, sham: np.ndarray) -> float:
    n1, n0 = len(ft), len(sham)
    s1, s0 = np.var(ft, ddof=1), np.var(sham, ddof=1)
    pooled = math.sqrt(((n1 - 1) * s1 + (n0 - 1) * s0) / max(1, n1 + n0 - 2))
    if pooled == 0:
        return 0.0
    d = (float(np.mean(ft)) - float(np.mean(sham))) / pooled
    correction = 1.0 - 3.0 / max(1.0, 4.0 * (n1 + n0) - 9.0)
    return d * correction


def permutation_p(ft: np.ndarray, sham: np.ndarray, rng: np.random.Generator) -> float:
    values = np.concatenate([ft, sham])
    n1 = len(ft)
    observed = abs(float(np.mean(ft) - np.mean(sham)))
    total_combinations = math.comb(len(values), n1)
    if total_combinations <= 50000:
        extreme = 0
        for idx in combinations(range(len(values)), n1):
            mask = np.zeros(len(values), dtype=bool); mask[list(idx)] = True
            diff = abs(float(values[mask].mean() - values[~mask].mean()))
            extreme += diff >= observed - 1e-12
        return extreme / total_combinations
    extreme = 0
    n_perm = 20000
    for _ in range(n_perm):
        perm = rng.permutation(values)
        extreme += abs(float(perm[:n1].mean() - perm[n1:].mean())) >= observed - 1e-12
    return (extreme + 1) / (n_perm + 1)


def bh(p: pd.Series) -> pd.Series:
    n = len(p)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=p.index)


def analyze_region(path: Path, region: str, modules: pd.DataFrame, rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    data = pd.read_csv(path, sep="\t", index_col=0)
    parsed = pd.DataFrame([parse_column(c) for c in data.columns], columns=["cell_type", "treatment", "replicate"], index=data.columns)
    records: list[dict] = []
    coverage: dict[str, dict] = {}
    lower_index = {str(g).lower(): str(g) for g in data.index}
    for process, module in modules.groupby("process", sort=False):
        genes = []
        signs = []
        for row in module.itertuples():
            actual = lower_index.get(str(row.gene).lower())
            if actual is not None:
                genes.append(actual); signs.append(float(row.sign))
        coverage[process] = {"requested": int(len(module)), "present": int(len(genes)), "genes": genes}
        if len(genes) < 3:
            continue
        expression = np.log1p(data.loc[genes].astype(float))
        signed = expression.mul(np.asarray(signs), axis=0)
        sample_scores = signed.mean(axis=0)
        for cell_type in parsed.cell_type.unique():
            cols = parsed.index[parsed.cell_type == cell_type]
            labels = parsed.loc[cols, "treatment"]
            ft = sample_scores.loc[cols[labels == "FT"]].to_numpy(dtype=float)
            sham = sample_scores.loc[cols[labels == "Sham"]].to_numpy(dtype=float)
            if len(ft) < 3 or len(sham) < 3:
                continue
            records.append({
                "region": region,
                "cell_type": cell_type,
                "process": process,
                "n_ft": len(ft),
                "n_sham": len(sham),
                "n_genes": len(genes),
                "mean_ft": float(np.mean(ft)),
                "mean_sham": float(np.mean(sham)),
                "difference": float(np.mean(ft) - np.mean(sham)),
                "hedges_g": standardized_effect(ft, sham),
                "p_permutation": permutation_p(ft, sham, rng),
                "genes_present": ";".join(genes),
            })
    result = pd.DataFrame(records)
    audit = {
        "region": region,
        "genes": int(data.shape[0]),
        "pseudobulk_profiles": int(data.shape[1]),
        "cell_types": int(parsed.cell_type.nunique()),
        "treatment_profiles": parsed.treatment.value_counts().to_dict(),
        "module_coverage": coverage,
    }
    return result, audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dorsal", required=True)
    parser.add_argument("--ventral", required=True)
    parser.add_argument("--modules", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=20260825)
    args = parser.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    modules = pd.read_csv(args.modules)
    rng = np.random.default_rng(args.seed)
    dorsal, d_audit = analyze_region(Path(args.dorsal), "dorsal_DG", modules, rng)
    ventral, v_audit = analyze_region(Path(args.ventral), "ventral_DG", modules, rng)
    all_results = pd.concat([dorsal, ventral], ignore_index=True)
    all_results["q_bh"] = bh(all_results.p_permutation)
    all_results.to_csv(outdir / "module_effects.csv", index=False)
    process = all_results.groupby("process", as_index=False).agg(
        n_tests=("hedges_g", "count"),
        median_abs_g=("hedges_g", lambda x: float(np.median(np.abs(x)))),
        max_abs_g=("hedges_g", lambda x: float(np.max(np.abs(x)))),
        median_g=("hedges_g", "median"),
        min_q=("q_bh", "min"),
        significant_tests=("q_bh", lambda x: int(np.sum(x < 0.10))),
    )
    process.to_csv(outdir / "process_omics_evidence.csv", index=False)
    audit = {
        "source": "Rayan et al. 2022 authors' deposited replicate-level TP10K pseudobulk matrices",
        "inference_unit": "pooled biological replicate within region and cell type",
        "hard_qc_repeated": False,
        "reason_qc_not_repeated": "analysis begins from post-QC deposited matrices; pre-QC barcodes and scDblFinder calls were not deposited in this processed object",
        "regions": [d_audit, v_audit],
    }
    (outdir / "transcriptomic_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(process.to_string(index=False))


if __name__ == "__main__":
    main()
