from pathlib import Path
import argparse
import itertools
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from analyze_public_transcriptomes import programs_by_species, read_delimited, gpl_mapping
from transcriptome_core import geo_metadata, series_expression, log_cpm, hedges_g


def score(matrix, left, right, mode):
    if mode == 'sample_rank':
        z = matrix.rank(axis=1, pct=True)
    else:
        z = (matrix - matrix.mean()) / matrix.std(ddof=1)
    a, b = z[left].mean(axis=1), z[right].mean(axis=1)
    if mode == 'unit_variance':
        a, b = a / a.std(ddof=1), b / b.std(ddof=1)
    return (a - b).to_numpy()


def test(values, labels):
    n = len(labels)
    masks = np.zeros((len(list(itertools.combinations(range(n), int(labels.sum())))), n))
    for i, chosen in enumerate(itertools.combinations(range(n), int(labels.sum()))):
        masks[i, list(chosen)] = 1
    weights = masks / labels.sum() - (1 - masks) / (~labels).sum()
    observed = values[labels].mean() - values[~labels].mean()
    null = weights @ values
    return float(np.mean(np.abs(null) >= abs(observed) - 1e-12)), hedges_g(values[labels], values[~labels])[0]


def audit(matrix, labels, programs, dataset, contrast, rng):
    matrix = matrix.loc[:, matrix.var() > 1e-12]
    measured = set(matrix.columns)
    left = sorted(measured & programs['integration_exclusive'])
    right = sorted(measured & programs['extent_exclusive'])
    rows = []
    for mode in ['unit_variance', 'gene_z', 'sample_rank']:
        values = score(matrix, left, right, mode)
        p, g = test(values, labels)
        loo = []
        for i in range(len(matrix)):
            keep = np.arange(len(matrix)) != i
            sub = matrix.iloc[keep]
            x = score(sub, left, right, mode)
            loo.append(test(x, labels[keep])[1])
        rows.append(dict(dataset=dataset, contrast=contrast, scoring=mode, p_value=p, hedges_g=g, loo_min_g=min(loo), loo_max_g=max(loo), loo_positive_fraction=float(np.mean(np.asarray(loo)>0)), integration_genes=len(left), extent_genes=len(right)))
    splits = []
    for iteration in range(200):
        a, b = rng.permutation(left), rng.permutation(right)
        effects = []
        for half in [0, 1]:
            x = score(matrix, list(a[half::2]), list(b[half::2]), 'gene_z')
            effects.append(hedges_g(x[labels], x[~labels])[0])
        splits.append(dict(dataset=dataset, contrast=contrast, split=iteration, first_half_g=effects[0], second_half_g=effects[1]))
    return rows, splits


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    root, out = Path(args.data), Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    programs = programs_by_species(root / 'gene_programs.tsv')['mouse']
    rng = np.random.default_rng(20260905)
    rows, splits = [], []
    for cell in ['GC', 'HC', 'MC']:
        raw = read_delimited(root / f'GSE309750_raw_count_{cell}_Samples.txt.gz').set_index('Symbol')
        counts = raw.groupby(level=0).sum().T.astype(float)
        labels = counts.index.str.contains('_FLX_')
        a, b = audit(log_cpm(counts), labels, programs, 'GSE309750', cell, rng)
        rows.extend(a)
        splits.extend(b)
    path = root / 'GSE43261_series_matrix.txt.gz'
    matrix = series_expression(path).astype(float)
    mapping = gpl_mapping(root / 'GPL1261.annot.gz')
    matrix = matrix[[x for x in matrix.columns if x in mapping]].rename(columns=mapping).T.groupby(level=0).median().T
    metadata = geo_metadata(path)
    for region in sorted(metadata.tissue.unique()):
        ids = metadata.index[(metadata.tissue == region) & metadata.response.isin(['Responder','Resistant'])]
        labels = (metadata.loc[ids,'response']=='Responder').to_numpy()
        a,b = audit(matrix.loc[ids], labels, programs, 'GSE43261', region, rng)
        rows.extend(a)
        splits.extend(b)
    pd.DataFrame(rows).to_csv(out / 'relative_score_robustness.csv', index=False)
    split_frame = pd.DataFrame(splits)
    split_frame.to_csv(out / 'relative_score_gene_splits.csv', index=False)
    summary = split_frame.assign(both_positive=lambda x:(x.first_half_g>0)&(x.second_half_g>0)).groupby(['dataset','contrast']).both_positive.mean()
    print(pd.DataFrame(rows).to_string(index=False))
    print(summary.to_string())


if __name__ == '__main__':
    main()
