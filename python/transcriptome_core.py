from __future__ import annotations

from pathlib import Path
import csv
import gzip
import io
import itertools
import re
import warnings
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    out = np.empty(len(values))
    out[order] = np.minimum(adjusted, 1.0)
    return out


def geo_metadata(path: Path) -> pd.DataFrame:
    entries = []
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            row = next(csv.reader([line.rstrip()], delimiter="\t"))
            key, values = row[0].removeprefix("!Sample_"), row[1:]
            if key in {"title", "geo_accession", "description"}:
                entries.append((key, values))
            elif key == "characteristics_ch1":
                parsed = [value.split(": ", 1) if ": " in value else ["characteristic", value] for value in values]
                if len({item[0] for item in parsed}) == 1:
                    entries.append((parsed[0][0], [item[1] for item in parsed]))
    ids = next(values for key, values in entries if key == "geo_accession")
    data = {key: values for key, values in entries if len(values) == len(ids)}
    return pd.DataFrame(data, index=ids)


def series_expression(path: Path) -> pd.DataFrame:
    rows = []
    active = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                active = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if active:
                rows.append(line)
    return pd.read_csv(io.StringIO("".join(rows)), sep="\t", index_col=0).T


def parse_gtf(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    genes, transcripts = {}, {}
    pattern = re.compile(r'(gene_id|transcript_id|gene_name) "([^"]+)"')
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "transcript"}:
                continue
            attributes = dict(pattern.findall(fields[8]))
            gene_id = attributes.get("gene_id", "").split(".")[0]
            name = attributes.get("gene_name", gene_id)
            if gene_id:
                genes[gene_id] = name
            transcript_id = attributes.get("transcript_id", "").split(".")[0]
            if transcript_id:
                transcripts[transcript_id] = name
    return genes, transcripts


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    libraries = counts.sum(axis=1).replace(0, np.nan)
    return np.log2(counts.div(libraries, axis=0) * 1_000_000 + 0.5)


def hedges_g(treated: np.ndarray, control: np.ndarray) -> tuple[float, float]:
    n1, n0 = len(treated), len(control)
    pooled = ((n1 - 1) * np.var(treated, ddof=1) + (n0 - 1) * np.var(control, ddof=1)) / max(n1 + n0 - 2, 1)
    if pooled <= 0:
        return 0.0, float("inf")
    correction = 1.0 - 3.0 / max(4.0 * (n1 + n0) - 9.0, 1.0)
    effect = correction * (np.mean(treated) - np.mean(control)) / np.sqrt(pooled)
    standard_error = np.sqrt((n1 + n0) / (n1 * n0) + effect * effect / (2.0 * max(n1 + n0 - 2, 1)))
    return float(effect), float(standard_error)


def program_tests(expression: pd.DataFrame, groups: pd.Series, treated: str, control: str, programs: dict[str, set[str]], dataset: str, contrast: str, context: str) -> pd.DataFrame:
    selected = groups[groups.isin([treated, control])].index
    matrix = expression.loc[selected].copy()
    labels = groups.loc[selected]
    matrix = matrix.loc[:, matrix.var(axis=0) > 0]
    standardized = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=1)
    score_columns, sizes = {}, {}
    measured = set(matrix.columns)
    for program, genes in programs.items():
        overlap = sorted(measured & genes)
        if len(overlap) >= 10:
            score_columns[program] = standardized[overlap].mean(axis=1)
            sizes[program] = len(overlap)
    scores = pd.DataFrame(score_columns, index=selected)
    if {"integration_composite", "extent_composite"}.issubset(scores.columns):
        integration = scores["integration_composite"]
        extent = scores["extent_composite"]
        integration_sd = integration.std(ddof=1)
        extent_sd = extent.std(ddof=1)
        if integration_sd > 0 and extent_sd > 0:
            scores["integration_bias"] = (integration - integration.mean()) / integration_sd - (extent - extent.mean()) / extent_sd
            sizes["integration_bias"] = len(measured & (programs["integration_composite"] | programs["extent_composite"]))
    if {"integration_exclusive", "extent_exclusive"}.issubset(scores.columns):
        integration = scores["integration_exclusive"]
        extent = scores["extent_exclusive"]
        integration_sd = integration.std(ddof=1)
        extent_sd = extent.std(ddof=1)
        if integration_sd > 0 and extent_sd > 0:
            scores["exclusive_integration_bias"] = (integration - integration.mean()) / integration_sd - (extent - extent.mean()) / extent_sd
            sizes["exclusive_integration_bias"] = len(measured & (programs["integration_exclusive"] | programs["extent_exclusive"]))
    n_treated = int((labels == treated).sum())
    positions = np.arange(len(selected))
    combinations = list(itertools.combinations(positions, n_treated))
    assignments = np.zeros((len(combinations), len(selected)), dtype=bool)
    for row, choice in enumerate(combinations):
        assignments[row, list(choice)] = True
    records = []
    observed_mask = (labels == treated).to_numpy()
    for program in scores.columns:
        values = scores[program].to_numpy()
        observed = values[observed_mask].mean() - values[~observed_mask].mean()
        null = np.array([values[mask].mean() - values[~mask].mean() for mask in assignments])
        p_value = float(np.mean(np.abs(null) >= abs(observed) - 1e-12))
        effect, standard_error = hedges_g(values[observed_mask], values[~observed_mask])
        records.append({
            "dataset": dataset,
            "contrast": contrast,
            "context": context,
            "treated": treated,
            "control": control,
            "program": program,
            "measured_genes": sizes[program],
            "n_treated": n_treated,
            "n_control": len(selected) - n_treated,
            "score_difference": observed,
            "hedges_g": effect,
            "standard_error": standard_error,
            "p_value": p_value,
            "permutations": len(combinations),
        })
    return pd.DataFrame(records)


def welch_table(expression: pd.DataFrame, groups: pd.Series, treated: str, control: str, dataset: str, contrast: str) -> pd.DataFrame:
    selected = groups[groups.isin([treated, control])].index
    matrix = expression.loc[selected]
    matrix = matrix.loc[:, matrix.var(axis=0) > 1e-12]
    left = matrix.loc[groups.loc[selected] == treated]
    right = matrix.loc[groups.loc[selected] == control]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        test = ttest_ind(left.to_numpy(), right.to_numpy(), axis=0, equal_var=False, nan_policy="omit")
    frame = pd.DataFrame({
        "dataset": dataset,
        "contrast": contrast,
        "gene": matrix.columns,
        "effect": left.mean(axis=0).to_numpy() - right.mean(axis=0).to_numpy(),
        "statistic": test.statistic,
        "p_value": test.pvalue,
    }).dropna(subset=["p_value"])
    frame["q_value"] = bh(frame.p_value.to_numpy())
    return frame
