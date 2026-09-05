from __future__ import annotations

from pathlib import Path
import argparse
import gzip
import io
import itertools
import json
import re
import tarfile
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from transcriptome_core import bh, geo_metadata, hedges_g, log_cpm, parse_gtf, program_tests, series_expression, welch_table


def programs_by_species(path: Path) -> dict[str, dict[str, set[str]]]:
    frame = pd.read_csv(path, sep="\t")
    out = {}
    for species, group in frame.groupby("species"):
        species_programs = {program: set(values.gene) for program, values in group.groupby("program")}
        species_programs["extent_exclusive"] = species_programs["extent_composite"] - species_programs["integration_composite"]
        species_programs["integration_exclusive"] = species_programs["integration_composite"] - species_programs["extent_composite"]
        out[species] = species_programs
    return out


def read_delimited(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rb") as handle:
        prefix = handle.read(2)
    encoding = "utf-16" if prefix in {b"\xff\xfe", b"\xfe\xff"} else "utf-8"
    return pd.read_csv(path, sep="\t", compression="gzip", encoding=encoding)


def deseq_table(counts: pd.DataFrame, groups: pd.Series, treated: str, control: str, dataset: str, contrast: str) -> pd.DataFrame:
    selected = groups[groups.isin([treated, control])].index
    matrix = counts.loc[selected].round().astype(int)
    matrix = matrix.loc[:, (matrix.sum(axis=0) >= 10) & ((matrix > 0).sum(axis=0) >= 2)]
    metadata = pd.DataFrame({"condition": groups.loc[selected]}, index=selected)
    dds = DeseqDataSet(counts=matrix, metadata=metadata, design="~condition", refit_cooks=True, quiet=True, n_cpus=1)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", treated, control], quiet=True, n_cpus=1)
    stats.summary()
    frame = stats.results_df.reset_index()
    frame = frame.rename(columns={frame.columns[0]: "gene", "log2FoldChange": "effect", "stat": "statistic", "pvalue": "p_value", "padj": "q_value"})
    frame.insert(0, "contrast", contrast)
    frame.insert(0, "dataset", dataset)
    return frame[["dataset", "contrast", "gene", "effect", "statistic", "p_value", "q_value"]]


def gpl_mapping(path: Path) -> dict[str, str]:
    rows = []
    active = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("ID\t"):
                active = True
            if active:
                rows.append(line)
    frame = pd.read_csv(io.StringIO("".join(rows)), sep="\t", dtype=str)
    frame = frame.loc[frame["Gene symbol"].notna() & ~frame["Gene symbol"].str.contains("///", regex=False) & (frame["Gene symbol"] != "---")]
    return dict(zip(frame.ID, frame["Gene symbol"]))


def quant_gene_matrix(path: Path, transcript_map: dict[str, str]) -> pd.DataFrame:
    samples = {}
    with tarfile.open(path) as archive:
        for member in archive.getmembers():
            if not member.name.endswith(".sf.gz"):
                continue
            raw = archive.extractfile(member).read()
            quant = pd.read_csv(gzip.GzipFile(fileobj=io.BytesIO(raw)), sep="\t")
            quant["gene"] = quant.Name.str.split(".").str[0].map(transcript_map)
            values = quant.dropna(subset=["gene"]).groupby("gene").TPM.sum()
            name = re.sub(r"^GSM\d+_", "", member.name).removesuffix(".sf.gz")
            samples[name] = values
    return pd.DataFrame(samples).T.fillna(0.0)


def random_effects(frame: pd.DataFrame) -> dict[str, float]:
    y = frame.hedges_g.to_numpy(float)
    v = np.square(frame.standard_error.to_numpy(float))
    fixed = 1.0 / v
    fixed_mean = np.sum(fixed * y) / np.sum(fixed)
    q = np.sum(fixed * np.square(y - fixed_mean))
    c = np.sum(fixed) - np.sum(np.square(fixed)) / np.sum(fixed)
    tau2 = max(0.0, (q - len(y) + 1) / c) if c > 0 else 0.0
    weights = 1.0 / (v + tau2)
    estimate = np.sum(weights * y) / np.sum(weights)
    se = np.sqrt(1.0 / np.sum(weights))
    return {"studies": int(len(y)), "estimate": float(estimate), "se": float(se), "ci_low": float(estimate - 1.96 * se), "ci_high": float(estimate + 1.96 * se), "tau2": float(tau2), "q": float(q)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--rayan-dorsal", required=True)
    parser.add_argument("--rayan-ventral", required=True)
    parser.add_argument("--programs", required=True)
    parser.add_argument("--mouse-gtf", required=True)
    parser.add_argument("--rat-gtf", required=True)
    parser.add_argument("--gpl1261", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    root = Path(args.data)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    programs = programs_by_species(Path(args.programs))
    mouse_genes, mouse_transcripts = parse_gtf(Path(args.mouse_gtf))
    rat_genes, _ = parse_gtf(Path(args.rat_gtf))
    tests, genes, qc = [], [], []

    for region, path in [("dorDG", Path(args.rayan_dorsal)), ("venDG", Path(args.rayan_ventral))]:
        raw = pd.read_csv(path, sep="\t")
        raw.index.name = "gene"
        expression = np.log2(raw.T.astype(float) + 0.5)
        metadata = pd.DataFrame(index=expression.index)
        parsed = metadata.index.to_series().str.extract(r"^(.*)_(FT|Sham)_(Rep\d+)$")
        metadata[["cell_type", "condition", "replicate_name"]] = parsed.to_numpy()
        qc.append({"dataset": "GSE197622", "context": region, "samples": len(expression), "genes": expression.shape[1], "median_library": float(raw.sum(axis=0).median())})
        for cell_type in sorted(metadata.cell_type.dropna().unique()):
            ids = metadata.index[metadata.cell_type == cell_type]
            contrast = f"{region}_{cell_type}_FT_vs_Sham"
            result = program_tests(expression.loc[ids], metadata.condition.loc[ids], "FT", "Sham", programs["rat"], "GSE197622", contrast, f"{region}:{cell_type}")
            result["role"] = "cell_type_map"
            tests.append(result)
            genes.append(welch_table(expression.loc[ids], metadata.condition.loc[ids], "FT", "Sham", "GSE197622", contrast).nsmallest(200, "p_value"))

    matrix_path = root / "GSE43261_series_matrix.txt.gz"
    expression = series_expression(matrix_path).astype(float)
    mapping = gpl_mapping(Path(args.gpl1261))
    keep = [probe for probe in expression.columns if probe in mapping]
    expression = expression[keep].rename(columns=mapping).T.groupby(level=0).median().T
    metadata = geo_metadata(matrix_path)
    qc.append({"dataset": "GSE43261", "context": "dentate_gyrus", "samples": len(expression), "genes": expression.shape[1], "median_library": float("nan")})
    for region in sorted(metadata.tissue.unique()):
        ids = metadata.index[metadata.tissue == region]
        short = "dorsal" if region.startswith("Dorsal") else "ventral"
        for treated, control, role in [("Responder", "Resistant", "response_validation"), ("Responder", "control", "response_reference"), ("Resistant", "control", "response_reference")]:
            contrast = f"{short}_{treated}_vs_{control}"
            result = program_tests(expression.loc[ids], metadata.response.loc[ids], treated, control, programs["mouse"], "GSE43261", contrast, short)
            result["role"] = role
            tests.append(result)
            if role == "response_validation":
                genes.append(welch_table(expression.loc[ids], metadata.response.loc[ids], treated, control, "GSE43261", contrast))

    for code, label in [("GC", "granule_cells"), ("HC", "hippocampal_circuit"), ("MC", "mossy_cells")]:
        path = root / f"GSE309750_raw_count_{code}_Samples.txt.gz"
        raw = read_delimited(path).set_index("Symbol")
        counts = raw.T.groupby(level=0).sum().astype(float)
        counts = counts.T.groupby(level=0).sum().T
        group = pd.Series(np.where(counts.index.str.contains("_FLX_"), "Fluoxetine", "Vehicle"), index=counts.index)
        expression = log_cpm(counts)
        contrast = f"{label}_Fluoxetine_vs_Vehicle"
        result = program_tests(expression, group, "Fluoxetine", "Vehicle", programs["mouse"], "GSE309750", contrast, label)
        result["role"] = "circuit_validation"
        tests.append(result)
        genes.append(deseq_table(counts, group, "Fluoxetine", "Vehicle", "GSE309750", contrast))
        qc.append({"dataset": "GSE309750", "context": label, "samples": len(counts), "genes": counts.shape[1], "median_library": float(counts.sum(axis=1).median())})

    counts = pd.read_csv(root / "GSE222756_merged_gene_counts_sorted.txt.gz", sep="\t").set_index("gene_names").T.astype(float)
    counts = counts.T.groupby(level=0).sum().T
    metadata = geo_metadata(root / "GSE222756_series_matrix.txt.gz")
    metadata.index = metadata.title
    expression = log_cpm(counts)
    qc.append({"dataset": "GSE222756", "context": "ten_regions", "samples": len(counts), "genes": counts.shape[1], "median_library": float(counts.sum(axis=1).median())})
    for region in sorted(metadata["brain region_name"].unique()):
        ids = metadata.index[metadata["brain region_name"] == region]
        for treated in ("FT", "Bupropion", "Desipramine"):
            contrast = f"{region}_{treated}_vs_Control"
            result = program_tests(expression.loc[ids], metadata.treatment.loc[ids], treated, "Control", programs["rat"], "GSE222756", contrast, region)
            result["role"] = "regional_atlas"
            tests.append(result)
            if treated == "FT" and region in {"dorDG", "venDG"}:
                genes.append(deseq_table(counts.loc[ids], metadata.treatment.loc[ids], treated, "Control", "GSE222756", contrast))

    raw = pd.read_excel(root / "GSE205325_CUS_Rats_Counts.xlsx")
    raw["gene"] = raw.Geneid.str.split(".").str[0].map(rat_genes)
    counts = raw.dropna(subset=["gene"]).drop(columns=["Geneid"]).groupby("gene").sum().T.astype(float)
    groups = pd.Series([re.match(r"[A-Z]+", name).group(0) for name in counts.index], index=counts.index)
    expression = log_cpm(counts)
    result = program_tests(expression, groups, "F", "S", programs["rat"], "GSE205325", "stress_Fluoxetine_vs_stress", "hippocampus")
    result["role"] = "disease_context"
    tests.append(result)
    genes.append(deseq_table(counts, groups, "F", "S", "GSE205325", "stress_Fluoxetine_vs_stress"))
    qc.append({"dataset": "GSE205325", "context": "hippocampus", "samples": len(counts), "genes": counts.shape[1], "median_library": float(counts.sum(axis=1).median())})

    tpm = quant_gene_matrix(root / "GSE292948_RAW.tar", mouse_transcripts)
    expression = np.log2(tpm + 0.5)
    metadata = pd.DataFrame(index=expression.index)
    parsed = metadata.index.to_series().str.extract(r"^(vitro|vivo)_(Fluoxetine|Sertraline|Citalopram|Control)_(\d+)$")
    metadata[["context", "drug", "replicate_name"]] = parsed.to_numpy()
    qc.append({"dataset": "GSE292948", "context": "cross_ssri", "samples": len(tpm), "genes": tpm.shape[1], "median_library": float(tpm.sum(axis=1).median())})
    for context in ("vitro", "vivo"):
        ids = metadata.index[metadata.context == context]
        for drug in ("Fluoxetine", "Sertraline", "Citalopram"):
            contrast = f"{context}_{drug}_vs_Control"
            result = program_tests(expression.loc[ids], metadata.drug.loc[ids], drug, "Control", programs["mouse"], "GSE292948", contrast, context)
            result["role"] = "cross_ssri"
            tests.append(result)
            genes.append(welch_table(expression.loc[ids], metadata.drug.loc[ids], drug, "Control", "GSE292948", contrast))

    modules = pd.concat(tests, ignore_index=True)
    modules["q_global"] = bh(modules.p_value.to_numpy())
    modules["q_dataset"] = modules.groupby("dataset").p_value.transform(lambda values: bh(values.to_numpy()))
    modules.to_csv(outdir / "program_effects.csv", index=False)
    gene_frame = pd.concat(genes, ignore_index=True)
    gene_frame["within_contrast_rank"] = gene_frame.groupby(["dataset", "contrast"]).p_value.rank(method="first")
    gene_frame = gene_frame.loc[(gene_frame.within_contrast_rank <= 250) | (gene_frame.q_value < 0.05)]
    gene_frame.to_csv(outdir / "gene_effects_top.csv", index=False)
    pd.DataFrame(qc).to_csv(outdir / "dataset_qc.csv", index=False)

    integration = modules.loc[modules.program == "integration_composite"].copy()
    bias = modules.loc[modules.program == "integration_bias"].copy()
    meta_rows = []
    for frame, suffix in [(integration, "integration"), (bias, "integration_bias")]:
        chosen = pd.concat([
            frame.loc[(frame.dataset == "GSE197622") & frame.context.str.contains("Granule")],
            frame.loc[(frame.dataset == "GSE309750") & (frame.context == "granule_cells")],
            frame.loc[(frame.dataset == "GSE222756") & (frame.control == "Control") & (frame.treated == "FT") & frame.context.isin(["dorDG", "venDG"])],
            frame.loc[frame.dataset == "GSE205325"],
            frame.loc[(frame.dataset == "GSE292948") & (frame.context == "vivo") & (frame.treated == "Fluoxetine")],
        ])
        collapsed = chosen.groupby("dataset", as_index=False).agg(hedges_g=("hedges_g", "mean"), standard_error=("standard_error", lambda x: float(np.sqrt(np.square(x).sum()) / len(x))))
        meta_rows.append({"analysis": f"fluoxetine_{suffix}_across_datasets", **random_effects(collapsed)})
        response = frame.loc[frame.role == "response_validation", ["hedges_g", "standard_error"]]
        meta_rows.append({"analysis": f"responder_vs_resistant_{suffix}", **random_effects(response)})
    pd.DataFrame(meta_rows).to_csv(outdir / "meta_analysis.csv", index=False)

    atlas = integration.loc[(integration.dataset == "GSE222756") & (integration.treated == "FT") & (integration.control == "Control")]
    effects = atlas.set_index("context").hedges_g
    regions = effects.index.tolist()
    observed = effects.loc[["dorDG", "venDG"]].mean() - effects.drop(["dorDG", "venDG"]).mean()
    null = []
    for selected in itertools.combinations(regions, 2):
        other = [region for region in regions if region not in selected]
        null.append(effects.loc[list(selected)].mean() - effects.loc[other].mean())
    cross = integration.loc[integration.role == "cross_ssri"]
    cross_pairs = modules.loc[modules.role == "cross_ssri"].pivot_table(index=["treated", "program"], columns="context", values="hedges_g").reset_index()
    cross_pairs["same_direction"] = np.sign(cross_pairs.vitro) == np.sign(cross_pairs.vivo)
    cross_pairs["context_difference"] = cross_pairs.vivo - cross_pairs.vitro
    cross_pairs.to_csv(outdir / "cross_ssri_context.csv", index=False)
    atlas_drugs = modules.loc[(modules.dataset == "GSE222756") & (modules.program.isin(["extent_composite", "integration_composite", "integration_bias", "exclusive_integration_bias"]))]
    atlas_wide = atlas_drugs.pivot_table(index=["context", "program"], columns="treated", values="hedges_g").reset_index()
    atlas_wide["three_drug_direction"] = np.where((atlas_wide[["FT", "Bupropion", "Desipramine"]] > 0).all(axis=1), "positive", np.where((atlas_wide[["FT", "Bupropion", "Desipramine"]] < 0).all(axis=1), "negative", "mixed"))
    atlas_wide.to_csv(outdir / "regional_drug_concordance.csv", index=False)
    summary = {
        "atlas_dg_specificity": {
            "difference": float(observed),
            "exact_p": float(np.mean(np.abs(null) >= abs(observed) - 1e-12)),
            "permutations": len(null),
        },
        "cross_ssri_integration": {
            "positive_contrasts": int((cross.hedges_g > 0).sum()),
            "contrasts": int(len(cross)),
            "median_effect": float(cross.hedges_g.median()),
        },
        "cross_ssri_context_concordance": {
            "same_direction": int(cross_pairs.same_direction.sum()),
            "comparisons": int(len(cross_pairs)),
        },
        "globally_significant_program_tests": int((modules.q_global < 0.05).sum()),
        "dataset_significant_program_tests": int((modules.q_dataset < 0.05).sum()),
        "total_program_tests": int(len(modules)),
    }
    (outdir / "transcriptome_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
