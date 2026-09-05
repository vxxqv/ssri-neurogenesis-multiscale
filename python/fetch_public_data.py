from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import io
import zipfile
import requests
import pandas as pd


FILES = {
    "GSE309750_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309750/matrix/GSE309750_series_matrix.txt.gz",
    "GSE309750_raw_count_GC_Samples.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309750/suppl/GSE309750_raw_count_GC_Samples.txt.gz",
    "GSE309750_raw_count_HC_Samples.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309750/suppl/GSE309750_raw_count_HC_Samples.txt.gz",
    "GSE309750_raw_count_MC_Samples.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309750/suppl/GSE309750_raw_count_MC_Samples.txt.gz",
    "GSE292948_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292948/matrix/GSE292948_series_matrix.txt.gz",
    "GSE292948_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292948/suppl/GSE292948_RAW.tar",
    "GSE222756_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222756/matrix/GSE222756_series_matrix.txt.gz",
    "GSE222756_merged_gene_counts_sorted.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222756/suppl/GSE222756_merged_gene_counts_sorted.txt.gz",
    "GSE205325_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205325/matrix/GSE205325_series_matrix.txt.gz",
    "GSE205325_CUS_Rats_Counts.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205325/suppl/GSE205325_CUS_Rats_Counts.xlsx",
    "GSE43261_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43261/matrix/GSE43261_series_matrix.txt.gz",
    "go-basic.obo": "http://purl.obolibrary.org/obo/go/go-basic.obo",
    "mgi.gaf.gz": "https://current.geneontology.org/annotations/mgi.gaf.gz",
    "rgd.gaf.gz": "https://current.geneontology.org/annotations/rgd.gaf.gz",
    "Mus_musculus.GRCm39.116.gtf.gz": "https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz",
    "Rattus_norvegicus.GRCr8.116.gtf.gz": "https://ftp.ensembl.org/pub/release-116/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.116.gtf.gz",
    "GPL1261.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1261/annot/GPL1261.annot.gz",
}

RAYAN_URL = "https://codeload.github.com/arulrayan/Integrative-multi-omics-landscape-of-fluoxetine-action-across-27-brain-regions/zip/c3ae3f6d2a98a81d6b0208f19b70fcc56dc50c5f"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, target: Path) -> None:
    if target.exists():
        return
    with requests.get(url, stream=True, timeout=180) as response:
        response.raise_for_status()
        with target.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                handle.write(chunk)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, url in FILES.items():
        target = outdir / name
        download(url, target)
        rows.append({"file": name, "source": url, "sha256": sha256(target), "bytes": target.stat().st_size})
    archive = outdir / "rayan2022.zip"
    download(RAYAN_URL, archive)
    with zipfile.ZipFile(archive) as source:
        for region in ("dorDG", "venDG"):
            suffix = f"single-cell DEGs/Single-cell collated tp10k_{region}.txt"
            member = next(name for name in source.namelist() if name.endswith(suffix))
            target = outdir / f"GSE197622_tp10k_{region}.txt"
            if not target.exists():
                target.write_bytes(source.read(member))
            rows.append({"file": target.name, "source": RAYAN_URL, "sha256": sha256(target), "bytes": target.stat().st_size})
    pd.DataFrame(rows).to_csv(outdir / "source_checksums.tsv", sep="\t", index=False)
    print(f"verified {len(rows)} public files")


if __name__ == "__main__":
    main()
