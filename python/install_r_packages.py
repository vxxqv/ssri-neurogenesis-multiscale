from __future__ import annotations

from pathlib import Path
import argparse
import gzip
import re
import subprocess
import requests


INDEX = "https://cloud.r-project.org/bin/windows/contrib/4.6"


def records(raw: str) -> dict[str, dict[str, str]]:
    parsed = {}
    for block in raw.split("\n\n"):
        fields = {}
        key = None
        for line in block.splitlines():
            if line.startswith(" ") and key:
                fields[key] += line.strip()
            elif ": " in line:
                key, value = line.split(": ", 1)
                fields[key] = value
        if "Package" in fields:
            parsed[fields["Package"]] = fields
    return parsed


def dependencies(package: str, index: dict[str, dict[str, str]]) -> list[str]:
    values = []
    for field in ("Depends", "Imports", "LinkingTo"):
        for item in index[package].get(field, "").split(","):
            name = re.sub(r"\s*\(.*\)", "", item).strip()
            if name and name != "R" and name in index:
                values.append(name)
    return values


def resolve(package: str, index: dict[str, dict[str, str]], ordered: list[str], active: set[str]) -> None:
    if package in ordered or package in active:
        return
    active.add(package)
    for item in dependencies(package, index):
        resolve(item, index, ordered, active)
    active.remove(package)
    ordered.append(package)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--library", required=True)
    parser.add_argument("--cache", required=True)
    args = parser.parse_args()
    library = Path(args.library).resolve()
    cache = Path(args.cache).resolve()
    library.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    check = f'.libPaths(c("{library.as_posix()}",.libPaths()));quit(status=if(requireNamespace("ggplot2",quietly=TRUE)) 0 else 1)'
    if subprocess.run([args.rscript, "-e", check]).returncode == 0:
        print("reused installed R packages")
        return
    response = requests.get(f"{INDEX}/PACKAGES.gz", timeout=120)
    response.raise_for_status()
    index = records(gzip.decompress(response.content).decode())
    ordered = []
    resolve("ggplot2", index, ordered, set())
    archives = []
    for package in ordered:
        version = index[package]["Version"]
        target = cache / f"{package}_{version}.zip"
        if not target.exists():
            item = requests.get(f"{INDEX}/{target.name}", timeout=120)
            item.raise_for_status()
            target.write_bytes(item.content)
        archives.append(target)
    paths = ",".join(f'"{path.as_posix()}"' for path in archives)
    expression = f'.libPaths(c("{library.as_posix()}",.libPaths()));install.packages(c({paths}),repos=NULL,type="win.binary",lib=.libPaths()[1])'
    subprocess.run([args.rscript, "-e", expression], check=True)
    print(f"installed {len(archives)} R packages")


if __name__ == "__main__":
    main()
