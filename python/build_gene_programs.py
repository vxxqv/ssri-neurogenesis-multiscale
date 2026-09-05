from __future__ import annotations

from pathlib import Path
import argparse
import gzip
import re
from collections import defaultdict
import pandas as pd


def ontology(path: Path) -> dict[str, set[str]]:
    parents = defaultdict(set)
    current = None
    obsolete = False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip()
            if line == "[Term]":
                current = None
                obsolete = False
            elif line.startswith("id: GO:"):
                current = line.split("id: ", 1)[1]
            elif line == "is_obsolete: true":
                obsolete = True
            elif current and not obsolete and line.startswith("is_a: GO:"):
                parents[current].add(line.split()[1])
            elif current and not obsolete and line.startswith("relationship: part_of GO:"):
                parents[current].add(line.split()[2])
    return parents


def ancestors(term: str, parents: dict[str, set[str]], memo: dict[str, set[str]]) -> set[str]:
    if term in memo:
        return memo[term]
    found = {term}
    for parent in parents.get(term, set()):
        found.update(ancestors(parent, parents, memo))
    memo[term] = found
    return found


def gaf_programs(path: Path, targets: dict[str, str], parents: dict[str, set[str]], species: str) -> list[dict[str, str]]:
    records = []
    memo = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 15 or "NOT" in fields[3].split("|") or fields[6] == "ND":
                continue
            symbol, term = fields[2], fields[4]
            term_ancestors = ancestors(term, parents, memo)
            for program, target in targets.items():
                if target in term_ancestors:
                    records.append({"species": species, "program": program, "gene": symbol, "source_term": term})
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--programs", required=True)
    parser.add_argument("--obo", required=True)
    parser.add_argument("--mouse-gaf", required=True)
    parser.add_argument("--rat-gaf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    programs = pd.read_csv(args.programs)
    targets = dict(zip(programs.program, programs.go_id))
    parents = ontology(Path(args.obo))
    records = gaf_programs(Path(args.mouse_gaf), targets, parents, "mouse")
    records.extend(gaf_programs(Path(args.rat_gaf), targets, parents, "rat"))
    frame = pd.DataFrame(records).drop_duplicates(["species", "program", "gene"])
    composites = []
    layers = dict(zip(programs.program, programs.layer))
    for species, group in frame.groupby("species"):
        for layer in ("extent", "integration"):
            genes = group.loc[group.program.map(layers) == layer, "gene"].unique()
            composites.extend({"species": species, "program": f"{layer}_composite", "gene": gene, "source_term": "composite"} for gene in genes)
    frame = pd.concat([frame, pd.DataFrame(composites)], ignore_index=True).drop_duplicates(["species", "program", "gene"])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.sort_values(["species", "program", "gene"]).to_csv(out, sep="\t", index=False)
    print(frame.groupby(["species", "program"]).size().to_string())


if __name__ == "__main__":
    main()
