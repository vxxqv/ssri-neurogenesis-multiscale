from __future__ import annotations

from pathlib import Path
import argparse
import concurrent.futures
import shutil
import subprocess
import pandas as pd


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def combine(paths: list[Path], target: Path) -> None:
    with target.open("wb") as output:
        for index, path in enumerate(paths):
            with path.open("rb") as source:
                if index:
                    source.readline()
                shutil.copyfileobj(source, output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulator", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--aggregate", required=True)
    parser.add_argument("--replicates", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    for path in workdir.glob("shard_*.*"):
        path.unlink()
    design = pd.read_csv(args.design, dtype={"model": str}, low_memory=False)
    design = design.loc[design.reps > 0].reset_index(drop=True)
    workers = min(args.workers, len(design))
    commands, aggregate_paths, replicate_paths = [], [], []
    for worker in range(workers):
        shard = design.iloc[worker::workers]
        design_path = workdir / f"shard_{worker:02d}.design.csv"
        aggregate_path = workdir / f"shard_{worker:02d}.aggregate.csv"
        replicate_path = workdir / f"shard_{worker:02d}.replicates.csv"
        shard.to_csv(design_path, index=False)
        commands.append([args.simulator, str(design_path), str(aggregate_path), str(replicate_path)])
        aggregate_paths.append(aggregate_path)
        replicate_paths.append(replicate_path)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(run, commands))
    combine(aggregate_paths, Path(args.aggregate))
    combine(replicate_paths, Path(args.replicates))
    print(f"simulated {len(design)} stochastic systems with {workers} workers")


if __name__ == "__main__":
    main()
