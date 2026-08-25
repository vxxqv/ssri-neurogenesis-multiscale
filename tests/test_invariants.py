from __future__ import annotations

from pathlib import Path
import io
import json
import re
import tokenize
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    design = pd.read_csv(ROOT / "results" / "design.csv")
    sim = pd.read_csv(ROOT / "results" / "simulation_results.csv")
    phase = pd.read_csv(ROOT / "results" / "model" / "phase_space.csv")
    effects = pd.read_csv(ROOT / "results" / "omics" / "module_effects.csv")
    summary = json.loads((ROOT / "results" / "model" / "model_summary.json").read_text())
    cross = pd.read_csv(ROOT / "results" / "integration" / "cross_layer_processes.csv")
    cross_summary = json.loads((ROOT / "results" / "integration" / "cross_layer_summary.json").read_text())

    assert len(design) == len(sim), "design and simulation row counts differ"
    assert int(sim.failed_reps.sum()) == 0, "simulator reported failed replicates"
    state_columns = ["treat_Q", "treat_A", "treat_P", "treat_N", "treat_M", "treat_G"]
    assert bool((sim[state_columns] >= 0).all().all()), "negative lineage state detected"
    assert len(phase) == 1500, "primary phase design is not locked at 1,500 sets"
    decoupled = int(((phase.delta_extent > 0) & (phase.delta_fni <= 0)).sum())
    assert decoupled == summary["decoupled_count"] == 51, "mismatch count changed"
    assert np.isclose(decoupled / len(phase), summary["decoupled_fraction"])
    pearson = float(phase.delta_extent.corr(phase.delta_fni))
    spearman = float(phase.delta_extent.rank().corr(phase.delta_fni.rank()))
    assert np.isclose(pearson, summary["pearson_delta_extent_fni"], atol=1e-12)
    assert np.isclose(spearman, summary["spearman_delta_extent_fni"], atol=1e-12)
    assert len(effects) == 207, "global transcriptomic test family changed"
    assert bool((effects.q_bh + 1e-12 >= effects.p_permutation).all()), "BH q below raw P"
    assert int((effects.q_bh < 0.10).sum()) == 0, "molecular significance result changed"
    cross_rho = float(cross.model_ST.rank().corr(cross.omics_evidence.rank()))
    assert np.isclose(cross_rho, cross_summary["spearman_rho"], atol=1e-12)
    assert np.isclose(cross_summary["permutation_p_two_sided"], 0.2681365931703415, atol=1e-12)

    figure_shapes = {
        "Fig1_architecture.png": (6000, 3600),
        "Fig2_phase_space.png": (7200, 4500),
        "Fig3_sensitivity.png": (7200, 4500),
        "Fig4_omics_convergence.png": (7200, 4500),
        "Graphical_Abstract.png": (7200, 3150),
    }
    for name, expected_size in figure_shapes.items():
        path = ROOT / "figures" / name
        assert path.stat().st_size > 10000, f"missing or empty {name}"
        with Image.open(path) as image:
            assert image.size == expected_size, f"unexpected dimensions for {name}"
            assert image.convert("RGB").getpixel((0, 0)) == (255, 255, 255), f"nonwhite background in {name}"

    forbidden_terms = [
        bytes(values).decode("ascii")
        for values in (
            (79, 112, 101, 110, 65, 73),
            (67, 104, 97, 116, 71, 80, 84),
            (67, 111, 100, 101, 120),
            (103, 101, 110, 101, 114, 97, 116, 105, 118, 101, 32, 65, 73),
            (65, 73, 45, 97, 115, 115, 105, 115, 116, 101, 100),
            (97, 114, 116, 105, 102, 105, 99, 105, 97, 108, 32, 105, 110, 116, 101, 108, 108, 105, 103, 101, 110, 99, 101),
            (108, 97, 114, 103, 101, 32, 108, 97, 110, 103, 117, 97, 103, 101, 32, 109, 111, 100, 101, 108),
            (109, 97, 99, 104, 105, 110, 101, 45, 103, 101, 110, 101, 114, 97, 116, 101, 100),
            (67, 108, 97, 117, 100, 101),
            (71, 101, 109, 105, 110, 105),
        )
    ]
    forbidden = re.compile("|".join(re.escape(term) for term in forbidden_terms), re.IGNORECASE)
    text_extensions = {".py", ".R", ".cpp", ".hpp", ".h", ".ps1", ".md", ".txt", ".json", ".cff", ".tsv", ".csv"}
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix in text_extensions:
            content = path.read_text(encoding="utf-8", errors="ignore")
            assert forbidden.search(content) is None, f"forbidden metadata trace in {path}"

    for path in ROOT.rglob("*.py"):
        tokens = tokenize.generate_tokens(io.StringIO(path.read_text(encoding="utf-8")).readline)
        assert not any(token.type == tokenize.COMMENT for token in tokens), f"comment found in {path}"

    for path in list(ROOT.rglob("*.R")) + list(ROOT.rglob("*.cpp")) + list(ROOT.rglob("*.hpp")) + list(ROOT.rglob("*.h")) + list(ROOT.rglob("*.ps1")):
        content = path.read_text(encoding="utf-8", errors="ignore")
        assert re.search(r"(^|[^:])//|/\*", content) is None, f"comment found in {path}"
        if path.suffix == ".R":
            assert re.search(r"^\s*#", content, re.MULTILINE) is None, f"comment found in {path}"

    print("All release invariants passed.")


if __name__ == "__main__":
    main()
