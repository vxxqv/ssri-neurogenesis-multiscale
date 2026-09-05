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
    model = ROOT / "results" / "model"
    transcript = ROOT / "results" / "transcriptomics"
    phase = pd.read_csv(model / "phase_space.csv")
    summary = json.loads((model / "model_summary.json").read_text(encoding="utf-8"))
    prediction = pd.read_csv(model / "predictive_information.csv").set_index("model")
    temporal = pd.read_csv(model / "temporal_summary.csv")
    numerical = json.loads((model / "numerical_validation.json").read_text(encoding="utf-8"))
    programs = pd.read_csv(transcript / "program_effects.csv")
    qc = pd.read_csv(transcript / "dataset_qc.csv")
    robustness = pd.read_csv(transcript / "relative_score_robustness.csv")
    splits = pd.read_csv(transcript / "relative_score_gene_splits.csv")

    assert len(phase) == 2048
    assert phase.n_rep.eq(16).all()
    assert np.isclose(phase.delta_extent.corr(phase.delta_fni), summary["pearson_extent_fni"], atol=1e-12)
    assert np.isclose(phase.delta_extent.rank().corr(phase.delta_fni.rank()), summary["spearman_extent_fni"], atol=1e-12)
    assert summary["robust_mismatch_count"] == 1
    assert np.isclose(summary["replicate_mismatch_fraction"], 0.040130615234375)
    assert prediction.loc["process_aware", "r2"] > prediction.loc["cell_composition", "r2"] > prediction.loc["extent_only", "r2"]
    assert prediction.loc["process_aware", "r2"] > 0.85
    assert temporal.loc[temporal.mean_mismatch_probability.idxmax(), "t_end"] == 14
    assert numerical["rk4_max_relative_error"] < 2e-7
    assert numerical["null_max_absolute_extent"] == 0
    assert numerical["null_max_absolute_fni"] == 0
    assert numerical["shapley_max_efficiency_error"] < 1e-10

    assert set(qc.dataset) == {"GSE197622", "GSE43261", "GSE309750", "GSE222756", "GSE205325", "GSE292948"}
    assert len(programs) == 966
    assert programs[["hedges_g", "p_value", "q_global", "q_dataset"]].notna().all().all()
    assert programs.p_value.between(0, 1).all()
    assert programs.q_global.between(0, 1).all()
    assert programs.q_dataset.between(0, 1).all()
    mossy = programs[(programs.dataset == "GSE309750") & (programs.contrast == "mossy_cells_Fluoxetine_vs_Vehicle")].set_index("program")
    assert mossy.loc["integration_bias", "hedges_g"] > 3.7
    assert mossy.loc["exclusive_integration_bias", "hedges_g"] > 3.7
    assert mossy.loc["integration_bias", "q_dataset"] < 0.01
    response = programs[(programs.dataset == "GSE43261") & (programs.contrast == "dorsal_Responder_vs_Resistant")].set_index("program")
    assert response.loc["exclusive_integration_bias", "hedges_g"] > 1.8
    wide = programs[programs.program.isin(["integration_bias", "exclusive_integration_bias"])].pivot_table(index=["dataset", "contrast"], columns="program", values="hedges_g").dropna()
    assert wide.integration_bias.rank().corr(wide.exclusive_integration_bias.rank()) > 0.98
    mossy_robustness = robustness[(robustness.dataset == "GSE309750") & (robustness.contrast == "MC")]
    assert len(mossy_robustness) == 3
    assert mossy_robustness.hedges_g.min() > 2.2
    assert mossy_robustness.loo_min_g.min() > 2.0
    mossy_splits = splits[(splits.dataset == "GSE309750") & (splits.contrast == "MC")]
    assert len(mossy_splits) == 200
    assert ((mossy_splits.first_half_g > 0) & (mossy_splits.second_half_g > 0)).all()

    names = ["Fig1_study_architecture", "Fig2_model_results", "Fig3_mechanism", "Fig4_transcriptomic_evidence"]
    for name in names:
        for suffix in (".png", ".tiff", ".pdf"):
            path = ROOT / "figures" / f"{name}{suffix}"
            assert path.stat().st_size > 5000
        with Image.open(ROOT / "figures" / f"{name}.png") as image:
            assert min(image.size) > 4000
            assert image.convert("RGB").getpixel((0, 0)) == (255, 255, 255)

    for path in ROOT.rglob("*.py"):
        if any(part in {".venv", ".r-lib", "external", "work"} for part in path.parts):
            continue
        tokens = tokenize.generate_tokens(io.StringIO(path.read_text(encoding="utf-8")).readline)
        assert not any(token.type == tokenize.COMMENT for token in tokens)
    for path in list(ROOT.rglob("*.R")) + list(ROOT.rglob("*.cpp")) + list(ROOT.rglob("*.ps1")):
        if any(part in {".venv", ".r-lib", "external", "work"} for part in path.parts):
            continue
        content = path.read_text(encoding="utf-8", errors="ignore")
        assert re.search(r"(^|[^:])//|/\*", content) is None
        if path.suffix == ".R":
            assert re.search(r"^\s*#", content, re.MULTILINE) is None

    print("All release invariants passed.")


if __name__ == "__main__":
    main()
