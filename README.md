# SSRI neurogenesis model

This project asks a simple question: does producing more new hippocampal cells always produce more useful neurogenesis? It combines an exact stochastic lineage model with six public mouse and rat transcriptomic datasets. No new animal or human experiment was performed.

## What the model does

Think of neurogenesis as an assembly line. Cells move from stem cells to integrated neurons. SSRIs can change six parts of that line: activation, proliferation, maturation, survival, integration, and efficacy. The model compares cell number with a functional index that also considers whether integrated cells work effectively.

## Main results

- Cell number predicted some functional change, but a process-aware model predicted much more in held-out simulations: R2 = 0.864 versus 0.470.
- A cell-gain without functional-gain pattern was most common at day 14 and became less common by day 56.
- In independent mossy-cell data, fluoxetine shifted gene activity toward integration rather than cell-production programs. The effect survived three scoring methods, every leave-one-out test, and 200 random gene splits.
- Responder data showed the same direction, but this comparison remains exploratory after the stricter multiple-testing correction.

These are computational and secondary-data results. They do not show clinical benefit or prove that the model is biologically complete.

## Run it

Install Python, R, and a C++17 compiler, then run:

```powershell
python -m pip install -r requirements.txt
.\scripts\run_analysis.ps1
```

The script downloads public data, runs the analyses, rebuilds the figures, and checks the results. See [RESULTS.md](RESULTS.md) for a short explanation and [source_data_manifest.tsv](source_data_manifest.tsv) for data credits.

Code is released under the MIT License. Public source data keep their original terms.
