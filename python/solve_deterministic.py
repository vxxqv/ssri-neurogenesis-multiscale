from __future__ import annotations

from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def derivative(time: float, state: np.ndarray, frame: pd.DataFrame, treated: bool) -> np.ndarray:
    q, a, p, n, m, g, fni = state.T
    tau = frame.tau.to_numpy()

    def multiplier(name: str, switch: str) -> np.ndarray:
        if not treated:
            return np.ones(len(frame))
        target = frame[name].to_numpy()
        active = frame[switch].to_numpy(bool)
        value = np.exp(np.log(np.maximum(target, 1e-9)) * (1.0 - np.exp(-time / np.maximum(tau, 1e-6))))
        return np.where(active, value, 1.0)

    ma = multiplier("tx_activation", "activation_on")
    mp = multiplier("tx_prolif", "proliferation_on")
    mm = multiplier("tx_maturation", "maturation_on")
    mi = multiplier("tx_integration", "integration_on")
    me = multiplier("tx_eff", "efficacy_on")
    ms = multiplier("tx_survival", "survival_on")
    qa = frame.k_qa.to_numpy() * ma * q
    aq = frame.k_aq.to_numpy() * a
    ap = frame.b_a.to_numpy() * mp * a
    pp = frame.b_p.to_numpy() * mp * p
    pn = frame.k_pn.to_numpy() * p
    nm = frame.k_nm.to_numpy() * mm * n
    mg = frame.k_mg.to_numpy() * mi * m
    return np.column_stack([
        -qa + aq,
        qa - aq - frame.d_a.to_numpy() * a,
        ap + pp - pn - frame.d_p.to_numpy() * p,
        pn - nm - frame.d_n.to_numpy() * ms * n,
        nm - mg - frame.d_m.to_numpy() * ms * m,
        mg - frame.d_g.to_numpy() * ms * g,
        mg * frame.eff_mean.to_numpy() * me - frame.d_g.to_numpy() * ms * fni,
    ])


def solve(frame: pd.DataFrame, treated: bool, step: float = 0.25) -> np.ndarray:
    state = np.column_stack([
        frame.Q0, frame.A0, frame.P0, frame.N0, frame.M0, frame.G0,
        frame.G0.to_numpy() * frame.eff_mean.to_numpy(),
    ]).astype(float)
    time = 0.0
    end = float(frame.t_end.iloc[0])
    while time < end - 1e-12:
        dt = min(step, end - time)
        k1 = derivative(time, state, frame, treated)
        k2 = derivative(time + dt / 2.0, state + dt * k1 / 2.0, frame, treated)
        k3 = derivative(time + dt / 2.0, state + dt * k2 / 2.0, frame, treated)
        k4 = derivative(time + dt, state + dt * k3, frame, treated)
        state += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        time += dt
    return state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--batch-size", type=int, default=4096)
    args = parser.parse_args()
    design = pd.read_csv(args.design, dtype={"model": str}, low_memory=False)
    records = []
    for end, time_group in design.groupby("t_end", sort=True):
        for start in range(0, len(time_group), args.batch_size):
            frame = time_group.iloc[start : start + args.batch_size].reset_index(drop=True)
            control = solve(frame, False)
            treatment = solve(frame, True)
            result = frame[["set_id", "base_id", "design_group", "model", "t_end"]].copy()
            result["control_extent"] = control[:, 4] + control[:, 5]
            result["treat_extent"] = treatment[:, 4] + treatment[:, 5]
            result["delta_extent"] = result.treat_extent - result.control_extent
            result["control_fni"] = control[:, 6]
            result["treat_fni"] = treatment[:, 6]
            result["delta_fni"] = result.treat_fni - result.control_fni
            for index, state in enumerate("QAPNMG"):
                result[f"control_{state}"] = control[:, index]
                result[f"treat_{state}"] = treatment[:, index]
            records.append(result)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(records, ignore_index=True).to_csv(out, index=False)
    print(f"solved {len(design)} deterministic systems")


if __name__ == "__main__":
    main()
