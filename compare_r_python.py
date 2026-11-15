"""Compare the R results (results/r/, from run_all.R) with the Python ones (results/, from run_all.py).

    python compare_r_python.py

Writes results/r_vs_python.md (the table in the README) and figures/r_vs_python.png.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PY, R = ROOT / "results", ROOT / "results" / "r"
ACCENT, GREY, INK2 = "#2a78d6", "#8a8984", "#52514e"


def maxdiff(a, b):
    return float(np.max(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def fmt(x):
    if x == 0:
        return "0"
    return f"{x:.1e}" if abs(x) < 0.01 else f"{x:.2f}"


def main():
    rows = []

    def add(what, py, r, diff, unit=""):
        rows.append((what, py, r, diff, unit))

    # Two-step factors and in-sample fit
    f_py = pd.read_csv(PY / "two_step_factors.csv", index_col=0)
    f_r = pd.read_csv(R / "two_step_factors.csv", index_col=0)
    add("Two-step factors, 348 months × 3", "", "", maxdiff(f_py, f_r), "%")
    fit_py = pd.read_csv(PY / "in_sample_rmse_bp.csv", index_col=0)
    fit_r = pd.read_csv(R / "in_sample_rmse_bp.csv", index_col=0)
    for col in ["Nelson-Siegel (two-step)", "Svensson (two-step)", "Nelson-Siegel (Kalman)"]:
        p, q = np.sqrt((fit_py[col] ** 2).mean()), np.sqrt((fit_r[col] ** 2).mean())
        add(f"In-sample RMSE, {col}", f"{p:.3f}", f"{q:.3f}", maxdiff(fit_py[col], fit_r[col]), "bp")

    # Kalman estimates
    k_py, k_r = json.load(open(PY / "kalman_params.json")), json.load(open(R / "kalman_params.json"))
    add("Kalman λ", f"{k_py['lambda']:.5f}", f"{k_r['lambda']:.5f}", abs(k_py["lambda"] - k_r["lambda"]))
    add("Kalman log-likelihood", f"{k_py['loglik']:.3f}", f"{k_r['loglik']:.3f}", abs(k_py["loglik"] - k_r["loglik"]))
    add("Kalman transition matrix A (9 values)", "", "", maxdiff(k_py["A"], k_r["A"]))
    add("Kalman factor means μ", "", "", maxdiff(k_py["mu"], k_r["mu"]), "%")
    add("Kalman measurement error sd (17 values)", "", "", maxdiff(k_py["measurement_sd_bp"], k_r["measurement_sd_bp"]), "bp")
    kf_py = pd.read_csv(PY / "kalman_factors.csv", index_col=0)
    kf_r = pd.read_csv(R / "kalman_factors.csv", index_col=0)
    add("Kalman filtered factors, 348 months × 3", "", "", maxdiff(kf_py, kf_r), "%")

    # Out-of-sample RMSE
    o_py = pd.read_csv(PY / "oos_rmse_bp.csv", index_col=[0, 1])
    o_r = pd.read_csv(R / "oos_rmse_bp.csv", index_col=[0, 1])
    o_r.columns = o_py.columns
    o_r = o_r.loc[o_py.index]
    kal = o_py.index.get_level_values(0).str.contains("Kalman")
    add("Out-of-sample RMSE, 5 two-step models and benchmarks (255 values)", "", "",
        maxdiff(o_py[~kal], o_r[~kal]), "bp")
    if kal.any():
        add("Out-of-sample RMSE, Kalman model (51 values)", "", "", maxdiff(o_py[kal], o_r[kal]), "bp")

    # Significance decisions
    d_py = pd.read_csv(PY / "dm_tests_vs_random_walk.csv")
    d_r = pd.read_csv(R / "dm_tests_vs_random_walk.csv")
    m = d_py.merge(d_r, on=["model", "horizon", "maturity_months"], suffixes=("_py", "_r"))
    flips = int(((m.p_value_py < 0.05) != (m.p_value_r < 0.05)).sum())
    add("Diebold-Mariano tests with a different 5% verdict", "", "", flips, f"of {len(m)}")

    # Runtime
    t_py, t_r = json.load(open(PY / "runtime.json")), json.load(open(R / "runtime.json"))
    add("Runtime, one Kalman estimation", f"{t_py['kalman_fit_seconds']:.0f} s", f"{t_r['kalman_fit_seconds']:.0f} s", None)
    add("Runtime, full pipeline", f"{t_py['total_seconds'] / 60:.1f} min", f"{t_r['total_seconds'] / 60:.1f} min", None)

    lines = [f"Python {t_py['python_version']} (NumPy, SciPy) vs {t_r['r_version']} (base R only), same data and settings.", "",
             "| Quantity | Python | R | Largest absolute difference |", "|---|---:|---:|---:|"]
    for what, p, r, diff, unit in rows:
        d = "" if diff is None else (f"{diff}" if isinstance(diff, int) else fmt(diff))
        lines.append(f"| {what} | {p} | {r} | {d} {unit if diff is not None else ''}".rstrip() + " |")
    (PY / "r_vs_python.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    # Figure: every out-of-sample RMSE, R against Python
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    lo, hi = 10, max(o_py.values.max(), o_r.values.max()) * 1.05
    ax.plot([lo, hi], [lo, hi], color=GREY, lw=1, ls="--")
    ax.scatter(o_py[~kal].values.ravel(), o_r[~kal].values.ravel(), s=14, color=GREY, alpha=0.7,
               label="Two-step models and benchmarks")
    if kal.any():
        ax.scatter(o_py[kal].values.ravel(), o_r[kal].values.ravel(), s=18, color=ACCENT, label="Kalman model")
    ax.set_xlabel("Python, out-of-sample RMSE (bp)")
    ax.set_ylabel("R, out-of-sample RMSE (bp)")
    ax.set_title("R and Python give the same forecast errors", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.grid(color="#ecebe8")
    fig.savefig(ROOT / "figures" / "r_vs_python.png", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
