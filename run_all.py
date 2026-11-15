"""Reproduce every result and figure of the README.

    python run_all.py           # full run, about 10 minutes (the Kalman backtest dominates)
    python run_all.py --fast    # skip the Kalman out-of-sample backtest (about 30 seconds)

Outputs go to results/ (CSV, JSON, Markdown tables) and figures/ (PNG).
"""

import argparse
import json
import platform
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dnsmodel import (MATURITIES, KalmanDNS, TwoStepDNS, fit_factors, load_yields,
                      nelson_siegel_loadings, svensson_loadings)
from dnsmodel.backtest import ar1_yields, diebold_mariano, expanding_window_errors, random_walk, rmse

ROOT = Path(__file__).resolve().parent
FIG, RES = ROOT / "figures", ROOT / "results"
LAMBDA = 0.0609            # Diebold & Li (2006)
LAMBDA2 = 2 * LAMBDA       # second Svensson decay, as in the original R code
HORIZONS = [1, 6, 12]
FIRST_ORIGIN = "1994-01-31"  # out-of-sample period starts here, as in Diebold & Li (2006)
SHOW = [3, 12, 36, 60, 120]  # maturities shown in the README tables

# Figure style: one accent colour for the main result, greys for the rest.
ACCENT, ACCENT_LIGHT = "#2a78d6", "#9cc2ee"
GREY, GREY_LIGHT, INK, INK2 = "#8a8984", "#c9c8c3", "#0b0b0b", "#52514e"
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": GREY_LIGHT,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": "#ecebe8", "grid.linewidth": 0.8, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "legend.frameon": False, "lines.linewidth": 2,
})


def mat_label(m):
    return f"{m}m" if m < 12 else f"{m // 12}y" if m % 12 == 0 else f"{m}m"


def in_sample(y: pd.DataFrame):
    Y = y.to_numpy()
    out = {}
    for name, Z in [("Nelson-Siegel", nelson_siegel_loadings(LAMBDA, MATURITIES)),
                    ("Svensson", svensson_loadings(LAMBDA, LAMBDA2, MATURITIES))]:
        betas, fitted = fit_factors(Y, Z)
        out[name] = {"betas": betas, "fitted": fitted, "Z": Z}
    return out


def fig_curves(y, fits):
    dates = ["1974-08-31", "1989-04-30", "1992-12-31", "1998-11-30"]
    grid = np.linspace(1, 120, 240)
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.2), sharey=False)
    for ax, d in zip(axes, dates):
        i = y.index.get_loc(pd.Timestamp(d))
        b = fits["Nelson-Siegel"]["betas"][i]
        ax.plot(grid, nelson_siegel_loadings(LAMBDA, grid) @ b, color=ACCENT, label="Nelson-Siegel fit")
        ax.scatter(MATURITIES, y.iloc[i], s=18, color=INK, zorder=3, label="Observed yields")
        ax.set_title(pd.Timestamp(d).strftime("%b %Y"), fontweight="normal")
        ax.set_xlabel("Maturity (months)")
    axes[0].set_ylabel("Yield (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=2, fontsize=8)
    fig.suptitle("Three factors reproduce very different curve shapes: inverted, humped, steep, flat",
                 x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "curve_fits.png")
    plt.close(fig)


def fig_factors(y, betas):
    Y = y
    level_p = Y[120]
    slope_p = Y[3] - Y[120]
    curv_p = 2 * Y[24] - Y[3] - Y[120]
    panels = [("Level  β1", betas[:, 0], level_p, "10y yield"),
              ("Slope  β2", betas[:, 1], slope_p, "3m − 10y spread"),
              ("Curvature  β3", betas[:, 2], curv_p, "2 × 2y − 3m − 10y")]
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    for ax, (title, b, proxy, plabel) in zip(axes, panels):
        ax.plot(y.index, proxy, color=GREY_LIGHT, lw=1.5, label=f"Data proxy: {plabel}")
        ax.plot(y.index, b, color=ACCENT, lw=1.5, label="Estimated factor")
        corr = np.corrcoef(b, proxy)[0, 1]
        ax.set_title(f"{title}   (correlation with proxy: {corr:.2f})", fontweight="normal", fontsize=10)
        ax.legend(loc="upper right", fontsize=8, ncol=2)
        ax.set_ylabel("%")
    fig.suptitle("The estimated factors are the level, slope and curvature of the curve",
                 x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "factors.png")
    plt.close(fig)


def fig_lambda(y, kalman_lambda):
    Y = y.to_numpy()
    grid = np.linspace(0.02, 0.2, 91)
    r = [np.sqrt(np.mean((Y - fit_factors(Y, nelson_siegel_loadings(l, MATURITIES))[1]) ** 2)) * 100 for l in grid]
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(grid, r, color=GREY)
    for lam, lab, c, off in [(LAMBDA, f"Diebold-Li λ = {LAMBDA}", ACCENT, (-10, 12))] + (
            [(kalman_lambda, f"Kalman MLE λ = {kalman_lambda:.4f}", INK, (8, -22))] if kalman_lambda else []):
        v = np.sqrt(np.mean((Y - fit_factors(Y, nelson_siegel_loadings(lam, MATURITIES))[1]) ** 2)) * 100
        ax.scatter([lam], [v], color=c, zorder=3)
        ax.annotate(lab, (lam, v), xytext=off, textcoords="offset points", color=c, fontsize=9,
                    ha="right" if off[0] < 0 else "left")
    ax.set_xlabel("λ (per month)")
    ax.set_ylabel("In-sample RMSE (bp)")
    ax.set_title("Fit quality is flat around the conventional λ")
    fig.savefig(FIG / "lambda_choice.png")
    plt.close(fig)


def fig_ratio(rm):
    fig, ax = plt.subplots(figsize=(8, 4))
    styles = {1: (ACCENT_LIGHT, 1.5), 6: (ACCENT, 2.5), 12: ("#1b4f8f", 1.5)}
    for h in HORIZONS:
        ratio = rm.loc[("DNS two-step (AR1)", h)] / rm.loc[("Random walk", h)]
        c, lw = styles[h]
        ax.plot(MATURITIES, ratio.values, marker="o", ms=4, color=c, lw=lw, label=f"{h}-month horizon")
    ax.axhline(1, color=GREY, ls="--", lw=1.2)
    ax.text(121, 1.005, "random walk", color=INK2, fontsize=8, ha="right", va="bottom")
    ax.set_xlabel("Maturity (months)")
    ax.set_ylabel("RMSE ratio, DNS / random walk")
    r6 = rm.loc[("DNS two-step (AR1)", 6)] / rm.loc[("Random walk", 6)]
    ax.set_title(f"DNS beats the random walk at short maturities "
                 f"(−{(1 - r6.iloc[0]) * 100:.0f}% error on the 3m yield at 6 months)")
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(FIG / "oos_rmse_ratio.png")
    plt.close(fig)


def main(fast: bool):
    global T_START
    T_START = time.perf_counter()
    FIG.mkdir(exist_ok=True)
    RES.mkdir(exist_ok=True)
    y = load_yields()
    Y = y.to_numpy()

    # ---- In-sample fit -------------------------------------------------------
    fits = in_sample(y)
    t0 = time.perf_counter()
    kalman = KalmanDNS(MATURITIES).fit(Y)
    kalman_seconds = time.perf_counter() - t0
    kp = kalman.params
    kz = nelson_siegel_loadings(kp["lambda"], MATURITIES)
    k_states = kalman.filtered_factors()
    k_fitted = k_states @ kz.T
    cols = ["level", "slope", "curvature"]
    pd.DataFrame(fits["Nelson-Siegel"]["betas"], index=y.index.strftime("%Y-%m-%d"), columns=cols).rename_axis("date") \
        .to_csv(RES / "two_step_factors.csv")
    pd.DataFrame(k_states, index=y.index.strftime("%Y-%m-%d"), columns=cols).rename_axis("date") \
        .to_csv(RES / "kalman_factors.csv")
    fit_tbl = pd.DataFrame({
        "Nelson-Siegel (two-step)": rmse(Y - fits["Nelson-Siegel"]["fitted"]) * 100,
        "Svensson (two-step)": rmse(Y - fits["Svensson"]["fitted"]) * 100,
        "Nelson-Siegel (Kalman)": rmse(Y - k_fitted) * 100,
    }, index=pd.Index(MATURITIES.astype(int), name="maturity_months"))
    fit_tbl.round(2).to_csv(RES / "in_sample_rmse_bp.csv")
    json.dump({"lambda": kp["lambda"], "mu": kp["mu"].tolist(), "A": kp["A"].tolist(),
               "Q": kp["Q"].round(5).tolist(), "measurement_sd_bp": (kp["meas_sd"] * 100).tolist(),
               "loglik": kalman.loglik_}, open(RES / "kalman_params.json", "w"), indent=2)

    fig_curves(y, fits)
    fig_factors(y, fits["Nelson-Siegel"]["betas"])
    fig_lambda(y, kp["lambda"])

    # ---- Out-of-sample backtest ---------------------------------------------
    first = y.index.get_loc(pd.Timestamp(FIRST_ORIGIN))
    Zns = nelson_siegel_loadings(LAMBDA, MATURITIES)
    Zsv = svensson_loadings(LAMBDA, LAMBDA2, MATURITIES)
    models = {
        "Random walk": random_walk,
        "AR(1) on each yield": ar1_yields,
        "DNS two-step (AR1)": TwoStepDNS(Zns, "ar1").forecast,
        "DNS two-step (VAR1)": TwoStepDNS(Zns, "var1").forecast,
        "Svensson two-step (AR1)": TwoStepDNS(Zsv, "ar1").forecast,
    }
    refit = {}
    if not fast:
        models["DNS Kalman (VAR1)"] = KalmanDNS(MATURITIES)
        refit["DNS Kalman (VAR1)"] = 12  # parameters re-estimated every 12 months
    t1 = time.perf_counter()
    errors = expanding_window_errors(Y, models, HORIZONS, first, refit)
    backtest_seconds = time.perf_counter() - t1
    n_origins = len(next(iter(errors.values()))[1])
    rows, dm_rows = [], []
    for name, eh in errors.items():
        for h, e in eh.items():
            rows.append(pd.Series(rmse(e) * 100, index=MATURITIES.astype(int), name=(name, h)))
            if name != "Random walk":
                for j, m in enumerate(MATURITIES.astype(int)):
                    stat, p = diebold_mariano(e[:, j], errors["Random walk"][h][:, j], h)
                    dm_rows.append({"model": name, "horizon": h, "maturity_months": m, "dm_stat": stat, "p_value": p})
    rm = pd.DataFrame(rows)
    rm.index = pd.MultiIndex.from_tuples(rm.index, names=["model", "horizon"])
    rm.round(2).to_csv(RES / "oos_rmse_bp.csv")
    dm = pd.DataFrame(dm_rows)
    dm.round(4).to_csv(RES / "dm_tests_vs_random_walk.csv", index=False)
    fig_ratio(rm)

    # ---- Markdown summary for the README ------------------------------------
    period = f"{y.index[first]:%b %Y} – {y.index[first + n_origins - 1]:%b %Y}"
    lines = [f"Out-of-sample RMSE in basis points, forecast origins {period} ({n_origins} months), "
             f"expanding window. **Bold** = best model; † = significantly different from the random walk "
             f"(Diebold-Mariano, 5%).", ""]
    for h in HORIZONS:
        lines += [f"**{h}-month horizon**", "",
                  "| Model | " + " | ".join(mat_label(m) for m in SHOW) + " |",
                  "|---|" + "---:|" * len(SHOW)]
        sub = rm.xs(h, level="horizon")[SHOW]
        for name, r in sub.iterrows():
            cells = []
            for m in SHOW:
                v = f"{r[m]:.1f}"
                if name != "Random walk":
                    p = dm.query("model == @name and horizon == @h and maturity_months == @m").p_value.iloc[0]
                    v += "†" if p < 0.05 else ""
                cells.append(f"**{v}**" if r[m] == sub[m].min() else v)
            lines.append(f"| {name} | " + " | ".join(cells) + " |")
        lines.append("")
    lines += ["In-sample RMSE in basis points, 1972–2000:", "",
              "| Model | " + " | ".join(mat_label(m) for m in SHOW) + " | All |", "|---|" + "---:|" * (len(SHOW) + 1)]
    for col in fit_tbl.columns:
        allm = np.sqrt((fit_tbl[col] ** 2).mean())
        lines.append(f"| {col} | " + " | ".join(f"{fit_tbl.loc[m, col]:.1f}" for m in SHOW) + f" | {allm:.1f} |")
    json.dump({"total_seconds": round(time.perf_counter() - T_START, 1), "kalman_fit_seconds": round(kalman_seconds, 1),
               "backtest_seconds": round(backtest_seconds, 1), "fast": fast, "python_version": platform.python_version()},
              open(RES / "runtime.json", "w"))
    (RES / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nKalman lambda = {kp['lambda']:.4f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="skip the Kalman out-of-sample backtest")
    main(ap.parse_args().fast)
