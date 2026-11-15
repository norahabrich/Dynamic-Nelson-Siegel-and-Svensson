"""Pseudo out-of-sample backtest with an expanding window, and the Diebold-Mariano test.

At each forecast origin t, every model is estimated on data up to t only, and
forecasts the whole curve h months ahead. Errors are realised minus forecast.
"""

import numpy as np
from scipy import stats

from .two_step import fit_ar1


def random_walk(history: np.ndarray, h: int) -> np.ndarray:
    """No-change forecast: the curve in h months equals today's curve."""
    return history[-1]


def ar1_yields(history: np.ndarray, h: int) -> np.ndarray:
    """Direct h-step AR(1) on each yield separately."""
    c, phi = fit_ar1(history, h)
    return c + phi * history[-1]


def expanding_window_errors(yields: np.ndarray, forecasters: dict, horizons, first_origin: int,
                            refit_every: dict | None = None) -> dict:
    """Forecast errors for each forecaster and horizon: errors[name][h] is (origins x maturities).

    `forecasters` maps a name to a function (history, h) -> forecast curve, or to an
    object with `.fit(history, warm_start)`, `.update(history)` and `.forecast(h)`.
    For the latter, `refit_every[name]` sets how many months pass between
    re-estimations of the parameters (factors are always filtered up to the origin).
    The same origins are used for every horizon, so that horizons are comparable.
    """
    yields = np.asarray(yields, dtype=float)
    horizons = list(horizons)
    origins = range(first_origin, len(yields) - max(horizons))
    errors = {name: {h: [] for h in horizons} for name in forecasters}
    for name, f in forecasters.items():
        last_fit = None
        step = (refit_every or {}).get(name, 1)
        for t in origins:
            history = yields[: t + 1]
            if not callable(f):
                if last_fit is None or t - last_fit >= step:
                    f.fit(history, warm_start=last_fit is not None)
                    last_fit = t
                else:
                    f.update(history)
            for h in horizons:
                fc = f(history, h) if callable(f) else f.forecast(h)
                errors[name][h].append(yields[t + h] - fc)
        errors[name] = {h: np.array(e) for h, e in errors[name].items()}
    return errors


def rmse(errors: np.ndarray) -> np.ndarray:
    return np.sqrt(np.mean(errors ** 2, axis=0))


def diebold_mariano(e1: np.ndarray, e2: np.ndarray, h: int) -> tuple[float, float]:
    """Diebold-Mariano (1995) test of equal squared-error loss, with the Harvey,
    Leybourne & Newbold (1997) small-sample correction and h-1 Newey-West lags.

    Returns (statistic, two-sided p-value). A negative statistic means model 1 is
    more accurate than model 2.
    """
    d = np.asarray(e1) ** 2 - np.asarray(e2) ** 2
    n = len(d)
    dbar = d.mean()
    dc = d - dbar
    var = dc @ dc / n
    for k in range(1, h):
        var += 2 * (1 - k / h) * (dc[k:] @ dc[:-k]) / n
    if var <= 0:
        return np.nan, np.nan
    dm = dbar / np.sqrt(var / n)
    dm *= np.sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
    return dm, 2 * stats.t.sf(abs(dm), df=n - 1)
