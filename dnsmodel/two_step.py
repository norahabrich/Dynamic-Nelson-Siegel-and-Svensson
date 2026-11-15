"""Two-step Dynamic Nelson-Siegel, as in Diebold & Li (2006).

Step 1: for each month, regress the cross-section of yields on the fixed
loadings (OLS) to get the factors beta_t.
Step 2: model the factors as AR(1) processes (one per factor, as in the paper)
or as a VAR(1), and map factor forecasts back to yields through the loadings.
"""

from dataclasses import dataclass

import numpy as np


def fit_factors(yields: np.ndarray, Z: np.ndarray):
    """OLS of each yield curve on the loadings. Returns (betas T x k, fitted yields T x N)."""
    yields = np.asarray(yields, dtype=float)
    betas = np.linalg.lstsq(Z, yields.T, rcond=None)[0].T
    return betas, betas @ Z.T


def fit_var1(x: np.ndarray, h: int = 1):
    """Direct h-step VAR(1) by OLS: x[t+h] = c + A x[t] + e. Returns (c, A, residual covariance)."""
    x = np.asarray(x, dtype=float)
    X = np.column_stack([np.ones(len(x) - h), x[:-h]])
    coef = np.linalg.lstsq(X, x[h:], rcond=None)[0]
    resid = x[h:] - X @ coef
    return coef[0], coef[1:].T, np.cov(resid, rowvar=False)


def fit_ar1(x: np.ndarray, h: int = 1):
    """Direct h-step AR(1) per column: x[t+h, j] = c_j + phi_j x[t, j]. Returns (c, phi)."""
    x = np.asarray(x, dtype=float)
    c, phi = np.empty(x.shape[1]), np.empty(x.shape[1])
    for j in range(x.shape[1]):
        X = np.column_stack([np.ones(len(x) - h), x[:-h, j]])
        c[j], phi[j] = np.linalg.lstsq(X, x[h:, j], rcond=None)[0]
    return c, phi


@dataclass
class TwoStepDNS:
    """Two-step DNS forecaster. `dynamics` is "ar1" (Diebold-Li) or "var1"."""

    Z: np.ndarray
    dynamics: str = "ar1"

    def forecast(self, yields_history: np.ndarray, h: int) -> np.ndarray:
        """Forecast the yield curve h months after the last row of `yields_history`."""
        betas, _ = fit_factors(yields_history, self.Z)
        if self.dynamics == "ar1":
            c, phi = fit_ar1(betas, h)
            beta_f = c + phi * betas[-1]
        elif self.dynamics == "var1":
            c, A, _ = fit_var1(betas, h)
            beta_f = c + A @ betas[-1]
        else:
            raise ValueError(f"unknown dynamics {self.dynamics!r}")
        return self.Z @ beta_f
