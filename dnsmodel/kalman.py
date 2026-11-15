"""One-step Dynamic Nelson-Siegel in state-space form, estimated by maximum likelihood.

    Measurement:  y_t = Z(lambda) beta_t + eps_t,          eps_t ~ N(0, H),  H diagonal
    Transition:   beta_t - mu = A (beta_{t-1} - mu) + eta_t, eta_t ~ N(0, Q)

This is the specification of Diebold, Rudebusch & Aruoba (2006). Unlike the
two-step approach, lambda, the factor dynamics and the measurement errors are
estimated jointly, and the factors are filtered rather than fitted month by month.
The two-step estimates are used as starting values, as the original R code
(`legacy/R/DNS-TS.R`) prepared but never used.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from .loadings import nelson_siegel_loadings
from .two_step import fit_factors, fit_var1

K = 3


def _unpack(theta, n):
    i = 0
    lam = np.exp(theta[i]); i += 1
    h = np.exp(2 * theta[i:i + n]); i += n
    A = theta[i:i + K * K].reshape(K, K); i += K * K
    mu = theta[i:i + K]; i += K
    L = np.zeros((K, K))
    L[np.tril_indices(K)] = theta[i:i + K * (K + 1) // 2]
    L[np.diag_indices(K)] = np.exp(L[np.diag_indices(K)])
    return lam, h, A, mu, L @ L.T


def _stationary_cov(A, Q):
    """Solve P = A P A' + Q (discrete Lyapunov equation); None if A is not stable."""
    if np.max(np.abs(np.linalg.eigvals(A))) >= 0.999:
        return None
    P = np.linalg.solve(np.eye(K * K) - np.kron(A, A), Q.reshape(-1)).reshape(K, K)
    return (P + P.T) / 2


def kalman_filter(y, Z, h, A, mu, Q):
    """Run the filter. Returns (log-likelihood, filtered states T x 3)."""
    T, n = y.shape
    P = _stationary_cov(A, Q)
    if P is None:
        return -np.inf, None
    b = mu.copy()
    loglik = 0.0
    states = np.empty((T, K))
    for t in range(T):
        # prediction (the first period starts from the unconditional moments)
        if t > 0:
            b = mu + A @ (b - mu)
            P = A @ P @ A.T + Q
        v = y[t] - Z @ b
        F = Z @ P @ Z.T + np.diag(h)
        try:
            C = np.linalg.cholesky(F)
        except np.linalg.LinAlgError:
            return -np.inf, None
        Finv_v = np.linalg.solve(C.T, np.linalg.solve(C, v))
        Finv_Z = np.linalg.solve(C.T, np.linalg.solve(C, Z))
        loglik -= 0.5 * (n * np.log(2 * np.pi) + 2 * np.log(np.diag(C)).sum() + v @ Finv_v)
        # update
        PZt = P @ Z.T
        b = b + PZt @ Finv_v
        P = P - PZt @ Finv_Z @ P
        states[t] = b
    return loglik, states


@dataclass
class KalmanDNS:
    """DNS estimated jointly by the Kalman filter and maximum likelihood."""

    maturities: np.ndarray
    lam0: float = 0.0609
    maxiter: int = 2000
    theta_: np.ndarray = field(default=None, repr=False)
    loglik_: float = None

    def _start(self, y):
        Z = nelson_siegel_loadings(self.lam0, self.maturities)
        betas, fitted = fit_factors(y, Z)
        c, A, S = fit_var1(betas, 1)
        mu = np.linalg.solve(np.eye(K) - A, c)
        L = np.linalg.cholesky(S)
        L[np.diag_indices(K)] = np.log(np.diag(L))
        sd = np.std(y - fitted, axis=0)
        return np.concatenate([[np.log(self.lam0)], np.log(sd), A.ravel(), mu, L[np.tril_indices(K)]])

    def _negloglik(self, theta, y):
        lam, h, A, mu, Q = _unpack(theta, y.shape[1])
        ll, _ = kalman_filter(y, nelson_siegel_loadings(lam, self.maturities), h, A, mu, Q)
        return 1e10 if not np.isfinite(ll) else -ll

    def fit(self, yields, warm_start: bool = False):
        y = np.asarray(yields, dtype=float)
        theta0 = self.theta_ if (warm_start and self.theta_ is not None) else self._start(y)
        res = minimize(self._negloglik, theta0, args=(y,), method="L-BFGS-B",
                       options={"maxiter": self.maxiter, "maxfun": 50 * self.maxiter})
        self.theta_, self.loglik_ = res.x, -res.fun
        self._y = y
        return self

    def update(self, yields):
        """Keep the estimated parameters and filter the factors on new data."""
        self._y = np.asarray(yields, dtype=float)
        return self

    @property
    def params(self):
        lam, h, A, mu, Q = _unpack(self.theta_, self._y.shape[1])
        return {"lambda": lam, "meas_sd": np.sqrt(h), "A": A, "mu": mu, "Q": Q}

    def filtered_factors(self):
        p = self.params
        Z = nelson_siegel_loadings(p["lambda"], self.maturities)
        return kalman_filter(self._y, Z, p["meas_sd"] ** 2, p["A"], p["mu"], p["Q"])[1]

    def forecast(self, h: int) -> np.ndarray:
        """Yield-curve forecast h months after the end of the estimation sample."""
        p = self.params
        b = self.filtered_factors()[-1]
        b_f = p["mu"] + np.linalg.matrix_power(p["A"], h) @ (b - p["mu"])
        return nelson_siegel_loadings(p["lambda"], self.maturities) @ b_f
