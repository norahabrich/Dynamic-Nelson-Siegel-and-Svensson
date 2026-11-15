"""Factor loadings of the Nelson-Siegel (1987) and Svensson (1994) curves.

    y(tau) = b1 + b2 * (1 - exp(-l tau)) / (l tau)
                + b3 * ((1 - exp(-l tau)) / (l tau) - exp(-l tau))
                [+ b4 * ((1 - exp(-l2 tau)) / (l2 tau) - exp(-l2 tau))]   (Svensson)

b1 is the level, b2 the slope and b3 the curvature factor. With tau in months,
Diebold & Li (2006) fix l = 0.0609, which puts the maximum of the curvature
loading at 30 months.
"""

import numpy as np


def _slope(lam: float, tau: np.ndarray) -> np.ndarray:
    x = lam * np.asarray(tau, dtype=float)
    return (1.0 - np.exp(-x)) / x


def _curvature(lam: float, tau: np.ndarray) -> np.ndarray:
    return _slope(lam, tau) - np.exp(-lam * np.asarray(tau, dtype=float))


def nelson_siegel_loadings(lam: float, tau) -> np.ndarray:
    """Matrix of shape (len(tau), 3): level, slope and curvature loadings."""
    tau = np.asarray(tau, dtype=float)
    return np.column_stack([np.ones_like(tau), _slope(lam, tau), _curvature(lam, tau)])


def svensson_loadings(lam1: float, lam2: float, tau) -> np.ndarray:
    """Matrix of shape (len(tau), 4): Nelson-Siegel loadings plus a second curvature."""
    tau = np.asarray(tau, dtype=float)
    return np.column_stack([nelson_siegel_loadings(lam1, tau), _curvature(lam2, tau)])
