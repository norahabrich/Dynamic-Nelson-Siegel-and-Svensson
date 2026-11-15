import numpy as np
import pytest

from dnsmodel import MATURITIES, fit_factors, fit_var1, load_yields, nelson_siegel_loadings, svensson_loadings
from dnsmodel.backtest import diebold_mariano, expanding_window_errors, random_walk
from dnsmodel.kalman import kalman_filter


def test_data_shape_and_dates():
    y = load_yields()
    assert y.shape == (348, 17)
    assert str(y.index[0].date()) == "1972-01-31" and str(y.index[-1].date()) == "2000-12-31"
    assert y.notna().all().all()


def test_nelson_siegel_loadings_properties():
    Z = nelson_siegel_loadings(0.0609, MATURITIES)
    assert Z.shape == (17, 3)
    assert np.allclose(Z[:, 0], 1)
    assert np.all(np.diff(Z[:, 1]) < 0)              # slope loading decays from 1 to 0
    tau = np.arange(1, 121)
    curv = nelson_siegel_loadings(0.0609, tau)[:, 2]
    assert tau[curv.argmax()] == 29                   # curvature peaks near 30 months (Diebold-Li)


def test_svensson_nests_nelson_siegel():
    S = svensson_loadings(0.0609, 0.1218, MATURITIES)
    assert np.allclose(S[:, :3], nelson_siegel_loadings(0.0609, MATURITIES))


def test_fit_factors_recovers_exact_betas():
    Z = nelson_siegel_loadings(0.0609, MATURITIES)
    true = np.array([[6.0, -2.0, 1.0], [8.0, 1.0, -0.5]])
    betas, fitted = fit_factors(true @ Z.T, Z)
    assert np.allclose(betas, true)
    assert np.allclose(fitted, true @ Z.T)


def test_var1_recovers_coefficients():
    rng = np.random.default_rng(0)
    A = np.array([[0.9, 0.05, 0.0], [0.0, 0.8, 0.1], [0.0, 0.0, 0.7]])
    x = np.zeros((5000, 3))
    for t in range(1, len(x)):
        x[t] = 0.1 + A @ x[t - 1] + rng.normal(scale=0.1, size=3)
    _, A_hat, _ = fit_var1(x)
    assert np.allclose(A_hat, A, atol=0.03)


def test_kalman_filter_tracks_simulated_factors():
    rng = np.random.default_rng(1)
    Z = nelson_siegel_loadings(0.0609, MATURITIES)
    A, mu, Q = np.diag([0.98, 0.95, 0.9]), np.array([7.0, -1.5, 0.0]), np.diag([0.1, 0.15, 0.3]) ** 2
    h = np.full(17, 0.05 ** 2)
    b = mu.copy()
    states, y = [], []
    for _ in range(300):
        b = mu + A @ (b - mu) + rng.multivariate_normal(np.zeros(3), Q)
        states.append(b)
        y.append(Z @ b + rng.normal(scale=0.05, size=17))
    ll, filt = kalman_filter(np.array(y), Z, h, A, mu, Q)
    assert np.isfinite(ll)
    assert np.sqrt(np.mean((filt[50:] - np.array(states)[50:]) ** 2)) < 0.1


def test_kalman_rejects_nonstationary_dynamics():
    Z = nelson_siegel_loadings(0.0609, MATURITIES)
    ll, _ = kalman_filter(np.zeros((5, 17)), Z, np.ones(17), np.eye(3) * 1.01, np.zeros(3), np.eye(3))
    assert ll == -np.inf


def test_diebold_mariano_direction():
    rng = np.random.default_rng(2)
    e_good, e_bad = rng.normal(scale=1, size=200), rng.normal(scale=2, size=200)
    stat, p = diebold_mariano(e_good, e_bad, h=1)
    assert stat < 0 and p < 0.01
    assert np.isnan(diebold_mariano(e_good, e_good, h=1)[0])


def test_backtest_uses_only_past_data():
    y = np.arange(40, dtype=float)[:, None] * np.ones((1, 3))
    errs = expanding_window_errors(y, {"rw": random_walk}, [1, 2], first_origin=10)
    # random walk on a series rising by 1 a month: the error equals the horizon
    assert np.allclose(errs["rw"][1], 1) and np.allclose(errs["rw"][2], 2)
    assert errs["rw"][1].shape == errs["rw"][2].shape == (28, 3)
