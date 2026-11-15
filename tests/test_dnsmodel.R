# Unit tests for the base-R implementation. Run from the repository root:
#   Rscript tests/test_dnsmodel.R
source("R/dnsmodel.R")
check <- function(name, ok) { if (!isTRUE(ok)) stop("FAILED: ", name); cat("ok  ", name, "\n") }

Y <- load_yields("data/bonds.csv")
check("data shape and dates", all(dim(Y) == c(348, 17)) && rownames(Y)[1] == "1972-01-31" && rownames(Y)[348] == "2000-12-31")

Z <- nelson_siegel_loadings(0.0609, MATURITIES)
check("level loading is 1", all(Z[, 1] == 1))
check("slope loading decays", all(diff(Z[, 2]) < 0))
check("curvature peaks near 30 months", which.max(nelson_siegel_loadings(0.0609, 1:120)[, 3]) == 29)
check("Svensson nests Nelson-Siegel", isTRUE(all.equal(svensson_loadings(0.0609, 0.1218, MATURITIES)[, 1:3], Z, check.attributes = FALSE)))

true_b <- rbind(c(6, -2, 1), c(8, 1, -0.5))
f <- fit_factors(true_b %*% t(Z), Z)
check("OLS recovers exact factors", max(abs(f$betas - true_b)) < 1e-10)

set.seed(0)
A <- matrix(c(0.9, 0.05, 0, 0, 0.8, 0.1, 0, 0, 0.7), 3, 3, byrow = TRUE)
x <- matrix(0, 5000, 3)
for (t in 2:5000) x[t, ] <- 0.1 + A %*% x[t - 1, ] + rnorm(3, sd = 0.1)
check("VAR(1) recovers coefficients", max(abs(fit_var1(x)$A - A)) < 0.03)

check("non-stationary dynamics rejected",
      kalman_filter(matrix(0, 5, 17), Z, rep(1, 17), diag(3) * 1.01, rep(0, 3), diag(3))$loglik == -Inf)

set.seed(2)
d <- diebold_mariano(rnorm(200, sd = 1), rnorm(200, sd = 2), 1)
check("Diebold-Mariano direction", d[["stat"]] < 0 && d[["p"]] < 0.01)

yl <- outer(0:39, rep(1, 3))
e <- expanding_window_errors(yl, list(rw = random_walk), c(1, 2), 11)
check("backtest uses only past data", all(e$rw[["1"]] == 1) && all(e$rw[["2"]] == 2) && nrow(e$rw[["1"]]) == 28)

# Cross-language check: R two-step factors equal the Python ones
ref <- read.csv("tests/reference_python_two_step_factors.csv")
fr <- fit_factors(Y, Z)$betas
check("same two-step factors as Python (max diff < 1e-9)", max(abs(fr - as.matrix(ref[, -1]))) < 1e-9)
cat("All R tests passed\n")
