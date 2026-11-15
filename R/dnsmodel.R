# Dynamic Nelson-Siegel and Svensson models in base R (no package needed).
# Mirror of the Python package dnsmodel/: same models, same conventions, so that
# results can be compared line by line (see compare_r_python.py).

MATURITIES <- c(3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120)

# ---- Data --------------------------------------------------------------------

load_yields <- function(path) {
  df <- read.csv(path, sep = ";", check.names = FALSE)
  Y <- as.matrix(df[, -1])
  ym <- strsplit(df$date, ":")
  first <- as.Date(sprintf("%s-%s-01", sapply(ym, `[`, 1), sapply(ym, `[`, 2)))
  dates <- seq(first[1], by = "month", length.out = length(first) + 1)[-1] - 1  # month ends
  colnames(Y) <- sub("^M", "", colnames(Y))
  rownames(Y) <- format(dates)
  stopifnot(all(as.numeric(colnames(Y)) == MATURITIES))
  Y
}

# ---- Loadings ----------------------------------------------------------------

ns_slope <- function(lam, tau) (1 - exp(-lam * tau)) / (lam * tau)
ns_curv  <- function(lam, tau) ns_slope(lam, tau) - exp(-lam * tau)

nelson_siegel_loadings <- function(lam, tau) cbind(1, ns_slope(lam, tau), ns_curv(lam, tau))
svensson_loadings <- function(lam1, lam2, tau) cbind(nelson_siegel_loadings(lam1, tau), ns_curv(lam2, tau))

# ---- Two-step DNS (Diebold & Li, 2006) -----------------------------------------

fit_factors <- function(Y, Z) {
  betas <- t(qr.solve(Z, t(Y)))          # OLS of each curve on the loadings
  list(betas = betas, fitted = betas %*% t(Z))
}

fit_var1 <- function(x, h = 1) {          # direct h-step VAR(1): x[t+h] = c + A x[t]
  n <- nrow(x)
  X <- cbind(1, x[1:(n - h), , drop = FALSE])
  coef <- qr.solve(X, x[(h + 1):n, , drop = FALSE])
  resid <- x[(h + 1):n, , drop = FALSE] - X %*% coef
  list(c = coef[1, ], A = t(coef[-1, , drop = FALSE]), S = cov(resid))
}

fit_ar1 <- function(x, h = 1) {           # direct h-step AR(1) on each column
  n <- nrow(x)
  co <- sapply(seq_len(ncol(x)), function(j) qr.solve(cbind(1, x[1:(n - h), j]), x[(h + 1):n, j]))
  list(c = co[1, ], phi = co[2, ])
}

two_step_forecast <- function(Z, dynamics = "ar1") {
  function(history, h) {
    b <- fit_factors(history, Z)$betas
    last <- b[nrow(b), ]
    bf <- if (dynamics == "ar1") { m <- fit_ar1(b, h); m$c + m$phi * last }
          else { m <- fit_var1(b, h); m$c + m$A %*% last }
    as.vector(Z %*% bf)
  }
}

random_walk <- function(history, h) history[nrow(history), ]

ar1_yields <- function(history, h) {
  m <- fit_ar1(history, h)
  m$c + m$phi * history[nrow(history), ]
}

# ---- One-step DNS with the Kalman filter (Diebold, Rudebusch & Aruoba, 2006) --

K <- 3

unpack <- function(theta, n) {
  i <- 1
  lam <- exp(theta[i]); i <- i + 1
  hv <- exp(2 * theta[i:(i + n - 1)]); i <- i + n
  A <- matrix(theta[i:(i + 8)], 3, 3, byrow = TRUE); i <- i + 9
  mu <- theta[i:(i + 2)]; i <- i + 3
  L <- matrix(0, 3, 3)
  idx <- rbind(c(1, 1), c(2, 1), c(2, 2), c(3, 1), c(3, 2), c(3, 3))  # same order as numpy tril_indices
  for (k in 1:6) L[idx[k, 1], idx[k, 2]] <- theta[i + k - 1]
  diag(L) <- exp(diag(L))
  list(lam = lam, h = hv, A = A, mu = mu, Q = L %*% t(L))
}

stationary_cov <- function(A, Q) {
  if (max(Mod(eigen(A, only.values = TRUE)$values)) >= 0.999) return(NULL)
  P <- matrix(solve(diag(9) - kronecker(A, A), as.vector(Q)), 3, 3)
  (P + t(P)) / 2
}

kalman_filter <- function(Y, Z, h, A, mu, Q) {
  n_t <- nrow(Y); n <- ncol(Y)
  P <- stationary_cov(A, Q)
  if (is.null(P)) return(list(loglik = -Inf, states = NULL))
  b <- mu
  loglik <- 0
  states <- matrix(NA_real_, n_t, K)
  H <- diag(h, n)
  for (t in seq_len(n_t)) {
    if (t > 1) {
      b <- mu + A %*% (b - mu)
      P <- A %*% P %*% t(A) + Q
    }
    v <- Y[t, ] - Z %*% b
    F <- Z %*% P %*% t(Z) + H
    C <- tryCatch(chol(F), error = function(e) NULL)   # upper triangular, F = C'C
    if (is.null(C)) return(list(loglik = -Inf, states = NULL))
    Finv_v <- backsolve(C, forwardsolve(t(C), v))
    Finv_Z <- backsolve(C, forwardsolve(t(C), Z))
    loglik <- loglik - 0.5 * (n * log(2 * pi) + 2 * sum(log(diag(C))) + sum(v * Finv_v))
    PZt <- P %*% t(Z)
    b <- b + PZt %*% Finv_v
    P <- P - PZt %*% Finv_Z %*% P
    states[t, ] <- b
  }
  list(loglik = loglik, states = states)
}

kalman_start <- function(Y, lam0 = 0.0609, tau = MATURITIES) {
  Z <- nelson_siegel_loadings(lam0, tau)
  fs <- fit_factors(Y, Z)
  v <- fit_var1(fs$betas, 1)
  mu <- solve(diag(3) - v$A, v$c)
  L <- t(chol(v$S))
  sd <- apply(Y - fs$fitted, 2, sd) * sqrt((nrow(Y) - 1) / nrow(Y))  # population sd, as numpy
  c(log(lam0), log(sd), as.vector(t(v$A)), mu,
    log(L[1, 1]), L[2, 1], log(L[2, 2]), L[3, 1], L[3, 2], log(L[3, 3]))
}

kalman_negloglik <- function(theta, Y, tau = MATURITIES) {
  p <- unpack(theta, ncol(Y))
  ll <- kalman_filter(Y, nelson_siegel_loadings(p$lam, tau), p$h, p$A, p$mu, p$Q)$loglik
  if (!is.finite(ll)) 1e10 else -ll
}

kalman_fit <- function(Y, theta0 = NULL, tau = MATURITIES, maxit = 2000) {
  if (is.null(theta0)) theta0 <- kalman_start(Y, tau = tau)
  res <- optim(theta0, kalman_negloglik, Y = Y, tau = tau, method = "L-BFGS-B",
               control = list(maxit = maxit))
  list(theta = res$par, loglik = -res$value, convergence = res$convergence, Y = Y, tau = tau)
}

kalman_params <- function(fit) unpack(fit$theta, ncol(fit$Y))

kalman_filtered <- function(fit, Y = fit$Y) {
  p <- kalman_params(fit)
  kalman_filter(Y, nelson_siegel_loadings(p$lam, fit$tau), p$h, p$A, p$mu, p$Q)$states
}

kalman_forecast <- function(fit, h, Y = fit$Y) {
  p <- kalman_params(fit)
  s <- kalman_filtered(fit, Y)
  b <- s[nrow(s), ]
  Ah <- diag(3); for (k in seq_len(h)) Ah <- Ah %*% p$A
  as.vector(nelson_siegel_loadings(p$lam, fit$tau) %*% (p$mu + Ah %*% (b - p$mu)))
}

# ---- Backtest ----------------------------------------------------------------

rmse <- function(e) sqrt(colMeans(e^2))

# forecasters: named list of functions (history, h) -> curve, or the string "kalman"
expanding_window_errors <- function(Y, forecasters, horizons, first_origin, kalman_refit = 12) {
  origins <- first_origin:(nrow(Y) - max(horizons))
  out <- list()
  for (name in names(forecasters)) {
    f <- forecasters[[name]]
    errs <- lapply(horizons, function(h) matrix(NA_real_, length(origins), ncol(Y)))
    fit <- NULL; last_fit <- NA
    for (k in seq_along(origins)) {
      t <- origins[k]
      history <- Y[1:t, , drop = FALSE]
      if (identical(f, "kalman") && (is.null(fit) || t - last_fit >= kalman_refit)) {
        fit <- kalman_fit(history, theta0 = if (is.null(fit)) NULL else fit$theta)
        last_fit <- t
      }
      for (j in seq_along(horizons)) {
        h <- horizons[j]
        fc <- if (identical(f, "kalman")) kalman_forecast(fit, h, history) else f(history, h)
        errs[[j]][k, ] <- Y[t + h, ] - fc
      }
    }
    names(errs) <- horizons
    out[[name]] <- errs
  }
  out
}

diebold_mariano <- function(e1, e2, h) {
  d <- e1^2 - e2^2
  n <- length(d)
  dc <- d - mean(d)
  v <- sum(dc^2) / n
  if (h > 1) for (k in 1:(h - 1)) v <- v + 2 * (1 - k / h) * sum(dc[(k + 1):n] * dc[1:(n - k)]) / n
  if (v <= 0) return(c(stat = NA, p = NA))
  dm <- mean(d) / sqrt(v / n) * sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
  c(stat = dm, p = 2 * pt(abs(dm), df = n - 1, lower.tail = FALSE))
}
