# Objective function: log(L) for independent censoring ---------
loglik_am <- function(params, Y, X, C, Z = NULL, data, negate = FALSE) {
  if (!is.null(Z)) {
    # Analysis model parameters P(Y|X) ----------------
    beta0 <- params[1]
    beta1 <- params[2]
    beta2 <- params[3]
    sigY <- params[4]
    # ---------------- Analysis model parameters P(Y|X)
    # Censored model parameters P(X) ------------------
    sigX <- params[7]
    # ------------------ Censored model parameters P(X)
    # Log-likelihood contribution of observed X -------
    obs <- !is.na(data[, X])
    muY <- beta0 + beta1 * data[obs, X] + beta2 * data[obs, Z]
    eY <- data[obs, Y] - muY
    pYgivX_obs <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    muX <- params[5] + params[6] * data[obs, Z]
    eX <- data[obs, X] - muX
    pX_obs <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    pYX_obs <- pYgivX_obs * pX_obs
    ll <- sum(log(pYX_obs))
    # ------- Log-likelihood contribution of observed X
    # Log-likelihood contribution of censored X -------
    n_cens <- sum(!obs)
    m <- length(unique(data[obs, X]))
    muY <- beta0 + beta1 * rep(x = unique(data[obs, X]), each = n_cens) + beta2 * rep(x = data[!obs, Z], times = m)
    eY <- rep(x = data[!obs, Y], times = m) - muY
    pYgivX_cens <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    muX <- params[5] + params[6] * rep(x = data[!obs, Z], times = m)
    eX <- rep(x = unique(data[obs, X]), each = n_cens) - muX
    pX_cens <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    pYX_cens <- rowsum(x = c(pYgivX_cens * pX_cens), group = rep(x = 1:n_cens, times = m))
    ll <- ll + sum(log(pYX_cens))
    # ------- Log-likelihood contribution of censored X
  } else {
    # Analysis model parameters P(Y|X) ----------------
    beta0 <- params[1]
    beta1 <- params[2]
    sigY <- params[3]
    # ---------------- Analysis model parameters P(Y|X)
    # Censored model parameters P(X) ------------------
    muX <- params[4]
    sigX <- params[5]
    # ------------------ Censored model parameters P(X)
    # Log-likelihood contribution of observed X -------
    obs <- data[, X] <= data[, C]
    muY <- beta0 + beta1 * data[obs, X]
    e <- data[obs, Y] - muY
    pYX_obs <- 1 / (2 * pi * sigY * sigX) * exp(- e ^ 2 / (2 * sigY ^ 2) - (data[obs, X] - muX) ^ 2 / (2 * sigX ^ 2))
    ll <- sum(log(pYX_obs))
    # ------- Log-likelihood contribution of observed X
    # Log-likelihood contribution of censored X -------
    Q <- sqrt((sigX ^ (-2)) + beta1 ^ 2 / sigY ^ 2)
    emu <- data[!obs, Y] - beta0 - beta1 * muX
    pYX_cens <- 1 / (2 * pi * (sigY ^ 2 + beta1 ^ 2 * sigX ^ 2)) *
      pnorm(Q * (muX + (beta1 * emu) / (sigY ^ 2 * Q ^ 2) - data[!obs, C])) *
      exp(- (emu ^ 2 / (2 * sigY ^ 2)) * (1 - (beta1 ^ 2 / (sigY ^ 2 * Q ^ 2))))
    ll <- ll + sum(log(pYX_cens))
    # ------- Log-likelihood contribution of censored X
  }
  if (negate) return(-ll) else return(ll)
}
