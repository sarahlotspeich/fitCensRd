#' Maximum likelihood estimator (MLE) for right-censored covariates in normal linear regression
#'
#' @param params0 Starting parameter values.
#' @param Y Name of column variable.
#' @param X Name of censored covariate variable.
#' @param C Name of observed version of \code{X}.
#' @param Z (Optional) name(s) of additional fully observed covariates. Default is \code{Z = NULL}
#' @param data A dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param steptol (Fed to \code{nlm()}) A positive scalar providing the minimum allowable relative step length. Default is \code{steptol = 1e-6}.
#' @param iterlim (Fed to \code{nlm()}) A positive integer specifying the maximum number of iterations to be performed before the program is terminated. Default is \code{iterlim = 100}.
#'
#' @return A list with the following two elements:
#' \item{coeff}{A dataframe containing coefficient and standard error estimates at convergence.}
#' \item{code}{An integer indicating why the optimization process terminated. See \code{?nlm} for details on values.}
#'
#' @export
#'
linreg_am <- function(params0, Y, X, C, Z = NULL, data, steptol = 1e-6, iterlim = 100) {
  suppressWarnings(
    ols <- nlm(f = loglik_am, p = params0, steptol = steptol, iterlim = iterlim, hessian = TRUE,
               Y = Y, X = X, C = C, Z = Z, data = data, negate = TRUE)
  )
  param_est <- ols$estimate
  if (!is.null(Z)) {
    param_se <- c(sqrt(diag(solve(ols$hessian)))[1:3], NA, sqrt(diag(solve(ols$hessian)))[5:6], NA)
    param_df <- data.frame(est = param_est, se = param_se)
    rownames(param_df) <- c("beta0", "beta1", "beta2", "sigmaY", "alpha0", "alpha1", "sigmaX")
  } else {
    param_se <- c(sqrt(diag(solve(ols$hessian)))[1:2], NA, sqrt(diag(solve(ols$hessian)))[4], NA)
    param_df <- data.frame(est = param_est, se = param_se)
    rownames(param_df) <- c("beta0", "beta1", "sigmaY", "alpha0", "sigmaX")
  }
  return(list(coeff = param_df, code = ols$code))
}
