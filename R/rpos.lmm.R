#' Title
#'
#' @param n a
#' @param p a
#' @param mu0 a
#' @param sig2mu a
#' @param a1 a
#' @param b1 a
#' @param a2 a
#' @param b2 a
#' @param chains a
#' @param iter a
#' @param nburn a
#'
#' @return a
#' @export
#'
#' @examples a
rpos.lmm <- function(n, p, mu0 = 0, sig2mu = 1E2, a1 = 2.1, b1 = 4,
                     a2 = 2.1, b2 = 4, chains = 2, iter = 1E2, nburn = 1E2) {
  subj <- rep(1:p, n)
  data_stan <- list(N = n*p, p = p, subj = subj, y = y, mu0 = mu0,
                    sig2mu = sig2mu, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
  fit <- rstan::sampling(stanmodels$lmm, data = data_stan, pars = "mu",
                         chains = chains, iter = iter, warmup = nburn)
  return(fit)
}
