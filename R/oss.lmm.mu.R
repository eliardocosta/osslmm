#' Title
#'
#' @param n a
#' @param c a
#' @param mu0 a
#' @param sig2mu a
#' @param a1 a
#' @param b1 a
#' @param a2 a
#' @param b2 a
#' @param lrep a
#' @param pmax a
#' @param pby a
#' @param prep a
#' @param chains a
#' @param iter a
#' @param plot a
#' @param ...
#'
#' @return a
#' @export
#'
#' @examples a
oss.lmm.mu <- function(n, c, mu0 = 0, sig2mu = 1E2, a1 = 2.1, b1 = 4,
                       a2 = 2.1, b2 = 4, lrep = 5E1, pmax = 5E2,  pby = 1E2,
                       prep = 3, nchains = 2, niter = 9E2, nburn = 3E2, nthin = 1,
                       ncores = 1, plot = TRUE, ...) {
  cl <- match.call()
  ps <- c(1:9, rep(seq(10, pmax, pby), each = prep))
  tc <- numeric(length(ps))
  m.rhat <- numeric(length(ps))

  pb <- txtProgressBar(min = 1, max = length(tc), style = 3)
  for(i in 1:length(tc)) {
    loss <- numeric(lrep)
    rhat <- numeric(lrep)
    subj <- rep(1:ps[i], n)
    for (j in 1:lrep) {
      mu <- rnorm(1, mu0, sqrt(sig2mu))
      sig2a <- 1/rgamma(1, shape = a1, rate = b1)
      sig2 <- 1/rgamma(1, shape = a2, rate = b2)

      y <- rlmm(n = n, p = ps[i], mu = mu, sig2 = sig2, sig2a = sig2a)
      data_stan <- list(N = n*ps[i], p = ps[i], subj = subj, y = y, mu0 = mu0,
                        sig2mu = sig2mu, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
      fit <- rstan::sampling(stanmodels$lmm, data = data_stan, pars = "mu",
                             chains = nchains, iter = niter, warmup = nburn,
                             thin = nthin, cores = ncores, refresh = 0)
      mu.pos <- as.vector(rstan::extract(fit, pars = "mu")$mu)
      loss[j] <- var(mu.pos) + c*n*ps[i]
      rhat[j] <- rstan::summary(fit, pars = "mu")$summary[10]
    }
    tc[i] <- mean(loss)
    m.rhat[i] <- max(rhat)
    setTxtProgressBar(pb, i)
  }
  Y <- log(tc - c*ps*n)
  mod <- stats::lm(Y ~ I(log(ps + 1)))
  E <- as.numeric(exp(mod$coef[1]))
  G <- as.numeric(-mod$coef[2])
  popt <- ceiling((E*G/(c*n))^(1/(G + 1)) - 1)
  if (plot == TRUE) {
    plot(ps, tc, xlim = c(0, pmax), ylim = c(min(tc) - 0.5, max(tc) + 0.5), xlab = "p", ylab = "TC(p)")
    curve.tc <- function(x) {c*x*n + E/(1 + x)^G}
    plot(function(x) curve.tc(x), 0, pmax, col = "blue", lty = 2, add = TRUE)
    graphics::abline(v = popt, col = "red", lty = 2)
  }
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nOptimal number of subjects:\n")
  cat("p  = ", popt, "\n", sep = "")
  cat("Rhat = ", round(mean(m.rhat), 4), " (", round(sd(m.rhat), 4), ")\n", sep = "")
}
