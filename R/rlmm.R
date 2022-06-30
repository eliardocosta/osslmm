#' Title
#'
#' @param n a
#' @param p a
#' @param mu a
#' @param sig2 a
#' @param sig2a a
#'
#' @return a
#' @export
#'
#' @examples a
rlmm <- function(n, p, mu, sig2, sig2a) {
  a <- rnorm(p, 0, sqrt(sig2a))
  subj <- rep(1:p, n)
  y <- numeric(n*p)
  for (i in 1:(n*p)) {
    y[i] <- mu + a[subj[i]] + rnorm(1, 0, sqrt(sig2))
  }
  return(y)
}
