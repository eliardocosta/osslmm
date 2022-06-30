data {
  int<lower=0> N; // total observations n*p
  int<lower=0> p; // number of subjects p;
  int<lower=1,upper=p> subj[N]; // index for subject predictor;
  vector[N] y; // response variable;
  real mu0; // hyperparameter for mu
  real<lower=0> sig2mu; // hyperparameter for mu
  real<lower=0> a1; // hyperparameter for sig2a
  real<lower=0> b1; // hyperparameter for sig2a
  real<lower=0> a2; // hyperparameter for sig2
  real<lower=0> b2; // hyperparameter for sig2
}
parameters {
  vector[p] a; // random effects;
  real mu;
  real<lower=0> sig2a; // variance of a
  real<lower=0> sig2; // variance of e
}
model {
  vector[N] eta;
  for (i in 1:N) {
    eta[i] = mu + sqrt(sig2a)*a[subj[i]];
  }
  mu ~ normal(mu0, sqrt(sig2mu));
  a ~ normal(0, 1);
  sig2a ~ inv_gamma(a1, b1);
  sig2 ~ inv_gamma(a2, b2);
  y ~ normal(eta, sqrt(sig2));
}
