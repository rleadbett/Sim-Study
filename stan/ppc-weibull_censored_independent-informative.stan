//
// This file is for prior predictive checks of;
// two-parameter-weibull_censored_independent-informative.stan
//

data {
  
  int<lower=0> N;
  
}
model {
  
}
generated quantities {
  
  real beta = gamma_rng(156.25, 120.19);
  real lambda = gamma_rng(6.25, 42.59);
  real log_lambda;
  real<lower = 0> eta;
  vector[N] y_sim;
  
  log_lambda = log(lambda);
  eta = exp((-1 / beta) * log_lambda);
  
  for (n in 1:N) y_sim[N] = weibull_rng(beta, eta);
  
}

