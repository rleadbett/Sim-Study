// Joint-miss-specified-informative-prior
//
// Stan model for censoring paper.
//
// This is a Bayesian model for right-censored lifetime data 
// which followes a single two-paramater Weibull distribution with 
// an informative joint prior on the shape and scale paramaters of 
// the weibull distribution. (Kaminskiy2017)
//
// The prior has been miss specified so that the effect of accuracy
// of the prior on the posterior inference can be assessed.
//

functions {
  real fn(real tCDF){
    return log(-log1m(tCDF));
  }
}

data {
	
	int N_obs;
	int N_cens;
	real lifetime_obs[N_obs];
	real lifetime_cens[N_cens];
  real t1;   // should be 4.2
	real t2;   // should be 15
	
}

parameters {

	real<lower = 0> t1CDF;
	real<lower = t1CDF> t2CDF;  // the CDF at t2 must be greater than at t1

}

transformed parameters {
  
	real<lower = 0> beta;
	real<lower = 0> eta;

  // calculate Weibull paramaters based on the
  // draws from the CDF at t1 and t2.
  beta = (fn(t2CDF) - fn(t1CDF)) / log(t2 / t1);
  eta = exp(log(t1) - (fn(t1CDF) / beta));
  
}

model{

// Likelihood
// non-censored portion
for(i in 1:N_obs){
	target += weibull_lpdf(lifetime_obs[i]|beta, eta);
}
// censored portion
for(j in 1:N_cens){
	target += weibull_lccdf(lifetime_cens[j]|beta, eta);
}
  
// Prior models
// The prior was constructed by simulateing 100 datasets of size 
// n = 100 from the true Weibull distribution and estimating the 
// paramaters via MLE and calculating to value of the estimated 
// CDF at t1 and t2 to get a distribution.
t1CDF ~ beta(70.36, (140.72 - 70.36));
t2CDF ~ beta(40, (41.67 - 40));

}
