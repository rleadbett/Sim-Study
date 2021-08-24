// Joint-miss-specified-informative-prior
//
// Stan model for censoring paper.
//
// This is a Bayesian model for right-censored lifetime data 
// which followes a single two-paramater Weibull distribution with 
// an informative joint prior on the shape and scale paramaters of 
// the weibull distribution. (Kaminskiy2017)
//
// The prior has been misspecified to be more difuse by placing only 
// 80% of the mas within the 95% centred interval of the "true prior".
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
  real t1;   // should be 3.82
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
// CDF at t1 and t2 to get a distribution. The standard deviation
// of both priors was then multiplied by 1.5.

t1CDF ~ beta(33.33, 33.07);
t2CDF ~ beta(84.75, 3.07);

}
