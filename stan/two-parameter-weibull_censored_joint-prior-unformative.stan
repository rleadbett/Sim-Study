// Joint-uninformative-prior
// 
// Stan model for censoring paper.
//
// This is a Bayesian model for right-censored lifetime data 
// which followes a single two-paramater Weibull distribution with 
// a joint-uninformative prior on the shape and scale paramaters. (Kaminskiy2017)
// 
// It is considered uninformative because a uniform(0, 1) is placed on the
// the CDF value at t1 and t2. (I am not to sure if this is truely uninformative)
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
	real<lower = t1CDF> t2CDF;

}

transformed parameters {
  
	real<lower = 0> beta;
	real<lower = 0> eta;


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
// beta(1, 1) is a uniform distiribution with support (0, 1)
t1CDF ~ beta(1, 1);
t2CDF ~ beta(1, 1);
}
