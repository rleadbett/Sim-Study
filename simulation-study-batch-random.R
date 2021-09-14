### Set up ### ----
# set wd
setwd("/home/rstudio/Sim-Study")
.libPaths("/home/rstudio/library")
# set seed
set.seed(64467364)

# set variables
batchNo <- 2
runsPerIter <- 100
nruns <- 500

### Librarys ### ----
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### Custom Functions ### ----
# To generate censored datasets
source("scripts//progDataGenerate_random.R")

# To pass data into Stan
source("scripts//stanPrep.R")

### Load in Stan Models ### ----

uninformativeJoint <- stan_model("stan//two-parameter-weibull_censored_joint-prior-unformative.stan",
                                 verbose = TRUE)
Sys.sleep(0.1)
informativeJoint <- stan_model("stan//two-parameter-weibull_censored_joint-prior-informative.stan",
                               verbose = TRUE)
Sys.sleep(0.1)
uninformativeInd <- stan_model("stan//two-parameter-weibull_censored_independed-uninformative.stan",
                               verbose = TRUE)
Sys.sleep(0.1)
informativeInd <- stan_model("stan//two-parameter-weibull_censored_independed-informative.stan",
                             verbose = TRUE)
Sys.sleep(0.1)
informativeJointMiss <- stan_model("stan//two-parameter-weibull_censored_joint-miss-specified-prior-informative.stan",
                                   verbose = TRUE)
Sys.sleep(0.1)
informativeJointMiss_sd1 <- stan_model("stan//two-parameter-weibull_censored_joint-miss-specified-prior-informative_var1.stan",
                                       verbose = TRUE)
Sys.sleep(0.1)
informativeJointMiss_sd2 <- stan_model("stan//two-parameter-weibull_censored_joint-miss-specified-prior-informative_var2.stan",
                                       verbose = TRUE)

### Result Storage ### ----
# Create a multi-dimensional list to hold the simulation results
# 5 models
models <- list(uninformativeInd = NULL, 
               informativeInd = NULL, 
               uninformativeJoint = NULL, 
               informativeJoint = NULL,
               missSpecifiedJoint = NULL, 
               missSpecifiedjoint_sd1 = NULL,
               missSpecifiedjoint_sd2 = NULL,
               missSpecifiedjoint_t1.1 = NULL,
               missSpecifiedjoint_t1.2 = NULL,
               missSpecifiedjoint_t0.9 = NULL)
# 5 levels of censoring
censoringLevels <- list(prop_0.00 = models, 
                        prop_0.25 = models,
                        prop_0.50 = models,
                        prop_0.75 = models, 
                        prop_0.90 = models)
# 1000 runs of the simulation
runsEmpty <- list()

for (i in 1:runsPerIter) {
  runsEmpty[[i]] <- censoringLevels
}

### Simulation ### ----

# Define the data generating constants
shape <- 1.15
scale <- 5.254
proportions <- c(0, 
                 0.25, 
                 0.5, 
                 0.75, 
                 0.90)

runs <- runsEmpty

if(batchNo == 1){m <- 1}
if(batchNo == 2){m <- 6}

n <- 1

# Loop through Simulation
for(i in 1:nruns) { #length(runs)
  
  print(str_c(i, "/", nruns))
  
  # generate the data sets
  dfTemp <- try(progDataGenerate_random(eta = scale, 
                                        beta = shape, 
                                        props = proportions), silent = TRUE)
  
  # re sample if the frame arrangements doesn't work properly
  while(class(dfTemp) == "try-error") {
    dfTemp <- try(progDataGenerate_TypeII(eta = scale, 
                                          beta = shape, 
                                          props = proportions), silent = TRUE)
  }
  
  j <- 1
  
  # fit all models to all data sets
  for(dfName in names(dfTemp)){
    # prep data
    inputData <- stanPrep(dfTemp[[dfName]])
    
    # fit the models
    m1 <- sampling(object = uninformativeInd, 
                   data = inputData, 
                   chains = 4, 
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    m2 <- sampling(object = informativeInd, 
                   data = inputData, 
                   chains = 4, 
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    
    inputData$t1 <- 3.82   # need to add some extra inputs for the joint models
    inputData$t2 <- 15
    
    m3 <- sampling(object = uninformativeJoint, 
                   data = inputData, 
                   chains = 4, 
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    
    m4 <- sampling(object = informativeJoint, 
                   data = inputData, 
                   chains = 4, 
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    
    m6 <- sampling(object = informativeJointMiss_sd1, 
                   data = inputData, 
                   chains = 4, 
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    
    m7 <- sampling(object = informativeJointMiss_sd2, 
                   data = inputData, 
                   chains = 4, 
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    
    # need to change the t1 for the miss specified model
    inputData$t1 <- 4.21
    
    m5 <- sampling(object = informativeJointMiss,
                   data = inputData,
                   chains = 4,
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    
    # for the 10% over estimated priors we need to reset
    inputData$t1 <- 3.82 * 1.1
    inputData$t2 <- 15 * 1.1
    
    m8 <- sampling(object = informativeJoint,
                   data = inputData,
                   chains = 4,
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    # for the 20% over estimated priors we need to reset
    inputData$t1 <- 3.82 * 1.2
    inputData$t2 <- 15 * 1.2
    
    m9 <- sampling(object = informativeJoint,
                   data = inputData,
                   chains = 4,
                   cores = 4,
                   iter = 1500,
                   warmup = 500, 
                   refresh = 0)
    
    # for the 10% under estimated priors we need to reset
    inputData$t1 <- 3.82 * 0.9
    inputData$t2 <- 15 * 0.9
    
    m10 <- sampling(object = informativeJoint,
                    data = inputData,
                    chains = 4,
                    cores = 4,
                    iter = 1500,
                    warmup = 500, 
                    refresh = 0)  
    
    # extract the posteriors and save in list object
    
    runs[[n]][[j]]$uninformativeInd <- extract(m1)[c(1, 4)]
    runs[[n]][[j]]$informativeInd <- extract(m2)[c(1, 4)]
    runs[[n]][[j]]$uninformativeJoint <- extract(m3)[c(3, 4)]
    runs[[n]][[j]]$informativeJoint <- extract(m4)[c(3, 4)]
    runs[[n]][[j]]$missSpecifiedJoint <- extract(m5)[c(3, 4)]
    runs[[n]][[j]]$missSpecifiedjoint_sd1 <- extract(m6)[c(3, 4)]
    runs[[n]][[j]]$missSpecifiedjoint_sd1 <- extract(m7)[c(3, 4)]
    runs[[n]][[j]]$missSpecifiedjoint_t1.1 <- extract(m8)[c(3, 4)]
    runs[[n]][[j]]$missSpecifiedjoint_t1.2 <- extract(m9)[c(3, 4)]
    runs[[n]][[j]]$missSpecifiedjoint_t0.9 <- extract(m10)[c(3, 4)]
    
    j <- j + 1
    
    
  }
  
  ### Save Output every X runs###
  if((i %% runsPerIter) == 0){
    save(runs, file = str_c("results//random-censoring//simulation_random_", m, ".RData"))
    m <- m + 1
    n <- 1
    print(str_c((i / runsPerIter), " out of ", (nruns / runsPerIter), " saves completed"))
  } else {n <- n + 1}
  
}