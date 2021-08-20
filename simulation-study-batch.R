### Set wd ### ----
setwd("/home/rstudio/Sim-Study")

runsPerIter <- 100

### Librarys ### ----
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### Custom Functions ### ----
# To generate censored datasets
source("scripts//progDataGenerate.R")
source("scripts//dataGenerate.R")

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


### Result Storage ### ----
# Create a multi-dimensional list to hold the simulation results
# 5 models
models <- list(uninformativeInd = NULL, 
               informativeInd = NULL, 
               uninformativeJoint = NULL, 
               informativeJoint = NULL,
               missSpecifiedJoint = NULL)
# 5 levels of censoring
censoringLevels <- list(prop_0.00 = models, 
                        prop_0.25 = models,
                        prop_0.50 = models,
                        prop_0.75 = models, 
                        prop_0.95 = models)
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
                 0.95)

runs <- runsEmpty
m <- 1
n <- 1

# Loop through Simulation
for(i in 1:10) { #length(runs)
  
  print(str_c(i, "/1000"))
  
  # generate the data sets
  dfTemp <- try(progDataGenerate(eta = scale, 
                                 beta = shape, 
                                 props = proportions), silent = TRUE)
  
  # re sample if the frame arrangements doesn't work properly
  while(class(dfTemp) == "try-error") {
    dfTemp <- try(progDataGenerate(eta = scale, 
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
    
    # need to change the t1 for the miss specified model
    inputData$t1 <- 4.21
    
    m5 <- sampling(object = informativeJointMiss,
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
    
    j <- j + 1
    
   
  }
  
  ### Save Output every X runs###
  if((i %% runsPerIter) == 0){
    save(runs, file = str_c("simulation_true-data-gen_", m, ".RData"))
    m <- m + 1
    n <- 1
    print(srt_c((i / runsPerIter), " out of ", (1000 / runsPerIter), " saves completed"))
  } else {n <- n + 1}
  
}