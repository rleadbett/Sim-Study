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

### Libraries ### ----
library(dplyr)
library(fitdistrplus)

### Custom Functions ### ----
# To generate censored datasets
source("scripts//progDataGenerate.R")
source("scripts//dataGenerate.R")

# To pass data into fitdist
toFitdist <- function(df){
  
  dfFitDist <- data.frame(left = df$Time, 
                          right = df$Time)
  
  dfFitDist$right[df$cens == 1] <- NA
  
  return(dfFitDist)
  
}


### Result Storage ### ----
# Create a multi-dimensional list to hold the simulation results
# 5 models
models <- list(maximumLikelihood = NULL)

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
    inputData <- toFitdist(as.data.frame(dfTemp[[dfName]]))
    
    # fit the models
    m1 <- fitdistcens(inputData, distr = "weibull")
    
    # extract the posteriors and save in list object
    
    runs[[n]][[j]]$maximumLikelihood <- m1
    
    j <- j + 1
    
    
  }
  
  ### Save Output every X runs###
  if((i %% runsPerIter) == 0){
    save(runs, file = str_c("ML_simulation_true-data-gen_", m, ".RData"))
    m <- m + 1
    n <- 1
    print(str_c((i / runsPerIter), " out of ", (nruns / runsPerIter), " saves completed"))
  } else {n <- n + 1}
  
}
