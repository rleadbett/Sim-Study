setwd("/home/rstudio/Sim-Study")
.libPaths("/home/rstudio/library")

library(dplyr)
library(stringr)

#### read in required data ----
base::load("results/true-censoring-mech/simulation_true-data-gen_1.RData") # renv masks the base r load function.
simStudyTrue <- runs
for(i in 2:10){
  fileStr <- str_c("results/true-censoring-mech/simulation_true-data-gen_", i, ".RData")
  simStudyTrue <- append(simStudyTrue, runs)
}

# window information
load("scripts/progDataGenerate-propData_shape-1.15_scale-5.254.RData")


#### preperation ----
windowSizes <- c(Inf)

for (i in c(0.25, 0.50, 0.75, 0.95)){
  pos <- which(censProp >= i) %>% max()
  windowSizes <- c(windowSizes, windows[pos]) 
}

windowSizes <- windowSizes / 365

# function to calculate the mean lifetime
meanLife <- function(shape, scale) scale * gamma(1 + (1 / shape))

# set up mesh for KS statistic
trueQs <- seq(0, 100, length.out = 1000)
mesh <- pweibull(trueQs, shape = 1.15, scale = 5.254)

windowQs <- c(NA, pweibull(windowSizes[-1], shape = 1.15, scale = 5.254))

# get mesh value where the window Q is reached
windowMeshCutoff <- c(length(mesh))
for (Q in windowQs[-1]){
  meshPos <- which(mesh <= Q) %>% max()
  windowMeshCutoff <- c(windowMeshCutoff, meshPos)
}


#### computation ----
start <- Sys.time()

# Function to apply operations over the whole of the simulation data

ApplyAll <- function(FUN, 
                     sims = 1:1000, 
                     props = names(simStudyTrue[[1]]),
                     models = names(simStudyTrue[[1]][[1]])) {
  lapply(sims, 
         function(x1) lapply(props, 
                             function(x2) lapply(models,
                                                 function(x3) FUN(simStudyTrue[[x1]][[x2]][[x3]])
                                                 )
                             )
         )
}

# Calculate the mean life posteriors

meanLifes <- ApplyAll(function(df) meanLife(shape = df$beta, scale = df$eta))

# Function to calculate the KS statistic

ksStat <- function(shape, scale, prop) {
  res <- mesh - pweibull(trueQs, shape, scale)
  
  cutOff <- windowMeshCutoff[which(names(simStudyTrue[[1]]) == prop)]
  
  pre <- max(abs(res[1:cutOff]))
  if(prop == "prop_0.00") {
    post <- NA
  } else {
    post <- max(abs(res[cutOff:1000]))
  }
  
  return(list(ksPre = pre, ksPost = post))
}

# Apply KS stat across the sim study
KSstatistics <- lapply(names(simStudyTrue[[1]]), 
                       function(prop) ApplyAll(function(df) ksStat(shape = df$beta, 
                                                                   scale = df$eta,
                                                                   prop = prop),
                                               props = c(prop)
                                               )
                       ) 
                       
# Add the mean and KS information back into original list
for(i in 1:1000) {
  for(prop in 1:5){
    for(model in 1:10){
      simStudyTrue[[i]][[prop]][[model]]$meanLife <- meanLifes[[i]][[prop]][[model]]
      simStudyTrue[[i]][[prop]][[model]]$KSPre <- KSstatistics[[prop]][[i]][[1]][[model]]$ksPre
      simStudyTrue[[i]][[prop]][[model]]$KSPost <- KSstatistics[[prop]][[i]][[1]][[model]]$ksPost
    }
  }
}

rm(meanLifes, KSstatistics, runs)
                       
end <- Sys.time()
print(end - start)

#### saving ----
save(simStudyTrue, file = "results/trueSimPostProcess.RData")                     
                      
                       
                       