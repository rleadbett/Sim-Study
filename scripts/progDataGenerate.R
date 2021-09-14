#' A new function which uses the `dataGenerate.R` to 
#' generate multiple datasets of progresively more sevier
#' censoring 
#' 
#' @param eta The scale parameter of the weibull distribution
#' @param beta The shape parameter of the weibull distribution
#' @param StartCensProp The stating proportion of censoring
#' @param EndCensProp The end proportion of censoring
#' @param Steps The number of steps between \code{StartCensProp} 
#' and \code{EndCensProp}
#' 
#' @return A list with `n = Steps` elements, each of which is a
#' data frame with progresively more censoring

source("scripts//dataGenerate.R")
library(stringr)
library(dplyr)

progDataGenerate <- function(eta, 
                             beta, 
                             StartCensProp = 0.5, 
                             EndCensProp = 0.95, 
                             props = c(0, 0.25, 0.5, 0.75, 0.95)){
  
  # Check to see if the function has been run before with these parameters
  cacheString <- str_c("scripts//progDataGenerate-propData_shape-", 
                       beta, 
                       "_scale-", 
                       eta, 
                       ".RData")
  
  if (!file.exists(cacheString)) {
    
    # Firstly test the different window sizes to get desires proportions
    windows <- seq(500, 28000, 100)
    
    censProp <- rep(0, length(windows))
    i <- 1
    
    for (win in windows){
      samples <- rep(0,5)
      
      if(win <10000){
        samplesPerFrame <- 5000
      } else {
        samplesPerFrame <- win
      }
      
      for (j in 1:10) {
        SimData <- try(dataGenerate(eta = eta, 
                                beta = beta,
                                duration = win,
                                frames = 50,
                                samplesPerFrame = samplesPerFrame), silent = T)
        while(class(SimData) == "try-error"){
          SimData <- try(dataGenerate(eta = eta, 
                                      beta = beta,
                                      duration = win,
                                      frames = 50,
                                      samplesPerFrame = samplesPerFrame), silent = T)
        }
        samples[j] <- sum(SimData$cens == 1) / nrow(SimData)
      }
      
      censProp[i] <- mean(samples)
      i <- i + 1
      
    }
    
    base::save(windows, censProp, file = cacheString)
    
  } else {
    
    load(cacheString)
    
  }
  
  # Get the window size required for the correct proportions of censoring
  windowSizeForCens <- rep(NA, length(props)) 
  
  for (i in 1:length(props)) {
    
    prop <- props[i]
    windowSizeForCens[i] <- windows[min(which(censProp <= prop))]
  }
  
  
  # Generate the data sets
  censoredDataSets <- as.list(rep(NA, length(props)))
  
  samples <- rweibull(50000, 
                      shape = beta, 
                      scale = eta)
  
  for(i in 1:length(props)) {
    if (is.na(windowSizeForCens[i])){      # for the case when we want no cens
      dfTemp <- data.frame(Frame = seq(1, 100), 
                           Time = samples[1:100],
                           Startcens = rep(0, 100), 
                           Endcens = rep(0, 100),
                           cens = rep(0, 100))
      censoredDataSets[[i]] <- dfTemp
    } else {
      dfTemp <- dataGenerate(samples = samples,
                             duration = windowSizeForCens[i],
                             frames = 100)
      
      lifesPerFrame <- dfTemp %>% group_by(Frame) %>% summarise(n = n())
      lifesPerFrame$n <- cumsum(lifesPerFrame$n)
      
      selectFrames <- lifesPerFrame$Frame[lifesPerFrame$n < 100]
      
      censoredDataSets[[i]] <- dfTemp[dfTemp$Frame %in% selectFrames, ]
    }
    
  }
  
  # Name the elements in the data set list
  dfNames <- str_c(rep("Prop", length(props)), props)
  names(censoredDataSets) <- dfNames
  
  return(censoredDataSets)
}

