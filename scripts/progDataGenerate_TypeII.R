library(stringr)

progDataGenerate_TypeII <- function(eta, 
                                    beta, 
                                    props = c(0, 0.25, 0.5, 0.75, 0.95)){
  sampleSize <- 100
  
  # draw samples
  draws <- rweibull(n = sampleSize, 
                    shape = beta, 
                    scale = eta)
  
  # make object to store the datasets
  censoredDataSets <- as.list(rep(NA, length(props)))
  names(censoredDataSets) <- str_c("Prop", props)
  
  # censor the data
  i <- 1
  for(prop in props){
    
    # Order the data
    drawsTemp <- draws[order(draws)]
    
    # Censor the correct proportion of rows
    censRow <- sampleSize - (prop * sampleSize)
    
    cens <- rep(0, sampleSize)
    
    if (prop > 0){
      
      drawsTemp[(censRow + 1):sampleSize] <- drawsTemp[censRow]
      cens[(censRow + 1):sampleSize] <- 1
      
    } 
    
    # organise into a data frame
    dfTemp <- data.frame(Time = drawsTemp, 
                         cens = cens)
    
    # save the data set in list object
    censoredDataSets[[i]] <- dfTemp
    
    i <- i + 1
  }
  
  return(censoredDataSets)
  
}