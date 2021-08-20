library(stringr)

progDataGenerate_random <- function(eta, 
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
    drawsTemp <- draws
    
    # randomely select the rows to censor
    censSamples <- sample(1:sampleSize, 
                          size = (sampleSize*prop), 
                          replace = FALSE)
    # randomly censor the selected rows
    drawsTemp[censSamples] <- runif(n = length(censSamples),
                                    min = 0,
                                    max = drawsTemp[censSamples])
    
    # organise into a data frame
    dfTemp <- data.frame(Time = drawsTemp, 
                         cens = rep(0, sampleSize))
    dfTemp$cens[censSamples] <- 1
    
    # save the data set in list object
    censoredDataSets[[i]] <- dfTemp
    
    i <- i + 1
  }
  
  return(censoredDataSets)
  
}
