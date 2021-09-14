# Define data generating function
dataGenerate <- function(eta = NULL, beta = NULL, 
                         samples = c(), 
                         frames = 50, 
                         startDay = 1000, 
                         duration = 2100,
                         samplesPerFrame = 1000) {
  # Takes input eta(scale) and beta(shape) and generates
  # lifetime data from Weibull distribution. These life-
  # times are then evenly assigned to 'frames' which 
  # represent frames of idlers along a section of a conve-
  # or. The lifetimes simulate the ongoing replacement of
  # idler frames on the conveyor.StartDay and duration
  # control the experimental time window.
  
  # convert the start day and duration into years
  
  startDay <-  startDay / 365
  duration <-  duration / 365
  
  # Each frame gets assigned 10 lifetimes.
  sampleSize <- frames * samplesPerFrame
  
  # Sample from the specified Weibull distribution
  if (length(samples) > 0) {
    lifetimeRaw <- samples
  } else {
    lifetimeRaw <- rweibull(sampleSize, 
                            shape = beta, 
                            scale = eta)
  }
  
  
  # Assign each frame 10 lifetimes.
  rawFrameLivesInd <- matrix(lifetimeRaw, nrow = frames, ncol = 20)
  
  # Apply cumulative sum to each row of the matrix 
  # (this is equivalent to time in days from start of simulation)
  rawFrameLivesCumulative <- c()
  
  for (i in 1:nrow(rawFrameLivesInd)) {
    
    rawFrameLivesCumulative <- rbind(rawFrameLivesCumulative, 
                                     cumsum(rawFrameLivesInd[i, ]))
    
  }
  
  # Catch error if the end Day + duration is longer than the minimum tenth frame failure date.
  if (startDay + duration > min(rawFrameLivesCumulative[, ncol(rawFrameLivesCumulative)])){
    stop('dataGenerate: duration exceeds population boundary')
  }
  
  # Subtract the start date from the cumulative values (shifting t = 0)
  t2CensLifetimes <- rawFrameLivesCumulative - startDay
  
  # Remove any negative values and any failures which occur after the 
  # end of the observation time (right censoring).
  t2CensLifetimes[(t2CensLifetimes < 0) | (t2CensLifetimes > duration)] <- 0
  
  # Rearrange the data into a long format similar to the industry data set of idler lifetimes.
  survSim <- c()
  
  for (i in 1:nrow(t2CensLifetimes)){
    
    # Get the lifetimes for each frame
    lifetimes <- t2CensLifetimes[i, ]
    
    # Select only those lifetimes which are greater than zero
    lifetimes <- lifetimes[lifetimes > 0]
    
    if (length(lifetimes) == 0) {   # If the lifetime spans the 
      # entire observation period
      
      # Then the lifetime equals the length of observation period with 
      # censoring at both beginning and end of life.
      lifeDF <- data.frame(Frame = i, 
                           Time = duration, 
                           Startcens = 1, 
                           Endcens = 1)
      
    } else {    # The case when there is at least one observed failure event
      
      # Add the last censored lifetime 
      lifetimes <- append(lifetimes, duration)
      
      # Reverse the cumulative sum so that we are left once
      # again with raw lifetimes.
      lifetimes <- c(lifetimes[1], diff(lifetimes))
      
      # Add frame number and start and end of live censoring 
      # indicator
      frame <- rep(i, length(lifetimes))
      startCens <- c(1, rep(0, length(lifetimes) - 1))
      endtCens <- c(rep(0, length(lifetimes) - 1), 1)
      
      # Create data frame
      lifeDF <- data.frame(Frame = frame, 
                           Time = lifetimes, 
                           Startcens = startCens, Endcens = endtCens)
    }
    
    # Append to the data frame `survSim`
    survSim <- rbind(survSim, lifeDF)
    
  }
  
  # Create a general censoring indicator which is irrespective of whether 
  # censoring occurs at the start or end censoring of the lifetime
  survSim$cens <- ifelse((survSim$Startcens == 1) | (survSim$Endcens == 1), 1, 0)
  
  # Round the Time to the nearest Day 
  #survSim$Time <- round(survSim$Time)
  
  # Return the simulated data set as the output
  return(survSim)
}
