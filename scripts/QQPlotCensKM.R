#' Generates a Q-Q plot for the censored data (Kaplan-Meier) and the true weibull distribution
#' 
#' @param df A data frame with one column `lifetime` containing 
#' the lifetimes and another column `cens` containing the censoring
#' indicator for each of the lifetimes.
#' @param shape The shape parameter of the true Weibull distribution
#' @param scale The scale parameter of the true Weibull distribution
#' 

QQPlotCensKM <- function(df, shape = 1.15, scale = 5.254){
  
  df <- df[order(df$Time), ]
  
  # Get the exposure times of failures
  ts <- unique(df$Time[df$cens == 0])
  ts <- ts[order(ts)]
  
  # Get the failures at each exposure time
  ds <- sapply(X = ts, FUN = function(x){length(which(df$Time == x))})
  
  # Get the at risk units at each exposure time
  ns <- sapply(X = ts, FUN = function(x){length(which(df$Time >= x))})
  
  # Calculate the Kaplan-Meier nonparametric approximation
  # of the survival function
  survEst <- cumprod(1 - (ds/ns))
  
  # Convert this to the CDF
  cdfEst <- 1 - survEst
  
  # Calculate the true quantiles at each exposure
  trueQuantiles <- pweibull(ts, shape = shape, scale = scale)
  
  # Generate Q-Q Plot
  plot(x = cdfEst, y = trueQuantiles,
       main = "Q-Q Plot (K-M)",
       xlab = "Data",
       ylab = "True Weibull",
       xlim = c(0, 1),
       ylim = c(0, 1),
       pch = 4)
  abline(a = c(0, 1),
         col = "red")
}