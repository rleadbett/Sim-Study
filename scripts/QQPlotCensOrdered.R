#' Generates a Q-Q plot for the censored data (using ordered CDF) and the true weibull distribution
#' 
#' @param df A data frame with one column `lifetime` containing 
#' the lifetimes and another column `cens` containing the censoring
#' indicator for each of the lifetimes.
#' @param shape The shape parameter of the true Weibull distribution
#' @param scale The scale parameter of the true Weibull distribution
#' 

QQPlotCensOrdered <- function(df, shape = 1.15, scale = 5.254){
  
  df <- df[order(df$Time), ]
  df$cdf <- seq(1, nrow(df)) / nrow(df)
  
  # Get the exposure times of failures
  ts <- df$Time[df$cens == 0]
  
  # Convert this to the CDF
  cdfEst <- df$cdf[df$cens == 0]
  
  # Calculate the true quantiles at each exposure
  trueQuantiles <- pweibull(ts, shape = shape, scale = scale)
  
  # Generate Q-Q Plot
  plot(x = cdfEst, y = trueQuantiles,
       main = "Q-Q Plot (Ordered)",
       xlab = "Data",
       ylab = "True Weibull",
       xlim = c(0, 1),
       ylim = c(0, 1),
       pch = 4)
  abline(a = c(0, 1),
         col = "red")
}