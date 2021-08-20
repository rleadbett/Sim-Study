#' This function creates a similar censoring CDF plot to fitdistrplus but takes a more conventional input type
#'
#' @Inputs
#' lifetimes A vector of lifetimes
#' censoring A binary vector; 1 if lifetime is censored and 0 otherwise
#' xlabel 
#' ylabel
#' title
#' 
#' @returns A plot of the lifetimes ordered by lifetime and divided into censored and non-censored observations
#' 

censPlot <- function(lifetimes, 
                     censoring, 
                     xlabel = "Lifetime", 
                     ylabel = "Cumulative Proportion", 
                     title = "Survival Plot"){
  
  if(length(lifetimes) != length(censoring)){stop("vectors lifetimes and censorning should be the same length")}
  
  cens <- lifetimes[censoring == 1]
  obs <- lifetimes[censoring == 0]
  
  cens <- cens[order(cens)]
  obs <- obs[order(obs)]
  
  xrange <- max(lifetimes) * 1.05
  yind <- (1:length(lifetimes)) / length(lifetimes)
  
  yindObs <- yind[1:length(obs)]
  yindCens <- yind[(length(obs) + 1): length(yind)]
  
  plot(x = c(), y = c(), xlim = c(0, xrange), ylim = c(0, 1), xlab = xlabel, ylab = ylabel, main = title)
  segments(x0 = cens, x1 = rep(xrange, length(yindCens)), y0 = yindCens, y1 = yindCens)
  points(x = obs, y = yindObs, pch = 4)
}