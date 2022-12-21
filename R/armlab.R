armlab <- function(sm, backtransf, big.mark = gs("big.mark")) {
  
  lab <- ""
  
  if (sm == "SMD")
    lab <- "Std. Mean"
  ##
  else if (sm %in% c("MD", "ROM"))
    lab <- "Mean"
  ##
  else if (sm %in% c("RD", "RR"))
    lab <- "Risk"
  ##
  else if (sm %in% c("OR", "DOR"))
    lab <- "Odds"
  ##
  else if (sm == "HR")
    lab <- "Hazard"
  ##
  else if (sm == "ASD")
    lab <- "AS"
  ##
  else if (sm %in% c("IRD", "IRR"))
    lab <- "Rate"
  ##
  else if (sm == "VE")
    lab <- "Effectiveness"
  ##
  if (!backtransf) {
    if (is.relative.effect(sm))
      lab <- paste("Log", lab)
    ##
    else if (sm == "VE")
      lab <- "Log Risk"
  }
  
  lab
}
