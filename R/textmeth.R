textmeth <- function(x, random = FALSE, print.tau2 = FALSE, print.tau = FALSE,
                     text.tau2 = gs("text.tau2"), text.tau = gs("text.tau"),
                     digits.tau2 = gs("digits.tau2"),
                     digits.tau = gs("digits.tau"),
                     print.I2 = FALSE, text.I2 = gs("text.I2"),
                     big.mark = gs("big.mark"), forest = FALSE) {
  text.details <- ""
  #
  if (print.tau2 & print.tau)
    print.tau <- FALSE
  #
  if (forest) {
    text.tau2 <- "tau^2"
    text.tau <- "tau"
    #
    text.I2 <- "I^2"
  }
  #
  if (inherits(x, c("netmeta", "netcomb"))) {
    if (inherits(x, "netmetabin")) {
      if (x$method == "MH") {
        text.details <- "- Mantel-Haenszel method\n"
        #
        if (x$cc.pooled & x$incr != 0)
          text.details <-
            paste0(text.details,
                   paste("- Continuity correction of", x$incr),
                   "\n")
      }
      else if (x$method == "NCH")
        text.details <-
          "- Based on the non-central hypergeometric distribution"
    }
    else {
      text.details <- "- Frequentist graph-theoretical approach\n"
      if (inherits(x, "netcomb"))
        text.details <-
          paste0(text.details, "- Component network meta-analysis\n")
      #
      if (random | print.tau2 | print.tau | print.I2) {
        if (!is.null(x$tau.preset)) {
          if (print.tau2) {
            tau2 <- x$tau.preset^2
            tau2 <- formatPT(tau2, lab = TRUE, labval = text.tau2,
                             digits = digits.tau2,
                             lab.NA = "NA", big.mark = big.mark)
            #
            text.details <-
              paste0(text.details,
                     "- Preset between-study variance: ", tau2, "\n")
          }
          else if (print.tau) {
            tau <- x$tau.preset
            tau <- formatPT(tau, lab = TRUE, labval = text.tau,
                            digits = digits.tau,
                            lab.NA = "NA", big.mark = big.mark)
            #
            text.details <-
              paste0(text.details,
                     "- Preset between-study standard deviation: ", tau, "\n")
          }
        }
        else if (random | print.tau2 | print.tau) {
          text.details <-
            paste0(text.details,
                   if (x$method.tau == "DL")
                     "- DerSimonian-Laird estimator"
                   #
                   else if (x$method.tau == "REML")
                     "- Restricted maximum-likelihood estimator"
                   ##
                   else if (x$method.tau == "ML")
                     "- Maximum-likelihood estimator")
          #
          text.details <-
            paste0(text.details, " for ",
                   if (print.tau) text.tau else text.tau2, "\n")
        }
      }
      #
      if (print.I2) {
        text.details <-
          paste0(text.details, "- Calculation of ", text.I2, " based on Q\n")
      }
    }
  }
  else if (inherits(x, "netimpact")) {
    text.details <- "- Frequentist graph-theoretical approach\n"
    #
    text.details <-
      paste0(text.details,
             if (x$method.tau == "DL")
               "- DerSimonian-Laird estimator"
             #
             else if (x$method.tau == "REML")
               "- Restricted maximum-likelihood estimator"
             ##
             else if (x$method.tau == "ML")
               "- Maximum-likelihood estimator")
    #
    text.details <-
      paste0(text.details, " for ",
             if (print.tau) text.tau else text.tau2, "\n")
    #
    if (print.I2) {
      text.details <-
        paste0(text.details, "- Calculation of ", text.I2, " based on Q\n")
    }
  }
  else if (inherits(x, "rankogram")) {
    if (random && !is.null(x$x)) {
      if (!is.null(x$x$tau.preset)) {
          tau2 <- formatPT(x$x$tau.preset^2,
                           lab = TRUE, labval = text.tau2,
                           digits = digits.tau2,
                           lab.NA = "NA", big.mark = big.mark)
          #
          text.details <-
            paste0(text.details,
                   "- Preset between-study variance: ", tau2, "\n")
      }
      else {
        text.details <-
          paste0(text.details,
                 if (x$x$method.tau == "DL")
                   "- DerSimonian-Laird estimator"
                 #
                 else if (x$x$method.tau == "REML")
                   "- Restricted maximum-likelihood estimator"
                 ##
                 else if (x$x$method.tau == "ML")
                   "- Maximum-likelihood estimator")
        #
        text.details <-
          paste0(text.details, " for ",
                 if (print.tau) text.tau else text.tau2, "\n")
      }
    }
  }
  else if (inherits(x, "subgroup.netmeta")) {
    if (random) {
      text.details <-
        paste0(text.details,
               if (x$method.tau == "DL")
                 "- DerSimonian-Laird estimator"
               #
               else if (x$method.tau == "REML")
                 "- Restricted maximum-likelihood estimator"
               ##
               else if (x$method.tau == "ML")
                 "- Maximum-likelihood estimator")
      #
      text.details <-
        paste0(text.details, " for ",
               if (print.tau) text.tau else text.tau2, "\n")
    }
  }
  #
  if (text.details != "")
    text.details <-
      paste0("\nDetails of network meta-analysis methods:\n",
             text.details)
  #
  text.details
}
