catmeth <- function(x, random, text.tau2, digits.tau2, big.mark) {
  text.details <- ""
  #
  if (!inherits(x, "netmetabin"))
    text.details <-
      "- Frequentist graph-theoretical approach\n"
  else {
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
  #
  if (random) {
    if (!is.null(x$tau.preset)) {
      tau2 <- x$tau.preset^2
      tau2 <- formatPT(tau2, lab = TRUE, labval = text.tau2,
                       digits = digits.tau2,
                       lab.NA = "NA", big.mark = big.mark)
      #
      text.details <- paste0(text.details,
                             "- Preset between-study variance: ", tau2, "\n")
    }
    else {
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
      text.details <- paste0(text.details, " for ", text.tau2, "\n")
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
