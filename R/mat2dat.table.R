mat2dat.table <- function(x, pooled = "common", dat.trts,
                          backtransf,
                          digits, digits.I2, digits.pval,
                          scientific.pval, zero.pval, JAMA.pval,
                          big.mark, text.NA,
                          writexl) {
  
  pooled <- setchar(pooled, c("common", "random"))
  ##
  name.n <- paste0("nma.", pooled)
  name.d <- paste0("direct.", pooled)
  name.i <- paste0("indirect.", pooled)
  
  
  res <-
    data.frame(treat1 = rep_len(NA, nrow(dat.trts)),
               treat2 = NA,
               k = NA, n = NA, I2 = NA,
               direct = NA, indirect = NA, nma = NA,
               p.value = NA,
               ##
               TE.direct = NA,
               seTE.direct = NA,
               lower.direct = NA,
               upper.direct = NA,
               TE.indirect = NA,
               seTE.indirect = NA,
               lower.indirect = NA,
               upper.indirect = NA,
               TE.nma = NA,
               lower.nma = NA,
               upper.nma = NA,
               ##
               stringsAsFactors = FALSE)
  ##
  for (i in seq_len(nrow(dat.trts))) {
    t1.i <- dat.trts$treat1[i]
    t2.i <- dat.trts$treat2[i]
    ##
    res$treat1[i] <- t1.i
    res$treat2[i] <- t2.i
    ##
    res$k[i] <- x$x$A.matrix[t1.i, t2.i]
    if (!is.null(x$x$n.matrix))
      res$n[i] <- x$x$n.matrix[t1.i, t2.i]
    res$I2[i] <- x$x$I2.direct[t1.i, t2.i]
    ##
    res$TE.direct[i] <- x[[name.d]]$TE[t1.i, t2.i]
    res$seTE.direct[i] <- x[[name.d]]$seTE[t1.i, t2.i]
    res$lower.direct[i] <- x[[name.d]]$lower[t1.i, t2.i]
    res$upper.direct[i] <- x[[name.d]]$upper[t1.i, t2.i]
    ##
    res$TE.indirect[i] <- x[[name.i]]$TE[t1.i, t2.i]
    res$seTE.indirect[i] <- x[[name.i]]$seTE[t1.i, t2.i]
    res$lower.indirect[i] <- x[[name.i]]$lower[t1.i, t2.i]
    res$upper.indirect[i] <- x[[name.i]]$upper[t1.i, t2.i]
    ##
    res$TE.nma[i] <- x[[name.n]]$TE[t1.i, t2.i]
    res$lower.nma[i] <- x[[name.n]]$lower[t1.i, t2.i]
    res$upper.nma[i] <- x[[name.n]]$upper[t1.i, t2.i]
  }
  ##
  res$p.value <-
    suppressWarnings(
      metagen(res$TE.direct - res$TE.indirect,
              sqrt(res$seTE.direct^2 + res$seTE.indirect^2),
              level = x$x$level.ma,
              method.tau = "DL", method.tau.ci = "")$pval)
  res$p.value <-
    formatPT(res$p.value, digits = digits.pval,
             scientific = scientific.pval,
             zero = zero.pval, JAMA = JAMA.pval,
             lab.NA = text.NA)
  ##
  if (backtransf & is.relative.effect(x$x$sm)) {
    res$TE.nma <- exp(res$TE.nma)
    res$lower.nma <- exp(res$lower.nma)
    res$upper.nma <- exp(res$upper.nma)
    ##
    res$TE.direct <- exp(res$TE.direct)
    res$lower.direct <- exp(res$lower.direct)
    res$upper.direct <- exp(res$upper.direct)
    ##
    res$TE.indirect <- exp(res$TE.indirect)
    res$lower.indirect <- exp(res$lower.indirect)
    res$upper.indirect <- exp(res$upper.indirect)
  }
  ##
  ## Round and round ...
  ##
  res$TE.nma <- round(res$TE.nma, digits = digits)
  res$lower.nma <- round(res$lower.nma, digits = digits)
  res$upper.nma <- round(res$upper.nma, digits = digits)
  ##
  res$TE.direct <- round(res$TE.direct, digits = digits)
  res$lower.direct <- round(res$lower.direct, digits = digits)
  res$upper.direct <- round(res$upper.direct, digits = digits)
  ##
  res$TE.indirect <- round(res$TE.indirect, digits = digits)
  res$lower.indirect <- round(res$lower.indirect, digits = digits)
  res$upper.indirect <- round(res$upper.indirect, digits = digits)
  ##
  ## Format results
  ##
  res$nma <-
    ifelse(is.na(res$TE.nma),
           text.NA,
           paste(formatN(res$TE.nma, digits = digits,
                         text.NA = text.NA, big.mark = big.mark),
                 formatCI(formatN(res$lower.nma,
                                  digits = digits,
                                  text.NA = text.NA, big.mark = big.mark),
                          formatN(res$upper.nma,
                                  digits = digits,
                                  text.NA = text.NA, big.mark = big.mark))))
  ##
  res$direct <-
    ifelse(is.na(res$TE.direct),
           text.NA,
           paste(formatN(res$TE.direct, digits = digits,
                         text.NA = text.NA, big.mark = big.mark),
                 formatCI(formatN(res$lower.direct,
                                  digits = digits,
                                  text.NA = text.NA, big.mark = big.mark),
                          formatN(res$upper.direct,
                                  digits = digits,
                                  text.NA = text.NA, big.mark = big.mark))))
  ##
  res$indirect <-
    ifelse(is.na(res$TE.indirect),
           text.NA,
           paste(formatN(res$TE.indirect, digits = digits,
                         text.NA = text.NA, big.mark = big.mark),
                 formatCI(formatN(res$lower.indirect,
                                  digits = digits,
                                  text.NA = text.NA, big.mark = big.mark),
                          formatN(res$upper.indirect,
                                  digits = digits,
                                  text.NA = text.NA, big.mark = big.mark))))
  ##
  res$I2 <- round(100 * res$I2, digits = digits.I2)
  res$I2 <-
    ifelse(is.na(res$I2),
           text.NA,
           paste0(formatN(res$I2, digits = digits.I2), "%"))
  res$direct[is.na(res$direct)] <- text.NA
  res$indirect[is.na(res$indirect)] <- text.NA
  res$p.value[is.na(res$p.value)] <- text.NA
  ##
  ## Drop unnecessary columns
  ##
  res$TE.direct <- res$seTE.direct <-
    res$lower.direct <- res$upper.direct <-
      res$TE.indirect <- res$seTE.indirect <-
        res$lower.indirect <- res$upper.indirect <-
          res$TE.nma <- res$lower.nma <- res$upper.nma <- NULL
  ##
  ## Change variable names
  ##
  nam <- c("Arm 1", "Arm 2", "k", "n", "I2",
           "Direct estimate", "Indirect estimate",
           "Network meta-analysis",
           if (writexl) "Incoherence P-value" else "Incoherence")
  names(res) <- nam
  ##
  if (!is.null(x$x$outcome.name)) {
    res$Outcome <- x$x$outcome.name
    res <- res[, c("Outcome", nam)]
  }
  ##
  res
}
