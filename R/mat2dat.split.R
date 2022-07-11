mat2dat.split <- function(x, pooled = "common", dat.trts) {

  pooled <- setchar(pooled, c("common", "random"))
  ##
  name.n <- paste0("nma.", pooled)
  name.d <- paste0("direct.", pooled)
  name.i <- paste0("indirect.", pooled)
  
  
  nma <- direct <- indirect <-
    data.frame(comparison = dat.trts$comparison,
               TE = NA, seTE = NA, lower = NA, upper = NA,
               statistic = NA, p = NA)
  ##
  direct$I2 <- direct$tau <- direct$tau2 <- direct$Q <-
    direct$n <- NA
  ##
  predict <- data.frame(comparison = dat.trts$comparison,
                        lower = NA, upper = NA)
  
  
  k <- nma$TE
  ##
  for (i in seq_len(nrow(dat.trts))) {
    t1.i <- dat.trts$treat1[i]
    t2.i <- dat.trts$treat2[i]
    ##
    k[i] <- x$x$A.matrix[t1.i, t2.i]
    ##
    nma$TE[i] <- x[[name.n]]$TE[t1.i, t2.i]
    nma$seTE[i] <- x[[name.n]]$seTE[t1.i, t2.i]
    nma$lower[i] <- x[[name.n]]$lower[t1.i, t2.i]
    nma$upper[i] <- x[[name.n]]$upper[t1.i, t2.i]
    nma$statistic[i] <- x[[name.n]]$statistic[t1.i, t2.i]
    nma$p[i] <- x[[name.n]]$p[t1.i, t2.i]
    ##
    direct$TE[i] <- x[[name.d]]$TE[t1.i, t2.i]
    direct$seTE[i] <- x[[name.d]]$seTE[t1.i, t2.i]
    direct$lower[i] <- x[[name.d]]$lower[t1.i, t2.i]
    direct$upper[i] <- x[[name.d]]$upper[t1.i, t2.i]
    direct$statistic[i] <- x[[name.d]]$statistic[t1.i, t2.i]
    direct$p[i] <- x[[name.d]]$p[t1.i, t2.i]
    ##
    if (!is.null(x$x$n.matrix))
      direct$n[i] <- x$x$n.matrix[t1.i, t2.i]
    direct$Q[i] <- x$x$Q.direct[t1.i, t2.i]
    direct$tau2[i] <- x$x$tau2.direct[t1.i, t2.i]
    direct$tau[i] <- x$x$tau.direct[t1.i, t2.i]
    direct$I2[i] <- x$x$I2.direct[t1.i, t2.i]
    ##
    indirect$TE[i] <- x[[name.i]]$TE[t1.i, t2.i]
    indirect$seTE[i] <- x[[name.i]]$seTE[t1.i, t2.i]
    indirect$lower[i] <- x[[name.i]]$lower[t1.i, t2.i]
    indirect$upper[i] <- x[[name.i]]$upper[t1.i, t2.i]
    indirect$statistic[i] <- x[[name.i]]$statistic[t1.i, t2.i]
    indirect$p[i] <- x[[name.i]]$p[t1.i, t2.i]
    ##
    if (pooled == "random") {
      predict$lower[i] <- x$x$lower.predict[t1.i, t2.i]
      predict$upper[i] <- x$x$upper.predict[t1.i, t2.i]
    }
  }
  ##
  m <-
    suppressWarnings(metagen(direct$TE - indirect$TE,
                             sqrt(direct$seTE^2 +
                                  indirect$seTE^2),
                             level = x$x$level.ma,
                             method.tau = "DL", method.tau.ci = ""))
  ##
  compare <-
    data.frame(comparison = dat.trts$comparison,
               TE = m$TE, seTE = m$seTE,
               lower = m$lower, upper = m$upper,
               statistic = m$statistic, p = m$pval,
               z = m$statistic,
               stringsAsFactors = FALSE)

  res <- list(k = k, nma = nma,
              direct = direct, indirect = indirect,
              compare = compare, predict = predict)
  ##
  res
}
