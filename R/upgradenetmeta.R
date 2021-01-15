upgradenetmeta <- function(x) {
  
  
  if (inherits(x, "netmeta")) {
    if (!is.null(x$version)) {
      version <- as.numeric(unlist(strsplit(x$version, "-")))
      major <- version[1]
      minor <- version[2]
      ##
      if (!((major == 0.9 & minor > 5) | major > 0.9)) {
        upgrade.096 <- TRUE
        upgrade.097 <- TRUE
      }
      if (!((major == 0.9 & minor > 6) | major > 0.9)) {
        upgrade.096 <- FALSE
        upgrade.097 <- TRUE
      }
      else {
        upgrade.096 <- FALSE
        upgrade.097 <- FALSE
      }
    }
    else {
      upgrade.096 <- TRUE
      upgrade.097 <- TRUE
    }

    
    if (upgrade.096) {
      x$prediction <- FALSE
      x$df.Q <- x$df
      ##
      x$d <- nma.krahn(x)$d
      if (is.null(x$d))
        x$d <- 1
      ##
      if (x$d > 1) {
        dd <- decomp.design(x)
        x$Q.heterogeneity <- dd$Q.decomp$Q[2]
        x$Q.inconsistency <- dd$Q.decomp$Q[3]
        ##
        x$df.Q.heterogeneity <- dd$Q.decomp$df[2]
        x$df.Q.inconsistency <- dd$Q.decomp$df[3]
        ##
        x$pval.Q.heterogeneity <- dd$Q.decomp$pval[2]
        x$pval.Q.inconsistency <- dd$Q.decomp$pval[3]
      }
      else {
        x$Q.heterogeneity <- NA
        x$Q.inconsistency <- NA
        ##
        x$df.Q.heterogeneity <- NA
        x$df.Q.inconsistency <- NA
        ##
        x$pval.Q.heterogeneity <- NA
        x$pval.Q.inconsistency <- NA
      }
      ##
      x$df <- NULL
      ##
      x$baseline.reference <- TRUE
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    
    
    if (upgrade.097)
      x$backtransf <- TRUE


    if (is.null(x$multiarm)) {
      if (any(x$narms > 2)) {
        tdata1 <- data.frame(studlab = x$studlab,
                             .order = seq(along = x$studlab))
        tdata2 <- data.frame(studlab = as.character(x$studies),
                             narms = x$narms)
        ##
        tdata12 <- merge(tdata1, tdata2,
                         by = "studlab", all.x = TRUE, all.y = FALSE,
                         sort = FALSE)
        tdata12 <- tdata12[order(tdata12$.order), ]
        x$n.arms <- tdata12$narms
        x$multiarm <- tdata12$narms > 2
      }
      else {
        x$n.arms <- rep(2, length(x$studlab))
        x$multiarm <- rep(FALSE, length(x$studlab))
      }
    }
  }
  
  
  ##
  ## Use statistic.* instead of zval.*
  ##
  if (is.null(x$version) || x$version != packageDescription("netmeta")$Version)
    if (inherits(x, "netmeta") & is.null(x$statistic.fixed)) {
      x$statistic.fixed <- x$zval.fixed
      x$statistic.random <- x$zval.random
      x$statistic.direct.fixed <- x$zval.direct.fixed
      x$statistic.direct.random <- x$zval.direct.random
      x$statistic.indirect.fixed <- x$zval.indirect.fixed
      x$statistic.indirect.random <- x$zval.indirect.random
      x$statistic.nma.fixed <- x$zval.nma.fixed
      x$statistic.nma.random <- x$zval.nma.random
      ##
      x$zval.fixed <- x$zval.random <-
        x$zval.nma.fixed <- x$zval.nma.random <-
          x$zval.direct.fixed <- x$zval.direct.random <-
            x$zval.indirect.fixed <- x$zval.indirect.random <- NULL
    }
    else if (inherits(x, c("discomb", "netcomb")) &
             is.null(x$statistic.fixed)) {
      x$statistic.fixed <- x$zval.fixed
      x$statistic.random <- x$zval.random
      x$statistic.nma.fixed <- x$zval.nma.fixed
      x$statistic.nma.random <- x$zval.nma.random
      x$statistic.cnma.fixed <- x$zval.cnma.fixed
      x$statistic.cnma.random <- x$zval.cnma.random
      x$statistic.Comb.fixed <- x$zval.Comb.fixed
      x$statistic.Comb.random <- x$zval.Comb.random
      x$statistic.Comp.fixed <- x$zval.Comp.fixed
      x$statistic.Comp.random <- x$zval.Comp.random
      ##
      x$zval.fixed <- x$zval.random <-
        x$zval.nma.fixed <- x$zval.nma.random <-
          x$zval.cnma.fixed <- x$zval.cnma.random <-
            x$zval.Comb.fixed <- x$zval.Comb.random <-
              x$zval.Comp.fixed <- x$zval.Comp.random <- NULL
    }
    else if (inherits(x, "netsplit") &
             is.null(x$fixed$statistic) & !is.null(x$fixed$z)) {
      nam <- names(x$fixed)
      names(x$fixed)[nam == "z"] <- "statistic"
      nam <- names(x$random)
      names(x$random)[nam == "z"] <- "statistic"
      ##
      nam <- names(x$direct.fixed)
      names(x$direct.fixed)[nam == "z"] <- "statistic"
      nam <- names(x$direct.random)
      names(x$direct.random)[nam == "z"] <- "statistic"
      ##
      nam <- names(x$indirect.fixed)
      names(x$indirect.fixed)[nam == "z"] <- "statistic"
      nam <- names(x$indirect.random)
      names(x$indirect.random)[nam == "z"] <- "statistic"
      ##
      nam <- names(x$compare.fixed)
      names(x$compare.fixed)[nam == "z"] <- "statistic"
      nam <- names(x$compare.random)
      names(x$compare.random)[nam == "z"] <- "statistic"
    }
    else if (inherits(x, "netrank"))
      x$x <- upgradenetmeta(x$x)
  
  
  x
}
