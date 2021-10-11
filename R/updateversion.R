updateversion <- function(x) {
  
  if (is.null(x$version)) {
    major <- 0
    minor <- 0
  }
  else {
    version <- as.numeric(unlist(strsplit(x$version, "-")))
    major <- version[1]
    minor <- version[2]
  }
  ##
  update.0.9.6 <- FALSE
  update.0.9.7 <- FALSE
  update.1.3.0 <- FALSE
  update.2.0.0 <- FALSE
  ##
  if (!((major == 0.9 & minor > 5) | major > 0.9)) {
    update.0.9.6 <- TRUE
    update.0.9.7 <- TRUE
  }
  if (!((major == 0.9 & minor > 6) | major > 0.9))
    update.0.9.7 <- TRUE
  ##
  if (major < 1.3)
    update.1.3.0 <- TRUE
  ##
  if (major < 2.0)
    update.2.0.0 <- TRUE
  
  
  ##
  ##  (1) Update netmeta object
  ##
  if (inherits(x, "netmeta")) {
    ##
    if (update.0.9.6) {
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
    ##
    if (update.0.9.7)
      x$backtransf <- TRUE
    ##
    if (update.1.3.0) {
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
      ##
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
        ##
        x$n.arms <- tdata12$narms
        x$multiarm <- tdata12$narms > 2
      }
      else {
        x$n.arms <- rep(2, length(x$studlab))
        x$multiarm <- rep(FALSE, length(x$studlab))
      }
    }
    ##
    if (update.2.0.0) {
      x$fixed <- x$comb.fixed
      x$random <- x$comb.random
      x$level.ma <- x$level.comb
      ##
      x$comb.fixed <- x$comb.random <- x$level.comb <- NULL
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (2) Update summary.netmeta / summary.netcomb object
  ##
  if (inherits(x, c("summary.netmeta", "summary.netcomb"))) {
    if (update.2.0.0) {
      x$level.ma <- x$level.comb
      x$x$fixed <- x$comb.fixed
      x$x$random <- x$comb.random
      ##
      x$comb.fixed <- x$comb.random <- x$level.comb <- NULL
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (3) Update netcomb object
  ##
  if (inherits(x, "netcomb") && !inherits(x, "discomb")) {
    if (update.1.3.0) {
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
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    if (update.2.0.0) {
      x$fixed <- x$comb.fixed
      x$random <- x$comb.random
      x$level.ma <- x$level.comb
      ##
      x$comb.fixed <- x$comb.random <- x$level.comb <- NULL
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (4) Update discomb object
  ##
  if (inherits(x, "discomb")) {
    if (update.1.3.0) {
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
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    if (update.2.0.0) {
      x$fixed <- x$comb.fixed
      x$random <- x$comb.random
      x$level.ma <- x$level.comb
      ##
      x$comb.fixed <- x$comb.random <- x$level.comb <- NULL
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (5) Update netsplit object
  ##
  if (inherits(x, "netsplit")) {
    if (update.1.3.0) {
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
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    if (update.2.0.0) {
      x$x$fixed <- x$comb.fixed
      x$x$random <- x$comb.random
      x$level.ma <- x$level.comb
      ##
      x$comb.fixed <- x$comb.random <- x$level.comb <- NULL
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (6) Update netrank object
  ##
  if (inherits(x, "netrank")) {
    if (update.2.0.0) {
      x$x <- updateversion(x$x)
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (7) Update rankogram object
  ##
  if (inherits(x, "rankogram")) {
    if (update.2.0.0) {
      if (is.null(x$cumulative.rankprob))
        x$cumulative.rankprob <- FALSE
      if (is.null(x$nchar.trts))
        x$nchar.trts <- 666
      ##
      x$fixed <- x$comb.fixed
      x$random <- x$comb.random
      ##
      x$comb.fixed <- x$comb.random <- NULL
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (8) Update netimpact object
  ##
  if (inherits(x, "netimpact")) {
    if (update.2.0.0) {
      x$x <- updateversion(x$x)
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (9) Update netbind object
  ##
  if (inherits(x, "netbind")) {
    if (update.2.0.0) {
      x$x$fixed <- x$comb.fixed
      x$x$random <- x$comb.random
      x$x$level.ma <- x$level.comb
      ##
      x$comb.fixed <- x$comb.random <- x$level.comb <- NULL
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (10) Update netposet object
  ##
  if (inherits(x, "netposet")) {
    if (update.2.0.0) {
      x$fixed <- x$comb.fixed
      x$random <- x$comb.random
      ##
      x$comb.fixed <- x$comb.random <- NULL
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    return(x)
  }
  
  
  ##
  ##  (11) Update netcontrib object
  ##
  if (inherits(x, "netcontrib")) {
    if (update.2.0.0) {
      x$x$fixed <- x$comb.fixed
      x$x$random <- x$comb.random
      ##
      x$comb.fixed <- x$comb.random <- NULL
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    return(x)
  }

  x
}
