updateversion <- function(x, verbose = FALSE) {
  
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
  ##
  if (!((major == 0.9 & minor > 5) | major > 0.9)) {
    update.0.9.6 <- TRUE
    update.0.9.7 <- TRUE
  }
  if (!((major == 0.9 & minor > 6) | major > 0.9))
    update.0.9.7 <- TRUE
  ##
  if (verbose) {
    if (update.0.9.6)
      message("Update to netmeta, version 0.9-6")
    if (update.0.9.7)
      message("Update to netmeta, version 0.9-7")
  }
  ##
  update.1.3.0 <- update_needed(x$version, 1, 3, verbose)
  update.2.0.0 <- update_needed(x$version, 2, 0, verbose)
  update.2.5.0 <- update_needed(x$version, 2, 5, verbose)
  update.2.8.0 <- update_needed(x$version, 2, 8, verbose)
  update.2.9.0 <- update_needed(x$version, 2, 9, verbose)
  
  
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
    if (update.2.5.0) {
      x$common <- x$fixed
      ##
      x$seTE.adj.common <- x$seTE.adj.fixed
      ##
      if (!inherits(x, "netmetabin")) {
        x$TE.nma.common <- x$TE.nma.fixed
        x$seTE.nma.common <- x$seTE.nma.fixed
        x$lower.nma.common <- x$lower.nma.fixed
        x$upper.nma.common <- x$upper.nma.fixed
        x$statistic.nma.common <- x$statistic.nma.fixed
        x$pval.nma.common <- x$pval.nma.fixed
        x$leverage.common <- x$leverage.fixed
        x$w.common <- x$w.fixed
        x$Q.common <- x$Q.fixed
      }
      ##
      x$TE.common <- x$TE.fixed
      x$seTE.common <- x$seTE.fixed
      x$lower.common <- x$lower.fixed
      x$upper.common <- x$upper.fixed
      x$statistic.common <- x$statistic.fixed
      x$pval.common <- x$pval.fixed
      ##
      x$prop.direct.common <- x$prop.direct.fixed
      ##
      x$TE.direct.common <- x$TE.direct.fixed
      x$seTE.direct.common <- x$seTE.direct.fixed
      x$lower.direct.common <- x$lower.direct.fixed
      x$upper.direct.common <- x$upper.direct.fixed
      x$statistic.direct.common <- x$statistic.direct.fixed
      x$pval.direct.common <- x$pval.direct.fixed
      ##
      x$TE.indirect.common <- x$TE.indirect.fixed
      x$seTE.indirect.common <- x$seTE.indirect.fixed
      x$lower.indirect.common <- x$lower.indirect.fixed
      x$upper.indirect.common <- x$upper.indirect.fixed
      x$statistic.indirect.common <- x$statistic.indirect.fixed
      x$pval.indirect.common <- x$pval.indirect.fixed
      ##
      if (!inherits(x, "netmetabin")) {
        x$L.matrix.common <- x$L.matrix.fixed
        x$Lplus.matrix.common <- x$Lplus.matrix.fixed
        x$H.matrix.common <- x$H.matrix.fixed
      }
      x$P.common <- x$P.fixed
      x$Cov.common <- x$Cov.fixed               
    }
    ##
    if (update.2.8.0)
      x$small.values <- setsv(x$small.values)
    ##
    if (update.2.9.0) {
      if (is.null(x$keepdata) || x$keepdata) {
        if (isCol(x$data, "subset"))
          sel.s <- x$data$subset
        else
          sel.s <- rep(TRUE, nrow(x$data))
        ##
        x$data$.seTE.adj.common <- NA
        x$data$.seTE.adj.random <- NA
        ##
        x$data$.seTE.adj.common[sel.s] <- x$seTE.adj.common
        x$data$.seTE.adj.random[sel.s] <- x$seTE.adj.random
      }
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
    if (update.2.5.0)
      x$x$common <- x$x$fixed
    ##
    if (update.2.8.0)
      x$x$small.values <- setsv(x$x$small.values)
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
    if (update.2.5.0) {
      x$common <- x$fixed
      ##
      x$seTE.adj.common <- x$seTE.adj.fixed
      x$TE.nma.common <- x$TE.nma.fixed
      x$seTE.nma.common <- x$seTE.nma.fixed
      x$lower.nma.common <- x$lower.nma.fixed
      x$upper.nma.common <- x$upper.nma.fixed
      x$statistic.nma.common <- x$statistic.nma.fixed
      x$pval.nma.common <- x$pval.nma.fixed
      x$TE.cnma.common <- x$TE.cnma.fixed
      x$seTE.cnma.common <- x$seTE.cnma.fixed
      x$lower.cnma.common <- x$lower.cnma.fixed
      x$upper.cnma.common <- x$upper.cnma.fixed
      x$statistic.cnma.common <- x$statistic.cnma.fixed
      x$pval.cnma.common <- x$pval.cnma.fixed
      ##
      x$TE.common <- x$TE.fixed
      x$seTE.common <- x$seTE.fixed
      x$lower.common <- x$lower.fixed
      x$upper.common <- x$upper.fixed
      x$statistic.common <- x$statistic.fixed
      x$pval.common <- x$pval.fixed
      ##
      x$Comp.common <- x$Comp.fixed
      x$seComp.common <- x$seComp.fixed
      x$lower.Comp.common <- x$lower.Comp.fixed
      x$upper.Comp.common <- x$upper.Comp.fixed
      x$statistic.Comp.common <- x$statistic.Comp.fixed
      x$pval.Comp.common <- x$pval.Comp.fixed
      ##
      x$Comb.common <- x$Comb.fixed
      x$seComb.common <- x$seComb.fixed
      x$lower.Comb.common <- x$lower.Comb.fixed
      x$upper.Comb.common <- x$upper.Comb.fixed
      x$statistic.Comb.common <- x$statistic.Comb.fixed
      x$pval.Comb.common <- x$pval.Comb.fixed
      ##
      x$L.matrix.common <- x$L.matrix.fixed
      x$Lplus.matrix.common <- x$Lplus.matrix.fixed
      x$H.matrix.common <- x$H.matrix.fixed        
    }
    ##
    if (update.2.8.0)
      x$x$small.values <- setsv(x$x$small.values)
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
    if (update.2.5.0) {
      x$common <- x$fixed
      ##
      x$seTE.adj.common <- x$seTE.adj.fixed
      x$TE.nma.common <- x$TE.nma.fixed
      x$seTE.nma.common <- x$seTE.nma.fixed
      x$lower.nma.common <- x$lower.nma.fixed
      x$upper.nma.common <- x$upper.nma.fixed
      x$statistic.nma.common <- x$statistic.nma.fixed
      x$pval.nma.common <- x$pval.nma.fixed
      x$TE.cnma.common <- x$TE.cnma.fixed
      x$seTE.cnma.common <- x$seTE.cnma.fixed
      x$lower.cnma.common <- x$lower.cnma.fixed
      x$upper.cnma.common <- x$upper.cnma.fixed
      x$statistic.cnma.common <- x$statistic.cnma.fixed
      x$pval.cnma.common <- x$pval.cnma.fixed
      ##
      x$TE.common <- x$TE.fixed
      x$seTE.common <- x$seTE.fixed
      x$lower.common <- x$lower.fixed
      x$upper.common <- x$upper.fixed
      x$statistic.common <- x$statistic.fixed
      x$pval.common <- x$pval.fixed
      ##
      x$Comp.common <- x$Comp.fixed
      x$seComp.common <- x$seComp.fixed
      x$lower.Comp.common <- x$lower.Comp.fixed
      x$upper.Comp.common <- x$upper.Comp.fixed
      x$statistic.Comp.common <- x$statistic.Comp.fixed
      x$pval.Comp.common <- x$pval.Comp.fixed
      ##
      x$Comb.common <- x$Comb.fixed
      x$seComb.common <- x$seComb.fixed
      x$lower.Comb.common <- x$lower.Comb.fixed
      x$upper.Comb.common <- x$upper.Comb.fixed
      x$statistic.Comb.common <- x$statistic.Comb.fixed
      x$pval.Comb.common <- x$pval.Comb.fixed
      ##
      x$L.matrix.common <- x$L.matrix.fixed
      x$Lplus.matrix.common <- x$Lplus.matrix.fixed
      x$H.matrix.common <- x$H.matrix.fixed        
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
    if (update.2.5.0) {
      x$common <- x$fixed
      ##
      x$prop.common <- x$prop.fixed
      x$direct.common <- x$direct.fixed
      x$indirect.common <- x$indirect.fixed
      x$compare.common <- x$compare.fixed
    }
    ##
    if (update.2.8.0)
      x$x$small.values <- setsv(x$x$small.values)
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
    if (update.2.5.0) {
      x$common <- x$fixed
      ##
      x$ranking.common <- x$ranking.fixed
      x$Pmatrix.common <- x$Pmatrix.fixed
    }
    ##
    if (update.2.8.0)
      x$small.values <- setsv(x$small.values)
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
    if (update.2.5.0) {
      x$common <- x$fixed
      ##
      x$ranking.common <- x$ranking.fixed
      x$ranking.matrix.common <- x$ranking.matrix.fixed
      x$cumrank.matrix.common <- x$cumrank.matrix.fixed
    }
    ##
    if (update.2.8.0)
      x$small.values <- setsv(x$small.values)
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
    if (update.2.5.0) {
      x$common <- x$fixed
      x$impact.common <- x$impact.fixed
      ##
      x$version <- packageDescription("netmeta")$Version
    }              
    ##
    if (update.2.8.0)
      x$x$small.values <- setsv(x$x$small.values)
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
    if (update.2.5.0) {
      x$common <- x$fixed
      x$x$common <- x$x$fixed
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
    if (update.2.5.0) {
      x$common <- x$fixed
      ##
      x$P.common <- x$P.fixed
      x$M0.common <- x$M0.fixed
      x$M.common <- x$M.fixed
      x$O.common <- x$O.fixed
      ##
      x$version <- packageDescription("netmeta")$Version
    }              
    ##
    if (update.2.8.0)
      x$small.values <- setsv(x$small.values)
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
    if (update.2.5.0) {
      res$common <- res$fixed
      ##
      if (!is.null(res$tictoc.fixed))
        res$tictoc.common <- res$tictoc.fixed
    }
    ##
    if (update.2.8.0)
      x$x$small.values <- setsv(x$x$small.values)
    ##
    return(x)
  }
  
  
  ##
  ##  (12) Update netcomparison object
  ##
  if (inherits(x, "netcomparison")) {
    if (update.2.5.0) {
      res$common <- res$fixed
      ##
      res$TE.common <- res$TE.fixed
      res$seTE.common <- res$seTE.fixed
      res$lower.common <- res$lower.fixed
      res$upper.common <- res$upper.fixed
      res$statistic.common <- res$statistic.fixed
      res$pval.common <- res$pval.fixed
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    if (update.2.8.0)
      x$x$small.values <- setsv(x$x$small.values)
    ##
    return(x)
  }
  
  
  ##
  ##  (13) Update netcomplex object
  ##
  if (inherits(x, "netcomplex")) {
    if (update.2.5.0) {
      res$common <- res$fixed
      ##
      res$Comb.common <- res$Comb.fixed
      res$seComb.common <- res$seComb.fixed
      res$lower.Comb.common <- res$lower.Comb.fixed
      res$upper.Comb.common <- res$upper.Comb.fixed
      res$statistic.Comb.common <- res$statistic.Comb.fixed
      res$pval.Comb.common <- res$pval.Comb.fixed
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    if (update.2.8.0)
      x$x$small.values <- setsv(x$x$small.values)
    ##
    return(x)
  }
  
  
  ##
  ##  (14) Update nettable object
  ##
  if (inherits(x, "nettable")) {
    if (update.2.5.0) {
      res$common <- res$fixed
      res$x$common <- res$x$fixed
      ##
      x$version <- packageDescription("netmeta")$Version
    }
    ##
    return(x)
  }
  
  x
}


update_needed <- function(version, major = 0, minor = 0,
                          verbose = FALSE) {
  if (is.null(version)) {
    version <- 0.1
    major.cur <- 0
    minor.cur <- 1
  }
  else {
    version <- unlist(strsplit(version, "-")[1])
    major.cur <-
      as.numeric(unlist(strsplit(version, ".", fixed = TRUE))[1])
    minor.cur <-
      as.numeric(unlist(strsplit(version, ".", fixed = TRUE))[2])
  }
  ##
  res <-
    ifelse(major.cur < major,
           TRUE, ifelse(major.cur > major,
                        FALSE, minor.cur < minor))
  if (res & verbose)
    message(paste0("Update to netmeta, version ", major, ".", minor))
  ##
  res
}
