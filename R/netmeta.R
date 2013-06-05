netmeta <- function(TE, seTE,
                    treat1, treat2,
                    studlab, data=NULL,
                    sm="",
                    level=0.95, level.comb=0.95,
                    comb.fixed=TRUE, comb.random=FALSE,
                    reference.group="",
                    all.treatments=NULL,
                    title="",
                    warn=TRUE
                    ){
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch TE, treat1, treat2, seTE, studlab from data:
  ##
  mf <- match.call()
  mf$sm <- mf$level <- mf$level.comb <- NULL
  mf$data <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)

  TE <- mf$TE
  treat1 <- mf$treat1
  treat2 <- mf$treat2
  seTE <- mf$seTE
  ##
  if (length(mf$studlab)!=0){
    if (is.factor(mf$studlab))
      studlab <- as.character(mf$studlab)
  }
  else
    studlab <- seq(along=TE)
  
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  
  ##
  ## Check for levels of confidence interval
  ##
  if (!is.numeric(level) | length(level)!=1)
    stop("parameter 'level' must be a numeric of length 1")
  if (level <= 0 | level >= 1)
    stop("parameter 'level': no valid level for confidence interval")
  
  ##
  ## Check for levels of confidence interval
  ##
  if (!is.numeric(level.comb) | length(level.comb)!=1)
    stop("parameter 'level.comb' must be a numeric of length 1")
  if (level.comb <= 0 | level.comb >= 1)
    stop("parameter 'level.comb': no valid level.comb for confidence interval")

  if (is.null(all.treatments)){
    if (reference.group=="")
      all.treatments <- TRUE
    else
      all.treatments <- FALSE
  }
      
  
  ## Generate ordered data set, with added numbers of arms per study
  ##
  p0 <- prepare(TE, seTE, treat1, treat2, studlab)
  ##
  ## Study overview
  ##
  tdata <- data.frame(studies=p0$studlab, narms=p0$narms)
  tdata <- unique(tdata[order(tdata$studies, tdata$narms),])
  studies <- tdata$studies
  narms <- tdata$narms                 
  ##
  ##
  ## Network meta-analysis based on prepared data set
  ##
  res <- network(p0$TE, sqrt(1/p0$w.fixed),
                 p0$treat1, p0$treat2,
                 p0$treat1.pos, p0$treat2.pos,
                 p0$narms, p0$studlab,
                 sm=sm,
                 level=level, level.comb=level.comb,
                 comb.fixed=comb.fixed, comb.random=comb.random)

  o <- order(p0$order)
  ##
  res$studlab <- res$studlab[o]
  res$treat1 <- res$treat1[o]
  res$treat2 <- res$treat2[o]
  ##
  res$TE <- res$TE[o]
  res$seTE <- res$seTE[o]
  ##
  res$TE.nma.fixed <- res$TE.nma.fixed[o]
  res$seTE.nma.fixed <- res$seTE.nma.fixed[o]
  res$lower.nma.fixed <- res$lower.nma.fixed[o]
  res$upper.nma.fixed <- res$upper.nma.fixed[o]
  ##
  res$leverage.fixed <- res$leverage.fixed[o]
  res$w.fixed <- res$w.fixed[o]
  res$w.random <- res$w.random[o]
  ##
  res$Q.fixed <- res$Q.fixed[o]
  ##
  res$treat1.pos <- res$treat1.pos[o]
  res$treat2.pos <- res$treat2.pos[o]
  ##
  res$studies <- studies
  res$narms <- narms
  ##
  res$G.matrix <- res$G.matrix[o,o]
  res$H.matrix <- res$H.matrix[o,o]
  
  res$all.treatments <- all.treatments
  res$reference.group <- reference.group
  ##
  res$title <- title
  ##
  res$warn <- warn
  res$version <- packageDescription("netmeta")$Version
  
  class(res) <- "netmeta"
  
  res
}
