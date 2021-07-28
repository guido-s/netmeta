createC <- function(x,
                    sep.comps = "+",
                    inactive = NULL,
                    ncol, ncomb = 1) {
  

  if (!missing(x)) {
    meta:::chkclass(x, c("netconnection", "netmeta"))
    ##
    ## Set of all treatments
    ##
    if (class(x) == "netconnection")
      trts <- rownames(x$D.matrix)
    else
      trts <- x$trts
    ##
    trts.all <- trts
    ##
    ## Set of all components
    ##
    components <- unique(sort(unlist(compsplit(trts, sep.comps))))
    ##
    ## Inactive treatment must be empty or one of the treatments
    ##
    inactive.given <- !is.null(inactive)
    ##
    if (inactive.given) {
      inactive <- meta:::setchar(inactive, components)
      ##
      trts <- trts[trts != inactive]
      ##
      if (length(trts) == 0)
        stop("All treatments equal to reference treatment ",
             "(argument 'inactive')",
             call. = FALSE)
      ##
      ## Remove inactive treatment from list of components
      ##
      components <- components[components != inactive]
    }
    
    
    ##
    ## Create list with all treatment components
    ##
    components.list <- compsplit(trts, sep.comps)
    ##
    ## Remove blanks (at start and end)
    ##
    components.list <- lapply(components.list, meta:::rmSpace)
    components.list <- lapply(components.list, meta:::rmSpace, end = TRUE)
    
    
    ##
    ## Create C matrix
    ##
    C <- matrix(NA, nrow = length(trts), ncol = length(components))
    ##
    ## Each list element
    ## - corresponds to a treatment combination
    ## - consists of a vector with treatment components
    ##
    equal <- function(x, pattern) x == pattern
    ##
    for (i in seq(along = components)) {
      ## Logical grep whether components[i] is part of a treatment combination
      C[, i] <- unlist(lapply(lapply(components.list, equal,
                                     pattern = components[i]),
                              sum)) > 0
    }
    ##
    ## Add row for reference group (and convert to matrix with numeric values)
    ##
    if (inactive.given && inactive %in% trts.all)
      C <- rbind(C, rep(0, length(components)))
    else
      C <- 1L * C
    ##
    ## Convert to data frame
    ##
    C <- data.frame(C)
    names(C) <- components
    ##
    if (inactive.given && inactive %in% trts.all)
      rownames(C) <- c(trts, inactive)
    else
      rownames(C) <- trts
  }
  else {
    if (missing(ncol))
      stop("Either argument 'x' or 'ncol' must be provided.",
           call. = TRUE)
    ##
    meta:::chknumeric(ncomb, min = 1, max = min(ncol, 5), length = 1)
    ##
    inactive <- NULL
    ##
    nrow <- choose(ncol, ncomb)
    C <- matrix(0, nrow = nrow, ncol = ncol)
    ##
    i <- 0
    ##
    if (ncomb == 1)
      diag(C) <- 1
    else if (ncomb == 2) {
      for (pos1.i in 1:(ncol - 1)) {
        for (pos2.i in (pos1.i + 1):ncol) {
          i <- i + 1
          C[i, pos1.i] <- C[i, pos2.i] <- 1
        }
      }
    }
    else if (ncomb == 3) {
      for (pos1.i in 1:(ncol - 2)) {
        for (pos2.i in (pos1.i + 1):(ncol - 1)) {
          for (pos3.i in (pos2.i + 1):ncol) {
            i <- i + 1
            C[i, pos1.i] <- C[i, pos2.i] <- C[i, pos3.i] <- 1
          }
        }
      }
    }
    else if (ncomb == 4) {
      for (pos1.i in 1:(ncol - 3)) {
        for (pos2.i in (pos1.i + 1):(ncol - 2)) {
          for (pos3.i in (pos2.i + 1):(ncol - 1)) {
            for (pos4.i in (pos3.i + 1):ncol) {
              i <- i + 1
              C[i, pos1.i] <- C[i, pos2.i] <- C[i, pos3.i] <- C[i, pos4.i] <- 1
            }
          }
        }
      }
    }
    else if (ncomb == 5) {
      for (pos1.i in 1:(ncol - 4)) {
        for (pos2.i in (pos1.i + 1):(ncol - 3)) {
          for (pos3.i in (pos2.i + 1):(ncol - 2)) {
            for (pos4.i in (pos3.i + 1):(ncol - 1)) {
              for (pos5.i in (pos4.i + 1):ncol) {
              i <- i + 1
              C[i, pos1.i] <- C[i, pos2.i] <- C[i, pos3.i] <-
                C[i, pos4.i] <- C[i, pos5.i] <- 1
              }
            }
          }
        }
      }
    }
  }
  
  
  attr(C, "inactive") <- inactive
  ##
  C
}
