createC <- function(x,
                    sep.comps = "+",
                    inactive = NULL,
                    ncol, ncomb = 1) {
  

  if (!missing(x)) {
    chkclass(x, c("netconnection", "netmeta"))
    ##
    ## Set of all treatments
    ##
    if (inherits(x, "netconnection"))
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
      inactive <- setchar(inactive, components)
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
    components.list <- lapply(components.list, rmSpace)
    components.list <- lapply(components.list, rmSpace, end = TRUE)
    
    
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
    chknumeric(ncomb, min = 1, max = ncol, length = 1)
    ##
    inactive <- NULL
    ##
    C <- createC.full(ncol, ncomb)
  }
  
  
  attr(C, "inactive") <- inactive
  ##
  C
}





createC.full <- function(n, k) {
  if (k > n)
    stop("Error: k may not exceed n")
  else if (k == 0)
    res <- as.matrix(t(rep(0, n)))
  else if (k == n)
    res <- as.matrix(t(rep(1, n)))
  else
    res <- rbind(cbind(rep(1, choose(n - 1, k - 1)),
                       createC.full(n - 1, k - 1)),
                 cbind(rep(0, choose(n - 1, k)),
                       createC.full(n - 1, k)))
  ##
  res
}
