createC <- function(x,
                    sep.components = "+",
                    inactive = NULL) {
  
  
  meta:::chkclass(x, c("netconnection", "netmeta"))
  ##
  ## Set of all treatments
  ##
  if (class(x) == "netconnection")
    trts <- rownames(x$D.matrix)
  else
    trts <- rownames(x$TE.fixed)
  ##
  ## Inactive treatment must be empty or one of the treatments
  ##
  inactive.given <- !is.null(inactive)
  ##
  if (inactive.given) {
    inactive <- meta:::setchar(inactive, trts)
    ##
    trts <- trts[trts != inactive]
    if (length(trts) == 0)
      stop("All treatments equal to reference treatment (argument 'inactive')",
           call. = FALSE)
  }
  
  
  ##
  ## Create list with all treatment components
  ##
  if (length(sep.components) != 1 || !is.character(sep.components) ||
      nchar(sep.components) != 1)
    stop(paste("Argument '", "sep.components",
               "' must be a single character.", sep = ""))
  ##
  components.list <- compsplit(trts, sep.components)
  ##
  ## Remove blanks (at start and end)
  ##
  components.list <- lapply(components.list, meta:::rmSpace)
  components.list <- lapply(components.list, meta:::rmSpace, end = TRUE)
  
  
  ##
  ## Determine treatment components of interest
  ##
  components <- unique(sort(unlist(components.list)))
  ##
  ## Remove inactive treatment from list of components
  ##
  if (inactive.given)
    components <- components[components != inactive]
  
  
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
  if (inactive.given)
    C <- rbind(C, rep(0, length(components)))
  else
    C <- 1L * C
  ##
  ## Convert to data frame
  ##
  C <- data.frame(C)
  names(C) <- components
  ##
  if (inactive.given)
    rownames(C) <- c(trts, inactive)
  else
    rownames(C) <- trts
  
  
  C
}
