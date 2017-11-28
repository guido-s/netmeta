createC <- function(x,
                    reference.group = x$reference.group,
                    sep.components = "+") {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  ## Set of all treatments
  ##
  trts <- rownames(x$TE.fixed)
  ##
  ## Reference group must be one of the treatments
  ##
  reference.group <- meta:::setchar(reference.group, trts)
  
  
  ##
  ## Remove reference group from treatment vector
  ##
  if (all(trts != reference.group))
    stop("No treatment equal to reference treatment (argument 'reference.group')",
         call. = FALSE)
  ##
  trts <- trts[trts != reference.group]
  if (length(trts) == 0)
    stop("All treatments equal to reference treatment (argument 'reference.group')",
         call. = FALSE)
  
  
  ##
  ## Create list with all treatment components
  ##
  if (length(sep.components) != 1 || !is.character(sep.components) ||
      nchar(sep.components) != 1)
    stop(paste("Argument '", "sep.components",
               "' must be a single character.", sep = ""))
  ##
  if (sep.components == "+")
    sep.components <- "\\+"
  else if (sep.components == ".")
    sep.components <- "\\."
  else if (sep.components == "&")
    sep.components <- "\\&"
  else if (sep.components == "$")
    sep.components <- "\\$"
  else if (sep.components == "#")
    sep.components <- "\\#"
  else if (sep.components == "|")
    sep.components <- "\\|"
  ##
  components.list <- strsplit(trts, sep.components)
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
  ## Remove reference treatment from list of components
  ##
  components <- components[components != reference.group]
  

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
  C <- rbind(C, rep(0, length(components)))
  ##
  ## Convert to data frame
  ##
  C <- data.frame(C)
  names(C) <- components
  rownames(C) <- c(trts, reference.group)
  
  
  C
}
