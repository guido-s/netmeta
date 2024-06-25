createC_full <- function(n, k) {
  if (k > n)
    stop("Error: k may not exceed n")
  else if (k == 0)
    res <- as.matrix(t(rep(0, n)))
  else if (k == n)
    res <- as.matrix(t(rep(1, n)))
  else
    res <- rbind(cbind(rep(1, choose(n - 1, k - 1)),
                       createC_full(n - 1, k - 1)),
                 cbind(rep(0, choose(n - 1, k)),
                       createC_full(n - 1, k)))
  #
  res
}


createC_trts_inactive <- function(trts, inactive = NULL, sep.comps) {
  trts.all <- trts
  #
  # Set of all components
  #
  components <- unique(sort(unlist(compsplit(trts, sep.comps))))
  #
  # Inactive treatment must be empty or one of the treatments
  #
  avail.inactive <- !is.null(inactive)
  #
  if (avail.inactive) {
    inactive <- setchar(inactive, components)
    #
    trts <- trts[trts != inactive]
    #
    if (length(trts) == 0)
      stop("All treatments equal to reference treatment ",
           "(argument 'inactive')",
           call. = FALSE)
    #
    # Remove inactive treatment from list of components
    #
    components <- components[components != inactive]
  }
  #
  # Create list with all treatment components
  #
  components.list <- compsplit(trts, sep.comps)
  #
  # Remove blanks (at start and end)
  #
  components.list <- lapply(components.list, rmSpace)
  components.list <- lapply(components.list, rmSpace, end = TRUE)
  #
  # Create C matrix (without interaction terms)
  #
  C.matrix <- matrix(NA, nrow = length(trts), ncol = length(components))
  #
  # Each list element
  # - corresponds to a treatment combination
  # - consists of a vector with treatment components
  #
  equal <- function(x, pattern) x == pattern
  #
  for (i in seq(along = components)) {
    # Logical grep whether components[i] is part of a treatment combination
    C.matrix[, i] <- unlist(lapply(lapply(components.list, equal,
                                          pattern = components[i]),
                                   sum)) > 0
  }
  #
  # Add row for reference group (and convert to matrix with numeric values)
  #
  if (avail.inactive && inactive %in% trts.all)
    C.matrix <- rbind(C.matrix, rep(0, length(components)))
  else
    C.matrix <- 1L * C.matrix
  #
  colnames(C.matrix) <- components
  #
  if (avail.inactive && inactive %in% trts.all)
    rownames(C.matrix) <- c(trts, inactive)
  else
    rownames(C.matrix) <- trts
  #
  attr(C.matrix, "inactive") <- inactive
  #
  C.matrix
}
