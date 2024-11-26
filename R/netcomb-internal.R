get_unident <- function(X.matrix, details.chkident) {
  comps.unident <- character(0)
  #
  Xplus <- ginv(X.matrix)
  colnames(Xplus) <- rownames(X.matrix)
  rownames(Xplus) <- colnames(X.matrix)
  e <- eigen(Xplus %*% X.matrix)$values
  E <- eigen(Xplus %*% X.matrix)$vectors
  rownames(E) <- rownames(Xplus %*% X.matrix)
  M <- as.matrix(E[, is.zero(e, n = 100)])
  #
  # Identify all unidentifiable components
  #
  if (dim(M)[2] > 0) {
    for (i in 1:dim(M)[2])
      comps.unident <-
        c(comps.unident, names(M[, i])[!is.zero(M[, i], n = 100)])
    #
    comps.unident <- unique(sort(comps.unident))
  }
  #
  # Identify components only appearing in combination
  #
  Xid <- X.matrix[, is_identical(X.matrix), drop = FALSE]
  #
  ncols <- ncol(Xid)
  list.id <- vector("list")
  j <- 0
  #
  while (ncols > 1) {
    sel.i <- vector()
    for (i in 1 + seq_len(ncols - 1)) {
      if (all(Xid[, 1] == Xid[, i]))
        sel.i <- c(sel.i, i)
    }
    #
    if (length(sel.i) > 0) {
      j <- j + 1
      list.id[[j]] <- colnames(Xid[, c(1, sel.i), drop = FALSE])
      Xid <- Xid[, -c(1, sel.i), drop = FALSE]
      ncols <- ncol(Xid)
    }
  }
  #
  if (length(list.id) > 0)
    warning(
      paste0("The following components are not uniquely identifiable as ",
             "they always appear together in a treatment: ",
             paste0("'", lapply(list.id, paste, collapse = " + "),
                    "'", collapse = "; "),
             ". The combined effect of these components is identifiable ",
             "by creating a combined component."),
      call. = FALSE)
  #
  # Warn about remaining unidentifiable components
  #
  if (dim(M)[2] > 0) {
    if (length(list.id) == 0)
      comps.unident.rest <- comps.unident
    else
      comps.unident.rest <-
        comps.unident[!(comps.unident %in% unique(unlist(list.id)))]
    #
    if (length(comps.unident.rest) > 0)
      warning(paste0("The following component",
                     if (length(comps.unident.rest) > 1)
                       "s are " else " is ",
                     if (length(list.id) > 0) "also ",
                     "not uniquely identifiable: ",
                     paste(paste0("'", comps.unident.rest, "'"),
                           collapse = ", ")),
              call. = FALSE)
    #
    if (details.chkident) {
      M[is.zero(M, n = 100)] <- 0
      prmatrix(M, quote = FALSE, right = TRUE)
    }
  }
  #
  comps.unident
}
