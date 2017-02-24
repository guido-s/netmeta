netposet <- function(..., outcomes, small.values) {
  
  
  args <- list(...)
  
  
  if (length(args) == 1) {
    if (!(is.matrix(args[[1]]) | is.data.frame(args[[1]])))
      stop("First argument must be a matrix or data frame if only a single argument is provided.",
           call. = FALSE)
    ##
    pscore.matrix <- args[[1]]
    ##
    if (ncol(pscore.matrix) == 1)
      stop("Minimum number of two outcomes (equal to number of columns of ranking matrix).",
           call. = FALSE)
    ##
    if (!missing(outcomes)) {
      if (length(outcomes) != ncol(pscore.matrix))
        stop("Number of outcomes provided by argument 'outcomes' is different from number of columns of ranking matrix.",
             call. = FALSE)
      colnames(pscore.matrix) <- outcomes
    }
    else
      outcomes <- colnames(pscore.matrix)
  }
  else {
    n.outcomes <- length(args)
    ##
    if (missing(outcomes)) {
      warning("Outcomes are labelled 'A' to '", LETTERS[n.outcomes],
              "' as argument 'outcomes' is missing.")
      outcomes <- LETTERS[1:n.outcomes]  
    }
    ##
    if (length(outcomes) != n.outcomes)
      stop("Number of outcomes provided by argument 'outcomes' and number of rankings differ.",
           call. = FALSE)
    ##
    pscore.list <- list()
    ##
    for (i in seq_along(args)) {
      ##
      args.i <- args[[i]]
      ##
      if (inherits(args.i, "netmeta")) {
        ##
        if (!missing(small.values))
          small.values <- meta:::setchar(small.values, c("good", "bad"))
        else {
          warning("R function netrank() called internally with argument small.values = 'good' as ")
          small.values <- rep("good", n.outcomes)
        }
        ##
        if (length(small.values) != n.outcomes)
          stop("Length of argument 'small.values' must be equal to number of outcomes.")
        ##
        pscore.list[[i]] <- netrank(args.i, small.values = small.values[i])$Pscore
      }
      else if (inherits(args.i, "netrank"))
        pscore.list[[i]] <- args.i$Pscore
      else if (is.vector(args.i))
        pscore.list[[i]] <- args.i
      else
        stop("Elements of argument '...' must be either a vector or of class 'netmeta' or 'netrank'.",
             call. = FALSE)
    }
    
    
    pscore.treatments <- lapply(pscore.list, names)
    n.treatments <- unlist(lapply(pscore.treatments, length))
    ##
    if (length(unique(n.treatments)) != 1) {
      sel.max <- seq_along(n.treatments)[n.treatments == max(n.treatments)][1]
      ##
      treatments <- pscore.treatments[[sel.max]]
      ##
      ## Act on rankings with missing treatments
      ##
      for (i in seq_along(n.treatments)[n.treatments < max(n.treatments)]) {
        treatments.i <- pscore.treatments[[i]]
        missing.i <- !(treatments %in% treatments.i)
        ##
        if (any(treatments.i != treatments[!missing.i]))
          stop("Treatment names of all rankings must be in same order.")
        ##
        pscore.i <- pscore.list[[i]]
        pscore.list[[i]] <- pscore.list[[sel.max]]
        pscore.list[[i]][treatments[!missing.i]] <- pscore.i
        pscore.list[[i]][treatments[missing.i]]  <- NA
        pscore.treatments[[i]] <- treatments
        }
    }
    else
      treatments <- pscore.treatments[[1]]
    ##
    for (i in seq_along(pscore.treatments))
      if (any(pscore.treatments[[i]] != treatments)) {
        if (all(sort(pscore.treatments[[i]]) == sort(treatments)))
          stop("Different order of treatments provided:\n ",
               paste(paste("'", treatments, "'", sep = ""),
                     collapse = " - "),
               "\n ",
               paste(paste("'", pscore.treatments[[i]], "'", sep = ""),
                     collapse = " - "),
               call. = FALSE)
        else
          stop("Different treatments provided:\n ",
               paste(paste("'", treatments, "'", sep = ""),
                     collapse = " - "),
               "\n ",
               paste(paste("'", pscore.treatments[[i]], "'", sep = ""),
                     collapse = " - "),
               call. = FALSE)
      }
    ##
    pscore.matrix <- matrix(unlist(pscore.list, use.names = FALSE),
                            ncol = length(outcomes), byrow = FALSE)
    rownames(pscore.matrix) <- treatments
    colnames(pscore.matrix) <- outcomes
  }
  
  
  n <- dim(pscore.matrix)[1]
  o <- dim(pscore.matrix)[2]
  ##
  mean.rank <- rowMeans(pscore.matrix)
  
  
  Pos <- M <- matrix(0, nrow = n, ncol = n)
  ##
  rownames(Pos) <- colnames(Pos) <- rownames(M) <- colnames(M) <- rownames(pscore.matrix)
  
  
  ## Pos[i, j] counts how many rankings judge treatment i superior to treatment j
  ## (number between 0 and o)
  ##
  for (i in 1:n)
    for (j in 1:n)
      for (k in 1:o)
        if (!is.na(pscore.matrix[i,k]) & !is.na(pscore.matrix[j,k]) &
            pscore.matrix[i, k] > pscore.matrix[j, k])
          Pos[i, j] <- Pos[i, j] + 1
  
  
  ## Calculate "full" Hasse matrix M
  ##
  M[Pos == o] <- 1
  
  
  ## Matrix with information about partial ordering
  ##
  PO <- M[1:2, ]
  PO[1, ] <- rowSums(M)
  PO[2, ] <- colSums(M)
  rownames(PO) <- c("inferior", "superior")
  
  
  ## Skipping each direct path where a path of length 2 is found
  ##
  M0 <- M - sign(M %*% M)
  
  
  res <- list(outcomes = outcomes,
              treatments = rownames(M0),
              M0 = M0, M = M, PO = PO,
              mean.rank = mean.rank,
              pscore.matrix = pscore.matrix,
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netposet"
  
  res
}
