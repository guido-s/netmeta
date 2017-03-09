netposet <- function(..., outcomes, treatments, small.values,
                     comb.fixed, comb.random) {
  
  
  args <- list(...)
  
  
  if (!missing(comb.fixed))
    meta:::chklogical(comb.fixed)
  if (!missing(comb.random))
    meta:::chklogical(comb.random)
  
  
  ##
  ## First argument is expected to be a ranking matrix if only a
  ## single argument is provided
  ##
  if (length(args) == 1) {
    ##
    ## (1) Check arguments
    ##
    if (!(is.matrix(args[[1]]) | is.data.frame(args[[1]])))
      stop("First argument must be a matrix or data frame if only a single argument is provided.",
           call. = FALSE)
    ##
    if (!missing(small.values))
      warning("Argument 'small.values' is ignored as first argument is a ranking matrix.")
    ##
    pscore.matrix <- args[[1]]
    ##
    if (any(pscore.matrix > 1) | any(pscore.matrix < 0))
      stop("All elements of ranking matrix must be between 0 and 1.")
    ##
    n.outcomes <- ncol(pscore.matrix)
    n.treatments <- nrow(pscore.matrix)
    ##
    if (n.outcomes == 1)
      stop("Minimum number of two outcomes (equal to number of columns of ranking matrix).",
           call. = FALSE)
    ##
    if (!missing(outcomes)) {
      if (length(outcomes) != n.outcomes)
        stop("Number of outcomes provided by argument 'outcomes' is different from number of columns of ranking matrix.",
             call. = FALSE)
    }
    else
      outcomes <- colnames(pscore.matrix)
    ##
    if (!missing(treatments)) {
      if (length(treatments) != n.treatments)
        stop("Number of treatments provided by argument 'treatments' is different from number of rows of ranking matrix.",
             call. = FALSE)
    }
    else
      treatments <- rownames(pscore.matrix)
    ##
    ## (2) P-Score matrices
    ##
    if (is.null(outcomes)) {
      warning("Outcomes are labelled 'A' to '", LETTERS[n.outcomes],
              "' as argument 'outcomes' is missing and column names are NULL.")
      outcomes <- LETTERS[1:n.outcomes]
    }
    ##
    if (is.null(treatments)) {
      warning("Treatments are labelled 'a' to '", letters[n.treatments],
              "' as row names are NULL.")
      treatments <- letters[1:n.treatments]
    }
    ##
    rownames(pscore.matrix) <- treatments
    colnames(pscore.matrix) <- outcomes
    ##
    pscore.matrix.fixed <- pscore.matrix
    pscore.matrix.random <- pscore.matrix
    ##
    if (missing(comb.fixed))
      comb.fixed <- FALSE
    if (missing(comb.random))
      comb.random <- FALSE
  }
  ##
  ## More than one element provided in '...' (netmeta or netrank objects)
  ##
  else {
    ##
    ## (1) Check arguments
    ##
    n.outcomes <- length(args)
    ##
    if (!missing(treatments))
      warning("Argument 'treatments' ignored (only needed if first argument is a ranking matrix).")
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
    any.netmeta <- any.netrank <- FALSE
    ##
    for (i in seq_along(args)) {
      if (inherits(args[[i]], "netmeta"))
        any.netmeta <- TRUE
      else if (inherits(args[[i]], "netrank"))
        any.netrank <- TRUE
      else
        stop("All elements of argument '...' must be of class 'netmeta' or 'netrank'.",
             call. = FALSE)
      ##
      if (!missing(small.values))
        small.values <- meta:::setchar(small.values, c("good", "bad"))
      else {
        if (any.netmeta)
          warning("R function netrank() called internally with argument small.values = 'good' as argument 'small.values' is missing.")
        small.values <- rep("good", n.outcomes)
      }
      ##
      if (length(small.values) != n.outcomes)
        stop("Length of argument 'small.values' must be equal to number of outcomes.",
             call. = TRUE)
    }
    ##
    ## (2) Extract P-Scores
    ##
    pscore.list.fixed <- pscore.list.random <- list()
    comb.fixeds <- comb.randoms <- rep_len(NA, length(args))
    ##
    for (i in seq_along(args)) {
      ##
      args.i <- args[[i]]
      ##
      if (inherits(args.i, "netmeta")) {
        pscore.list.fixed[[i]] <- netrank(args.i, small.values = small.values[i])$Pscore.fixed
        pscore.list.random[[i]] <- netrank(args.i, small.values = small.values[i])$Pscore.random
        comb.fixeds[i] <- args.i$comb.fixed
        comb.randoms[i] <- args.i$comb.random
      }
      else if (inherits(args.i, "netrank")) {
        pscore.list.fixed[[i]] <- args.i$Pscore.fixed
        pscore.list.random[[i]] <- args.i$Pscore.random
        comb.fixeds[i] <- args.i$x$comb.fixed
        comb.randoms[i] <- args.i$x$comb.random
      }
    }
    
    
    pscore.treatments <- lapply(pscore.list.fixed, names)
    n.treatments <- unlist(lapply(pscore.treatments, length))
    ##
    if (length(unique(n.treatments)) != 1) {
      sel.max <- seq_along(n.treatments)[n.treatments == max(n.treatments)][1]
      ##
      treatments <- pscore.treatments[[sel.max]]
      ##
      ## Act on rankings with missing treatments
      ##
      for (j in seq_along(n.treatments)[n.treatments < max(n.treatments)]) {
        treatments.j <- pscore.treatments[[j]]
        missing.j <- !(treatments %in% treatments.j)
        ##
        if (any(treatments.j != treatments[!missing.j]))
          stop("Treatment names of all rankings must be in same order.")
        ##
        pscore.j <- pscore.list.fixed[[j]]
        pscore.list.fixed[[j]] <- pscore.list.fixed[[sel.max]]
        pscore.list.fixed[[j]][treatments[!missing.j]] <- pscore.j
        pscore.list.fixed[[j]][treatments[missing.j]]  <- NA
        ##
        pscore.j <- pscore.list.random[[j]]
        pscore.list.random[[j]] <- pscore.list.random[[sel.max]]
        pscore.list.random[[j]][treatments[!missing.j]] <- pscore.j
        pscore.list.random[[j]][treatments[missing.j]]  <- NA
        ##
        pscore.treatments[[j]] <- treatments
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
    pscore.matrix.fixed <- matrix(unlist(pscore.list.fixed,
                                         use.names = FALSE),
                                  ncol = length(outcomes), byrow = FALSE)
    pscore.matrix.random <- matrix(unlist(pscore.list.random,
                                          use.names = FALSE),
                                   ncol = length(outcomes), byrow = FALSE)
    rownames(pscore.matrix.fixed) <- treatments
    rownames(pscore.matrix.random) <- treatments
    colnames(pscore.matrix.fixed) <- outcomes
    colnames(pscore.matrix.random) <- outcomes
    ##
    text.netmeta.netrank <- if (any.netmeta) "netmeta" else ""
    text.netmeta.netrank <- if (any.netmeta & any.netrank) "netmeta and netrank"
    text.netmeta.netrank <- if (!any.netmeta & any.netrank) "netrank"
    text.netmeta.netrank <- paste(text.netmeta.netrank, "objects.")
    ##
    if (missing(comb.fixed)) {
      comb.fixed <- unique(comb.fixeds)
      if (length(comb.fixed) != 1) {
        warning("Argument 'comb.fixed' set to TRUE as different values are available in",
                text.netmeta.netrank)
        comb.fixed <- TRUE
      }
    }
    ##
    if (missing(comb.random)) {
      comb.random <- unique(comb.randoms)
      if (length(comb.random) != 1) {
        warning("Argument 'comb.random' set to TRUE as different values are available in",
                text.netmeta.netrank)
        comb.random <- TRUE
      }
    }
  } 
  
  
  n <- dim(pscore.matrix.fixed)[1]
  o <- dim(pscore.matrix.fixed)[2]
  
  
  Pos.fixed <- M.fixed <- matrix(0, nrow = n, ncol = n)
  ##
  rownames(Pos.fixed) <- colnames(Pos.fixed) <- rownames(pscore.matrix.fixed)
  rownames(M.fixed) <- colnames(M.fixed) <- rownames(pscore.matrix.fixed)
  ##
  Pos.random <- M.random <- matrix(0, nrow = n, ncol = n)
  ##
  rownames(Pos.random) <- colnames(Pos.random) <- rownames(pscore.matrix.random)
  rownames(M.random) <- colnames(M.random) <- rownames(pscore.matrix.random)
  
  
  ## Pos[i, j] counts how many rankings judge treatment i superior to treatment j
  ## (number between 0 and o)
  ##
  for (i in 1:n)
    for (j in 1:n)
      for (k in 1:o)
        if (!is.na(pscore.matrix.fixed[i,k]) &
            !is.na(pscore.matrix.fixed[j,k]) &
            pscore.matrix.fixed[i, k] > pscore.matrix.fixed[j, k])
          Pos.fixed[i, j] <- Pos.fixed[i, j] + 1
  ##
  for (i in 1:n)
    for (j in 1:n)
      for (k in 1:o)
        if (!is.na(pscore.matrix.random[i,k]) &
            !is.na(pscore.matrix.random[j,k]) &
            pscore.matrix.random[i, k] > pscore.matrix.random[j, k])
          Pos.random[i, j] <- Pos.random[i, j] + 1
  
  
  ## Calculate "full" Hasse matrix M
  ##
  M.fixed[Pos.fixed == o] <- 1
  ##
  M.random[Pos.random == o] <- 1
  
  
  ## Matrix with information about partial ordering
  ##
  PO.fixed <- M.fixed[1:2, ]
  PO.fixed[1, ] <- rowSums(M.fixed)
  PO.fixed[2, ] <- colSums(M.fixed)
  rownames(PO.fixed) <- c("inferior", "superior")
  ##
  PO.random <- M.random[1:2, ]
  PO.random[1, ] <- rowSums(M.random)
  PO.random[2, ] <- colSums(M.random)
  rownames(PO.random) <- c("inferior", "superior")
  
  
  ## Skipping each direct path where a path of length 2 is found
  ##
  M0.fixed <- M.fixed - sign(M.fixed %*% M.fixed)
  ##
  M0.random <- M.random - sign(M.random %*% M.random)
  
  
  if (missing(small.values))
    small.values <- NULL
  
  
  res <- list(P.fixed = pscore.matrix.fixed,
              M0.fixed = M0.fixed,
              M.fixed = M.fixed,
              O.fixed = PO.fixed,
              P.random = pscore.matrix.random,
              M0.random = M0.random,
              M.random = M.random,
              O.random = PO.random,
              small.values = small.values,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              call = match.call(),
              version = packageDescription("netmeta")$Version)
  ##
  class(res) <- "netposet"
  
  res
}
