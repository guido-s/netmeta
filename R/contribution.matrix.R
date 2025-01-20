contribution.matrix <- function(x, method, model, hatmatrix.F1000,
                                verbose = FALSE) {
  if (verbose)
    cat("Calculate network contributions (", model, " effects model):\n",
        sep = "")
  ##
  if (method == "randomwalk")
    return(contribution.matrix.davies(x, model, verbose = verbose))
  else if (method == "shortestpath")
    return(contribution.matrix.tpapak(x, model, hatmatrix.F1000,
                                      verbose = verbose))
  else if (method == "cccp")
    return(contribution.matrix.ruecker.cccp(x, model, verbose = verbose))
  else if (method == "pseudoinverse")
    return(contribution.matrix.ruecker.pseudoinv(x, model, verbose = verbose))
  else
    return(NULL)
}





contribution.matrix.tpapak <- function(x, model, hatmatrix.F1000,
                                       verbose = FALSE) {
  
  chkclass(x, "netmeta")
  model <- setchar(model, c("common", "random"))
  chklogical(hatmatrix.F1000)
  chklogical(verbose)
  ##
  old <- hatmatrix.F1000
  
  
  ##
  ## Auxiliary R functions
  ##
  split <- function(dir)
    strsplit(dir, x$sep.trts)
  ##
  ## comparisonToEdge <- function (comp) unlist (split(comp))
  ##
  setWeights <- function(g, comparison, conMat)
    igraph::set.edge.attribute(g, "weight",
                               value = rep(1, dims[2]))
  ##
  getFlow <- function(g, edge)
    igraph::E(g)[edge]$flow
  ##
  sv <- function (comparison)
    split(comparison)[[1]][1][1]
  ##
  tv <- function (comparison)
    split(comparison)[[1]][2][1]
  ##
  initRowGraph <- function(x, comparison) {
    dedgeList <-
      lapply(seq_len(length(directs)),
             function(comp) {
               if (x[comparison, comp] > 0)
                 return(c(sv(directs[comp]), tv(directs[comp])))
               else
                 return(c(tv(directs[comp]), sv(directs[comp])))
             })
    ##
    dedgeList <- matrix(unlist(dedgeList), ncol = 2, byrow = TRUE)
    ##
    flows <- abs(x[comparison, ])
    dg <- igraph::graph_from_edgelist(dedgeList, directed = TRUE)
    igraph::E(dg)[]$weight <- rep(0, dims[2])
    igraph::E(dg)[]$flow <- abs(x[comparison, ])
    igraph::V(dg)[]$label <- igraph::V(dg)[]$name
    resg <-
      igraph::set.edge.attribute(dg, 'label', value = 1:igraph::gsize(dg))
    ##
    resg
  }
  ##
  reducePath <- function(g, comparison, spl) {
    pl <- length(spl[[1]])
    splE <- lapply(spl[[1]], function(e) {
      return (igraph::E(g)[e[]])
    })
    flow <- min(unlist(lapply(splE,
                              function(e) {
                                return(e$flow[])
                              })))
    gg <- Reduce(function(g, e) {
      elabel <- e$label
      pfl <- e$flow[]
      g <- igraph::set.edge.attribute(g, "flow", e, pfl-flow)
      cw <- e$weight[] + (flow[1] / pl) 
      weights[comparison, elabel] <<- cw
      return(igraph::set.edge.attribute(g, "weight", e, cw))},
      splE, g)
    emptyEdges <- Reduce(function(removedEdges, e) {
      e <- igraph::E(gg)[e[]]
      if(e$flow[[1]][[1]] == 0)
        removedEdges <- c(removedEdges, e)
      return(removedEdges)}, splE, c())
    ##
    igraph::delete_edges(gg, emptyEdges)
  }
  ##
  reduceGraph <- function(g, comparison, verbose, is.tictoc) {
    if (verbose)
      cat("- ", comparison, "\n", sep = "")
    ##
    if (verbose & is.tictoc)
      tictoc::tic()
    ##
    getshortest <- function (g, comparison) {
      getShortest <- function() {
        return(igraph::get.shortest.paths(g,
                                          sv(comparison),
                                          tv(comparison),
                                          mode = "out",
                                          output = "epath",
                                          weights = NA)$epath)
      }
      ##
      res <- suppressWarnings(getShortest())
      return(res)
    }
    spath <- getshortest(g, comparison)
    while (length(unlist(spath)) > 0) {
      g <- reducePath(g, comparison, spath)
      spath <- getshortest(g, comparison)
    }
    ##
    if (verbose & is.tictoc)
      tictoc::toc()
    ##
    g
  }
  
  
  if (old)
    H <- hatmatrix.F1000(x, model)
  else
    H <- hatmatrix.aggr(x, model, "long")
  ##
  directs <- colnames(H)
  
  
  directlist <- unlist(lapply(lapply(directs, split), unlist))
  edgeList <- matrix(directlist, ncol = 2, byrow = TRUE)
  ##
  g <- igraph::graph_from_edgelist(edgeList , directed = FALSE)
  g <- igraph::set.vertex.attribute(g, 'label', value = igraph::V(g))
  g <- igraph::set.edge.attribute(g, 'label', value = igraph::E(g))
  
  
  dims <- dim(H)
  ##
  contribution <- rep(0, dims[2])
  names(contribution) <- seq_len(dims[2])
  ##
  weights <- matrix(rep(0, dims[2] * dims[1]),
                    nrow = dims[1], ncol = dims[2], byrow = TRUE)
  rownames(weights) <- rownames(H)
  colnames(weights) <- seq_len(dims[2])
  
  
  is.tictoc <- is_installed_package("tictoc", stop = FALSE)
  ## rows of comparison matrix
  comparisons <- unlist(lapply(rownames(H), unlist))
  ##
  lapply(comparisons,
         function (comp)
           reduceGraph(initRowGraph(H, comp), comp, verbose, is.tictoc))
  
  
  colnames(weights) <- directs
  ##
  weights[is_zero(weights)] <- 0
  ##
  attr(weights, "model") <- model
  attr(weights, "hatmatrix.F1000") <- old
  ##
  res <- list(weights = weights)
  res
}





contribution.matrix.davies <- function(x, model, verbose = FALSE) {
  
  chkclass(x, "netmeta")
  model <- setchar(model, c("common", "random"))
  chklogical(verbose)
  ##
  is.tictoc <- is_installed_package("tictoc", stop = FALSE)
  
  
  ##
  ## 'Full' aggregated hat matrix
  ##
  H.full <- hatmatrix.aggr(x, model, type = "full")
  ## Total number of possible pairwise comparisons
  comps <- rownames(H.full)
  n.comps <- length(comps)
  
  
  ##
  ## Create prop contribution matrix (square matrix)
  ##
  weights <- matrix(0, nrow = n.comps, ncol = n.comps)
  rownames(weights) <- colnames(weights) <- comps
  
  
  ##
  ## Measure times
  ##
  if (is.tictoc) {
    tictoc <- rep(NA, n.comps)
    names(tictoc) <- comps
  }
  
  
  ##
  ## Cycle through comparisons
  ##
  n <- x$n
  r <- 0
  ##
  for (t1 in seq_len(n - 1)) {
    for (t2 in (t1 + 1):n) {
      r <- r + 1
      ##
      if (is.tictoc)
        tictoc::tic()
      ##
      if (verbose)
        cat("- ",
            paste(x$trts[t1], x$trts[t2], sep = x$sep.trts),
            " (", r, "/", n.comps, ")\n",
            sep = "")
      ##
      ## For each row, t1 is the source and t2 is the sink
      ##
      ## Create a transition matrix (P, n x n) from t1 (source) to t2
      ## (sink) using H.full
      ##
      ## Create a dummy n x n matrix Q (this is a step in getting P)
      ##
      ## Assign Q[i, j] equal to the element of H.full in row r and
      ## column representing comparison ij this defines only the upper
      ## half of the matrix Q
      ##
      Q <- matrix(0, nrow = n, ncol = n)
      idx <- 0
      for (i in seq_len(n - 1))
        for (j in (i + 1):n) {
          idx <- idx + 1
          Q[i, j] <- H.full[r, idx]
        }
      ##
      ## Now use H_ij = -H_ji to define the lower half of Q
      ##
      for (i in seq_len(n))
        for (j in seq_len(n))
          if (i != j)
            Q[j, i] <- -Q[i, j]
      ##
      ## RW can only move in direction of flow therefore if H_ij < 0
      ## then set Q[i, j] = 0
      ##
      for (i in seq_len(n))
        for (j in seq_len(n))
          if (i != j & Q[i, j] <= 0)
            Q[i, j] <- 0.0
      ##
      ## Create the transition matrix by normalising the values in
      ## each row of Q
      ##
      P <- ginv(diag(rowSums(Q))) %*% Q
      P[t2, ] <- rep(0, n)
      P[t2, t2] <- 1
      ##
      P[is_zero(P, n = 1000)] <- 0
      ##
      ## Find all possible paths
      ##
      ## Define an igraph using P as adjacency matrix. The graph is
      ## directed (and acyclic) and weighted
      ##
      Pgraph <-
        igraph::graph_from_adjacency_matrix(P, "directed",
                                            weighted = TRUE, diag = FALSE)
      ##
      ## Now find all possible paths from source to sink
      ##
      ## Simple paths means no vertex is visited more than once this
      ## is true for us as the graph is acyclic
      ##
      all.paths <- igraph::all_simple_paths(Pgraph, t1, t2, mode = "out")
      ##
      ## Calculate streams
      ##
      R <- matrix(0, nrow = n, ncol = n)
      colnames(R) <- rownames(R) <- x$seq
      ##
      for (p in seq_len(length(all.paths))) {
        path.p <- all.paths[[p]]
        L.p <- length(path.p) - 1
        prob.p <- 1 / L.p
        ##
        ## Prob of taking a path.p = product of transition
        ## probabilities in each edge along the path
        ##
        for (i in seq_len(L.p))
          prob.p <- prob.p * P[path.p[i], path.p[i + 1]]
        ##
        for (i in seq_len(L.p))
          R[path.p[i], path.p[i + 1]] <-
            R[path.p[i], path.p[i + 1]] + prob.p
        ##
        for (i in seq_len(n - 1))
          for (j in (i + 1):n)
            R[i, j] <- R[j, i] <- max(R[i, j], R[j, i])
      }
      ##
      weights[r, ] <- uppertri(R)
      ##
      if (is.tictoc) {
        tictoc.r <- tictoc::toc(func.toc = NULL)
        tictoc[r] <- as.numeric(tictoc.r$toc) - as.numeric(tictoc.r$tic)
        ##
        if (verbose)
          cat(round(tictoc[r], 3), "sec elapsed\n")
      }
    }
  }
  ##
  weights[is_zero(weights)] <- 0
  ##
  weights <- weights[, apply(weights, 2, sum) > 0, drop = FALSE]
  ##
  attr(weights, "model") <- model
  
  
  res <- list(weights = weights)
  ##
  if (is.tictoc)
    res$tictoc <- tictoc
  ##
  res
}





contribution.matrix.ruecker.cccp <- function (x, model, verbose = FALSE) {
  
  chkclass(x, "netmeta")
  ##
  model <- setchar(model, c("common", "random"))
  chklogical(verbose)
  ##
  is.tictoc <- is_installed_package("tictoc", stop = FALSE)
  
  
  H.full <- hatmatrix.aggr(x, model, type = "full")
  ## Total number of possible pairwise comparisons
  comps <- rownames(H.full)
  n.comps <- nrow(H.full)
  ##
  n <- x$n
  
  
  ##
  ## Create prop contribution matrix (square matrix)
  ##
  weights <- matrix(0, nrow = n.comps, ncol = n.comps)
  rownames(weights) <- colnames(weights) <- rownames(H.full)
  
  
  ##
  ## Create phi vector
  ##
  phi <- pl <- vector("list", n.comps)
  names(phi) <- names(pl) <- comps
  
  
  ##
  ## List with Z matrices
  ##
  zlist <- vector("list", choose(n,2))
  
  
  ##
  ## Measure times
  ##
  if (is.tictoc) {
    tictoc <- rep(NA, n.comps)
    names(tictoc) <- comps
  }
  
  
  ##
  ## Cycle through comparisons
  ##
  r <- 0
  ##
  for (t1 in seq_len(n - 1)) {
    for (t2 in (t1 + 1):n) {
      r <- r + 1
      ##
      if (is.tictoc)
        tictoc::tic()
      ##
      if (verbose) {
        cat("*** ",
            paste(x$trts[t1], x$trts[t2], sep = x$sep.trts),
            " (", r, " / ", n.comps, ") ***\n",
            sep = "")
      }
      ##
      Q <- matrix(0, nrow = n, ncol = n)
      idx <- 0
      for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
          idx <- idx + 1
          Q[i, j] <- H.full[r, idx]
        }
      }
      ##
      for (i in seq_len(n))
        for (j in seq_len(n))
          if (i != j)
            Q[j, i] <- -Q[i, j]
      ##
      for (i in seq_len(n))
        for (j in seq_len(n))
          if (i != j & Q[i, j] <= 0)
            Q[i, j] <- 0
      ##
      P <- ginv(diag(rowSums(Q))) %*% Q # Transition matrix
      P[t2, ] <- rep(0, n)
      P[t2, t2] <- 1
      P[is_zero(P, n = 1000)] <- 0
      Pgraph <-
        igraph::graph_from_adjacency_matrix(
                  P, "directed", weighted = TRUE, diag = FALSE)
      all.paths <- igraph::all_simple_paths(Pgraph, t1, t2, mode = "out")
      ##
      Z <- matrix(0, nrow = length(all.paths), ncol = dim(H.full)[2])
      rownames(Z) <- all.paths
      colnames(Z) <- colnames(H.full)
      ##
      for (p in 1:length(all.paths)) {
        l <- length(all.paths[[p]])
        for (i in 1:(l-1)) {
          s <- 1
          v1 <- all.paths[[p]][i]
          v2 <- all.paths[[p]][i+1]
          if (v2 < v1) {
            s <- -1
            v1 <- all.paths[[p]][i+1]
            v2 <- all.paths[[p]][i]
          }
          Z[p, n * (v1 - 1) - choose(v1 + 1, 2) + v2] <- s
        }
      }
      ##
      zlist[[r]] <- Z
      ##
      ## There is an exact least squares (L2) solution for 
      ## phi %*% Z = H.full[r, ]
      ## given by
      ## phi0 = H.full[r, ] %*% ginv(Z)
      ##
      ## The L1 solution is given by this idea (h = H.full[r, ]):
      ## If phi %*% Z = h is solvable, the full set of solutions is given by
      ##    phi = h %*% ginv(Z) + w %*% (I - Z %*% ginv(Z))
      ##        = phi0 + w %*% (I - Z %*% ginv(Z))
      ## Minimising |phi| means finding w such that 
      ##       |phi0 + w %*% (I - Z %*% ginv(Z))| is minimised
      ## Thus the problem is to minimise |phi0 + w %*% A| where
      ##       A = (I - Z %*% ginv(Z)) and phi0 = h %*% ginv(Z)
      ## If a solution w is found, the solution for phi is given by 
      ##       phi = phi0 + w %*% A
      ##
      ## Find LS solution phi0.r
      ##
      path.length <- rowSums(abs(Z))
      phi0.r <- H.full[r, ] %*% ginv(Z)
      ##
      ## Find L1 solution
      ##
      A <- diag(p) - Z %*% ginv(Z)
      ##
      ## Minimise |wA + phi0.r|
      ##
      invisible(capture.output(mn <-
                                 cccp::l1(A, -phi0.r,
                                          optctrl = cccp::ctrl(trace = FALSE)),
                               type = "message"))
      w <- cccp::getx(mn)[1:p]
      phi.r <- phi0.r + w %*% A
      ##
      weights[r, ] <- (phi.r / path.length) %*% abs(Z)
      phi[[r]] <- as.vector(phi.r)
      pl[[r]] <- path.length
      ##
      if (is.tictoc) {
        tictoc.r <- tictoc::toc(func.toc = NULL)
        tictoc[r] <- as.numeric(tictoc.r$toc) - as.numeric(tictoc.r$tic)
        ##
        if (verbose)
          cat(round(tictoc[r], 3), "sec elapsed\n")
      }
    }
  }
  ##
  w <- weights
  weights[is_zero(weights)] <- 0
  weights <- weights[, apply(weights, 2, sum) > 0, drop = FALSE]
  attr(weights, "model") <- model
  
  
  res <- list(weights = weights, phi = phi, w = w, pl = pl, zlist = zlist)
  ##
  if (is.tictoc) 
    res$tictoc <- tictoc
  ##
  res
}


contribution.matrix.ruecker.pseudoinv <- function (x, model, verbose = FALSE) {
  
  chkclass(x, "netmeta")
  model <- setchar(model, c("common", "random"))
  chklogical(verbose)
  ##
  is.tictoc <- is_installed_package("tictoc", stop = FALSE)
  
  
  ##
  ## 'Full' aggregated hat matrix
  ##
  H.full <- hatmatrix.aggr(x, model, type = "full")
  ## Total number of possible pairwise comparisons
  comps <- rownames(H.full)
  n.comps <- length(comps)
  n <- x$n
  
  
  ##
  ## Create prop contribution matrix (square matrix)
  ##
  weights <- matrix(0, nrow = n.comps, ncol = n.comps)
  rownames(weights) <- colnames(weights) <- comps


  ##
  ## List with phis (path weights) and pl (path lengths)
  ##
  phi <- pl <- vector("list", n.comps)
  names(phi) <- names(pl) <- comps
  
  ##                                       ####
  ## List with Z matrices                  ####
  ##                                       ####
  zlist <- vector("list", choose(n,2))     ####
  
  
  ##
  ## Measure times
  ##
  if (is.tictoc) {
    tictoc <- rep(NA, n.comps)
    names(tictoc) <- comps
  }
  
  
  ##
  ## Cycle through comparisons
  ##
  r <- 0
  ##
  for (t1 in seq_len(n - 1)) {
    for (t2 in (t1 + 1):n) {
      r <- r + 1
      ##
      if (is.tictoc)
        tictoc::tic()
      ##
      if (verbose)
        cat("- ",
            paste(x$trts[t1], x$trts[t2], sep = x$sep.trts),
            " (", r, "/", n.comps, ")\n",
            sep = "")
      ##
      ## For each row, t1 is the source and t2 is the sink
      ##
      ## Create a transition matrix (P, n x n) from t1 (source) to t2
      ## (sink) using H.full
      ##
      ## Create a dummy n x n matrix Q (this is a step in getting P)
      ##
      ## Assign Q[i, j] equal to the element of H.full in row r and
      ## column representing comparison ij this defines only the upper
      ## half of the matrix Q
      ##
      Q <- matrix(0, nrow = n, ncol = n)
      idx <- 0
      for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
          idx <- idx + 1
          Q[i, j] <- H.full[r, idx]
        }
      }
      ##
      ## Now use H_ij = -H_ji to define the lower half of Q
      ##
      for (i in seq_len(n))
        for (j in seq_len(n))
          if (i != j)
            Q[j, i] <- -Q[i, j]
      ##
      ## RW can only move in direction of flow therefore if H_ij < 0
      ## then set Q[i, j] = 0
      ##
      for (i in seq_len(n))
        for (j in seq_len(n))
          if (i != j & Q[i, j] <= 0)
            Q[i, j] <- 0.0
      ##
      ## Create the transition matrix by normalising the values in
      ## each row of Q
      ##
      P <- ginv(diag(rowSums(Q))) %*% Q
      P[t2, ] <- rep(0, n)
      P[t2, t2] <- 1
      ##
      P[is_zero(P, n = 1000)] <- 0
      ##
      ## Find all possible paths
      ##
      ## Define an igraph using P as adjacency matrix. The graph is
      ## directed (and acyclic) and weighted
      ##
      Pgraph <-
        igraph::graph_from_adjacency_matrix(P, "directed", weighted = TRUE,
                                            diag = FALSE)
      ##
      ## Now find all possible paths from source to sink
      ##
      ## Simple paths means no vertex is visited more than once 
      ## This is true for us as the graph is acyclic
      ##
      all.paths <- igraph::all_simple_paths(Pgraph, t1, t2, mode = "out")
      ##
      ## Calculate streams
      ##
      Z <- matrix(0, nrow = length(all.paths), ncol = dim(H.full)[2])
      rownames(Z) <- all.paths
      colnames(Z) <- colnames(H.full)
      ##
      ## The solution is unique <=> Z %*% ginv(Z) is an identity
      ## matrix, see round(Z %*% ginv(Z), 8).
      ##
      for (p in 1:length(all.paths)) {
        l <- length(all.paths[[p]])
        for (i in 1:(l-1)) {
          s <- 1
          v1 <- all.paths[[p]][i]
          v2 <- all.paths[[p]][i + 1]
          if (v2 < v1) {
            s <- -1
            v1 <- all.paths[[p]][i + 1]
            v2 <- all.paths[[p]][i]
          }
          Z[p, n * (v1 - 1) - choose(v1 + 1, 2) + v2] <- s
        }
      }
      zlist[[r]] <- Z               ####
      ##
      ## Calculate the least squares solution of
      ## phi %*% Z = H.full[r, ]
      ##
      path.length <- rowSums(abs(Z))
      phi.r <- H.full[r, ] %*% ginv(Z)
      weights[r, ] <- (phi.r / path.length) %*% abs(Z)
      phi[[r]] <- as.vector(phi.r)
      pl[[r]] <- path.length
      ##
      if (is.tictoc) {
        tictoc.r <- tictoc::toc(func.toc = NULL)
        tictoc[r] <- as.numeric(tictoc.r$toc) - as.numeric(tictoc.r$tic)
        ##
        if (verbose)
          cat(round(tictoc[r], 3), "sec elapsed\n")
      }
    }
  }
  ##  
  w <- weights
  weights[is_zero(weights)] <- 0
  weights <- weights[, apply(weights, 2, sum) > 0, drop = FALSE]
  attr(weights, "model") <- model
  
  
  res <- list(weights = weights, phi = phi, w = w, pl = pl, zlist = zlist)
  ##
  if (is.tictoc)
    res$tictoc <- tictoc
  ##
  res
}
