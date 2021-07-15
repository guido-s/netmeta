contribution.matrix <- function(x, method, model, hatmatrix.F1000) {
  if (method == "randomwalk")
    return(contribution.matrix.davies(x, model))
  else if (method == "stream")
    return(contribution.matrix.tpapak(x, model, hatmatrix.F1000))
  else
    return(NULL)
}





contribution.matrix.tpapak <- function(x, model, hatmatrix.F1000) {
  
  meta:::chkclass(x, "netmeta")
  model <- meta:::setchar(model, c("fixed", "random"))
  meta:::chklogical(hatmatrix.F1000)
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
      lapply(1:length(directs),
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
    dg <- igraph::set.edge.attribute(dg, 'label', value = igraph::E(dg))
    ##
    dg
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
  reduceGraph <- function (g, comparison) {
    getshortest <- function (g, compariston) {
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
  
  
  ## rows of comparison matrix
  comparisons <- unlist(lapply(rownames(H), unlist))
  ##
  lapply(comparisons,
         function (comp)
           reduceGraph(initRowGraph(H, comp), comp))
  
  
  colnames(weights) <- directs
  ##
  weights[is.zero(weights)] <- 0
  ##
  attr(weights, "model") <- model
  attr(weights, "hatmatrix.F1000") <- old
  ##
  weights
}





contribution.matrix.davies <- function(x, model) {
  
  meta:::chkclass(x, "netmeta")
  model <- meta:::setchar(model, c("fixed", "random"))
  
  
  ##
  ## Hat matrix
  ##
  ## 'Full' aggregated hat matrix
  ##
  H.full <- hatmatrix.aggr(x, model, type = "full")
  ## Total number of possible pairwise comparisons
  n.comps <- nrow(H.full)
  ## Create prop contribution matrix (square matrix)
  weights <- matrix(0, nrow = n.comps, ncol = n.comps)
  ##
  ## Treatment labels
  ##
  label <- vector("character", x$n) 
  for (t in 1:x$n) {
    foundt <- 0
    ## Look for treatment 't' in vector treat1.pos
    for (i in seq_along(x$treat1)) {
      if (x$treat1.pos[i] == t) {
        label[t] <- x$treat1[i] 
        foundt <- foundt + 1
      }
    }
    ## If treatment t is not in treat1 list look in treat2.pos
    if (foundt == 0){
      for (i in 1:length(x$treat2)){
        if (x$treat2.pos[i] == t){
          label[t] <- x$treat2[i]
          foundt <- foundt + 1
        }
      }
    }
  }
  
  ## Create labels for rows and columns of the contribution matrix
  label.names <- vector("character", n.comps)
  count <- 0
  for (t1 in 1:(x$n - 1)) {
    for (t2 in (t1 + 1):x$n) {
      count <- count + 1
      ## Label row and column by the treatment comparison
      label.names[count] <- paste(label[t1], label[t2], sep = x$sep.trts)
    }
  }
  ##
  rownames(weights) <- label.names
  colnames(weights) <- label.names
  ##
  ## Cycle through comparisons
  r <- 0
  for (t1 in 1:(x$n - 1)) {
    for (t2 in (t1 + 1):x$n) {
      r <- r + 1 
      ## For each row, t1 is the source and t2 is the sink
      ##
      ## Teransition matrix
      ##
      ## Create a transition matrix from t1 (source) to t2 (sink)
      ## using H.full
      ##
      ## Transition matrix (P) is n x n
      ##
      ## Create a dummy n x n matrix Q (this is a step in getting P)
      ##
      ## Assign Q[i,j] equal to the element of H.full in row r and
      ## column representing comparison ij this defines only the upper
      ## half of the matrix Q
      Q <- matrix(0, nrow = x$n, ncol = x$n)
      idx <- 0
      for (i in 1:(x$n - 1)) {
        for (j in (i + 1):x$n) {
          idx <- idx + 1
          Q[i,j] <- H.full[r, idx]
        }
      }
      ## Now use H_ij = -H_ji to define the lower half of Q
      for (i in 1:x$n) {
        for (j in 1:x$n) {
          if (i != j) {
            Q[j, i] <- -Q[i, j]
          }
        }
      }
      ## RW can only move in direction of flow therefore if H_ij < 0
      ## then set Q[i, j] = 0
      for (i in 1:x$n) {
        for (j in 1:x$n) {
          if (i != j) {
            if (Q[i, j] <= 0) {
              Q[i, j] <- 0.0
            }
          }
        }
      }
      ## Create the transition matrix by normalising the values in
      ## each row of Q
      P <- matrix(0, nrow = x$n, ncol = x$n)
      for (i in 1:x$n) {
        sum_row <- 0.0
        for (j in 1:x$n)
          sum_row <- sum_row + Q[i, j]
        ##
        if (sum_row != 0)
          for (j in 1:x$n)
            P[i, j] <- Q[i, j] / sum_row
      }
      ## Edit row corresponding to sink (t2)
      for (j in 1:x$n) {
        if (j == t2) {
          ## Once the walker reaches the sink it remains there forever
          P[t2, j] <- 1.0
        }
        else {
          ## No transitions can occur from the sink to another node
          P[t2, j] <- 0.0
        }
      }
      ##
      ## Find all possible paths
      ##
      ## Define an igraph using P as an adjacency matrix the graph is
      ## directed (and acyclic) and weighted
      Pgraph <-
        igraph::graph_from_adjacency_matrix(P,"directed",
                                            weighted = TRUE, diag = FALSE)
      ## Now find all possible paths from source to sink simple paths
      ## means no vertex is visited more than once this is true for us
      ## as the graph is acyclic
      all.paths <- igraph::all_simple_paths(Pgraph, t1, t2, mode = "out")
      ##
      ## Calculate streams
      ##
      ## Evidence flow through path = probability of a walker taking
      ## that path vector containing the probability of taking each
      ## path
      prob.paths <- vector("numeric", length(all.paths))
      ##
      for (i in 1:length(all.paths)) {
        prob <- 1.0
        for (j in 1:(length(all.paths[[i]]) - 1)) {
          ## Prob of taking a path = product of transition
          ## probabilities in each edge along the path
          prob <- prob * P[all.paths[[i]][j], all.paths[[i]][j + 1]]
        }
        prob.paths[i] <- prob
      }
      ##
      ## Calculate proportion contributions
      ##
      ## Contribution of an edge = sum over flow of evidence in each
      ## path...  containing that edge divided by the length of the
      ## path vector containing prop contirbutions for row r
      prob.edge <- vector("numeric", n.comps)
      ##
      ind <- 1 #index for path number
      for (i in 1:(x$n - 1)) {
        for (j in (i + 1):x$n) {
          for (k in 1:length(all.paths)) {
            ## Does path k contain the edge ij (or ji)?
            edgeij <- 0
            for (l in 1:(length(all.paths[[k]]) - 1)) {
              if (all.paths[[k]][l] == i & all.paths[[k]][l + 1] == j)
                edgeij <- 1
              else if (all.paths[[k]][l] == j & all.paths[[k]][l + 1] == i)
                edgeij <- 1
            }
            ## If yes: p_ab = sum_{paths i containing ab}
            ## phi_i/length(pi_i)
            if (edgeij == 1) {
              path.length <- length(all.paths[[k]]) - 1
              ## Add phi_i / length(pi_i) to p_ab
              prob.edge[ind] <- prob.edge[ind] + (prob.paths[k] / path.length)
            }
          }
          ## Assign this contribution to the contribution matrix
          weights[r, ind] <- prob.edge[ind]
          ##
          ind <- ind + 1
        }
      }
    }
  }
  
  
  weights[is.zero(weights)] <- 0
  ##
  weights <- weights[, apply(weights, 2, sum) > 0, drop = FALSE]
  ##
  attr(weights, "model") <- model
  ##
  weights
}
