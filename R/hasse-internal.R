# Copied from the R package hasseDiagram, version 0.2.0
# Original author: Krzysztof Ciomek <k.ciomek@gmail.com>
# Source: https://github.com/kciomek/hasseDiagram
# License: MIT
# 
# R functions are included in accordance with the terms of the
# MIT license and have been modified to work with the Hasse matrix
# created with netposet().
# 
# Arguments:
# - M          (logical Hasse matrix)
# - parameters (list with named elements)
# 
# Named elements in parameters:
# - newpage   (whether to call \code{grid.newpage()} before drawing)
# - shape     (shape of diagram nodes)
# - col.lines (edge colour)
# - col.nodes (node colour)

hasse_hasseDiagram <- function(M, parameters = list()) {
  
  stopifnot(is.matrix(M))
  stopifnot(nrow(M) > 0)
  stopifnot(nrow(M) == ncol(M))
  #
  stopifnot(is.list(parameters))
  
  # Set defaults
  #
  parameters$newpage <- replaceNULL(parameters$newpage, TRUE)
  #
  parameters$cluster <- TRUE
  parameters$clusterMerge <- FALSE
  parameters$clusterNonAdjacent <- FALSE
  parameters$arrow <- "forward"
  #
  parameters$shape <- replaceNULL(parameters$shape, "roundrect")
  parameters$col.lines <- replaceNULL(parameters$col.lines, "black")
  parameters$col.nodes <- replaceNULL(parameters$col.nodes, "black")
  parameters$lwd <- replaceNULL(parameters$lwd, 1)
  #
  # Node margins as a list with four numerical items:
  # - "tb" and "rl" for top-bottom and right-left margin,
  # - "otb" and "orl" for outer margin when multiple labels are present.
  #
  parameters$margin <- list()
  parameters$margin$rl <- parameters$margin$tb <- 0.125
  parameters$margin$orl <- parameters$margin$otb <- 0.08
  
  # Convert labels to list with named elements
  #
  labels <- as.list(rownames(M))
  names(labels) <- rownames(M)
  
  # Remove self-loops 
  #
  nrNodes <- nrow(M)
  #
  for (i in seq_len(nrNodes)) {
    M[i, i] <- FALSE
  }
  
  # Cluster
  #
  groups <- extractGroups(M, parameters$clusterNonAdjacent)
  toRemove <- c()
  
  for (group in groups) {
    for (i in group) {
      for (j in group) {
        M[i, j] <- FALSE
      }
    }
    
    if (parameters$cluster) {
      first <- group[1]
      rest <- group[-1]
      
      rownames(M)[first] <-
        colnames(M)[first] <-
        names(labels)[first] <- paste(rownames(M)[group], collapse = "")
      
      toRemove <- c(toRemove, rest)
      labels[[first]] <- c(unlist(labels[group]))
    }
  }
  
  if (!is.null(toRemove)) {
    M <- M[-toRemove, -toRemove]
    labels <- labels[-toRemove]
  }
  #
  nrNodes <- nrow(M)
  
  # Detect cycles
  #
  tmpM <- M
  toVisit <- which(sapply(1:nrow(M),
                          function(x) {length(which(tmpM[, x])) == 0}) == TRUE)
  
  while (length(toVisit) > 0) {
    n <- toVisit[1]
    toVisit <- toVisit[-1]
    
    for (m in which(tmpM[n, ] == TRUE)) {
      tmpM[n, m] <- FALSE
      
      if (length(which(tmpM[, m])) == 0) {
        toVisit <- c(toVisit, m)
      }
    }
  }
  
  notRemovedEdges <- which (tmpM == TRUE, arr.ind = TRUE)
  
  if (nrow(notRemovedEdges) > 0) {
    stop(paste("Cycle detected. Check edges: ",
               paste(sapply(seq_len(nrow(notRemovedEdges)),
                            function(x) {
                              paste(rownames(notRemovedEdges)[
                                notRemovedEdges[x, ]], collapse = "-")} ),
                     collapse = ", "),
               sep = ""))
  }
  
  # Perform transitive reduction
  #
  for (source in seq_len(nrNodes)) {
    stack <- which(M[source, ])
    visited <- rep(F, nrNodes)
    visited[stack] <- T
    
    while (length(stack) > 0) {
      element <- stack[1]
      stack <- stack[-1]
      
      children <- which(M[element, ])
      for (child in children) {
        M[source, child] = FALSE
        if (!visited[child]) {
          stack <- c(child, stack)
        }
      }
    }
  }
  
  # Calculate node levels
  #
  ranks <- rep(1, nrNodes)
  queue <- which(sapply(1:nrow(M),
                        function(x) {length(which(M[, x])) == 0}) == TRUE)
  distances <- rep(1, length(queue))
  
  while (length(queue) > 0) {
    element <- queue[1]
    queue <- queue[-1]
    dist <- distances[1]
    distances <- distances[-1]
    children <- which(M[element, ])
    
    for (i in seq_len(length(children))) {
      idx <- which(queue == children[i])
      
      if (length(idx) == 0) {
        ranks[children[i]] <- dist + 1
        queue <- c(queue, children[i])
        distances <- c(distances, dist + 1)
      } else {
        distances[idx] <- max(distances[idx], dist + 1)
        ranks[children[i]] <- max(ranks[children[i]], dist + 1)
      }
    }
  }
  
  # Build graph
  #
  graph <- as(graphAM(M, "directed"), "graphNEL")
  
  nAttrs <- list()
  nAttrs$width <- sapply(labels, function(x) { nWi(x, parameters) })
  nAttrs$height <- sapply(labels, function(x) { nHi(x, parameters) })
  nAttrs$fixedsize <- rep(TRUE, nrNodes)
  nAttrs <- lapply(nAttrs, function(x) { names(x) <- rownames(M); x})
  
  subGList <- list()
  
  for (i in seq_len(max(ranks))) {
    subGList[[length(subGList) + 1]] <- list(graph = subGraph(rownames(M)[which(ranks == i)], graph),
                                             cluster = FALSE)
  }
  
  ragraph <- agopen(graph,
                    name = "graph",
                    subGList = subGList,
                    attrs = list(node = list(shape = "box"),
                                 graph = list(rank = "same", rankdir = "TB")),
                    nodeAttrs = nAttrs)
  # Draw graph
  #
  if (parameters$newpage) {
    grid.newpage()
  }
  hGrob <- hasseGrob(ragraph, labels, parameters)
  grid.draw(hGrob)
  #return (hGrob)
}

isGroup <- function(data, i, j, groupNonAdjacent) {
  if ((data[i, j] == TRUE && data[j, i] == TRUE) || groupNonAdjacent == TRUE) {
    iParents <- data[, i]
    jParents <- data[, j]
    iChildren <- data[i, ]
    jChildren <- data[j, ]
    
    iParents[j] <- FALSE
    jParents[i] <- FALSE
    iChildren[j] <- FALSE
    jChildren[i] <- FALSE
    
    if (all(iParents == jParents) && all(iChildren == jChildren)) {
      return (TRUE)
    }
  }
  
  return (FALSE)
}

extractGroups <- function(data, groupNonAdjacent) {
  result <- list()
  itemGroup <- seq_len(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(nrow(data))) {
      if (isGroup(data, i, j, groupNonAdjacent)) {
        iGroup <- which(itemGroup == itemGroup[i])
        mergable <- TRUE
        
        for (k in iGroup) {
          if (k != i) {
            if (!isGroup(data, j, k, groupNonAdjacent)) {
              mergable <- FALSE
              break
            }
          }
        }
        
        if (mergable) {
          itemGroup[j] <- itemGroup[i]
        }
      }
    }
  }
  
  for (g in unique(itemGroup)) {
    items <- which(itemGroup == g)
    if (length(items) > 1) {
      result[[length(result) + 1]] <- items
    }
  }
  
  return (result)
}

# Node height by labels (in inches)
#
nHi <- function(labels, parameters) {
  result <- unit(1, "lines") + unit(parameters$margin$tb * 2, "inch")
  if (length(labels) > 1 && parameters$clusterMerge == FALSE)
    result <- result + unit(parameters$margin$otb * 2, "inch")
  
  return (convertY(result, "inches", TRUE))
}

# Node width by labels (in inches)
#
nWi <- function(labels, parameters) {
  result <- unit(0, "inch")
  for (label in labels)
    result <-
      result + stringWidth(label) + unit(parameters$margin$rl * 2, "inch")
  if (length(labels) > 1 && parameters$clusterMerge == FALSE)
    result <-
      result + (length(labels) + 1) * unit(parameters$margin$orl, "inch")
  
  return (convertX(result, "inches", TRUE))
}

drawNode <- function(x, y, width, height, labels, parameters, isInner=FALSE) {
  vp <- viewport(x,
                 y,
                 width,
                 height,
                 xscale = c(0, nWi(labels, parameters)),
                 yscale = c(0, nHi(labels, parameters)))
  pushViewport(vp)
  
  if (parameters$shape != "none" &&
      (isInner == FALSE || parameters$clusterMerge == FALSE)) {
    gp <- gpar(col = parameters$col.nodes, lwd = parameters$lwd)
    
    if (parameters$shape == "rect")
      grid.rect(gp = gp)
    else if (parameters$shape == "roundrect")
      grid.roundrect(gp = gp)
    else
      stop(paste("Unsupported node shape '", parameters$shape, "'.", sep = ""))
  }
  
  grid.clip()
  
  if (length(labels) == 1) {
    cex <- min(1.0 / (convertWidth(stringWidth(labels) +
                                     unit(parameters$margin$rl * 2, "inch"), 
                                   "npc", TRUE)),
               1.0 / (convertHeight(unit(1, "lines") +
                                      unit(parameters$margin$tb * 2, "inch"),
                                    "npc", TRUE)))
    
    grid.text(labels[1], gp = gpar(cex = cex))
  }
  else {
    xCenter <- unit(ifelse(parameters$clusterMerge == FALSE,
                           parameters$margin$orl, 0.0), "native")
    yCenter <- unit(0.5, "npc")
    
    for (label in labels) {
      drawNode(xCenter + unit(nWi(label, parameters), "native") * 0.5,
               yCenter,
               unit(nWi(label, parameters), "native"),
               unit(nHi(label, parameters), "native"),
               label,
               parameters,
               TRUE)
      xCenter <- xCenter + unit(nWi(label, parameters), "native") 
      if (parameters$clusterMerge == FALSE)
        xCenter <- xCenter + unit(parameters$margin$orl, "native")
    }
  }
  
  popViewport()
}


hasseGrob <- function(graph, labels, parameters) {
  grob(graph = graph, labels = labels,
       parameters = parameters, cl = "hasseGrob")
}


#' @method drawDetails hasseGrob
#' @export

drawDetails.hasseGrob <- function(x, recording) {
  g <- x$graph
  ur <- upRight(boundBox(g))
  bl <- botLeft(boundBox(g))
  
  vp <- viewport(width = unit(0.96, "npc"),
                 height = unit(0.96, "npc"),
                 xscale = c(getX(bl), getX(ur)),
                 yscale = c(getY(bl), getY(ur)))
  
  pushViewport(vp)
    
  # Draw edges before nodes
  #
  dir <- x$parameters$arrow
  gp <- gpar(col = x$parameters$col.lines, lwd = x$parameters$lwd)
  
  for (edge in AgEdge(g)) {
    nrLines <- length(edge@splines)
    
    for (i in seq_len(nrLines)) {
      arrow <- NULL
      arrowEnds <- NULL
      
      if (dir == "forward" && i == nrLines) {
        arrowEnds = "last"
      }
      else if (dir == "backward" && i == 1) {
        arrowEnds = "first"
      }
      else if (dir == "both") {
        if (nrLines == 1)
          arrowEnds = "both"
        else if (i == 1)
          arrowEnds = "first"
        else if (i == nrLines)
          arrowEnds = "last"
      }
      
      if (!is.null(arrowEnds)) {
        arrow <- arrow(angle = 30,
                       length = min(unit(4, "mm"), unit(0.02, "npc")),
                       ends = arrowEnds,
                       type = "open")
      }
      
      bp <- bezierPoints(edge@splines[[i]])
      grid.lines(bp[, 1], bp[, 2], default.units = "native",
                 arrow = arrow, gp = gp)
    }
  }
  
  # Draw nodes
  #
  for (agNode in AgNode(g)) {
    center <- getNodeCenter(agNode)
    centerX <- unit(getX(center), "native")
    centerY <- unit(getY(center), "native")
    width <- unit(getNodeRW(agNode) + getNodeLW(agNode), "native")
    height <- unit(getNodeHeight(agNode), "native")
    
    drawNode(centerX, centerY, width, height,
             unlist(x$labels[agNode@name]), x$parameters)
  }
    
  popViewport()
}
