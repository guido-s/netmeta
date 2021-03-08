#' Contribution matrix in network meta-analysis
#' 
#' @description
#' This function generates the contribution of direct comparisons to
#' every network treatment comparison as a different row
#' 
#' @aliases netcontrib print.netcontrib
#' 
#' @param x An object of class \code{netmeta}.
#' 
#' @details
#' In network meta-analysis, it is important to assess the influence of the limitations or 
#' other characteristics of individual studies on the estimates obtained from the network. 
#' The contribution matrix, shows how much each direct treatment effect 
#' contributes to each treatment effect estimate from network meta-analysis, is crucial in this context. 
#' We use ideas from graph theory to derive the proportion that is contributed by each direct treatment effect.
#' We start with the ‘projection’ matrix in a two-step network meta-analysis model, 
#' called the H matrix, which is analogous to the hat matrix in a linear regression model.
#' We develop a method to translate H entries to proportion contributions based on 
#' the observation that the rows of H can be interpreted as flow networks,
#' where a stream is defined as the composition of a path and its associated flow.
#' We present an algorithm that identifies the flow of evidence in each path and decomposes it into direct comparisons.
#' 
#' @return
#'  An object of class \code{netcontrib} with corresponding \code{print}
#' function. The object is a list containing the following components:
#' \item{fixed}{Numeric matrix of percentage contributions of direct comparisons
#' for each network comparison for the fixed effects model.}
#' \item{random}{Numeric matrix of percentage contributions of direct comparisons
#' for each network comparison for the random effects model.}
#' \item{x}{the \code{netmeta} object defined above.}
#' 
#' with the contribution matrices for fixed and random NMA.
#' Each matrix has the percentage contributions of each direct comparisons as columns 
#' for each network comparison, direct or indirect as rows.
#' 
#' @author Thodoris Papakonstantinou \email{dev@tpapak.com}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @references
#' Papakonstantinou, T., Nikolakopoulou, A., Rücker, G., Chaimani, A.,
#' Schwarzer, G., Egger, M., Salanti, G. (2018)
#' Estimating the contribution of studies in network meta-analysis: paths, flows and streams
#' \emph{F1000Research},
#' 
#' @keywords contribution
#' 
#' @examples
#' \dontrun{
#' # Use the Woods dataset
#' #
#' data("Woods2010")
#' p1 <- pairwise(treatment, event = r, n = N,
#'                studlab = author, data = Woods2010, sm = "OR")
#' 
#' net1 <- netmeta(p1)
#' cm = netcontrib(net1)
#' cm$random
#' 
#' }
#' 
#' @rdname netcontrib
#' @export netcontrib


netcontrib = function(x){
  if (!meta:::is.installed.package("igraph", stop = FALSE)) {
    warning(paste("Package 'igraph' missing.",
                  "\n  ",
                  "Please use the following R command for installation:",
                  "\n  ",
                  "install.packages(\"igraph\")",
                  sep = ""),
            call. = FALSE)
    return(invisible(NULL))
  }
  res <- list( fixed = contribution.matrix(x,"fixed")
             , random = contribution.matrix(x,"random")
             , x = x)
  class(res) <- "netcontrib"
  ##
  res
}

#' @rdname netcontrib
#' @param x An ojbect of class \code{netcontrib}
#' @param comb.fixed A logical indicating whether a league table
#'   should be printed for the fixed effects (common effects) network
#'   meta-analysis.
#' @param comb.random A logical indicating whether a league table
#'   should be printed for the random effects network meta-analysis.
#' @param digits number of rounding digits
#' @param \dots Additional arguments (ignored at the moment).
#' @method print netcontrib
#' @export
#' @export print.netcontrib


print.netcontrib <- function(x,
                            comb.fixed = x$x$comb.fixed,
                            comb.random = x$x$comb.random,
                            digits = 2,
                            ...) {
  
  meta:::chkclass(x, "netcontrib")
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chknumeric(digits, length = 1)

   matitle(x$x)
  if(comb.fixed){
    cat("Contribution matrix fixed model")
    cat("\n")
    prmatrix(round(x$fixed, digits))
    cat("\n")
  }
  if(comb.random){
    cat("Contribution matrix random model")
    cat("\n")
    prmatrix(round(x$random, digits))
  }
  invisible(NULL)
}


#' Internal function for \code{netcontrib}
#' @param x An object of class \code{netmeta}.
#' @param model "fixed" or "random" depending on the model of NMA
#' @return the contribution matrix
contribution.matrix <- function(x, model){
  contribution.hatmatrix <- function(metaNetw, model){
  
    #H matrix
    if(model=="fixed"){
      krahn = nma.krahn(metaNetw,tau.preset = 0)
    }
    
    if(model=="random"){
      krahn=nma.krahn(metaNetw,tau.preset = metaNetw$tau)
    }
    
    X.full=krahn$X.full
    direct=krahn$direct
    X=krahn$X.full[direct$comparison,,drop=FALSE]
    Vd=diag(direct$seTE^2,nrow=length(direct$seTE),ncol=length(direct$seTE))
    H <- X.full %*% solve(t(X) %*% solve(Vd) %*% X) %*% t(X) %*% solve(Vd)
    colnames(H)<-rownames(X)
    
    return(list( colNames=colnames(H)
               , rowNames=rownames(H)
               , H=H
               )
    )
  }
  
  hatmatrix = contribution.hatmatrix(x, model)

  directs <- hatmatrix$colNames

  hatMatrix <- hatmatrix$H
  
  rownames(hatMatrix) <- hatmatrix$rowNames

  split <- function (dir) {strsplit(dir,":")}

  dims <- dim(hatMatrix)

#rows of comparison matrix 
  comparisons <- unlist(lapply(rownames(hatMatrix),unlist))

  comparisonToEdge <- function (comp) unlist (split(comp))

  directlist <- unlist(lapply(lapply(directs,split),unlist))

  edgeList <- matrix( directlist, ncol = 2, byrow = TRUE)

  g <- igraph::graph_from_edgelist(edgeList , directed=FALSE)
  g <- igraph::set.vertex.attribute(g,'label',value = igraph::V(g))
  g <- igraph::set.edge.attribute(g,'label',value = igraph::E(g))

  setWeights <- function (g,comparison,conMat) {
    igraph::set.edge.attribute(g,"weight",value=rep(1,dims[2]))
  }


  getFlow <- function(g,edge) {return(igraph::E(g)[edge]$flow)}

  sv <- function (comparison) {split(comparison)[[1]][1][1]}

  tv <- function (comparison) {split(comparison)[[1]][2][1]}

  initRowGraph <- function(comparison) {
    dedgeList <- lapply(1:length(directs),function(comp) {
       if(hatMatrix[comparison,comp]>0){
         return (c(sv(directs[comp]),tv(directs[comp])))
       }else{
         return (c(tv(directs[comp]),sv(directs[comp])))
       }
    })
    dedgeList <- matrix( unlist(dedgeList), ncol = 2, byrow = TRUE)
    flows<-abs(hatMatrix[comparison,])
    dg <- igraph::graph_from_edgelist(dedgeList , directed = TRUE)
    igraph::E(dg)[]$weight <- rep(0,dims[2])
    igraph::E(dg)[]$flow <- abs(hatMatrix[comparison,])
    igraph::V(dg)[]$label <- igraph::V(dg)[]$name
    dg <- igraph::set.edge.attribute(dg,'label',value = igraph::E(dg))
    return(dg)
  }

  contribution = rep(0,dims[2])
  names(contribution) <- c(1:dims[2])

  weights <- matrix (
                     rep(0,dims[2]*dims[1]),
                     nrow=dims[1],
                     ncol=dims[2],
                     byrow=TRUE
                 )
  rownames(weights) <- rownames(hatMatrix)
  colnames(weights) <- c(1:dims[2])

  reducePath <- function (g,comparison,spl) {
    pl <- length(spl[[1]])
    splE <- lapply(spl[[1]], function(e){
       return (igraph::E(g)[e[]])
    })
    flow <- min(unlist(lapply(splE, function(e){
      return(e$flow[])
    })))
    gg <- Reduce(function(g, e){
      elabel <- e$label
      pfl <- e$flow[]
      g <- igraph::set.edge.attribute(g,"flow",e, pfl-flow)
      cw <-  e$weight[] + (flow[1]/pl) 
      weights[comparison,elabel] <<- cw
      return(igraph::set.edge.attribute(g,"weight",e, cw))
    },splE, g)
    emptyEdges <- Reduce(function(removedEdges, e){
      e <- igraph::E(gg)[e[]]
      if(e$flow[[1]][[1]]==0){
        removedEdges <- c(removedEdges, e)
      }
      return(removedEdges)
    },splE, c())
   return(igraph::delete_edges(gg, emptyEdges))
  }
  
 reduceGraph <- function (g,comparison) {
    getshortest <- function (g,compariston) {
      getShortest <- function() {
        return(igraph::get.shortest.paths(g
                                          ,sv(comparison)
                                          ,tv(comparison)
                                          ,mode="out"
                                          ,output="epath"
                                          ,weights=NA)$epath)
      }
      res <- suppressWarnings(getShortest())
      return(res)
      }
    spath <- getshortest(g,comparison)
    while(length(unlist(spath))>0){
      g <- reducePath(g,comparison,spath)
      spath <- getshortest(g,comparison)
    }
    return(g)
  }

  lapply (comparisons, function (comp) {
    reduceGraph (initRowGraph(comp), comp)
  })

  colnames(weights) <- directs
  weights <- 100 * weights
  totalSums <- colSums(weights)
  totalTotal <- sum(totalSums)
  totalWeights <- unlist(lapply(totalSums,function(comp){
                           100 * comp/ totalTotal
  }))
  return(weights)
}
