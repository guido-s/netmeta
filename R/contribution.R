#' Contribution matrix in network meta-analysis
#' 
#' @description
#' This function generates the contribution of direct comparisons to
#' every network treatment comparison as a different row
#' 
#' @param x An object of class \code{netmeta}.
#' @param model "fixed" or "random" nma.
#' 
#' @details
#' 
#' @return
#' A matrix with the percentage contributions of each direct comparisons as columns 
#' for each network comparison, direct or indirect as rows
#' 
#' @author Thodoris Papakonstantinou \email{hi@tpapak.com}
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
#' # Use Woods dataset
#' #
#' data("Woods2010")
#' p1 <- pairwise(treatment, event = r, n = N,
#'                studlab = author, data = Woods2010, sm = "OR")
#' 
#' net1 <- netmeta(p1)
#' cm = contribution.matrix(net1, model="fixed")
#' 
#' }
#' 
#' @export contribution.matrix

contribution.matrix <- function(x, model="random"){
  
  contribution.hatmatrix <- function(metaNetw, model){
  
    #H matrix
    if(model=="fixed"){
      krahn = netmeta:::nma.krahn(metaNetw,tau.preset = 0)
    }
    
    if(model=="random"){
      krahn=netmeta:::nma.krahn(metaNetw,tau.preset = metaNetw$tau)
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
  
  hatmatrix = contribution.hatmatrix(metaNetw, model)

  directs <- hatmatrix$colNames

  hatMatrix <- hatmatrix$H
  
  rownames(hatMatrix) <- hatmatrix$rowNames

  split <- function (dir) {strsplit(dir,":")}

  dims <- dim(hatMatrix)

#rows of comparison matrix 
  comparisons <- unlist(lapply(rownames(hatMatrix),unlist))

  comparisonToEdge <- function (comp) unlist (split(comp))

  directlist <- unlist(lapply(lapply(directs,split),unlist))

  edgeList <- matrix( directlist, nc = 2, byrow = TRUE)

  g <- graph_from_edgelist(edgeList , directed=FALSE)
  g <- set.vertex.attribute(g,'label',value = V(g))
  g <- set.edge.attribute(g,'label',value = E(g))

  setWeights <- function (g,comparison,conMat) {
    set.edge.attribute(g,"weight",value=rep(1,dims[2]))
  }


  getFlow <- function(g,edge) {return(E(g)[edge]$flow)}

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
    dedgeList <- matrix( unlist(dedgeList), nc = 2, byrow = TRUE)
    flows<-abs(hatMatrix[comparison,])
    dg <- graph_from_edgelist(dedgeList , directed = TRUE)
    E(dg)[]$weight <- rep(0,dims[2])
    E(dg)[]$flow <- abs(hatMatrix[comparison,])
    V(dg)[]$label <- V(dg)[]$name
    dg <- set.edge.attribute(dg,'label',value = E(dg))
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
       return (E(g)[e[]])
    })
    flow <- min(unlist(lapply(splE, function(e){
      return(e$flow[])
    })))
    gg <- Reduce(function(g, e){
      elabel <- e$label
      pfl <- e$flow[]
      g <- set.edge.attribute(g,"flow",e, pfl-flow)
      cw <-  e$weight[] + (flow[1]/pl) 
      weights[comparison,elabel] <<- cw
      return(set.edge.attribute(g,"weight",e, cw))
    },splE, g)
    emptyEdges <- Reduce(function(removedEdges, e){
      e <- E(gg)[e[]]
      if(e$flow[[1]][[1]]==0){
        removedEdges <- c(removedEdges, e)
      }
      return(removedEdges)
    },splE, c())
   return(delete_edges(gg, emptyEdges))
  }
  
 reduceGraph <- function (g,comparison) {
    getshortest <- function (g,compariston) {
      return(tryCatch({
        return(get.shortest.paths(g,sv(comparison),tv(comparison),mode="out",output="epath",weights=NA)$epath)
      }, error = function(e) {
        print(paste('error:', e))
      }, warning = function(w) {return({})}
      )
      )}
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
  totalSums <-colSums(weights)
  totalTotal <- sum(totalSums)
  totalWeights <- unlist(lapply(totalSums,function(comp){
                           100 * comp/ totalTotal
  }))
  return(weights)
}