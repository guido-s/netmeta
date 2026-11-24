#' Moore-Penrose Pseudoinverse of a Matrix
#' 
#' @description
#' Calculates the Moore-Penrose pseudoinverse of a square matrix \strong{X}.
#' 
#' @param X A square matrix.
#' 
#' @details
#' Since version 3.9-0, \bold{netmeta} uses
#' \code{\link[BigDataStatMeth]{pseudoinverse}} from R package
#' \bold{BigDataStatMeth} to calculate Moore-Penrose pseudoinverses.
#' Older versions used \code{invmat()}.
#'
#' For the Laplacian matrix \strong{L} of a connected network,
#' \code{invmat} calculates its Moore-Penrose pseudoinverse as
#' (Rücker, 2012):
#'
#' \strong{L\eqn{^+} = (L - J / \emph{n})\eqn{^{-1}} + J / \emph{n}},
#'
#' where \strong{J} is the \emph{n}x\emph{n} matrix of ones.
#' 
#' @return
#' The Moore-Penrose pseudoinverse for matrix \strong{X}.
#' 
#' @author
#' Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido Schwarzer
#'   \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{solve}}
#' 
#' @references
#' Rücker G (2012):
#' Network meta-analysis, electrical networks and graph theory.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 312--24
#' 
#' @examples
#' \donttest{
#' data(smokingcessation)
#' 
#' # Only consider first five studies (to reduce runtime of example)
#' #
#' smoke5 <- head(smokingcessation, n = 5)
#' 
#' # Transform data from arm-based format to contrast-based format
#' #
#' pw <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smoke5, sm = "OR")
#' 
#' nma <- netmeta(pw)
#' 
#' L1 <- nma$L.matrix.common
#' L2 <- invmat(nma$Lplus.matrix.common)
#' all.equal(round(L1, 10), round(L2, 10))
#' }
#' 
#' @export invmat

invmat <- function(X) {
  n <- nrow(X)
  m <- ncol(X)
  ##
  if (n != m)
    stop("Argument 'X' must be a square matrix", call. = FALSE)
  ##
  J <- matrix(1, nrow = n, ncol = n)
  ##
  res <- solve(X - J / n) + J / n
  ##
  res
}
