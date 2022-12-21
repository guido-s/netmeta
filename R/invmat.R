#' Moore-Penrose Pseudoinverse of a Matrix
#' 
#' @description
#' Calculates the Moore-Penrose pseudoinverse of a square matrix
#' \strong{X}.
#' 
#' @param X A square matrix.
#' 
#' @details
#' This function is used by default in R package \strong{netmeta} to
#' calculate the Moore-Penrose pseudoinverse \strong{L\eqn{^+}} of the
#' Laplacian matrix \strong{L} (Rücker, 2012):
#' 
#' \strong{L\eqn{^+} = (X - J / \emph{n})\eqn{^{-1}} + J / \emph{n}}
#' with identity matrix \strong{J} of dimension \emph{n}x\emph{n}.
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
#' data(smokingcessation)
#' 
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1)
#' 
#' invmat(net1$L.matrix.common)
#' 
#' \dontrun{
#' data(Senn2013)
#' 
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013)
#' L1 <- net2$L.matrix.common
#' L2 <- invmat(net2$Lplus.matrix.common)
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
