#' Abbreviate treatment names
#' 
#' @description
#' Auxiliary function to create uniquely abbreviated treatment names.
#' 
#' @details
#' This auxiliary function can be used to create uniquely abbreviated
#' treatment names (and is used internally in several R functions for
#' this purpose).
#' 
#' Initially, to construct uniquely abbreviated treatment names,
#' \code{\link{substring}} is used to extract the first
#' \code{nchar.trts} characters. If these abbreviated treatment names
#' are not unique, \code{\link{abbreviate}} with argument
#' \code{minlength = nchar.trts} is used.
#' 
#' @param x A vector with treatment names or a matrix with treatment
#'   names as row and / or column names.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param row A logical indicating whether row or column names should
#'   be used (only considered if argument \code{x} is a matrix).
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{print.netmeta}},
#'   \code{\link{print.summary.netmeta}}
#' 
#' @examples
#' data(Senn2013)
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2,
#'                 studlab, data = Senn2013)
#' 
#' # Use matrix with fixed effects estimates to create unique
#' # treatment names (with four characters)
#' #
#' treats(net1$TE.fixed, nchar.trts = 4)
#' 
#' # With two characters
#' #
#' treats(net1$TE.fixed, nchar.trts = 2)
#' 
#' # With one character
#' #
#' treats(net1$TE.fixed, nchar.trts = 1)
#' 
#' @export treats


treats <- function(x, nchar.trts = 8, row = TRUE) {
  
  meta:::chknumeric(nchar.trts, min = 1, single = TRUE)
  meta:::chklogical(row)
  
  if (is.matrix(x)) {
    if (row)
      trts <- rownames(x)
    else
      trts <- colnames(x)
  }
  else
    trts <- x
  ##
  ## Default: first 'nchar.trts' character of treatment names
  ##
  res <- substring(trts, 1, nchar.trts)
  ##
  ## Use abbreviated treatment names if necessary
  ##
  if (length(unique(res)) != length(unique(trts)))
    res <- as.character(abbreviate(trts, nchar.trts))
  ##
  res
}
