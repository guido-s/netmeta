#' Abbreviate treatment names
#' 
#' @description
#' Auxiliary functions to create uniquely abbreviated treatment names.
#' 
#' @details
#' These auxiliary functions can be used to create uniquely
#' abbreviated treatment names (and are used internally in several R
#' functions for this purpose).
#' 
#' In order to construct uniquely abbreviated treatment names,
#' \code{treats} uses \code{\link{substring}} to extract the first
#' \code{nchar.trts} characters. If these abbreviated treatment names
#' are not unique, \code{\link{abbreviate}} with argument
#' \code{minlength = nchar.trts} is used.
#'
#' In order to construct comparisons with uniquely abbreviated
#' treatment names, \code{comps} calls \code{treats} internally.
#' 
#' @param x A vector with treatment or comparison names or a matrix
#'   with treatment or comparison names as row and / or column names.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param row A logical indicating whether row or column names should
#'   be used (only considered if argument \code{x} is a matrix).
#' @param trts A character vector with treatment names.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{print.netmeta}},
#'   \code{\link{print.summary.netmeta}}
#' 
#' @examples
#' data(Senn2013)
#' #
#' net1 <- netmeta(TE, seTE, treat1.long, treat2.long,
#'                 studlab, data = Senn2013)
#' 
#' # Full treatment names
#' #
#' net1$trts
#' 
#' # Treatment names with four characters
#' #
#' treats(net1$trts, nchar.trts = 4)
#' 
#' # With two characters
#' #
#' treats(net1$trts, nchar.trts = 2)
#' 
#' # With one character (if possible)
#' #
#' treats(net1$trts, nchar.trts = 1)
#'
#' # Full comparison names
#' #
#' net1$comparisons
#' 
#' # Abbreviated comparison names
#' #
#' with(net1, comps(comparisons, trts, sep.trts, nchar = 4))
#' 
#' @export treats


treats <- function(x, nchar.trts = 8, row = TRUE) {
  
  meta:::chknumeric(nchar.trts, min = 1, length = 1)
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





#' @rdname treats
#' 
#' @export comps


comps <- function(x, trts, sep.trts, nchar.trts = 8, row = TRUE) {
  
  meta:::chknumeric(nchar.trts, min = 1, length = 1)
  meta:::chklogical(row)
  meta:::chkchar(sep.trts)
  ##
  trts.abbr <- treats(trts, nchar.trts)  
  
  charfac <- function(x, ...)
    as.character(factor(x, ...))
  
  if (is.matrix(x)) {
    if (row)
      comps <- rownames(x)
    else
      comps <- colnames(x)
  }
  else
    comps <- x
  
  trts.list <- compsplit(comps, sep.trts)
  trts.list.c <- lapply(trts.list, charfac, levels = trts, labels = trts.abbr)
  ##
  res <- sapply(trts.list.c, paste, collapse = sep.trts)

  res
}
