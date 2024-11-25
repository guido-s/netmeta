#' Merge pairwise object with additional data
#' 
#' @description
#' Merge pairwise object with additional data, e.g., information on
#' network connectivity.
#' 
#' @param x An object of class \code{pairwise}.
#' @param y A data frame or an object of class \code{netconnection}.
#' @param all.x A logical indicating whether to keep all observations from the
#'   pairwise object, i.e., also include observations not belonging to a
#'   subnetwork due to missing estimates or standard errors.
#' @param \dots Other arguments (passed on to \code{merge}).
#' 
#' @return An object of class \code{pairwise}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link[meta]{pairwise}}, \code{\link{netconnection}}
#' 
#' @examples
#' data(Woods2010)
#' 
#' # Transform data from long arm-based format to contrast-based
#' # format Argument 'sm' has to be used for odds ratio as summary
#' # measure; by default the risk ratio is used in the metabin
#' # function called internally.
#' #
#' p1 <- pairwise(treatment, event = r, n = N,
#'   studlab = author, data = Woods2010, sm = "OR")
#' head(p1)
#' 
#' # Add information on network connectivity
#' nc1 <- netconnection(p1)
#' p1nc1 <- merge(p1, nc1)
#' head(p1nc1)
#'
#' @method merge pairwise
#' @export


merge.pairwise <- function(x, y, all.x = TRUE, ...) {
  chkclass(x, "pairwise")
  #
  x$..order <- seq_len(nrow(x))
  xdat <- as.data.frame(x)
  #
  attribs <- attributes(x)
  attribs[["names"]] <- NULL
  attribs[["row.names"]] <- NULL
  #
  # Get rid of warning "no visible binding for global variable"
  #
  subnet <- design <- NULL
  #
  if (inherits(y, "netconnection")) {
    if (isCol(xdat, "subnet"))
      xdat <- rename(xdat, subnet.orig = subnet)
    #
    if (isCol(xdat, "design"))
      xdat <- rename(xdat, design.orig = design)
    #
    ydat <- as.data.frame(y)
    #
    res <- merge(xdat, ydat, by = c("studlab", "treat1", "treat2"),
                 all.x = all.x, ...)
  }
  else if (inherits(y, "data.frame"))
    res <- merge(xdat, y, all.x = all.x, ...)
  else
    stop("Argument 'y' must be a data frame or of class 'netconnection'.",
         call. = FALSE)
  #
  res <- res[order(res$..order), ]
  rownames(res) <- seq_len(nrow(res))
  res$..order <- NULL
  #
  for (i in names(attribs))
    attr(res, i) <- attribs[[i]]
  #
  res
}
