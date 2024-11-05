#' Auxiliary functions for component network meta-analysis
#' 
#' @description
#' Auxiliary functions to (i) create a combination / C matrix with information
#' on treatment combinations and interaction terms and (ii) a vector with all
#' combinations in a component network meta-analysis.
#' 
#' @param x A \code{\link{netcomb}}, \code{\link{netmeta}} or
#'   \code{\link{netconnection}} object, a matrix or the number of components.
#' @param inactive A character string defining the inactive treatment
#'   component (see \code{\link{netcomb}}).
#' @param comb.ia A character vector specifying treatment combinations which
#'   will be considered as interactions.
#' @param sep.comps A single character to define separator between
#'   treatment components.
#' @param sep.ia A single character to define separator for interactions.
#' @param n A single number specifying the number of components in combinations.
#' @param \dots Additional arguments.
#'
#' @return
#' R function \code{createC} returns a combination matrix / C matrix with
#' studies in rows and component and interaction terms in columns.
#' 
#' R function \code{combinations} returns a character vector with combinations.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de},
#'   Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}}
#' 
#' @examples
#' data(Linde2016)
#'
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#'
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = face, ref = "placebo", sm = "OR", common = FALSE)
#'
#' # Additive component network meta-analysis (with placebo as inactive
#' # treatment)
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#'
#' # Available combinations in CNMA
#' combinations(nc1)
#'
#' # Create C matrix with all available interactions, i.e., one interaction
#' # for each combination
#' createC(nc1)
#' 
#' # Run interaction CNMA model with all available interactions
#' # (same result as standard NMA)
#' netcomb(net1, C.matrix = createC(nc1))
#' 
#' \dontrun{
#' # Conduct random effects network meta-analysis on full dataset
#' #
#' net2 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = Linde2016, ref = "placebo", sm = "OR", common = FALSE)
#'
#' # Additive component network meta-analysis (with placebo as inactive
#' # treatment)
#' #
#' nc2 <- netcomb(net2, inactive = "placebo")
#'
#' # Available combinations in CNMA
#' combinations(nc2)
#'
#' # Create C matrix with all available interactions, i.e., one interaction
#' # for each combination
#' head(createC(nc2))
#' 
#' # Run interaction CNMA model with all available interactions
#' # (same result as standard NMA)
#' netcomb(net2, C.matrix = createC(nc2))
#' }
#' 
#' @export createC


createC <- function(x, ...)
  UseMethod("createC")


#' @rdname createC
#' @method createC matrix
#' @export

createC.matrix <- function(x, comb.ia, inactive = NULL,
                           sep.comps = gs("sep.comps"),
                           sep.ia = gs("sep.ia"), ...) {
  if (!is.matrix(x))
    stop("Argument 'x' must be a matrix.", call. = FALSE)
  #
  comps.all <- colnames(x)
  comps <- comps.all[!grepl(sep.ia, comps.all, fixed = TRUE)]
  #
  chkchar(comb.ia, length = 1)
  #
  if (!is.null(inactive)) {
    inactive <- setchar(inactive, comps.all, stop.at.error = FALSE)
    comps.all <- comps.all[!(comps.all %in% inactive)]
    comps <- comps[!(comps %in% inactive)]
    x <- x[, comps.all, drop = FALSE]
  }
  #
  missing.sep.ia <- missing(sep.ia)
  chkchar(sep.ia, nchar = 0:1, length = 1)
  if (sep.comps == sep.ia)
    stop("Input for arguments 'sep.comps' and 'sep.ia' must be different.",
         call. = FALSE)
  sep.ia <- setsep(comps, sep.ia, missing = missing.sep.ia)
  #
  mat.int <- x
  #
  for (i in length(comb.ia)) {
    comb.ia.i <- unlist(compsplit(comb.ia[i], sep.comps))
    #
    if (length(comb.ia.i) == 1)
      stop("Entry ", i, " of argument 'comb.ia' contains a single component.")
    #
    if (any(grepl(sep.ia, comb.ia.i, fixed = TRUE)))
      stop("Input for argument 'comb.ia' contains interaction separator '",
           sep.ia, "'.")
    #
    comb.ia.i <- setchar(comb.ia.i, comps, name = "comb.ia",
                         text = paste0("must be a combination of ",
                                       paste0("'", comps, "'",
                                              collapse = ", "),
                                       "; separated by '", sep.comps,
                                       "' (argument 'sep.comps')"))
    #
    if (length(comb.ia.i) != length(unique(comb.ia.i)))
      stop("Identical components provided in argument 'comb.ia'.")
    #
    mat.int.i <- x[, comb.ia.i[1], drop = FALSE]
    #
    for (i in (1 + seq_len(length(comb.ia.i) - 1)))
      mat.int.i <- mat.int.i * x[, comb.ia.i[i], drop = FALSE]
    #
    colnames(mat.int.i) <- paste(comb.ia.i, collapse = paste0(" ", sep.ia, " "))
    #
    if (colnames(mat.int.i) %in% comps.all) {
      warning("Interaction '", colnames(mat.int.i),
              "' is already part of the C matrix.",
              call. = FALSE)
      mat.int.i <- NULL
    }
    #
    if (!is.null(mat.int.i)) {
      i <- 1
      while (i <= ncol(x)) {
        if (all(unname(x[, i]) == unname(mat.int.i[, ncol(mat.int.i)]))) {
          warning("Information on interaction '", colnames(mat.int.i),
                  "' is already available in column '",
                  colnames(x)[i], "'.")
          mat.int.i <- NULL
          i <- ncol(x) + 1
        }
        else
          i <- i + 1
      }
    }
    #
    if (!is.null(mat.int.i)) {
      if (apply(abs(mat.int.i), 2, sum) == 0) {
        warning("Interaction '",
                colnames(mat.int.i),
                "' is ignored as it is inestimable ",
                "(column with zeros in C matrix).",
                call. = FALSE)
        mat.int.i <- NULL
      }
      #
      mat.int <- cbind(mat.int, mat.int.i)
    }
  }
  #
  attr(mat.int, "inactive") <- inactive
  attr(mat.int, "sep.comps") <- sep.comps
  attr(mat.int, "sep.ia") <- sep.ia
  #
  mat.int
}


#' @rdname createC
#' @method createC netcomb
#' @export

createC.netcomb <- function(x, comb.ia = NULL, inactive = NULL,
                            sep.ia = x$sep.ia,
                            ...) {
  chkclass(x, "netcomb")
  x <- updateversion(x)
  #
  C.matrix <- x$C.matrix
  #
  comps.all <- colnames(C.matrix)
  comps <- comps.all[!grepl(sep.ia, comps.all, fixed = TRUE)]
  #
  if (is.null(comb.ia))
    comb.ia <- combinations(x)
  #
  if (!is.null(inactive) & !is.null(x$inactive))
    stop("Inactive treatment was already defined in component NMA provided in ",
         "argument 'x'.",
         call. = FALSE)
  #
  missing.sep.ia <- missing(sep.ia)
  chkchar(sep.ia, nchar = 0:1, length = 1)
  if (x$sep.comps == sep.ia)
    stop("Input for argument 'sep.ia' identical to separator for components",
         call. = FALSE)
  sep.ia <- setsep(comps, sep.ia, missing = missing.sep.ia)
  #
  for (i in comb.ia)
    C.matrix <- createC(C.matrix, i, sep.comps = x$sep.comps, sep.ia = sep.ia,
                        inactive = inactive)
  #
  C.matrix
}


#' @rdname createC
#' @method createC netmeta
#' @export

createC.netmeta <- function(x, inactive = NULL,
                            sep.comps = gs("sep.comps"), ...) {
  chkclass(x, "netmeta")
  #
  chkchar(sep.comps, length = 1)
  #
  C.matrix <-
    createC_trts_inactive(x$trts, inactive = inactive, sep.comps = sep.comps)
  #
  C.matrix
}


#' @rdname createC
#' @method createC netconnection
#' @export

createC.netconnection <- function(x, inactive = NULL,
                                  sep.comps = gs("sep.comps"), ...) {
  chkclass(x, "netconnection")
  #
  chkchar(sep.comps, length = 1)
  #
  C.matrix <-
    createC_trts_inactive(rownames(x$D.matrix),
                          inactive = inactive, sep.comps = sep.comps)
  #
  C.matrix
}


#' @rdname createC
#' @method createC default
#' @export

createC.default <- function(x, n = 1, ...) {
  if (missing(x))
    stop("Argument 'x' must be provided.", call. = TRUE)
  #
  chknumeric(x, min = 1, length = 1)
  chknumeric(n, min = 1, max = x, length = 1)
  #
  C.matrix <- createC_full(x, n)
  #
  C.matrix
}



#' @rdname createC
#' @export combinations

combinations <- function(x, n = NULL) {
  chkclass(x, "netcomb")
  #
  comps <- x$comps
  trts <- x$trts
  #
  combs <- trts[!(trts %in% comps)]
  #
  if (!is.null(x$inactive))
    combs <- combs[!(combs %in% x$inactive)]
  
  if (!is.null(n)) {
    if (!is_wholenumber(n))
      stop("Argument 'n' must be a whole number.")
    sel <-
      unlist(lapply(strsplit(combs, x$sep.comps, fixed = TRUE), length)) == n
    combs <- combs[sel]
  }
  #
  combs
}
