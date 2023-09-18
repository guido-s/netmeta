#' Table with network meta-analysis results
#' 
#' @description
#' Construct a table with network, direct and indirect estimates from
#' one or more network meta-analyses.
#' 
#' @aliases nettable print.nettable
#' 
#' @param \dots Any number of network meta-analysis objects or a
#'   single list with network meta-analyses.
#' @param x An object of class \code{nettable}.
#' @param name An optional character vector providing descriptive
#'   names for network meta-analysis objects.
#' @param method A character string indicating which method to split
#'   direct and indirect evidence is to be used. Either
#'   \code{"Back-calculation"}, \code{"Edge-splitting"} or
#'   \code{"SIDDE"}, can be abbreviated. See Details.
#' @param order A optional character or numerical vector specifying
#'   the order of treatments in comparisons.
#' @param common A logical indicating whether table for the common
#'   effects network meta-analysis should be printed.
#' @param random A logical indicating whether table for the random
#'   effects network meta-analysis should be printed.
#' @param upper A logical indicating whether treatment comparisons
#'   should be selected from the lower or upper triangle of the
#'   treatment effect matrices (see list elements \code{TE.common} and
#'   \code{TE.random} in the \code{netmeta} object). Ignored if
#'   argument \code{order} is provided.
#' @param reference.group Reference treatment. Ignored if argument
#'   \code{order} is provided.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment or vice versa. This argument is only
#'   considered if \code{reference.group} is not equal to \code{""}
#'   and argument\code{order} is not provided.
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. For example, if \code{backtransf =
#'   TRUE}, results for \code{sm = "OR"} are printed as odds ratios
#'   rather than log odds ratios.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of test of agreement between direct and indirect evidence, see
#'   \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param big.mark A character used as thousands separator.
#' @param text.NA A character string specifying text printed for
#'   missing values.
#' @param bracket A character with bracket symbol to print lower
#'   confidence interval: "[", "(", "\{", "".
#' @param separator A character string with information on separator
#'   between lower and upper confidence interval.
#' @param lower.blank A logical indicating whether blanks between left
#'   bracket and lower confidence limit should be printed.
#' @param upper.blank A logical indicating whether blanks between
#'   separator and upper confidence limit should be printed.
#' @param tol.direct A numeric defining the maximum deviation of the
#'   direct evidence proportion from 0 or 1 to classify a comparison
#'   as providing only indirect or direct evidence, respectively.
#' @param writexl A logical indicating whether an Excel file should be
#'   created (R package \bold{writexl} must be available).
#' @param path A character string specifying the filename of the Excel
#'   file.
#' @param overwrite A logical indicating whether an existing Excel
#'   file should be overwritten.
#' @param warn A logical indicating whether warnings should be
#'   printed.
#' @param verbose A logical indicating whether progress information
#'   should be printed.
#' 
#' @details
#' Construct a table with network, direct and indirect estimates from
#' one or more network meta-analyses. The table looks very similar to
#' the statistical part of a GRADE table for a network meta-analysis
#' (Puhan et al., 2014).
#'
#' By default, an R object with the network tables is
#' generated. Alternatively, an Excel file is created if argument
#' \code{writexl = TRUE}.
#' 
#' Three methods to derive indirect estimates are available:
#' \itemize{
#' \item Separate Indirect from Direct Evidence (SIDE) using a
#'   back-calculation method (\code{method = "Back-calculation"})
#'   based on the \emph{direct evidence proportion} to calculate the
#'   indirect evidence (König et al., 2013);
#' \item Separate Indirect from Direct Evidence (SIDE) using
#'   node-splitting method in Dias et al. (2010) (\code{method =
#'   "Edge-splitting"});
#' \item Separate Indirect from Direct Design Evidence (SIDDE) as
#'   described in Efthimiou et al. (2019).
#' }
#' 
#' Note, for the back-calculation method, indirect treatment estimates
#' are already calculated in \code{\link{netmeta}} and this function
#' combines and prints these estimates in a user-friendly
#' way. Furthermore, this method is not available for the
#' Mantel-Haenszel and non-central hypergeometric distribution
#' approach implemented in \code{\link{netmetabin}}.
#'
#' Dias et al. (2010) used the term "node-splitting" method, however,
#' the method actually does not split nodes, i.e., treatments, but
#' edges, i.e., comparisons. Accordingly, we use the term
#' "side-splitting" method.
#' 
#' For the random-effects model, the direct treatment estimates are
#' based on the common between-study variance \eqn{\tau^2} from the
#' network meta-analysis, i.e. the square of list element
#' \code{x$tau}.
#'
#' The SIDDE approach can be compute-intensive in large
#' networks. Crude information on the computation progress is printed
#' for SIDDE if argument \code{verbose} is \code{TRUE}.
#'
#' @return
#' An object of class \code{nettable} with corresponding \code{print}
#' function if argument \code{writexl = FALSE}. The object is a list
#' containing the network tables in list elements 'common' and
#' 'random'. An Excel file is created if \code{writexl = TRUE}. In
#' this case, \code{NULL} is returned in R.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netsplit}}, \code{\link{netmeta}},
#'   \code{\link{netmetabin}}, \code{\link{netmeasures}}
#' 
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE (2010):
#' Checking consistency in mixed treatment comparison meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{29}, 932--44
#' 
#' Efthimiou O, Rücker G, Schwarzer G, Higgins J, Egger M, Salanti G
#' (2019):
#' A Mantel-Haenszel model for network meta-analysis of rare events.
#' \emph{Statistics in Medicine},
#' \bold{38}, 2992--3012
#' 
#' König J, Krahn U, Binder H (2013):
#' Visualizing the flow of evidence in network meta-analysis and
#' characterizing mixed treatment comparisons.
#' \emph{Statistics in Medicine},
#' \bold{32}, 5414--29
#' 
#' Puhan MA, Schünemann HJ, Murad MH, et al. (2014):
#' A GRADE working group approach for rating the quality of treatment
#' effect estimates from network meta-analysis.
#' \emph{British Medical Journal},
#' \bold{349}, g5630
#' 
#' @examples
#' data(Woods2010)
#' #
#' p1 <- pairwise(treatment, event = r, n = N,
#'   studlab = author, data = Woods2010, sm = "OR")
#' #
#' net1 <- netmeta(p1)
#' #
#' nt1 <- nettable(net1, digits = 2)
#' nt1
#' print(nt1, common = FALSE)
#' print(nt1, random = FALSE)
#' 
#' \dontrun{
#' # Create a CSV file with network table from random effects model
#' #
#' table1 <- nettable(net1, digits = 2, bracket = "(", separator = " to ")
#' #
#' write.table(table1$random, file = "table1-random.csv",
#'   row.names = FALSE, col.names = TRUE, sep = ",")
#' #
#' # Create Excel files with network tables
#' # (if R package writexl is available)
#' #
#' nettable(net1, digits = 2, bracket = "(", separator = " to ",
#'          path = tempfile(fileext = ".xlsx"))
#' }
#' 
#' @rdname nettable
#' @export nettable


nettable <- function(...,
                     name = NULL,
                     method = NULL,
                     ##
                     order = NULL, common, random,
                     ##
                     upper = TRUE,
                     reference.group = NULL,
                     baseline.reference = NULL,
                     ##
                     backtransf = NULL,
                     ##
                     digits = gs("digits"),
                     digits.I2 = gs("digits.I2"),
                     digits.pval = gs("digits.pval"),
                     ##
                     scientific.pval = gs("scientific.pval"),
                     zero.pval = gs("zero.pval"),
                     JAMA.pval = gs("JAMA.pval"),
                     ##
                     big.mark = gs("big.mark"),
                     text.NA = ".",
                     ##
                     bracket = gs("CIbracket"),
                     separator = gs("CIseparator"),
                     lower.blank = gs("CIlower.blank"),
                     upper.blank = gs("CIupper.blank"),
                     ##
                     tol.direct = 0.0005,
                     ##
                     writexl = !missing(path),
                     path = "nettable.xlsx",
                     overwrite = FALSE,
                     ##
                     warn = FALSE,
                     verbose = FALSE) {
  
  
  ##
  ##
  ## (1) Extract list elements and basic checks
  ##
  ##
  is.nma <- function(x)
    inherits(x, "netmeta")
  ##
  missing.common <- missing(common)
  missing.random <- missing(random)
  ##
  args <- list(...)
  ##
  n.netmeta <- length(args)
  n.i <- seq_len(n.netmeta)
  ##
  if (length(args) == 1) {
    if (!is.list(args[[1]]))
      stop("All elements of argument '...' must be of class 'netmeta'.",
           call. = FALSE)
    ##
    if (!is.nma(args[[1]])) {
      n.netmeta <- length(args[[1]])
      n.i <- seq_len(n.netmeta)
      ##
      args2 <- list()
      fix <- ran <- rep_len(NA, n.netmeta)
      for (i in n.i) {
        args2[[i]] <- args[[1]][[i]]
        fix[i] <- args[[1]][[i]]$common
        ran[i] <- args[[1]][[i]]$random
      }
      if (missing.common)
        common <- any(fix, na.rm = TRUE)
      if (missing.random)
        random <- any(ran, na.rm = TRUE)
      ##
      args <- args2
    }
    else {
      if (missing.common)
        common <- args[[1]]$common
      if (is.null(common))
        common <- args[[1]]$fixed
      if (is.null(common))
        common <- args[[1]]$comb.fixed
      if (missing.random)
        random <- args[[1]]$random
      if (is.null(random))
        random <- args[[1]]$comb.random
    }
  }
  ##
  for (i in n.i) {
    if (!is.nma(args[[i]]))
      stop("All elements of argument '...' must be of class 'netmeta'.",
           call. = FALSE)
    ##
    args[[i]] <- updateversion(args[[i]])
  }
  ##
  levs <- numeric(0)
  for (i in n.i)
    levs[i] <- args[[i]]$level.ma
  ##
  if (length(unique(levs)) != 1)
    stop("Different confidence levels used in network meta-analyses ",
         "(see list element 'level.ma').",
         call. = FALSE)
  ##
  if (n.netmeta > 1 & (missing.common | missing.random)) {
    fix <- ran <- rep_len(NA, n.netmeta)
    for (i in n.i) {
      if (is.nma(args[[i]])) {
        fix[i] <- args[[i]]$common
        ran[i] <- args[[i]]$random
      }
    }
    if (missing.common)
      common <- any(fix, na.rm = TRUE)
    if (missing.random)
      random <- any(ran, na.rm = TRUE)
  }
  ##
  sms <- character(0)
  for (i in n.i)
    sms[i] <- args[[i]]$sm
  ##
  if (length(unique(sms)) == 1)
    sms <- unique(sms)
  ##
  backtransfs <- logical(0)
  for (i in n.i)
    backtransfs[i] <- args[[i]]$backtransf
  ##
  if (length(unique(backtransfs)) == 1)
    backtransfs <- unique(backtransfs)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  if (!is.null(name))
    chklength(name, n.netmeta,
              paste("Argument 'name' must be of same length as",
                    "number of network meta-analyses"))
  else if (n.netmeta > 1)
    name <- paste("Outcome", n.i)
  ##
  if (!is.null(name))
    for (i in n.i)
      args[[i]]$outcome.name <- name[i]
  ##
  if (!missing(method))
    method <- setchar(method, c("Back-calculation", "Edge-splitting", "SIDDE"))
  ##
  chklogical(common)
  chklogical(random)
  ##
  chklogical(upper)
  if (!is.null(baseline.reference))
    chklogical(baseline.reference)
  ##
  chknumeric(tol.direct, min = 0, length = 1)
  if (!is.null(backtransf))
    chklogical(backtransf)
  ##
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  ##
  chklogical(writexl)
  chkchar(path, length = 1)
  chklogical(overwrite)
  ##
  chklogical(warn)
  chklogical(verbose)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  ##
  chkchar(text.NA)
  chkchar(big.mark)
  ##
  bracket.old <- gs("CIbracket")
  separator.old <- gs("CIseparator")
  lower.blank.old <- gs("CIlower.blank")
  upper.blank.old <- gs("CIupper.blank")
  ##
  cilayout(bracket, separator, lower.blank, upper.blank)
  on.exit(cilayout(bracket.old, separator.old,
                   lower.blank.old, upper.blank.old))
  
  
  ##
  ##
  ## (3) Generate network tables
  ##
  ##
  vars <- c(if (!is.null(name)) "Outcome",
            "Arm 1", "Arm 2", "k", "n", "I2",
            "Direct estimate", "Indirect estimate",
            "Network meta-analysis", "Incoherence")
  ##
  table.common <- table.random <-
    data.frame(matrix(nrow = 0, ncol = length(vars)))
  colnames(table.common) <- colnames(table.random) <- vars
  ##
  for (i in n.i) {
    table.i <- 
      nettable_internal(args[[i]], method,
                        upper, reference.group, baseline.reference,
                        order, tol.direct, backtransf,
                        digits, digits.I2, digits.pval,
                        scientific.pval, zero.pval, JAMA.pval,
                        big.mark, text.NA,
                        bracket, separator, lower.blank, upper.blank,
                        writexl,
                        warn, verbose)
    ##
    table.common <- rbind(table.common, table.i$common)
    table.random <- rbind(table.random, table.i$random)
  }
  ##
  if (all(is.na(table.common$n)))
    table.common$n <- NULL
  else
    table.common$n <- formatN(round(table.common$n), digits = 0,
                              text.NA = text.NA, big.mark = big.mark)
  ##
  if (all(is.na(table.random$n)))
    table.random$n <- NULL
  else
    table.random$n <- formatN(round(table.random$n), digits = 0,
                              text.NA = text.NA, big.mark = big.mark)
  
  
  ##
  ##
  ## (4) Save Excel file
  ##
  ##
  if (writexl) {
    if (!(common | random)) {
      warning("Excel file not generated as neither ",
              "argument 'common' nor 'random' is TRUE.")
      return(invisible(NULL))
    }
    ##
    if (!is.installed.package("writexl", stop = FALSE))
      stop(paste0("Package 'writexl' missing.",
                  "\n  ",
                  "Please use the following R command for installation:",
                  "\n  install.packages(\"writexl\")"),
           call. = FALSE)
    ##
    if (file.exists(path) & !overwrite)
      warning("File '", path, "' exists. ",
              "Use argument 'overwrite = TRUE' to overwrite file.",
              call. = FALSE)
    else {
      if (common & random)
        xlsx <- list(common = table.common, random = table.random)
      else if (common)
        xlsx <- list(common = table.common)
      else
        xlsx <- list(random = table.random)
      ##
      writexl::write_xlsx(xlsx, path = path, col_names = TRUE)
      message(paste0("Network table", if (common & random) "s",
                     " saved in file '", path, "'."))
    }
    ##
    return(invisible(NULL))
  }
  
  
  ##
  ##
  ## (5) Return network tables
  ##
  ##
  res <- list(common = table.common,
              random = table.random,
              ##
              upper = upper,
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              order = order,
              ##
              tol.direct = tol.direct,
              ##
              digits = digits,
              digits.I2 = digits.I2,
              digits.pval = digits.pval,
              ##
              scientific.pval = scientific.pval,
              zero.pval = zero.pval,
              JAMA.pval = JAMA.pval,
              ##
              big.mark = big.mark,
              text.NA = text.NA,
              ##
              bracket = bracket,
              separator = gs("CIseparator"),
              lower.blank = separator,
              upper.blank = upper.blank,
              ##
              x = list(common = common, random = random),
              ##
              backtransf = backtransfs,
              sm = sms,
              level.ma = unique(levs),
              ##
              version = packageDescription("netmeta")$Version
              )
  ##
  ## Backward compatibility
  ##
  res$fixed <- res$common
  res$x$fixed <- res$x$common
  ##
  class(res) <- "nettable"
  
  res
}





#' @rdname nettable
#' @method print nettable
#' @export


print.nettable <- function(x, common = x$x$common, random = x$x$random, ...) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkclass(x, "nettable")
  x <- updateversion(x)
  ##  
  ## All individual results in a single row - be on the save side:
  ##
  oldopts <- options(width = 200)
  on.exit(options(oldopts))
  ##
  args <- list(...)
  common <- deprecated(common, missing(common), args, "fixed", FALSE)
  chklogical(common)
  chklogical(random)
  
  
  ##
  ##
  ## (2) Print network table for common effects model
  ##
  ##
  if (common) {
    cat(paste0("Network table (", gs("text.w.common"),
               ") effects model):\n"))
    ##
    if (isCol(x$common, "Outcome")) {
      outcomes <- unique(x$common$Outcome)
      n.netmeta <- length(outcomes)
      backtransf <- x$backtransf
      if (length(backtransf) != n.netmeta)
        backtransf <- rep(backtransf, n.netmeta)
      for (i in seq_len(n.netmeta)) {
        mat.i <- x$common[x$common$Outcome == outcomes[i],
                          names(x$common) != "Outcome"]
        outcome.txt <- paste0("\nOutcome: ", outcomes[i])
        if (n.netmeta > 1 & length(x$sm) != 1)
          outcome.txt <-
            paste0(outcome.txt," (sm = '",
                   if (is.relative.effect(x$sm[i]) & !backtransf[i]) "log",
                   x$sm[i], "')")
        cat(paste0(outcome.txt, "\n"))
        prmatrix(mat.i,
                 quote = FALSE, right = TRUE,
                 rowlab = rep("", nrow(mat.i)), ...)
      }
    }
    else {
      cat("\n")
      prmatrix(x$common, quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(x$common)), ...)
    }
    if (random)
      cat("\n")
  }
  
  
  ##
  ##
  ## (3) Print network table for random effects model
  ##
  ##
  if (random) {
    cat(paste0("Network table (", gs("text.w.random"),
               ") effects model):\n"))
    ##
    if (isCol(x$random, "Outcome")) {
      outcomes <- unique(x$random$Outcome)
      n.netmeta <- length(outcomes)
      backtransf <- x$backtransf
      if (length(backtransf) != n.netmeta)
        backtransf <- rep(backtransf, n.netmeta)
      for (i in seq_len(n.netmeta)) {
        mat.i <- x$random[x$random$Outcome == outcomes[i],
                          names(x$random) != "Outcome"]
        outcome.txt <- paste0("\nOutcome: ", outcomes[i])
        if (n.netmeta > 1 & length(x$sm) != 1)
          outcome.txt <-
            paste0(outcome.txt," (sm = '",
                   if (is.relative.effect(x$sm[i]) & !backtransf[i]) "log",
                   x$sm[i], "')")
        cat(paste0(outcome.txt, "\n"))
        prmatrix(mat.i,
                 quote = FALSE, right = TRUE,
                 rowlab = rep("", nrow(mat.i)), ...)
      }
    }
    else {
      cat("\n")
      prmatrix(x$random, quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(x$random)), ...)
    }
  }
  
  
  invisible(NULL)
}
