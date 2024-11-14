#' Calculate rankogram from samples
#'
#' @description
#' This function calculates the probabilities of each treatment being
#' at each possible rank and the SUCRAs (Surface Under the Cumulative
#' RAnking curve) from a sample of treatment estimates in network
#' meta-analysis.
#'
#' @param x A matrix or data frame with samples.
#' @param pooled A character string indicating whether samples come from
#'   a common (\code{"common"}) or random effects model (\code{"random"}).
#'   Can be abbreviated.
#' @param small.values An optional character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated.
#' @param cumulative.rankprob A logical indicating whether cumulative
#'   ranking probabilities should be printed.
#' @param keep.samples A logical indicating whether to keep the generated
#'   samples.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' We derive a matrix showing the probability of each treatment being
#' at each possible rank. To this aim, we use samples and
#' summarise them using the ranking metric SUCRAs (Surface Under
#' the Cumulative RAnking curve).
#'
#' @return
#' An object of class \code{rankogram} with corresponding \code{print}
#' and \code{plot} function. The object is a list containing the
#' following components:
#' \item{ranking.matrix.common}{Numeric matrix giving the probability
#'   of each treatment being at each possible rank for the common
#'   effects model.}
#' \item{ranking.common}{SUCRA values for the common effects model.}
#' \item{ranking.matrix.random}{Numeric matrix giving the probability
#'   of each treatment being at each possible rank for the random
#'   effects model.}
#' \item{ranking.random}{SUCRA values for the random effects model.}
#' \item{cumrank.matrix.common}{Numeric matrix giving the cumulative
#'   ranking probability of each treatment for the
#'   common effects model.}
#' \item{cumrank.matrix.random}{Numeric matrix giving the cumulative
#'   ranking probability of each treatment for the random effects
#'   model.}
#' \item{nsim, common, random}{As defined above},
#' \item{small.values, x}{As defined above},
#'
#' @author Theodoros Papakonstantinou \email{dev@@tpapak.com}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{netmeta}}, \code{\link{netrank}},
#'   \code{\link{plot.rankogram}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#'
#' @examples
#' data(Woods2010)
#' p1 <- pairwise(treatment, event = r, n = N, studlab = author,
#'                data = Woods2010, sm = "OR")
#' net1 <- netmeta(p1, small.values = "desirable")
#'
#' ran1 <- rankogram(net1, nsim = 100, common = FALSE,
#'   keep.samples = TRUE)
#' ran1
#' 
#' rankogram(ran1$samples.random,
#'   pooled = "random", small.values = "d")
#'
#' @rdname rankogram.default
#' @method rankogram default
#' @export


rankogram.default <- function(x, pooled,
                              small.values = "",
                              cumulative.rankprob = FALSE,
                              keep.samples = FALSE,
                              nchar.trts = gs("nchar.trts"),
                              ...) {
  
  #
  #
  # (1) Check for matrix or data frame
  #
  #
  
  if (!inherits(x, c("matrix", "data.frame")))
    stop("Argument 'x' must be a matrix or data frame.",
         call. = FALSE)
  #
  if (any(!is.numeric(x)))
    stop("Input for argument 'x' must contain numeric values.",
         call. = FALSE)

  
  #
  #
  # (2) Check other arguments
  #
  #
  
  pooled <- setchar(pooled, c("common", "random"))
  small.values <- setsv(small.values)
  chklogical(cumulative.rankprob)
  chklogical(keep.samples)
  #
  if (is.null(nchar.trts))
    nchar.trts <- 666
  else
    chknumeric(nchar.trts, length = 1)
  #
  common <- pooled == "common"
  random <- !common
  
  
  #
  #
  # (3) Resampling to calculate ranking probabilites and SUCRAs
  #
  #
  
  sucras.common  <- ranking.matrix.common  <- cumrank.matrix.common  <- NULL
  sucras.random <- ranking.matrix.random <- rank.cum.random <- NULL
  #
  if (common) {
    res.f <- rankings(x)
    #
    sucras.common <- res.f$sucras
    ranking.matrix.common <- res.f$rankogram
    cumrank.matrix.common <- res.f$cumrank
    #
    samples.common <- x
    nsim <- res.f$nsim
  }
  #
  if (random) {
    res.r <- rankings(x)
    #
    sucras.random <- res.r$sucras
    ranking.matrix.random <- res.r$rankogram
    rank.cum.random <- res.r$cumrank
    #
    samples.random <- x
    nsim <- res.r$nsim
  }
  
  
  #
  #
  # (4) Create rankogram object
  #
  #
  
  res <- list(ranking.common = sucras.common,
              ranking.matrix.common = ranking.matrix.common,
              cumrank.matrix.common = cumrank.matrix.common,
              samples.common =
                if (common & keep.samples) samples.common else NULL,
              #
              ranking.random = sucras.random,
              ranking.matrix.random = ranking.matrix.random,
              cumrank.matrix.random = rank.cum.random,
              samples.random =
                if (random & keep.samples) samples.random else NULL,
              #
              nsim = nsim,
              #
              common = common,
              random = random,
              small.values = small.values,
              cumulative.rankprob = cumulative.rankprob,
              #
              nchar.trts = nchar.trts,
              #
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  #
  # Backward compatibility
  #
  res$fixed <- common
  #
  res$ranking.fixed <- sucras.common
  res$ranking.matrix.fixed <- ranking.matrix.common
  res$cumrank.matrix.fixed <- cumrank.matrix.common
  #
  class(res) <- "rankogram"
  
  res
}
