#' Calculate rankogram from treatment effect samples
#'
#' @description
#' This function calculates the probabilities of each treatment being
#' at each possible rank and the SUCRAs (Surface Under the Cumulative
#' RAnking curve) from a sample of treatment estimates in network
#' meta-analysis.
#'
#' @param x A matrix or data frame with treatment effects in columns and
#'  samples in rows.
#' @param pooled A character string indicating whether samples come from
#'   a common (\code{"common"}), random effects (\code{"random"}), or
#'   \code{"unspecified"} model, can be abbreviated.
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
#' The matrix / data frame in argument \code{x} must contain the sampled
#' effects for each treatment.
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
#'   \code{\link{plot.rankogram}},
#'   \code{\link[metadat]{dat.woods2010}},
#'   \code{\link[metadat]{dat.linde2015}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#'
#' @examples
#' p1 <- pairwise(treatment, event = r, n = N, studlab = author,
#'   data = dat.woods2010, sm = "OR")
#' net1 <- netmeta(p1, small.values = "desirable")
#'
#' set.seed(1909) # get reproducible results
#' ran1 <- rankogram(net1, nsim = 100, common = FALSE,
#'   keep.samples = TRUE)
#' ran1
#' 
#' rankogram(ran1$samples.random, pooled = "random")
#' 
#' \dontrun{
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' #
#' net2 <- netmeta(p2, common = FALSE, ref = "Placebo", small = "undesirable")
#' 
#' ran2 <- rankogram(net2, nsim = 100, common = FALSE,
#'   keep.samples = TRUE)
#' ran2
#' 
#' # Wrong ranking due to using the default,
#' # i.e., argument 'small.values = "desirable".
#' rankogram(ran2$samples.random, pooled = "random")
#' # Correct ranking 
#' rankogram(ran2$samples.random, pooled = "random",
#'   small.values = "undesirable")
#' }
#'
#' @rdname rankogram.default
#' @method rankogram default
#' @export

rankogram.default <- function(x, pooled = "unspecified",
                              small.values = "desirable",
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
  if (is.data.frame(x))
    x <- as.matrix(x)
  #
  if (any(!is.numeric(x)))
    stop("Input for argument 'x' must contain numeric values.",
         call. = FALSE)
  
  
  #
  #
  # (2) Check other arguments
  #
  #
  
  pooled <- setchar(pooled, c("common", "random", "unspecified"))
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
  
  if (small.values == "undesirable")
    x <- -x
  #
  sucras.common  <- ranking.matrix.common  <- cumrank.matrix.common  <- NULL
  sucras.random <- ranking.matrix.random <- rank.cum.random <- NULL
  #
  if (common) {
    res.c <- rankings(x)
    #
    sucras.common <- res.c$sucras
    ranking.matrix.common <- res.c$rankogram
    cumrank.matrix.common <- res.c$cumrank
    #
    samples.common <- x
    nsim <- res.c$nsim
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
              pooled = pooled,
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
