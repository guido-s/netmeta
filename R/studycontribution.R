# Gives contributions of individual studies for a single comparison
# arguments hatmatrix list from getHatMatrix and comparison the string
# identifying the hatmatrix row separated by ":"

getStudyContribution <- function(x, contribution.Matrix, comparison, model) {
  contribution <- contribution.Matrix[comparison, ]
  #
  if (model == "random")
    tau <- x$tau
  else
    tau <- 0
  # main data frame
  dfr <-
    data.frame(studlab = x$studlab,
               treat1 = x$treat1, treat2 = x$treat2,
               seTE.adj = x$seTE.adj)
  dfr$comp <- paste(dfr$treat1, dfr$treat2, sep = ":")
  dfr$w.adj <- 1 / (dfr$seTE.adj^2 + tau^2)
  
  studyContribution <- function(direct){
    aux <- dfr[dfr$comp == direct, ]
    normfac <- sum(aux[, "w.adj"])
    per <- contribution[direct]
    studycontr <-
      lapply(row.names(aux),
             function(row){
               w <- aux[row, "w.adj"]
               out <- data.frame(study = "", contribution = 0, comparison = "")
               out$study <- aux[row,"studlab"]
               out$contribution <- w * per / normfac
               out$comparison <- direct
               #
               out
             })
    #
    studycontr
  }
  
  outlist <- Reduce(
    function(acc, col) c(acc, studyContribution(col)),
  names(contribution),
  list(),
  accumulate = FALSE)
  
  SCR <- data.frame(
    study = sapply(outlist, function(r) r$study),
    contribution = sapply(outlist, function(r) r$contribution),
    comparison = sapply(outlist, function(r) r$comparison))
  #
  studyRow <-
    aggregate(SCR$contribution, by = list(study = SCR$study), FUN = sum)
  #
  studyRow <- cbind(comparison, studyRow) %>% rename(contribution = x)
  #
  studyRow
}
