ranksampling <- function(x, nsim,
                         pooled = "random", small.values = "desirable",
                         keep.samples = FALSE) {
  
  chknumeric(nsim, min = 1, length = 1)
  pooled <- setchar(pooled, c("common", "random"))
  small.values <- setsv(small.values)
  chklogical(keep.samples)
  
  
  #
  # Generate samples
  #
  
  samples <- samples_netmeta(x, nsim, pooled)
  
  
  #
  # Calculate rankings
  #
  
  if (small.values == "desirable")
    rankings <- rankings(samples$samples)
  else
    rankings <- rankings(-samples$samples)
  
  
  #
  # Return results
  #
  
  res <- list(sucras = rankings$sucras,
              rankogram = rankings$rankogram,
              cumrank = rankings$cumrank,
              #
              nsim = nsim,
              pooled = pooled,
              small.values = small.values,
              keep.samples = keep.samples,
              #
              compMatrix = samples$compMatrix)
  #
  if (keep.samples)
    res$samples <- samples$samples
  res
}
