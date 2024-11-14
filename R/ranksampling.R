ranksampling <- function(x, nsim,
                         pooled = "random", small.values = "desirable",
                         keep.samples = FALSE) {
  
  #
  # Generate samples
  #
  
  samples <- samples_netmeta(x, nsim, pooled, small.values)
  
  
  #
  # Calculate rankings
  #
  
  rankings <- rankings(samples$samples)
  
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
