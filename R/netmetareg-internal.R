nmr_results <- function(x,
                        covar = x$.netmeta$covar,
                        covar.name = x$.netmeta$covar.name) {
  
  # Link abbreviated treatment names with original treatment names
  #
  trts <- x$.netmeta$trts
  trts.abbr <- x$.netmeta$trts.abbr
  #
  assumption <- x$.netmeta$assumption
  reference.group <- x$.netmeta$reference.group
  
  # Get rid of warning "no visible binding for global variable"
  #
  treat <- treat_temp <- NULL
  
  dat <- data.frame(treat = "",
                    coef = x$b, se = x$se,
                    lower = x$ci.lb, upper = x$ci.ub,
                    z = x$zval, pval = x$pval)
  #
  # Remove back ticks which were created in the manual model matrix
  #
  rownames(dat) <- gsub("^`|`$", "", rownames(dat))
  #
  dat$treat <- rnam <- rownames(dat)
  #
  # Data set with main effects
  #
  sel_d <- !grepl(":", rnam)
  dat_d <- dat[sel_d, , drop = FALSE]
  dat_d$type <- "d"
  #
  rownames(dat_d) <- paste0("d[", rnam[sel_d], "]")
  #
  if (is.null(covar)) {
    return(dat_d)
  }
  else {
    is_num <- is.numeric(covar)
  }
  #
  # Add variables relevant for interaction terms
  #
  dat_d$covar <- NA
  #
  if (!is_num) {
    dat_d$cov_lvl <- NA
    dat_d$cov_ref <- NA
  }
  #
  # Data set with interaction effects
  #
  dat_beta <- dat[!sel_d, , drop = FALSE]
  #
  # Split 'rownames' using strsplit()
  #
  split_beta <- strsplit(rownames(dat_beta), "[:]")
  #
  trts_beta <- sapply(split_beta, first)
  lvls_beta <- sapply(split_beta, second)
  #
  dat_beta$type <- "beta"
  dat_beta$treat <- trts_beta
  dat_beta$covar <- covar.name
  #
  if (is_num) {
    row.names(dat_beta) <- paste0("beta[", rnam[!sel_d], "]")
  }  else {
    # Add original treatment names to merge results
    # (for the common assumption only the covariates observed for the
    #  reference treatment are relevant)
    #
    dat_beta %<>%
      mutate(treat_temp = ifelse(treat == "nonref", reference.group, treat),
             treat_long =
               as.character(factor(treat_temp,
                                   levels = trts.abbr, labels = trts))) %>%
      select(-treat_temp)
    #
    # Get observed covariate levels for each treatment
    #
    dat_cov <- x$.netmeta$x$data[, c("treat1", "treat2", covar.name)]
    #
    # Convert to long format
    #
    unique_lvls <- rbind(setNames(dat_cov[, c("treat1", covar.name)],
                                  c("treat_long", "lvl")),
                         setNames(dat_cov[, c("treat2", covar.name)],
                                  c("treat_long", "lvl")))
    #
    unique_lvls <- unique(unique_lvls)
    # Add treatment abbreviations
    unique_lvls <- merge(unique_lvls,
                         data.frame(treat_long = trts,
                                    treat = trts.abbr),
                         by = "treat_long", all.x = TRUE)
    #
    if (assumption == "common") {
      unique_lvls$treat <-
        ifelse(unique_lvls$treat == reference.group,
               "nonref",
               unique_lvls$treat)
    }
    #
    # Capture estimated factor levels
    #
    dat_beta$lvl <- gsub(covar.name, "", lvls_beta)
    #
    # Get the non-included covariate levels aka covariate references
    #
    ref_lvls <-
      merge(unique_lvls, dat_beta,
            by = c("treat_long", "lvl", "treat"),
            all.x = TRUE)
    #
    ref_lvls <- ref_lvls[is.na(ref_lvls$coef), c("treat_long", "lvl", "treat")]
    colnames(ref_lvls)[colnames(ref_lvls) == "lvl"] <- "cov_ref"
    #
    # Should be 1 to many. If there is an issue then it is likely in rma.mv()
    #
    dat_beta <-
      merge(dat_beta, ref_lvls,
            by = c("treat_long", "treat"), all.x = TRUE)
    colnames(dat_beta)[colnames(dat_beta) == "lvl"] <- "cov_lvl"
    #
    # Remove variable treat_long which we only used for merging
    #
    dat_beta <- dat_beta[, names(dat_beta) != "treat_long"]
    #
    rnam_beta <- paste0("beta[", rnam[!sel_d], "]")
    #
    for (i in seq_along(dat_beta$cov_lvl)) {
      rnam_beta[i] <- gsub(dat_beta$cov_lvl[i], "", rnam_beta[i])
    }
    #
    row.names(dat_beta) <- rnam_beta
  }
  #
  res <- rbind(dat_d, dat_beta)
  #
  attr(res, "reference.group") <- reference.group
  attr(res, "covar.name") <- covar.name
  #
  res <- res[, c("type", "treat",
                 if (!is_num) "cov_lvl", if (!is_num) "cov_ref",
                 "coef", "se", "lower", "upper", "z", "pval")]
  #
  res
}

nmr_full_results <- function(x) {
  
  # Extract NMR attributes/assumptions
  #
  assumption <- x$.netmeta$assumption
  reference.group <- x$.netmeta$reference.group
  consistency <- x$.netmeta$consistency
  covar <- x$.netmeta$covar
  # The issues with factors with more than two levels are the same for
  # character covariates with more than two levels
  #
  covar_is_factor <- is.factor(covar) | is.character(covar)
  #
  sep.trts <- x$.netmeta$x$sep.trts
  #
  # Extract estimates from rma.mv()
  #
  results <- x$results
  
  # Get rid of warning "no visible binding for global variable"
  #
  comparison <- cov <- cov.default <- interim.x <- interim.y <-
    ref.x <- ref.y <- treat <- type <- value <- treat1 <- treat2 <- NULL
  
  #
  # Generate matrix combinations for consistent parameters
  #
  sel.d <- results$type == "d"
  mat_d <- outer(c(0, results$coef[sel.d]), c(0, results$coef[sel.d]), "-")
  rownames(mat_d) <- colnames(mat_d) <- c(reference.group, results$treat[sel.d])
  #
  Cov <- vcov(x)
  #
  Cov_d <- Cov[sel.d, sel.d, drop = FALSE]
  Cov_d <- cbind(0, Cov_d)
  Cov_d <- rbind(0, Cov_d)
  rownames(Cov_d)[1] <- colnames(Cov_d)[1] <- reference.group
  #
  mat_se.d <- sqrt(outer(diag(Cov_d), diag(Cov_d), "+") - 2 * Cov_d)
  #
  # Convert to data frame
  #
  dat_d <- dat_se.d <- NULL
  #
  for (i in colnames(mat_d)) {
    # Data set with estimates
    #
    dat_d.i <- as.data.frame(mat_d)[i]
    names(dat_d.i) <- "d"
    #
    dat_d.i %<>%
      mutate(comparison = paste(rownames(dat_d.i), i, sep = sep.trts)) %>%
      select("comparison", "d")
    #
    dat_d <- rbind(dat_d, dat_d.i)
    #
    # Data set with standard errors
    #
    dat_se.d.i <- as.data.frame(mat_se.d)[i]
    names(dat_se.d.i) <- "se.d"
    #
    dat_se.d.i %<>%
      mutate(comparison = paste(rownames(dat_se.d.i), i, sep = sep.trts)) %>%
      select("comparison", "se.d")
    #
    dat_se.d <- rbind(dat_se.d, dat_se.d.i)
  }
  #
  rownames(dat_d) <- seq_len(nrow(dat_d))
  rownames(dat_se.d) <- seq_len(nrow(dat_se.d))
  
  #
  # Covariance for interaction terms
  #
  if (is.null(covar) | (covar_is_factor & length(unique(covar)) > 2)) {
    # We don't have a solution for factor / categorical covariates with more
    # than two levels. We lose the information on comparator / reference group
    # assignment which will later be important for more than two levels
    #
    return(merge(dat_d, dat_se.d, by = "comparison", all.x = TRUE))
  }
  else if (consistency) {
    if (assumption == "common") {
      # The common assumption has one beta per treatment
      #
      mat_beta <-
        outer(c(0, rep(results$coef[!sel.d], sum(sel.d))),
              c(0, rep(results$coef[!sel.d], sum(sel.d))), "-")
      rownames(mat_beta) <- colnames(mat_beta) <-
        c(reference.group, results$treat[sel.d])
      #
      # Single beta coefficient, i.e., no covariances
      #
      var_beta <- Cov[!sel.d, !sel.d]
      #
      mat_se.beta <-
        sqrt(outer(c(0, rep(var_beta, sum(sel.d))),
                   c(0, rep(var_beta, sum(sel.d))),
                   FUN = function(x, y) abs(x - y)))
      rownames(mat_se.beta) <- colnames(mat_se.beta) <-
        c(reference.group, results$treat[sel.d])
      #
      # Covariances
      #
      Cov_d_beta <- Cov[!sel.d, sel.d]
      #
      all_se.beta <- dat_d %>% select(comparison) %>% mutate(cov.default = 0)
      #
      dat_cov_d_beta <-
        data.frame(comparison =
                     paste(names(Cov_d_beta), reference.group, sep = sep.trts),
                   cov = Cov_d_beta) %>%
        # Include reverse comparisons by reversing the label; values are
        # equivalent.
        #
        bind_rows(data.frame(comparison =
                               paste(reference.group, names(Cov_d_beta),
                                     sep = sep.trts),
                             cov = Cov_d_beta))
      #
      dat_cov_d_beta <-
        merge(all_se.beta, dat_cov_d_beta, by = "comparison", all.x = TRUE) %>%
        mutate(cov = if_else(is.na(cov), cov.default, cov)) %>%
        select(-cov.default)
    }
    else if (assumption == "independent") {
      mat_beta <-
        outer(c(0, results$coef[!sel.d]), c(0, results$coef[!sel.d]), "-")
      rownames(mat_beta) <- colnames(mat_beta) <-
        c(reference.group, results$treat[!sel.d])
      #
      Cov_beta <- Cov[!sel.d, !sel.d, drop = FALSE]
      Cov_beta <- cbind(0, Cov_beta)
      Cov_beta <- rbind(0, Cov_beta)
      rownames(Cov_beta)[1] <- colnames(Cov_beta)[1] <- reference.group
      #
      mat_se.beta <-
        sqrt(outer(diag(Cov_beta), diag(Cov_beta), "+") - 2 * Cov_beta)
      rownames(mat_se.beta) <- colnames(mat_se.beta) <-
        c(reference.group, results$treat[!sel.d])
      #
      # Covariances
      #
      # Make the interaction vs treatment matrix a proper square matrix,
      # interactions in rows, treats in cols
      #
      Cov_d_beta <- Cov[!sel.d, sel.d, drop = FALSE]
      #
      # Drop interaction part from row names
      #
      rownames(Cov_d_beta) <- gsub(":.*|`", "", rownames(Cov_d_beta))
      #
      # Only keep treatments with interaction term
      #
      Cov_d_beta <- Cov_d_beta[, rownames(Cov_d_beta), drop = FALSE]
      #
      dat_cov_d_beta <- as.data.frame(Cov_d_beta)
      dat_cov_d_beta <- reshape(dat_cov_d_beta, 
                                idvar = "interaction", 
                                ids = rownames(dat_cov_d_beta), 
                                direction = "long", 
                                varying = list(names(dat_cov_d_beta)), 
                                v.names = "value", 
                                timevar = "treat", 
                                times = names(dat_cov_d_beta)) %>%
        mutate(type = if_else(treat == interaction, "ref", "interim"))
      #
      interim <- dat_cov_d_beta %>% select(-type) %>%
        rename(interim = value)
      
      ref <- dat_cov_d_beta %>% filter(type == "ref") %>%
        select(-treat, -type) %>%
        rename(ref = value)
      #
      comb1 <- merge(ref, interim, by = "interaction")
      comb1 <- merge(comb1, interim,
                     by.x = c("treat", "interaction"),
                     by.y = c("interaction", "treat"))
      comb1 <- merge(comb1, ref,
                     by.x = "treat", by.y = "interaction")
      #
      comb1 %<>%
        mutate(cov =
                 if_else(interaction == treat,
                         ref.x,
                         ref.x + ref.y - interim.x - interim.y))
      comb1 %<>%
        mutate(comparison =
                 if_else(interaction == treat,
                         paste(interaction, reference.group, sep = sep.trts),
                         paste(interaction, treat, sep = sep.trts))) %>%
        bind_rows(
          comb1 %<>%
            mutate(comparison =
                     if_else(interaction == treat,
                             paste(reference.group, interaction,
                                   sep = sep.trts),
                             paste(treat, interaction, sep = sep.trts)))) %>%
        select(comparison, cov)
      #
      dat_cov_d_beta <- unique(comb1)
    }
    #
    # Transform from matrix to data frame in long format
    #
    dat_beta <- dat_se.beta <- NULL
    #
    for (i in colnames(mat_beta)) {
      dat_beta.i <- as.data.frame(mat_beta)[i]
      names(dat_beta.i) <- "beta"
      dat_beta.i$comparison <- paste(rownames(dat_beta.i), i, sep = sep.trts)
      #
      dat_beta <- rbind(dat_beta, dat_beta.i)
      #
      dat_se.beta.i <- as.data.frame(mat_se.beta)[i]
      names(dat_se.beta.i) <- "se.beta"
      dat_se.beta.i$comparison <-
        paste(rownames(dat_se.beta.i), i, sep = sep.trts)
      dat_se.beta <- rbind(dat_se.beta, dat_se.beta.i)
    }
  }
  #
  # Combine elements for summary and all combinations
  #
  res <- merge(dat_d, dat_se.d, by = "comparison", all.x = TRUE)
  res <- merge(res, dat_beta, by = "comparison", all.x = TRUE)
  res <- merge(res, dat_se.beta, by = "comparison", all.x = TRUE)
  res <- merge(res, dat_cov_d_beta, by = "comparison", all.x = TRUE)
  #
  # Drop comparisons with themselves
  #
  split_comps <- strsplit(res$comparison, "[:]")
  res$treat1 <- sapply(split_comps, first)
  res$treat2 <- sapply(split_comps, second)
  res %<>% filter(treat1 != treat2) %>% filter(!duplicated(comparison))
  #
  res
}
