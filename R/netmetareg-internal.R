# utilities for extracting NMR results #########
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
  }
  else {
    # Add original treatment names to merge results
    # (for the common assumption only the covariates observed for the
    #  reference treatment are relevant)
    #
    dat_beta %<>%
      mutate(treat_temp = ifelse(treat == "nonref", reference.group, treat),
             treat_long =
               as.character(
                 factor(treat_temp, levels = trts.abbr, labels = trts))) %>%
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
  # Extract NMR attributes / assumptions
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
  }  else if (consistency) {
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
  } else if(consistency==FALSE){
    if(assumption=="common"){# common
      mat_beta<-matrix(nrow = length(x$results[x$results$type =="d","treat"])+1,
                       ncol=length(x$results[x$results$type =="d","treat"])+1,
                       data=x$results[x$results$type =="beta","coef"])
      diag(mat_beta)<-0
      dimnames(mat_beta)<-list(c(reference.group,x$results[x$results$type =="d","treat"]),
                               c(reference.group,x$results[x$results$type =="d","treat"]))
      
      mat_se.beta<-matrix(nrow = length(x$results[x$results$type =="d","treat"])+1,
                          ncol=length(x$results[x$results$type =="d","treat"])+1,
                          data=x$results[x$results$type =="beta","se"])
      diag(mat_se.beta)<-0
      dimnames(mat_se.beta)<-list(c(reference.group,x$results[x$results$type =="d","treat"]),
                                  c(reference.group,x$results[x$results$type =="d","treat"]))
      
      # covariances ###############
      temp_vcov<-Cov_d_beta <- Cov[!sel.d, sel.d]
      # rownames(temp_vcov)<-gsub(":.*|`","",rownames(temp_vcov))
      temp_vcov[reference.group]<-0
      mat_vcov<-as.numeric(temp_vcov)
      names(mat_vcov)<-names(temp_vcov)
      mat_vcov<-outer(mat_vcov,mat_vcov,"-")
      mat_vcov<-mat_vcov[order(rownames(mat_vcov)),order(colnames(mat_vcov))]
      
      mat_vcov<-as.data.frame(mat_vcov)
      
      
    }else if(assumption=="independent"){
      mat_beta<-x$results[x$results$type =="beta",c("coef","treat")]
      mat_se.beta<-x$results[x$results$type =="beta",c("se","treat")]
      # #reverse direction in case reverse direction occurs in original data
      mat_beta<-rbind.data.frame(mat_beta,
                                 data.frame("coef"=mat_beta$coef,# reversal depends on interpretation of ivals. Which is why users need the ival warning. Treatment direction is reversible
                                            "treat"=paste(gsub(".*_vs_","",mat_beta$treat), # without the vs this can get messy
                                                          gsub("_vs_.*","",mat_beta$treat),
                                                          sep="_vs_")))
      mat_beta$treat<-gsub(".comp_", "", mat_beta$treat) # since i added comp
      mat_beta$treat<-gsub("_vs_", ":", mat_beta$treat) #standardize with treatment effect matrices
      mat_se.beta<-rbind.data.frame(mat_se.beta)
      mat_se.beta<-rbind.data.frame(mat_se.beta,
                                    data.frame("se"=mat_se.beta$se,
                                               "treat"=paste(gsub(".*_vs_","",mat_se.beta$treat), # without the vs this can get messy
                                                             gsub("_vs_.*","",mat_se.beta$treat),
                                                             sep="_vs_")))
      mat_se.beta$treat<-gsub(".comp_", "", mat_se.beta$treat) # since i added comp
      mat_se.beta$treat<-gsub("_vs_", ":", mat_se.beta$treat) #standardize with treatment effect matrices
      names(mat_beta)<-c("beta","comparison")
      names(mat_se.beta)<-c("se.beta","comparison")
      
      # it is already a joinable list
      dat_beta<-mat_beta
      dat_se.beta<-mat_se.beta
      
      # covariances ###############
      temp_vcov<-Cov_beta <- Cov[!sel.d, sel.d]
      rownames(temp_vcov)<-gsub("\\.comp_|:.*|`","",rownames(temp_vcov))
      temp_vcov_interactions<-rownames(temp_vcov)
      temp_vcov<-as.data.frame(temp_vcov)
      temp_vcov[colnames(temp_vcov)[!colnames(temp_vcov)%in%temp_vcov_interactions],]<-0 # the dataframe appends rows whereas the matrix will not
      
      temp_vcov <- cbind(0, temp_vcov)
      temp_vcov <- rbind(0, temp_vcov)
      rownames(temp_vcov)[1] <- reference.group
      colnames(temp_vcov)[1] <- reference.group
      
      temp_vcov<-temp_vcov[order(rownames(temp_vcov)),order(colnames(temp_vcov))]
      mat_vcov<-temp_vcov
      
      mat_vcov<-as.data.frame(mat_vcov)
      mat_vcov<-reshape(mat_vcov,
                        idvar = "interaction",
                        ids=rownames(mat_vcov),
                        direction = "long", 
                        varying = list(names(mat_vcov)),
                        v.names = "value", 
                        timevar = "treat", 
                        times = names(mat_vcov))
      mat_vcov$int_treat1<-gsub(".*_vs_","",mat_vcov$interaction)
      mat_vcov$int_treat2<-gsub("_vs_.*","",mat_vcov$interaction)
      
      mat_vcov$type<-ifelse(mat_vcov$treat==mat_vcov$int_treat1&mat_vcov$treat==mat_vcov$int_treat2,"ref","interim")
      
      interim<-mat_vcov[mat_vcov$type=="interim"&
                          (mat_vcov$treat==mat_vcov$int_treat1|
                             mat_vcov$treat==mat_vcov$int_treat2),
                        !names(mat_vcov)%in%c("type")]
      names(interim)<-gsub("value","interim",names(interim))
      
      ref<-mat_vcov[mat_vcov$type=="ref",!names(mat_vcov)%in%c("treat","type","int_treat1", "int_treat2")]
      names(ref)<-gsub("value","ref",names(ref))
      
      comb1<-merge(interim,ref, by.x="int_treat2",by.y="interaction")
      comb1<-merge(comb1,ref, by.x="int_treat1",by.y="interaction")
      comb1<-merge(comb1,interim, by.x=c("interaction","int_treat1","int_treat2"),by.y=c("interaction","int_treat1","int_treat2"))
      comb1<-comb1[comb1$treat.x==comb1$int_treat2,]# only in the same order as the covariate
      
      comb1$cov<-comb1$ref.x+comb1$ref.y+comb1$interim.x-comb1$interim.y
      comb1$.comp<-paste(comb1$treat.x,comb1$treat.y,sep=":")
      combrev<-comb1
      combrev$.comp<-paste(combrev$treat.y,combrev$treat.x,sep=":")#reverse, does not change covariance
      comb1<-rbind(comb1,combrev)
      dat_cov_d_beta<-unique(comb1[,c(".comp","cov")])
      names(dat_cov_d_beta)<-c("comparison", "cov")
    }
    
    
    if (assumption == "common" | consistency) {
      dat_beta <- dat_se.beta<- dat_cov_d_beta<-NULL
    }
    #
    for (i in colnames(mat_beta)) {
      if (assumption == "common" | consistency) {
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
        #
        dat_cov_d_beta.i <- as.data.frame(mat_vcov)[i]
        names(dat_cov_d_beta.i) <- "cov"
        dat_cov_d_beta.i$comparison <-
          paste(rownames(dat_cov_d_beta.i), i, sep = sep.trts)
        dat_cov_d_beta <- rbind(dat_cov_d_beta, dat_cov_d_beta.i)
      }
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

fix_interaction_names <- function(names_vec, covar_name) {
  sapply(names_vec, function(nm) {
    parts <- strsplit(nm, ":", fixed = TRUE)[[1]] 
    if (length(parts) != 2) return(nm)
    if (parts[1] == covar_name) {
      paste(parts[2],covar_name, sep = ":") # If the name already has two parts, put the covar second
    } else {
      nm
    }
  })
}

# Compute covariance of any linear combination a' Sigma a
cov_combination <- function(a_coef, b_coef, temp_Sigma) {
  a <- b <- setNames(rep(0, nrow(temp_Sigma)), rownames(temp_Sigma))
  if (any(!names(a_coef) %in% names(a))) stop("invalid names")
  if (any(!names(b_coef) %in% names(b))) stop("invalid names")
  a[names(a_coef)] <- a_coef
  b[names(b_coef)] <- b_coef
  a <- as.matrix(a)
  b <- as.matrix(b)
  as.numeric(t(a) %*% temp_Sigma %*% b)
}

var_combination <- function(coef, temp_Sigma) {
  cov_combination(a_coef=coef, b_coef=coef, temp_Sigma)
}

compare_names <- function(x, reorder = FALSE) {
  if (!is.list(x) || is.data.frame(x)) {
    stop("Input must be a list")
  }
  if (length(x) < 2) return(x)
  
  ref_names <- names(x[[1]])
  if (is.null(ref_names) || length(ref_names) == 0) {
    stop("First vector has no names")
  }
  ref_names <- ref_names[!is.na(ref_names) & ref_names != ""]
  ref_names_sorted <- sort(ref_names)
  ref_len <- length(x[[1]])
  
  all_match <- TRUE
  mismatch_info <- character()
  needs_reorder <- FALSE
  
  for (i in seq_along(x)) {
    curr_names <- names(x[[i]])
    if (is.null(curr_names) || length(curr_names) == 0) {
      mismatch_info <- c(mismatch_info, paste0("Vector ", i, " has no names"))
      all_match <- FALSE
      next
    }
    curr_names <- curr_names[!is.na(curr_names) & curr_names != ""]
    curr_names_sorted <- sort(curr_names)
    
    if (!identical(ref_names_sorted, curr_names_sorted)) {
      missing <- setdiff(ref_names_sorted, curr_names_sorted)
      extra <- setdiff(curr_names_sorted, ref_names_sorted)
      if (length(missing) > 0) {
        mismatch_info <- c(mismatch_info, 
                           paste0("Vector ", i, " (", names(x)[i], ") missing: ", paste(missing, collapse = ", ")))
      }
      if (length(extra) > 0) {
        mismatch_info <- c(mismatch_info, 
                           paste0("Vector ", i, " (", names(x)[i], ") extra: ", paste(extra, collapse = ", ")))
      }
      all_match <- FALSE
    }
    
    curr_len <- length(x[[i]])
    if (curr_len != ref_len) {
      mismatch_info <- c(mismatch_info, 
                         paste0("Vector ", i, " (", names(x)[i], ") length: ", curr_len, " (expected: ", ref_len, ")"))
      all_match <- FALSE
    }
    
    if (!identical(ref_names, curr_names) && identical(ref_names_sorted, curr_names_sorted)) {
      needs_reorder <- TRUE
    }
  }
  
  if (!all_match) {
    stop("Name mismatch detected:\n", paste(mismatch_info, collapse = "\n"))
  }
  
  if (reorder && needs_reorder) {
    for (i in seq_along(x)) {
      curr_names <- names(x[[i]])
      if (!identical(ref_names, curr_names)) {
        order_idx <- match(ref_names, curr_names)
        x[[i]] <- x[[i]][order_idx]
      }
    }
    return(invisible(x))
  }
  
  return(x)
}

clean_var_names <- function(nms, covar.name = NULL) {
  if (is.null(covar.name)) {
    covar.name <- mynmri$.netmeta$covar.name
  }
  pattern <- paste0("`|:", covar.name, "|", covar.name, ":")
  gsub(pattern, "", nms, fixed = TRUE)
}

get_blup <- function(netmetareg_obj) {
  all_trts<-netmetareg_obj$.netmeta$trts[netmetareg_obj$.netmeta$trts!=netmetareg_obj$.netmeta$reference.group]
  all_trts<-make.names(all_trts)
  reference.group<-netmetareg_obj$.netmeta$reference.group
  reference.group<-make.names(reference.group)
  
  int_fixed <- netmetareg_obj$beta[grepl(":",rownames(netmetareg_obj$beta))==TRUE,]#fixed effect of interaction
  names(int_fixed)<-gsub(paste0("`|:",netmetareg_obj$.netmeta$covar.name),"",fix_interaction_names(names(int_fixed), covar_name = netmetareg_obj$.netmeta$covar.name))
  
  # todo check if X matrix has interactions for every treatment if the assumption is independent
  Z_slope<-netmetareg_obj$Z_slope
  covar.name<-netmetareg_obj$.netmeta$covar.name
  
  W <- chol2inv(chol(netmetareg_obj$M))
  stXWX <- chol2inv(chol(as.matrix(t(netmetareg_obj$X) %*% W %*% netmetareg_obj$X)))
  Hmat <- netmetareg_obj$X %*% stXWX %*% crossprod(netmetareg_obj$X, W)
  I <- diag(netmetareg_obj$k)
  
  crit <- qnorm(1-netmetareg_obj$level/2) 
  
  # if mynmre$tau2B=0 then it is the same as a CCIE model
  if(is.null(Z_slope)==FALSE & netmetareg_obj$tau2B!=0){
    withZ_slope<-TRUE
    colnames(Z_slope)<-fix_interaction_names(colnames(Z_slope), covar_name = netmetareg_obj$.netmeta$covar.name)
    D<-(netmetareg_obj$tau2B * diag(ncol(Z_slope)))
    rownames(D)<-colnames(D)<-colnames(Z_slope)
    
  } else if(netmetareg_obj$.netmeta$consistency==TRUE){
    withZ_slope<-FALSE
    Z_slope<-matrix(0,nrow = nrow(netmetareg_obj$X), ncol=length(netmetareg_obj$.netmeta$trts)-1)
    D<-(0 * diag(length(netmetareg_obj$.netmeta$trts)-1))
    rownames(D)<-colnames(D)<-colnames(Z_slope)<-netmetareg_obj$.netmeta$trts[netmetareg_obj$.netmeta$trts!=netmetareg_obj$.netmeta$reference.group]
  }else if(netmetareg_obj$.netmeta$consistency==FALSE){
    withZ_slope<-FALSE
    Z_slope<-matrix(0,nrow = nrow(netmetareg_obj$X), ncol=length(netmetareg_obj$.netmeta$x$comparisons))
    D<-(0 * diag(length(netmetareg_obj$.netmeta$x$comparisons)))
    rownames(D)<-colnames(D)<-colnames(Z_slope)<-paste0(".comp_",gsub(":","_vs_",netmetareg_obj$.netmeta$x$comparisons),":Ival")
  }
  
  DZtW<-D%*%t(Z_slope)%*%  W
  BLUP<-(DZtW %*% (netmetareg_obj$yi-predict(netmetareg_obj)$pred))[,1]
  var_BLUP <- D - (DZtW %*% (I - Hmat) %*% Z_slope %*% D)
  se_BLUP <- sqrt(diag(var_BLUP))
  pi.lb <- c(BLUP - crit * se_BLUP)
  pi.ub <- c(BLUP + crit * se_BLUP)
  
  #BLUP+BLUE of interaction
  if(netmetareg_obj$.netmeta$assumption=="independent"){
    names(int_fixed)<-make.names(names(int_fixed))
    names(BLUP)<-make.names(names(BLUP))
    temp_BLUPE_int <- merge(as.data.frame(BLUP), as.data.frame(int_fixed),
                            by = "row.names",        # use names as key
                            all = TRUE,              # keep all names
                            suffixes = c("BLUP", "int_fixed"))
  } else {
    # otherwise the interaction applies to all BLUPs
    BLUP_df<-as.data.frame(BLUP)
    BLUP_df$Row.names<-rownames(BLUP_df)
    temp_BLUPE_int <- merge(BLUP_df, as.data.frame(int_fixed),
                            all = TRUE,              # keep all names
                            suffixes = c("BLUP", "int_fixed"))
  }
  temp_BLUPE_int$BLUPE_int <- rowSums(temp_BLUPE_int[ , c("BLUP", "int_fixed")])
  BLUPE_int<-temp_BLUPE_int$BLUPE_int
  names(BLUPE_int)<-temp_BLUPE_int$Row.names
  
  
  # trying to get vcov of combined ranef and pred # exactshould be in Henderson (1975).
  if(withZ_slope){
    C_xx<-crossprod(netmetareg_obj$X, W %*% netmetareg_obj$X)          # X'WX
    C_xu <- crossprod(netmetareg_obj$X, W %*% Z_slope)          # X'WZ
    C_ux <- crossprod(Z_slope, W %*% netmetareg_obj$X)          # Z'WX
    C_uu <- crossprod(Z_slope, W %*% Z_slope) + solve(D)  # Z'WZ + G^{-1}
    C <- rbind(cbind(C_xx, C_xu),
               cbind(C_ux, C_uu))
    V_joint <- solve(C) # joint variance–covariance matrix of \((\hat\beta,\hat u)\) is the inverse of this coefficient matrix (up to scale):
    
  } else {
    V_joint <- netmetareg_obj$vb
  }
  # add missing interactions if any failed # Consistency!
  if(netmetareg_obj$.netmeta$assumption=="independent" & netmetareg_obj$.netmeta$consistency==TRUE){
    design.cols<-paste0("`",all_trts[all_trts!=reference.group],":",netmetareg_obj$.netmeta$covar.name,"`")
    missing.design <- setdiff(design.cols, colnames(V_joint))
    V_joint <- cbind(V_joint, matrix(0, nrow = nrow(V_joint), ncol = length(missing.design), dimnames = list(NULL, missing.design)))
    V_joint <- rbind(V_joint, matrix(0, ncol = ncol(V_joint), nrow = length(missing.design), dimnames = list(missing.design,NULL)))
    V_joint<-V_joint[rev(order(colnames(V_joint))),rev(order(colnames(V_joint)))]
    int_fixed[gsub(":.*","",missing.design)]<-0
    n_FE<-(length(netmetareg_obj$.netmeta$trts)-1)*2
  } else if(netmetareg_obj$.netmeta$assumption=="independent" & netmetareg_obj$.netmeta$consistency==FALSE){
    design.cols<-paste0("`.comp_",gsub(":","_vs_",netmetareg_obj$.netmeta$x$comparisons),":Ival:",covar.name,"`")
    missing.design <- setdiff(design.cols, colnames(V_joint))
    V_joint <- cbind(V_joint, matrix(0, nrow = nrow(V_joint), ncol = length(missing.design), dimnames = list(NULL, missing.design)))
    V_joint <- rbind(V_joint, matrix(0, ncol = ncol(V_joint), nrow = length(missing.design), dimnames = list(missing.design,NULL)))
    V_joint<-V_joint[rev(order(colnames(V_joint))),rev(order(colnames(V_joint)))]
    int_fixed[gsub(":.*","",missing.design)]<-NA
    n_FE<-(length(netmetareg_obj$.netmeta$trts)-1)*2
  } else {
    n_FE<-length(netmetareg_obj$beta)
  }
  
  # Extract blocks
  ntrt<-length(netmetareg_obj$.netmeta$trts)
  q<-ncol(V_joint) # total variance terms
  
  se_TE_fixed <- sqrt(diag(V_joint[1:(ntrt-1), 1:(ntrt-1), drop=FALSE])) #SE for fixed effect of TE
  se_int_fixed <- sqrt(diag(as.matrix(V_joint[ntrt:n_FE, ntrt:n_FE]))) #SE for fixed effect of interaction
  if(ntrt!=n_FE) {
    cov_TE_int_fixed <-diag(V_joint[1:(ntrt-1), ntrt:n_FE])
    names(cov_TE_int_fixed) <-rownames(V_joint[1:(ntrt-1), ntrt:n_FE])
  } else {
    cov_TE_int_fixed <- V_joint[1:(ntrt-1), ntrt:n_FE]} #Covariance of all fixed effects. TE and interaction fixed effect
  
  if(withZ_slope){
    se_int_ranef <- sqrt(diag(V_joint[(n_FE+1):q, (n_FE+1):q])) #SEs for ranef of interaction slopes
    cov_int_ranef_fixed <- V_joint[n_FE,(n_FE+1):q] #Covariance of fixed interaction effects and treatment specific random interaction effects
    cov_TE_ranef <- diag(V_joint[(n_FE+1):q,1:(ntrt-1), drop=FALSE]) #Covariance of random interaction effects and ALL fixed TE effects
    names(cov_TE_ranef)<-colnames(V_joint[(n_FE+1):q,1:(ntrt-1), drop=FALSE])
    names_nu <- setNames(rep(0, nrow(V_joint)), rownames(V_joint))
    
    if(netmetareg_obj$.netmeta$consistency){
      var_nu <- sapply(rownames(netmetareg_obj$beta)[1:ntrt-1], function(x) {
        names_nu[rownames(netmetareg_obj$beta)[ntrt:n_FE]] <- 1 # FE of interaction
        names_nu[paste0(x,":",covar.name)] <- 1 # REof interaction
        var_combination(names_nu, temp_Sigma=V_joint)
      })
    } else{
      
      var_nu <- sapply(colnames(Z_slope), function(x) {
        names_nu[rownames(netmetareg_obj$beta)[ntrt:n_FE]] <- 1 # FE of interaction
        names_nu[x] <- 1 # REof interaction
        var_combination(names_nu, temp_Sigma=V_joint)
      })
    }
  } else {
    # there are no random effects so there are no covariances
    names_nu <- setNames(rep(0, nrow(V_joint)), rownames(V_joint))
    # var_nu will have either length= 1 or length = number of interactions
    if(netmetareg_obj$.netmeta$consistency){
      se_int_ranef<- rep(0, length(netmetareg_obj$.netmeta$trts)-1)
      cov_int_ranef_fixed<- rep(0, length(netmetareg_obj$.netmeta$trts)-1)
      cov_TE_ranef<- rep(0, length(netmetareg_obj$.netmeta$trts)-1)
      names(se_int_ranef)<-names(cov_int_ranef_fixed)<-names(cov_TE_ranef)<-netmetareg_obj$.netmeta$trts[netmetareg_obj$.netmeta$trts!=netmetareg_obj$.netmeta$reference.group]
      
      var_nu <- sapply(rownames(V_joint)[ntrt:nrow(V_joint)], function(x) {
        # there are no ranefs
        names_nu[x] <- 1 # FE of interaction
        var_combination(names_nu, temp_Sigma=V_joint)
      })
      
      if (length(var_nu)==1){
        var_nu<-rep(var_nu,ntrt-1)
        names(var_nu)<-netmetareg_obj$.netmeta$trts[netmetareg_obj$.netmeta$trts!=netmetareg_obj$.netmeta$reference.group]
      }
      
    } else{
      se_int_ranef<- rep(0, length(netmetareg_obj$.netmeta$x$comparisons))
      cov_int_ranef_fixed<- rep(0, length(netmetareg_obj$.netmeta$x$comparisons))
      cov_TE_ranef<- rep(0, length(netmetareg_obj$.netmeta$x$comparisons))
      names(se_int_ranef)<-names(cov_int_ranef_fixed)<-names(cov_TE_ranef)<-paste0("`.comp_",gsub(":","_vs_",netmetareg_obj$.netmeta$x$comparisons),":Ival:",covar.name,"`")
      
      
      var_nu <- sapply(rownames(V_joint)[ntrt:nrow(V_joint)], function(x) {
        # there are no ranefs
        names_nu[x] <- 1 # FE of interaction
        var_combination(names_nu, temp_Sigma=V_joint)
      })
    }
    
    if (length(var_nu)==1){
      var_nu<-rep(var_nu,length(netmetareg_obj$.netmeta$x$comparisons))
      names(var_nu)<-paste0("`.comp_",gsub(":","_vs_",netmetareg_obj$.netmeta$x$comparisons),":Ival:",covar.name,"`")
    }
  }
  
  
  if (length(int_fixed)==1){
    int_fixed<-rep(int_fixed,ntrt-1)
    se_int_fixed<-rep(se_int_fixed,ntrt-1)
    names(int_fixed)<-names(se_int_fixed)<-netmetareg_obj$.netmeta$trts[netmetareg_obj$.netmeta$trts!=netmetareg_obj$.netmeta$reference.group]
  }
  
  
  
  res_list <- list(
    nu = BLUPE_int,
    var_nu = var_nu,
    BLUP = BLUP,
    se_BLUP = se_BLUP,#SEs for ranef of interaction slopes
    pi.lb = pi.lb,
    pi.ub = pi.ub,
    se_int_ranef=se_int_ranef,
    cov_int_ranef_fixed = cov_int_ranef_fixed
  )
  
  res_list <- lapply(res_list, function(vec) {
    names(vec) <- gsub(paste0("`|:",netmetareg_obj$.netmeta$covar.name,"|",netmetareg_obj$.netmeta$covar.name,":"),"",names(vec))
    vec
  })
  
  res_list<-lapply(res_list,function(x) {
    names(x)<-make.names(names(x))
    return(x)})
  res_list_order<-compare_names(res_list, reorder = TRUE)
  
  pred<-as.data.frame(do.call(cbind, res_list_order))
  pred$se_nu <- sqrt(pred$var_nu)
  pred$treat1<-names(res_list_order$nu)
  if(netmetareg_obj$.netmeta$consistency==TRUE){
    pred$treat2<-netmetareg_obj$.netmeta$reference.group
  } else {
    pred$treat2<-gsub("_vs_.*|:Ival|.comp_","",pred$treat1)
    pred$treat1<-gsub(".*_vs_|:Ival|.comp_","",pred$treat1)
  }
  # align with treat transformation in netmetareg
  pred$treat1 <- make.names(pred$treat1)
  pred$treat2 <- make.names(pred$treat2)
  
  return(list(pred = pred, V_joint = V_joint))
}

nmr_full_results_new <- function(netmetareg_obj) {
  
  reference.group<-make.names(netmetareg_obj$.netmeta$reference.group)
  covar.name<-make.names(netmetareg_obj$.netmeta$covar.name)
  
  Z_slope<-netmetareg_obj$Z_slope
  withZ_slope<-ifelse(is.null(Z_slope),FALSE,TRUE)
  
  # obtain per‑treatment totals and their SEs
  nu_res <- get_blup(netmetareg_obj)
  nu_df <- as.data.frame(nu_res$pred)
  rownames(nu_df)<-nu_df$treat1
  mu_df <- cbind.data.frame(netmetareg_obj$beta, netmetareg_obj$se)
  mu_df <- mu_df[grepl(":",rownames(mu_df))==FALSE,]
  rownames(mu_df) <-gsub(paste0("`|:",netmetareg_obj$.netmeta$covar.name,"|",netmetareg_obj$.netmeta$covar.name,":"),"",rownames(mu_df))
  V_joint <- nu_res$V_joint
  # add a null row for reference
  nu_df[reference.group,]<-NA
  nu_df[reference.group, sapply(nu_df, is.numeric)] <- lapply(nu_df[reference.group, sapply(nu_df, is.numeric)], function(x) ifelse(is.na(x), 0,x))
  nu_df[reference.group,c("treat1","treat2")]<-reference.group
  mu_df[reference.group,]<-c(0,0)
  names(mu_df)<-c("mu", "se_mu")
  
  if(withZ_slope){
    V_joint<-rbind(V_joint,rep(0,ncol(V_joint)),rep(0,ncol(V_joint)))
    V_joint<-cbind(V_joint,rep(0,nrow(V_joint)),rep(0,nrow(V_joint)))
    rownames(V_joint)<-c(rownames(nu_res$V_joint),reference.group, paste0( netmetareg_obj$.netmeta$covar.name,":",reference.group))
    colnames(V_joint)<-c(colnames(nu_res$V_joint),reference.group, paste0( netmetareg_obj$.netmeta$covar.name,":",reference.group))
  } else {
    V_joint<-rbind(V_joint,rep(0,ncol(V_joint)))
    V_joint<-cbind(V_joint,rep(0,nrow(V_joint)))
    rownames(V_joint)<-c(rownames(nu_res$V_joint),reference.group)
    colnames(V_joint)<-c(colnames(nu_res$V_joint),reference.group)
  }
  colnames(V_joint) <- gsub("`","", colnames(V_joint))
  rownames(V_joint) <- gsub("`","", rownames(V_joint))
  
  # ordered pairs (treat1, treat2)
  res_beta <- data.frame()
  if(netmetareg_obj$.netmeta$consistency==TRUE){
    for (i in seq_len(nrow(nu_df))) {
      for (j in seq_len(nrow(nu_df))) {
        nu_diff <- nu_df$nu[i] - nu_df$nu[j]
        
        names_nu <- setNames(rep(0, nrow(V_joint)), rownames(V_joint))
        # For potential SEM
        if(nu_df$treat1[i]==reference.group){
          names_nu[paste0("nonref:", netmetareg_obj$.netmeta$covar.name)] <- -1 
        } else if(nu_df$treat1[j]==reference.group){
          names_nu[paste0("nonref:", netmetareg_obj$.netmeta$covar.name)] <- 1 
        }else {
          names_nu[paste0("nonref:", netmetareg_obj$.netmeta$covar.name)] <- 0 # because they have the same interaction fixed effect
        }
        # for ICIE or ECIE
        if(i!=j){
          names_nu[paste0(rownames(nu_df)[i],":",covar.name)] <- 1
          names_nu[paste0(rownames(nu_df)[j],":",covar.name)] <- -1
        }
        # only keep those that match the model
        names_nu<-names_nu[names(names_nu)[names(names_nu) %in%  rownames(V_joint)]]
        var_nu_diff<-var_combination(names_nu, temp_Sigma = V_joint)
        se_nu_diff <- sqrt(var_nu_diff)
        res_beta <- rbind(res_beta, data.frame(
          treat1 = nu_df$treat1[i],
          treat2 = nu_df$treat1[j],
          beta = nu_diff,
          se.beta = se_nu_diff,
          var.beta = var_nu_diff
          
        ))
      }
    }
  } else {
    res_beta <- nu_df%>%
      select(treat1,treat2, beta=nu, se.beta=se_nu, var.beta=var_nu)
  }
  res_d <- data.frame()
  for (i in seq_len(nrow(mu_df))) {
    for (j in seq_len(nrow(mu_df))) {
      mu_diff <- mu_df$mu[i] - mu_df$mu[j]
      names_mu <- setNames(rep(0, nrow(V_joint)), rownames(V_joint))
      names_a <- setNames(rep(0, nrow(V_joint)), rownames(V_joint))
      names_b <- setNames(rep(0, nrow(V_joint)), rownames(V_joint))
      if(i!=j){
        # treatments in a (d)
        names_a[rownames(mu_df)[i]] <- names_mu[rownames(mu_df)[i]] <- 1
        names_a[rownames(mu_df)[j]] <- names_mu[rownames(mu_df)[j]] <- -1
        
        # interactions in b (nu)
        # For potential SEM
        if(nu_df$treat1[i]==reference.group){
          names_b[paste0("nonref:", netmetareg_obj$.netmeta$covar.name)] <- -1 
        } else if(nu_df$treat1[j]==reference.group){
          names_b[paste0("nonref:", netmetareg_obj$.netmeta$covar.name)] <- 1 
        } else {
          names_b[paste0("nonref:", netmetareg_obj$.netmeta$covar.name)] <- 0 # because they have the same interaction fixed effect
        }
        # for ICIE or ECIE
        names_b[paste0(rownames(mu_df)[i],":xcov")] <- 1
        names_b[paste0(rownames(mu_df)[j],":xcov")] <- -1
      }
      # only keep those that match the model
      names_b<-names_b[names(names_b)[names(names_b) %in%  rownames(V_joint)]]
      
      var_mu_diff<-var_combination(names_mu, temp_Sigma = V_joint)
      se_mu_diff <- sqrt(var_mu_diff)
      
      cov_theta_diff<-cov_combination(a_coef = names_a,b_coef = names_b, temp_Sigma = V_joint)
      
      res_d <- rbind(res_d, data.frame(
        treat1 = rownames(mu_df)[i],
        treat2 = rownames(mu_df)[j],
        d = mu_diff,
        se.d= se_mu_diff,
        var.d = var_mu_diff,
        cov = cov_theta_diff 
      ))
    }
  }
  res_d%>%
    mutate(comparison=paste0(treat1,":",treat2))%>%
    merge(res_beta, by=c("treat1", "treat2"), all.x = TRUE)
}


# utilities for implementing NMR ######################
mkIval <- function(x) {
  # Create directionality parameter as treatment order
  #
  ref <- sort(unique(c(x$treat1, x$treat2)))[1]
  #
  x$Ival1 <- ifelse(x$treat1 == ref, 0, 1)
  x$Ival2 <- ifelse(x$treat2 == ref, 0, 1)
  #
  return(x)
}

any_invalid_I <- function(x) {
  # Get rid of warning "no visible binding for global variable"
  #
  treat1 <- treat2 <- treat1_new <- treat2_new <- NULL
  #
  # Check Ival with ability to specify two custom variables
  #
  edges2 <- x %>%
    mutate(w = case_when(Ival_diff > 0 ~ -Ival_diff, TRUE ~ Ival_diff),
           treat1_new =
             case_when(Ival_diff %in% c(0, 1) ~ treat2, TRUE ~ treat1),
           treat2_new =
             case_when(Ival_diff %in% c(0, 1) ~ treat1, TRUE ~ treat2)
           ) %>%
    select(-treat1, -treat2) %>%
    rename(treat1 = treat1_new, treat2 = treat2_new)
  #
  treats <- unique(c(edges2$treat1, edges2$treat2))
  n <- length(treats)
  #
  dist <- rep(0, n)
  names(dist) <- treats
  #
  # For each treatment
  #
  for (i in seq_len(n)) {
    updated <- FALSE
    #
    # For each comparison
    #
    for (e in seq_len(nrow(edges2))) {
      # Distance for comparator treatment larger than the distance for
      # reference + the length
      if (dist[edges2$treat1[e]] > dist[edges2$treat2[e]] + edges2$w[e]) {
        dist[edges2$treat1[e]] <- dist[edges2$treat2[e]] + edges2$w[e]
        updated <- TRUE
      }
    }
    if (!updated)
      return(FALSE)
  }
  #
  return(TRUE)
}

check_Ival <- function(dat) {
  # Get rid of warning "no visible binding for global variable"
  #
  studlab <- treat <- treat1 <- treat2 <- Ival1 <- Ival2 <- NULL
  #
  # Check for invalid values
  #
  if (sum(!dat$Ival1 %in% c(-1, 0, 1)) > 0)
    stop("Error in directionality assignment. ",
         "Ival1 only accepts values of -1, 0, or 1.")
  #
  if (sum(!dat$Ival2 %in%c(-1, 0, 1)) > 0)
    stop("Error in directionality assignment. ",
         "Ival2 only accepts values of -1, 0, or 1.")
  #
  multiple_Ival <- dat %>%
    select(studlab, treat = treat1, Ival = Ival1) %>%
    bind_rows(dat %>% select(studlab, treat = treat2, Ival = Ival2)) %>%
    distinct() %>%
    count(studlab, treat)
  #
  longtrt_Ival <-
    unique(
      rbind(setNames(dat[c("studlab", "treat1", "Ival1")],
                     c("studlab", "treat", "Ival")),
            setNames(dat[c("studlab", "treat2", "Ival2")],
                     c("studlab", "treat", "Ival"))
            )
      )
  #
  ct_unique_Ival <-
    as.data.frame(table(longtrt_Ival$studlab, longtrt_Ival$treat))
  names(ct_unique_Ival) <- c("studlab", "treat", "n")
  #
  # Multiple treatment specific Ivalues within a study
  #
  if (sum(ct_unique_Ival$n > 1) > 0) 
    stop(paste("Error in directionality assignment. ",
               "Multiple directions defined for identical treatments in the ",
               "following studies: ",
               paste(unique(ct_unique_Ival[ct_unique_Ival$n > 1, "studlab"]),
                     collapse = ", ")))
  #
  # Inconsistent Ivalues between arms within a study
  #
  dat$Ival_diff <- dat$Ival1 - dat$Ival2
  Ival_consistency_list <-
    lapply(split(dat, dat$studlab), function(x) any_invalid_I(x))
  #
  if (length(Ival_consistency_list[Ival_consistency_list==TRUE]) > 0)
    stop(paste("Error in directionality assignment. ",
               "Inconsistent directions defined in the following studies: ",
               paste(names(Ival_consistency_list[Ival_consistency_list==TRUE]),
                     collapse = ", ")))
  #
  if (sum(!dat$Ival_diff %in% c(-1, 0, 1)) > 0)
    stop("Error in directionality assignment. ",
         "At least one study has a difference between directionality values ",
         "which is not in the accepted values of -1, 0, 1.")
  #
  NULL
}

#alternatively, could harcode these so that this syntax is less dependent on metafor updates
expose_metafor_helpers <- function() {
  # List every internal function that the original rma.mv/.ll.rma.mv code uses.
  internal_funs <- c(
    ".set.digits",".getx",".process.G.afterrmna",".process.G.aftersub",".anyNAv",
    ".chkpd", ".make.unique", ".set.btt", ".expand1",".chkopt", ".chkconv",".con.E",
    ".chkvccon", ".chksubset", ".getsubset",".getfromenv"
  )
  for (fn in internal_funs) {
    obj <- getFromNamespace(fn, "metafor")
    assign(fn, obj, envir = .GlobalEnv)
  }
  invisible(TRUE)
}
expose_metafor_helpers()

rma.mv.exch<-function (yi, V, W, mods, data, slab, subset, random, struct = "CS", 
                       Z_slope=NULL, tau2B,
                       intercept = TRUE, method = "REML", test = "z", dfs = "residual", 
                       level = 95, btt, R, Rscale = "cor", sigma2, tau2, rho, gamma2, 
                       phi, cvvc = FALSE, sparse = FALSE, digits, 
                       control, ...) 
{
  # set the level
  level<-ifelse(level >= 1, (100 - level)/100, level) # in case it's on a 100 pt scale
  level<-ifelse(level > 0.5, 1 - level, level) 
  
  verbose= FALSE
  if (!is.element(method, c("FE", "EE", "CE", "ML", "REML"))) 
    stop("Unknown 'method' specified.")
  if (any(!is.element(struct, c("CS", "HCS", "UN", "AR", "HAR", 
                                "CAR", "ID", "DIAG", "SPEXP", "SPGAU", "SPLIN", "SPRAT", 
                                "SPSPH", "GEN", "GDIAG")))) 
    stop("Unknown 'struct' specified.")
  if (length(struct) == 1L) 
    struct <- c(struct, struct)
  na.act <- getOption("na.action")
  on.exit(options(na.action = na.act), add = TRUE)
  if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
                            "na.pass"))) 
    stop("Unknown 'na.action' specified under options().")
  if (missing(random)) 
    random <- NULL
  if (missing(R)) 
    R <- NULL
  if (missing(sigma2)) 
    sigma2 <- NULL
  if (missing(tau2)) 
    tau2 <- NULL
  if (missing(tau2B)) 
    tau2B <- NULL
  if (missing(rho)) 
    rho <- NULL
  if (missing(gamma2)) 
    gamma2 <- NULL
  if (missing(phi)) 
    phi <- NULL
  if (missing(control)) 
    control <- list()
  if (missing(digits)) {
    digits <- .set.digits(dmiss = TRUE)
  }
  else {
    digits <- .set.digits(digits, dmiss = FALSE)
  }
  time.start <- proc.time()
  args <- list(...)
  if (is.null(args$lambda)) {
    lambda <- 0
  }
  else {
    lambda <- args$lambda
  }
  if (isFALSE(args$tdist)) 
    test <- "z"
  if (isTRUE(args$tdist)) 
    test <- "t"
  test <- tolower(test)
  if (!is.element(test, c("z", "t", "knha", "hksj", "adhoc"))) 
    stop("Invalid option selected for 'test' argument.")
  if (test == "hksj") 
    test <- "knha"
  if (is.character(dfs)) 
    dfs <- match.arg(dfs, c("residual", "contain"))
  if (test == "z") {
    if (is.numeric(dfs)) {
      test <- "t"
    }
    else {
      if (dfs == "contain") 
        test <- "t"
    }
  }
  if (is.character(Rscale)) 
    Rscale <- match.arg(Rscale, c("none", "cor", "cor0", 
                                  "cov0"))
  if (is.logical(Rscale)) 
    Rscale <- ifelse(Rscale, "cor", "none")
  if (is.numeric(Rscale)) {
    Rscale <- round(Rscale)
    if (Rscale > 3 | Rscale < 0) 
      stop("Unknown 'Rscale' value specified.")
    Rscale <- switch(as.character(Rscale), `0` = "none", 
                     `1` = "cor", `2` = "cor0", `3` = "cov0")
  }
  if (is.null(args$dist)) {
    args$dist <- list("euclidean", "euclidean")
  }
  else {
    if (is.data.frame(args$dist) || is.matrix(args$dist)) 
      args$dist <- list(args$dist)
    if (!inherits(args$dist, "list")) 
      args$dist <- as.list(args$dist)
    if (length(args$dist) == 1L) 
      args$dist <- c(args$dist, args$dist)
    dist.methods <- c("euclidean", "maximum", "manhattan", 
                      "gcd")
    for (j in 1:2) {
      if (is.data.frame(args$dist[[j]])) 
        args$dist[[j]] <- as.matrix(args$dist[[j]])
      if (!is.function(args$dist[[j]]) && !is.matrix(args$dist[[j]])) {
        args$dist[[j]] <- charmatch(args$dist[[j]], dist.methods, 
                                   nomatch = 0)
        if (args$dist[[j]] == 0) {
          stop("Argument 'dist' must be one of 'euclidean', 'maximum', 'manhattan', or 'gcd'.")
        }
        else {
          args$dist[[j]] <- dist.methods[args$dist[[j]]]
        }
      }
    }
    if (any(args$dist == "gcd")) {
      if (!requireNamespace("sp", quietly = TRUE)) 
        stop("Please install the 'sp' package to compute great-circle distances.")
    }
  }
  if (is.null(args$vccon)) {
    vccon <- NULL
  }
  else {
    vccon <- args$vccon
    sigma2 <- .chkvccon(vccon$sigma2, sigma2)
    tau2 <- .chkvccon(vccon$tau2, tau2)
    rho <- .chkvccon(vccon$rho, rho)
    gamma2 <- .chkvccon(vccon$gamma2, gamma2)
    phi <- .chkvccon(vccon$phi, phi)
  }
  formula.yi <- NULL
  formula.mods <- NULL
  con <- list(verbose = FALSE, optimizer = "nlminb", optmethod = "BFGS", 
              parallel = list(), cl = NULL, ncpus = 1L, REMLf = TRUE, 
              evtol = 1e-07, nearpd = FALSE, hessianCtrl = list(r = 8), 
              hesstol = .Machine$double.eps^0.5, hesspack = "numDeriv", 
              check.k.gtr.1 = TRUE)
  con.pos <- pmatch(names(control), names(con))
  con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
  if (missing(data)) 
    data <- NULL
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  }
  else {
    if (!is.data.frame(data)) 
      data <- data.frame(data)
  }
  mf <- match.call()
  yi <- .getx("yi", mf = mf, data = data)
  V <- .getx("V", mf = mf, data = data)
  W <- .getx("W", mf = mf, data = data)
  ni <- .getx("ni", mf = mf, data = data)
  slab <- .getx("slab", mf = mf, data = data)
  subset <- .getx("subset", mf = mf, data = data)
  mods <- .getx("mods", mf = mf, data = data)
  if (inherits(yi, "formula")) {
    formula.yi <- yi
    formula.mods <- formula.yi[-2]
    options(na.action = "na.pass")
    mods <- model.matrix(yi, data = data)
    attr(mods, "assign") <- NULL
    attr(mods, "contrasts") <- NULL
    yi <- model.response(model.frame(yi, data = data))
    options(na.action = na.act)
    names(yi) <- NULL
    intercept <- FALSE
  }
  if (is.data.frame(yi)) {
    if (ncol(yi) == 1L) {
      yi <- yi[[1]]
    }
    else {
      stop("The object/variable specified for the 'yi' argument is a data frame with multiple columns.")
    }
  }
  if (is.matrix(yi)) {
    if (nrow(yi) == 1L || ncol(yi) == 1L) {
      yi <- as.vector(yi)
    }
    else {
      stop("The object/variable specified for the 'yi' argument is a matrix with multiple rows/columns.")
    }
  }
  if (inherits(yi, "array")) 
    stop("The object/variable specified for the 'yi' argument is an array.")
  if (!is.numeric(yi)) 
    stop("The object/variable specified for the 'yi' argument is not numeric.")
  k <- length(yi)
  k.all <- k
  measure <- "GEN"
  if (!is.null(attr(yi, "measure"))) 
    measure <- attr(yi, "measure")
  attr(yi, "measure") <- measure
  if (is.null(V)) 
    stop("Must specify the 'V' argument.")
  if (identical(V, utils::vi)) 
    stop("Variable specified for 'V' argument cannot be found.")
  if (is.list(V) && !is.data.frame(V)) {
    V <- lapply(V, as.matrix)
    if (any(!sapply(V, function(X){NROW(X) == NCOL(X)}))) 
      stop("All list elements in 'V' must be square matrices.")
    if (sparse) {
      V <- bdiag(V)
    }
    else {
      V <- bldiag(V)
    }
  }
  if ((is.vector(V) && length(V) == 1L && V == 0) || (is.vector(V) && 
                                                       length(V) == k && !anyNA(V) && all(V == 0))) {
    V0 <- TRUE
  }
  else {
    V0 <- FALSE
  }
  if (V0 || is.vector(V) || nrow(V) == 1L || ncol(V) == 1L) {
    if (sparse) {
      V <- Diagonal(k, as.vector(V))
    }
    else {
      V <- diag(as.vector(V), nrow = k, ncol = k)
    }
  }
  if (is.data.frame(V)) 
    V <- as.matrix(V)
  if (!is.null(dimnames(V))) 
    V <- unname(V)
  if (!V0 && !NROW(V) == NCOL(V)) 
    stop(("'V' must be a square matrix."))
  if (!V0 && !isSymmetric(V)) 
    stop(("'V' must be a symmetric matrix."))
  if (nrow(V) != k) 
    stop((paste0("Length of 'yi' (", k, ") and the length/dimensions of 'V' (", 
                            nrow(V), ") are not the same.")))
  if (sparse && inherits(V, "matrix")) 
    V <- Matrix(V, sparse = TRUE)
  if (inherits(V, "matrix") && !is.numeric(V)) 
    stop(("The object/variable specified for the 'V' argument is not numeric."))
  if (!is.null(W)) {
    if (is.vector(W) || nrow(W) == 1L || ncol(W) == 1L) {
      W <- as.vector(W)
      W <- .expand1(W, k)
      A <- diag(W, nrow = length(W), ncol = length(W))
    }
    else {
      A <- W
    }
    if (is.data.frame(A)) 
      A <- as.matrix(A)
    if (!is.null(dimnames(A))) 
      A <- unname(A)
    if (!NROW(A) == NCOL(A)) 
      stop("'W' must be a square matrix.")
    if (!isSymmetric(A)) 
      stop("'W' must be a symmetric matrix.")
    if (nrow(A) != k) 
      stop(paste0("Length of 'yi' (", k, ") and length/dimensions of 'W' (", 
                              nrow(A), ") are not the same."))
    if (sparse && inherits(A, "matrix")) 
      A <- Matrix(A, sparse = TRUE)
    if (inherits(A, "matrix") && !is.numeric(A)) 
      stop("The object/variable specified for the 'W' argument is not numeric.")
  }
  else {
    A <- NULL
  }
  if (is.null(ni)) 
    ni <- attr(yi, "ni")
  if (!is.null(ni) && length(ni) != k) 
    ni <- NULL
  if (inherits(mods, "formula")) {
    formula.mods <- mods
    if (isTRUE(all.equal(formula.mods, ~1))) {
      mods <- matrix(1, nrow = k, ncol = 1)
      intercept <- FALSE
    }
    else {
      options(na.action = "na.pass")
      mods <- model.matrix(mods, data = data)
      attr(mods, "assign") <- NULL
      attr(mods, "contrasts") <- NULL
      options(na.action = na.act)
      intercept <- FALSE
    }
  }
  if (is.vector(mods)) 
    mods <- cbind(mods)
  if (is.data.frame(mods)) 
    mods <- as.matrix(mods)
  if (is.character(mods)) 
    stop("Model matrix contains character variables.")
  if (!is.null(mods) && nrow(mods) != k) 
    stop(paste0("Number of rows in the model matrix (", 
                            nrow(mods), ") do not match the length of the outcome vector (", 
                            k, ")."))
  if (!is.element(method, c("FE", "EE", "CE")) && !is.null(random)) {
    if (!is.list(random)) 
      random <- list(random)
    if (any(sapply(random, function(x) !inherits(x, "formula")))) 
      stop("All elements of 'random' must be formulas.")
    has.vbar <- sapply(random, function(f) grepl("|", paste0(f, 
                                                             collapse = ""), fixed = TRUE))
    if (any(!has.vbar)) 
      stop("All formulas in 'random' must contain a grouping variable after the | symbol.")
    has.dollar <- sapply(random, function(f) grepl("$", 
                                                   paste0(f, collapse = ""), fixed = TRUE))
    if (any(has.dollar)) 
      stop("Cannot use '$' notation in formulas in the 'random' argument (use the 'data' argument instead).")
    has.colon <- sapply(random, function(f) grepl(":", paste0(f, 
                                                              collapse = ""), fixed = TRUE))
    if (any(has.colon)) 
      stop("Cannot use ':' notation in formulas in the 'random' argument (use 'interaction()' instead).")
    has.in <- sapply(random, function(f) grepl("%in%", paste0(f, 
                                                              collapse = ""), fixed = TRUE))
    if (any(has.in)) 
      stop("Cannot use '%in%' notation in formulas in the 'random' argument (use 'interaction()' instead).")
    has.dblvbar <- sapply(random, function(f) grepl("||", 
                                                    paste0(f, collapse = ""), fixed = TRUE))
    random <- lapply(random, function(f) {
      if (grepl("||", paste0(f, collapse = ""), fixed = TRUE)) {
        f <- paste0(f, collapse = "")
        f <- gsub("||", "|", f, fixed = TRUE)
        f <- as.formula(f)
      }
      return(f)
    })
    formulas <- list(NULL, NULL)
    split.formulas <- sapply(random, function(f) strsplit(paste0(f, 
                                                                 collapse = ""), " | ", fixed = TRUE))
    is.inner.outer <- sapply(split.formulas, function(f) f[1] != 
                               "~1")
    if (sum(is.inner.outer) > 2) 
      stop("Only up to two '~ inner | outer' formulas allowed in the 'random' argument.")
    if (any(is.inner.outer)) 
      formulas[[1]] <- random[is.inner.outer][1][[1]]
    if (sum(is.inner.outer) == 2) 
      formulas[[2]] <- random[is.inner.outer][2][[1]]
    has.slash <- sapply(random, function(f) grepl("/", paste0(f, 
                                                              collapse = ""), fixed = TRUE))
    if (any(is.inner.outer & has.slash)) 
      stop("Cannot use '~ inner | outer1/outer2' type terms in the 'random' argument.")
    random.plus <- lapply(random, function(f) formula(sub("\\|", 
                                                          "+", paste0(f, collapse = ""))))
    mf.r <- list()
    io <- 0
    for (j in seq_along(is.inner.outer)) {
      if (is.inner.outer[j]) {
        io <- io + 1
        if (is.element(struct[io], c("GEN", "GDIAG"))) {
          f.inner <- as.formula(strsplit(paste(random[[j]], 
                                               collapse = ""), " | ", fixed = TRUE)[[1]][1])
          f.outer <- as.formula(paste("~", strsplit(paste(random[[j]], 
                                                          collapse = ""), " | ", fixed = TRUE)[[1]][2]))
          options(na.action = "na.pass")
          X.inner <- model.matrix(f.inner, data = data)
          options(na.action = na.act)
          is.int <- apply(X.inner, 2, function(x) all(abs(x - 1) < 1e-08))
          colnames(X.inner)[is.int] <- "intrcpt"
          mf.r[[j]] <- cbind(X.inner, model.frame(f.outer, 
                                                  data = data, na.action = na.pass))
          if (has.dblvbar[j]) 
            struct[io] <- "GDIAG"
        }
        else {
          mf.r[[j]] <- model.frame(random.plus[[j]], 
                                   data = data, na.action = na.pass)
        }
      }
      else {
        mf.r[[j]] <- model.frame(random.plus[[j]], data = data, 
                                 na.action = na.pass)
      }
    }
    mf.r.ncols <- sapply(mf.r, ncol)
    for (j in seq_along(has.slash)) {
      if (!has.slash[j]) 
        next
      for (p in mf.r.ncols[j]:1) {
        mf.r[[j]][, p] <- interaction(mf.r[[j]][1:p], 
                                      drop = TRUE, lex.order = TRUE, sep = "/")
        colnames(mf.r[[j]])[p] <- paste(colnames(mf.r[[j]])[1:p], 
                                        collapse = "/")
      }
    }
    if (any(has.slash)) {
      if (length(mf.r) == 1L) {
        mf.r <- lapply(seq(ncol(mf.r[[1]])), function(x) mf.r[[1]][x])
      }
      else {
        mf.r <- unlist(mapply(function(mf, sl) if (sl) 
          lapply(seq(mf), function(x) mf[x])
          else list(mf), mf.r, has.slash, SIMPLIFY = FALSE), 
          recursive = FALSE, use.names = FALSE)
      }
      mf.r.ncols <- sapply(mf.r, ncol)
    }
    mf.s <- mf.r[which(mf.r.ncols == 1)]
    mf.g <- mf.r[[which(mf.r.ncols >= 2)[1]]]
    mf.h <- mf.r[[which(mf.r.ncols >= 2)[2]]]
    if (length(mf.s) == 0L) 
      mf.s <- NULL
    withS <- !is.null(mf.s)
    withG <- !is.null(mf.g)
    withH <- !is.null(mf.h)
    withZslope <- !is.null(Z_slope)
    mf.r.nrows <- sapply(mf.r, nrow)
    if (any(mf.r.nrows != k)) 
      stop("Length of the variables specified via the 'random' argument does not match the length of the data.")
    mf.r <- lapply(random.plus, get_all_vars, data = data)
  }
  else {
    formulas <- list(NULL, NULL)
    mf.r <- NULL
    mf.s <- NULL
    mf.g <- NULL
    mf.h <- NULL
    withS <- FALSE
    withG <- FALSE
    withH <- FALSE
  }
  mf.r <- unname(mf.r)
  if (!withG && "struct" %in% names(mf)) 
    warning(("Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded."), 
            call. = FALSE)
  if (is.element(method, c("FE", "EE", "CE")) && "random" %in% 
      names(mf)) 
    warning((paste0("The 'random' argument is disregaded when method=\"", 
                                  method, "\".")), call. = FALSE)
  ids <- seq_len(k)
  if (is.null(slab)) {
    slab <- attr(yi, "slab")
    if (!is.null(slab) && length(slab) != k) 
      slab <- NULL
  }
  if (is.null(slab)) {
    slab.null <- TRUE
    slab <- ids
  }
  else {
    if (anyNA(slab)) 
      stop("NAs in study labels.")
    if (length(slab) != k) 
      stop(paste0("Length of the 'slab' argument (", 
                              length(slab), ") does not correspond to the size of the dataset (", 
                              k, ")."))
    if (is.factor(slab)) 
      slab <- as.character(slab)
    slab.null <- FALSE
  }
  if (!is.null(subset)) {
    subset <- .chksubset(subset, k)
    yi <- .getsubset(yi, subset)
    V <- .getsubset(V, subset, col = TRUE)
    A <- .getsubset(A, subset, col = TRUE)
    ni <- .getsubset(ni, subset)
    mods <- .getsubset(mods, subset)
    slab <- .getsubset(slab, subset)
    mf.r <- lapply(mf.r, .getsubset, subset)
    mf.s <- lapply(mf.s, .getsubset, subset)
    mf.g <- .getsubset(mf.g, subset)
    mf.h <- .getsubset(mf.h, subset)
    ids <- .getsubset(ids, subset)
    k <- length(yi)
    attr(yi, "measure") <- measure
    attr(yi, "ni") <- ni
  }
  if (anyDuplicated(slab)) 
    slab <- .make.unique(slab)
  attr(yi, "slab") <- slab
  vi <- diag(V)
  yi.f <- yi
  vi.f <- vi
  V.f <- V
  W.f <- A
  ni.f <- ni
  mods.f <- mods
  k.f <- k
  if (withS) {
    s.names <- sapply(mf.s, names)
    mf.s <- lapply(mf.s, function(x) factor(x[[1]]))
    if (any(sapply(mf.s, anyNA))) 
      stop("No NAs allowed in variables specified in the 'random' argument.")
    sigma2s <- length(mf.s)
    if (is.null(sigma2)) 
      sigma2 <- rep(NA_real_, sigma2s)
    sigma2 <- .expand1(sigma2, sigma2s)
    if (length(sigma2) != sigma2s) 
      stop(paste0("Length of the 'sigma2' argument (", 
                              length(sigma2), ") does not match the actual number of variance components (", 
                              sigma2s, ")."))
    if (any(sigma2 < 0, na.rm = TRUE)) 
      stop("Specified value(s) of 'sigma2' must be non-negative.")
    s.nlevels <- sapply(mf.s, nlevels)
    s.levels <- lapply(mf.s, levels)
    if (is.null(R)) {
      withR <- FALSE
      Rfix <- rep(FALSE, sigma2s)
    }
    else {
      withR <- TRUE
      if (is.data.frame(R) || !is.list(R)) 
        R <- list(R)
      if (is.null(names(R)) || any(nchar(names(R)) == 
                                   0L)) 
        stop(("Argument 'R' must be a *named* list."))
      R <- R[!sapply(R, is.null)]
      R <- lapply(R, as.matrix)
      R <- R[s.names]
      names(R) <- s.names
      Rfix <- !sapply(R, is.null)
      if (any(Rfix)) {
        if (any(!sapply(R[Rfix], function(X){NROW(X) == NCOL(X)}))) 
          stop(("Elements of 'R' must be square matrices."))
        if (any(!sapply(R[Rfix], function(x) isSymmetric(unname(x))))) 
          stop(("Elements of 'R' must be symmetric matrices."))
        for (j in seq_along(R)) {
          if (!Rfix[j]) 
            next
          R[[j]] <- symmpart(R[[j]])
          if (is.null(rownames(R[[j]]))) 
            rownames(R[[j]]) <- colnames(R[[j]])
          if (is.null(colnames(R[[j]]))) 
            colnames(R[[j]]) <- rownames(R[[j]])
          if (is.null(colnames(R[[j]]))) 
            stop(("Elements of 'R' must have dimension names."))
        }
        R[Rfix] <- lapply(R[Rfix], function(x) x[!duplicated(rownames(x)), 
                                                 !duplicated(colnames(x)), drop = FALSE])
        if (any(sapply(R[Rfix], function(x) length(colnames(x)) != 
                       length(unique(colnames(x)))))) 
          stop(("Each element of 'R' must have unique dimension names."))
        for (j in seq_along(R)) {
          if (!Rfix[j]) 
            next
          if (anyNA(R[[j]])) 
            stop(("No missing values allowed in matrices specified via 'R'."))
          if (any(!is.element(s.levels[[j]], colnames(R[[j]])))) 
            stop((paste0("There are levels in '", 
                                    s.names[j], "' for which there are no matching rows/columns in the corresponding 'R' matrix.")))
          if (any(!is.element(colnames(R[[j]]), s.levels[[j]]))) 
            warning((paste0("There are rows/columns in the 'R' matrix for '", 
                                          s.names[j], "' for which there are no data.")), 
                    call. = FALSE)
        }
      }
      else {
        warning(("Argument 'R' specified, but list name(s) not in 'random'."), 
                call. = FALSE)
        withR <- FALSE
        Rfix <- rep(FALSE, sigma2s)
        R <- NULL
      }
    }
  }
  else {
    sigma2s <- 1
    sigma2 <- 0
    s.nlevels <- NULL
    s.levels <- NULL
    s.names <- NULL
    withR <- FALSE
    Rfix <- FALSE
    R <- NULL
  }
  s.nlevels.f <- s.nlevels
  s.levels.f <- s.levels
  if (withG) {
    tmp <- .process.G.aftersub(mf.g, struct[1], formulas[[1]], 
                               tau2, rho, isG = TRUE, k, sparse, verbose=FALSE)
    mf.g <- tmp$mf.g
    g.names <- tmp$g.names
    g.nlevels <- tmp$g.nlevels
    g.levels <- tmp$g.levels
    g.values <- tmp$g.values
    tau2s <- tmp$tau2s
    rhos <- tmp$rhos
    tau2 <- tmp$tau2
    rho <- tmp$rho
    Z.G1 <- tmp$Z.G1
    Z.G2 <- tmp$Z.G2
  }
  else {
    tau2s <- 1
    rhos <- 1
    tau2 <- 0
    rho <- 0
    Z.G1 <- NULL
    Z.G2 <- NULL
    g.nlevels <- NULL
    g.levels <- NULL
    g.values <- NULL
    g.names <- NULL
  }
  mf.g.f <- mf.g
  if (withH) {
    tmp <- .process.G.aftersub(mf.h, struct[2], formulas[[2]], 
                               gamma2, phi, isG = FALSE, k, sparse, verbose=FALSE)
    mf.h <- tmp$mf.g
    h.names <- tmp$g.names
    h.nlevels <- tmp$g.nlevels
    h.levels <- tmp$g.levels
    h.values <- tmp$g.values
    gamma2s <- tmp$tau2s
    phis <- tmp$rhos
    gamma2 <- tmp$tau2
    phi <- tmp$rho
    Z.H1 <- tmp$Z.G1
    Z.H2 <- tmp$Z.G2
  }
  else {
    gamma2s <- 1
    phis <- 1
    gamma2 <- 0
    phi <- 0
    Z.H1 <- NULL
    Z.H2 <- NULL
    h.nlevels <- NULL
    h.levels <- NULL
    h.values <- NULL
    h.names <- NULL
  }
  mf.h.f <- mf.h
  if (withZslope) {
    # only allowing 1 heterogeneity for 1 interaction
    tau2Bs <- 1
  } else {
    Z_slope<-NULL
    tau2Bs <- 1
    tau2B <- 0
  }
  has.na <- is.na(yi) | (if (is.null(mods)) 
    FALSE
    else apply(is.na(mods), 1, any)) | (if (V0) 
      FALSE
      else .anyNAv(V)) | (if (is.null(A)) 
        FALSE
        else apply(is.na(A), 1, any))
  not.na <- !has.na
  if (any(has.na)) {
    if (na.act == "na.omit" || na.act == "na.exclude" || 
        na.act == "na.pass") {
      yi <- yi[not.na]
      V <- V[not.na, not.na, drop = FALSE]
      A <- A[not.na, not.na, drop = FALSE]
      vi <- vi[not.na]
      ni <- ni[not.na]
      mods <- mods[not.na, , drop = FALSE]
      mf.r <- lapply(mf.r, function(x) x[not.na, , drop = FALSE])
      mf.s <- lapply(mf.s, function(x) x[not.na])
      mf.g <- mf.g[not.na, , drop = FALSE]
      mf.h <- mf.h[not.na, , drop = FALSE]
      if (is.element(struct[1], c("SPEXP", "SPGAU", "SPLIN", 
                                  "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD"))) {
        Z.G1 <- Z.G1[not.na, not.na, drop = FALSE]
      }
      else {
        Z.G1 <- Z.G1[not.na, , drop = FALSE]
      }
      Z.G2 <- Z.G2[not.na, , drop = FALSE]
      if (is.element(struct[2], c("SPEXP", "SPGAU", "SPLIN", 
                                  "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD"))) {
        Z.H1 <- Z.H1[not.na, not.na, drop = FALSE]
      }
      else {
        Z.H1 <- Z.H1[not.na, , drop = FALSE]
      }
      Z.H2 <- Z.H2[not.na, , drop = FALSE]
      k <- length(yi)
      warning((paste(sum(has.na), ifelse(sum(has.na) > 
                                                         1, "rows", "row"), "with NAs omitted from model fitting.")), 
              call. = FALSE)
      attr(yi, "measure") <- measure
      attr(yi, "ni") <- ni
    }
    if (na.act == "na.fail") 
      stop(("Missing values in data."))
  }
  if (k <= 1) 
    stop(("Processing terminated since k <= 1."))
  if (any(vi <= 0)) {
    allvipos <- FALSE
    if (!V0) 
      warning(("There are outcomes with non-positive sampling variances."), 
              call. = FALSE)
    vi.neg <- vi < 0
    if (any(vi.neg)) {
      V[vi.neg, ] <- 0
      V[, vi.neg] <- 0
      vi[vi.neg] <- 0
      warning(("Negative sampling variances constrained to 0."), 
              call. = FALSE)
    }
  }
  else {
    allvipos <- TRUE
  }
  if (!V0 && !.chkpd(V)) 
    warning(("'V' appears to be not positive definite."), 
            call. = FALSE)
  vimaxmin <- max(vi)/min(vi)
  if (is.finite(vimaxmin) && vimaxmin >= 1e+07) 
    warning(("Ratio of largest to smallest sampling variance extremely large. May not be able to obtain stable results."), 
            call. = FALSE)
  if (is.null(mods) && !intercept) {
    warning(("Must either include an intercept and/or moderators in model.\nCoerced intercept into the model."), 
            call. = FALSE)
    intercept <- TRUE
  }
  if (!is.null(mods) && ncol(mods) == 0L) {
    warning(("Cannot fit model with an empty model matrix. Coerced intercept into the model."), 
            call. = FALSE)
    intercept <- TRUE
  }
  if (intercept) {
    X <- cbind(intrcpt = rep(1, k), mods)
    X.f <- cbind(intrcpt = rep(1, k.f), mods.f)
  }
  else {
    X <- mods
    X.f <- mods.f
  }
  tmp <- try(lm(yi ~ X - 1), silent = TRUE)
  if (inherits(tmp, "try-error")) {
    stop(("Error in check for redundant predictors."))
  }
  else {
    coef.na <- is.na(coef(tmp))
    if (any(coef.na)) {
      warning(("Redundant predictors dropped from the model."), 
              call. = FALSE)
      X <- X[, !coef.na, drop = FALSE]
      X.f <- X.f[, !coef.na, drop = FALSE]
    }
  }
  is.int <- apply(X, 2, function(x) all(abs(x - 1) < 1e-08))
  if (any(is.int)) {
    int.incl <- TRUE
    int.indx <- which(is.int, arr.ind = TRUE)
    X <- cbind(intrcpt = 1, X[, -int.indx, drop = FALSE])
    X.f <- cbind(intrcpt = 1, X.f[, -int.indx, drop = FALSE])
    intercept <- TRUE
  }
  else {
    int.incl <- FALSE
  }
  if (!.chkpd(crossprod(X), tol = con$evtol)) 
    stop(("Model matrix not of full rank. Cannot fit model."))
  p <- NCOL(X)
  colnames(X) <- colnames(X.f) <- .make.unique(colnames(X))
  if ((p == 1L) && all(abs(X - 1) < 1e-08)) {
    int.only <- TRUE
  }
  else {
    int.only <- FALSE
  }
  btt <- .set.btt(btt, p, int.incl, colnames(X))
  m <- length(btt)
  optbeta <- FALSE
  if (is.null(args$beta)) {
    beta.arg <- rep(NA_real_, p)
    beta.est <- rep(TRUE, p)
  }
  else {
    beta.arg <- args$beta
    if (length(beta.arg) != p) 
      stop((paste0("Length of the 'beta' argument (", 
                              length(beta.arg), ") does not match the actual number of fixed effects (", 
                              p, ").")))
    beta.est <- is.na(beta.arg)
  }
  if (withS) {
    mf.s <- lapply(mf.s, factor)
    s.nlevels <- sapply(mf.s, nlevels)
    s.levels <- lapply(mf.s, levels)
    if (any(is.na(sigma2) & s.nlevels == 1) && con$check.k.gtr.1) {
      sigma2[is.na(sigma2) & s.nlevels == 1] <- 0
      warning(("Single-level factor(s) found in 'random' argument. Corresponding 'sigma2' value(s) fixed to 0."), 
              call. = FALSE)
    }
    Z.S <- vector(mode = "list", length = sigma2s)
    for (j in seq_len(sigma2s)) {
      if (s.nlevels[j] == 1) {
        Z.S[[j]] <- cbind(rep(1, k))
      }
      else {
        if (sparse) {
          Z.S[[j]] <- sparse.model.matrix(~mf.s[[j]] - 
                                            1)
        }
        else {
          Z.S[[j]] <- model.matrix(~mf.s[[j]] - 1)
        }
      }
      attr(Z.S[[j]], "assign") <- NULL
      attr(Z.S[[j]], "contrasts") <- NULL
    }
  }
  else {
    Z.S <- NULL
  }
  if (withR) {
    for (j in seq_along(R)) {
      if (!Rfix[j]) 
        next
      R[[j]] <- R[[j]][s.levels[[j]], s.levels[[j]]]
    }
    if (Rscale == "cor" || Rscale == "cor0") {
      R[Rfix] <- lapply(R[Rfix], function(x) {
        if (any(diag(x) <= 0)) 
          stop(("Cannot use Rscale=\"cor\" or Rscale=\"cor0\" with non-positive values on the diagonal of an 'R' matrix."))
        tmp <- cov2cor(x)
        if (any(abs(tmp) > 1)) 
          warning(("Some values are larger than +-1 in an 'R' matrix after cov2cor() (see 'Rscale' argument)."), 
                  call. = FALSE)
        return(tmp)
      })
    }
    if (Rscale == "cor0") 
      R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x))/(1 - 
                                                             min(x)))
    if (Rscale == "cov0") 
      R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x)))
  }
  if (withS) {
    D.S <- vector(mode = "list", length = sigma2s)
    for (j in seq_len(sigma2s)) {
      if (Rfix[j]) {
        if (sparse) {
          D.S[[j]] <- Z.S[[j]] %*% Matrix(R[[j]], sparse = TRUE) %*% 
            t(Z.S[[j]])
        }
        else {
          D.S[[j]] <- Z.S[[j]] %*% R[[j]] %*% t(Z.S[[j]])
        }
      }
      else {
        D.S[[j]] <- tcrossprod(Z.S[[j]])
      }
    }
  }
  else {
    D.S <- NULL
  }
  if (withG) {
    tmp <- .process.G.afterrmna(mf.g, g.nlevels, g.levels, 
                                g.values, struct[1], formulas[[1]], tau2, rho, Z.G1, 
                                Z.G2, isG = TRUE, sparse, args$dist[[1]], con$check.k.gtr.1, 
                                verbose=FALSE)
    mf.g <- tmp$mf.g
    g.nlevels <- tmp$g.nlevels
    g.nlevels.f <- tmp$g.nlevels.f
    g.levels <- tmp$g.levels
    g.levels.f <- tmp$g.levels.f
    g.levels.r <- tmp$g.levels.r
    g.levels.k <- tmp$g.levels.k
    g.levels.comb.k <- tmp$g.levels.comb.k
    tau2 <- tmp$tau2
    rho <- tmp$rho
    G <- tmp$G
    g.Dmat <- tmp$Dmat
    g.rho.init <- tmp$rho.init
  }
  else {
    g.nlevels.f <- NULL
    g.levels.f <- NULL
    g.levels.r <- NULL
    g.levels.k <- NULL
    g.levels.comb.k <- NULL
    G <- NULL
    g.Dmat <- NULL
    g.rho.init <- NULL
  }
  if (withH) {
    tmp <- .process.G.afterrmna(mf.h, h.nlevels, h.levels, 
                                h.values, struct[2], formulas[[2]], gamma2, phi, 
                                Z.H1, Z.H2, isG = FALSE, sparse, args$dist[[2]], 
                                con$check.k.gtr.1, verbose=FALSE)
    mf.h <- tmp$mf.g
    h.nlevels <- tmp$g.nlevels
    h.nlevels.f <- tmp$g.nlevels.f
    h.levels <- tmp$g.levels
    h.levels.f <- tmp$g.levels.f
    h.levels.r <- tmp$g.levels.r
    h.levels.k <- tmp$g.levels.k
    h.levels.comb.k <- tmp$g.levels.comb.k
    gamma2 <- tmp$tau2
    phi <- tmp$rho
    H <- tmp$G
    h.Dmat <- tmp$Dmat
    h.phi.init <- tmp$rho.init
  }
  else {
    h.nlevels.f <- NULL
    h.levels.f <- NULL
    h.levels.r <- NULL
    h.levels.k <- NULL
    h.levels.comb.k <- NULL
    H <- NULL
    h.Dmat <- NULL
    h.phi.init <- NULL
  }
  Y <- as.matrix(yi)
  QE <- NA_real_
  if (!V0) {
      U <- try(suppressWarnings(chol(chol2inv(chol(V)))), 
               silent = TRUE)
  }
  if (V0 || inherits(U, "try-error") || any(is.infinite(U))) {
    total <- sigma(lm(Y ~ X - 1))^2
    if (is.na(total)) 
      total <- var(as.vector(Y))/100
  }
  else {
    sX <- U %*% X
    sY <- U %*% Y
    beta.FE <- try(solve(crossprod(sX), crossprod(sX, sY)), 
                   silent = TRUE)
    if (inherits(beta.FE, "try-error")) {
      total <- var(as.vector(Y))
    }
    else {
      total <- max(0.001 * (sigma2s + tau2s + gamma2s), 
                   var(as.vector(Y) - as.vector(X %*% beta.FE)) - 
                     1/mean(1/diag(V)))
      QE <- sum(as.vector(sY - sX %*% beta.FE)^2)
    }
  }
  sigma2.init <- rep(total/(sigma2s + tau2s + gamma2s+tau2Bs), sigma2s)
  tau2.init <- rep(total/(sigma2s + tau2s + gamma2s+tau2Bs), tau2s)
  tau2B.init <- rep(total/(sigma2s + tau2s + gamma2s+tau2Bs), tau2s)
  gamma2.init <- rep(total/(sigma2s + tau2s + gamma2s+tau2Bs), gamma2s)
  if (is.null(g.rho.init)) {
    rho.init <- rep(0.5, rhos)
  }
  else {
    rho.init <- g.rho.init
  }
  if (is.null(h.phi.init)) {
    phi.init <- rep(0.5, phis)
  }
  else {
    phi.init <- h.phi.init
  }
  con <- c(con, list(sigma2.init = sigma2.init, tau2.init = tau2.init, tau2B.init = tau2B.init,
                     rho.init = rho.init, gamma2.init = gamma2.init, phi.init = phi.init, 
                     cholesky = ifelse(is.element(struct, c("UN", "UNR", 
                                                            "GEN")), TRUE, FALSE)))
  con.pos <- pmatch(names(control), names(con))
  con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]
  if (isTRUE(args$restart)) {
    okrestart <- TRUE
    if (withS && (is.null(.getfromenv("rma.mv", "sigma2")) || 
                  length(.getfromenv("rma.mv", "sigma2")) != sigma2s)) 
      okrestart <- FALSE
    if (withG && (is.null(.getfromenv("rma.mv", "tau2")) || 
                  length(.getfromenv("rma.mv", "tau2")) != tau2s)) 
      okrestart <- FALSE
    if (withG && (is.null(.getfromenv("rma.mv", "rho")) || 
                  length(.getfromenv("rma.mv", "rho")) != rhos)) 
      okrestart <- FALSE
    if (withH && (is.null(.getfromenv("rma.mv", "gamma2")) || 
                  length(.getfromenv("rma.mv", "gamma2")) != gamma2s)) 
      okrestart <- FALSE
    if (withH && (is.null(.getfromenv("rma.mv", "phi")) || 
                  length(.getfromenv("rma.mv", "phi")) != phis)) 
      okrestart <- FALSE
    if (!okrestart) 
      stop((paste0("Restarting for a different model than the initial one.")))
    con$sigma2.init <- .getfromenv("rma.mv", "sigma2", default = con$sigma2.init)
    con$tau2.init <- .getfromenv("rma.mv", "tau2", default = con$tau2.init)
    con$rho.init <- .getfromenv("rma.mv", "rho", default = con$rho.init)
    con$gamma2.init <- .getfromenv("rma.mv", "gamma2", default = con$gamma2.init)
    con$phi.init <- .getfromenv("rma.mv", "phi", default = con$phi.init)
  }
  if (anyNA(con$sigma2.init)) 
    stop((paste0("No missing values allowed in 'sigma2.init'.")))
  if (anyNA(con$tau2.init)) 
    stop((paste0("No missing values allowed in 'tau2.init'.")))
  if (anyNA(con$rho.init)) 
    stop((paste0("No missing values allowed in 'rho.init'.")))
  if (anyNA(con$gamma2.init)) 
    stop((paste0("No missing values allowed in 'gamma2.init'.")))
  if (anyNA(con$phi.init)) 
    stop((paste0("No missing values allowed in 'phi.init'.")))
  con$sigma2.init <- .expand1(con$sigma2.init, sigma2s)
  con$tau2.init <- .expand1(con$tau2.init, tau2s)
  con$rho.init <- .expand1(con$rho.init, rhos)
  con$gamma2.init <- .expand1(con$gamma2.init, gamma2s)
  con$phi.init <- .expand1(con$phi.init, phis)
  if (withS && any(con$sigma2.init <= 0)) 
    stop(("Value(s) of 'sigma2.init' must be > 0"))
  if (withG && any(con$tau2.init <= 0)) 
    stop(("Value(s) of 'tau2.init' must be > 0."))
  if (withG && struct[1] == "CAR" && (con$rho.init <= 0 | 
                                      con$rho.init >= 1)) 
    stop(("Value(s) of 'rho.init' must be in (0,1)."))
  if (withG && is.element(struct[1], c("SPEXP", "SPGAU", "SPLIN", 
                                       "SPRAT", "SPSPH")) && any(con$rho.init <= 0)) 
    stop(("Value(s) of 'rho.init' must be > 0."))
  if (withG && is.element(struct[1], c("PHYPL", "PHYPD")) && 
      con$rho.init < 0) 
    stop(("Value(s) of 'rho.init' must be in >= 0."))
  if (withG && !is.element(struct[1], c("CAR", "SPEXP", "SPGAU", 
                                        "SPLIN", "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD")) && 
      any(con$rho.init <= -1 | con$rho.init >= 1)) 
    stop(("Value(s) of 'rho.init' must be in (-1,1)."))
  if (withH && any(con$gamma2.init <= 0)) 
    stop(("Value(s) of 'gamma2.init' must be > 0."))
  if (withH && struct[2] == "CAR" && (con$phi.init <= 0 | 
                                      con$phi.init >= 1)) 
    stop(("Value(s) of 'phi.init' must be in (0,1)."))
  if (withH && is.element(struct[2], c("SPEXP", "SPGAU", "SPLIN", 
                                       "SPRAT", "SPSPH")) && any(con$phi.init <= 0)) 
    stop(("Value(s) of 'phi.init' must be > 0."))
  if (withH && is.element(struct[2], c("PHYPL", "PHYPD")) && 
      con$phi.init < 0) 
    stop(("Value(s) of 'phi.init' must be in >= 0."))
  if (withH && !is.element(struct[2], c("CAR", "SPEXP", "SPGAU", 
                                        "SPLIN", "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD")) && 
      any(con$phi.init <= -1 | con$phi.init >= 1)) 
    stop(("Value(s) of 'phi.init' must be in (-1,1)."))
  con$cholesky <- .expand1(con$cholesky, 2L)
  if (!withG) 
    con$cholesky[1] <- FALSE
  if (con$cholesky[1] && !is.element(struct[1], c("UN", "UNR", 
                                                  "GEN"))) 
    con$cholesky[1] <- FALSE
  if (!withH) 
    con$cholesky[2] <- FALSE
  if (con$cholesky[2] && !is.element(struct[2], c("UN", "UNR", 
                                                  "GEN"))) 
    con$cholesky[2] <- FALSE
  sigma2.init <- con$sigma2.init
  tau2.init <- con$tau2.init
  rho.init <- con$rho.init
  gamma2.init <- con$gamma2.init
  phi.init <- con$phi.init
  con$sigma2.init <- log(sigma2.init)
  if (con$cholesky[1]) {
    if (struct[1] == "UNR") {
      G <- .con.vcov.UNR(tau2.init, rho.init)
    }
    else {
      G <- .con.vcov.UN(tau2.init, rho.init)
    }
    G <- try(chol(G), silent = TRUE)
    if (inherits(G, "try-error") || anyNA(G)) 
      stop(("Cannot take Choleski decomposition of initial 'G' matrix."))
    if (struct[1] == "UNR") {
      con$tau2.init <- log(tau2.init)
    }
    else {
      con$tau2.init <- diag(G)
      con$rho.init <- G[lower.tri(G)]
    }
    if (length(con$rho.init) == 0L) 
      con$rho.init <- 0
  }
  else {
    con$tau2.init <- log(tau2.init)
    if (struct[1] == "CAR") 
      con$rho.init <- qlogis(rho.init)
    if (is.element(struct[1], c("SPEXP", "SPGAU", "SPLIN", 
                                "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD"))) 
      con$rho.init <- log(rho.init)
    if (!is.element(struct[1], c("CAR", "SPEXP", "SPGAU", 
                                 "SPLIN", "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD"))) 
      con$rho.init <- atanh(rho.init)
  }
  if (con$cholesky[2]) {
    H <- .con.vcov.UN(gamma2.init, phi.init)
    H <- try(chol(H), silent = TRUE)
    if (inherits(H, "try-error") || anyNA(H)) 
      stop(("Cannot take Choleski decomposition of initial 'H' matrix."))
    con$gamma2.init <- diag(H)
    con$phi.init <- H[lower.tri(H)]
    if (length(con$phi.init) == 0L) 
      con$phi.init <- 0
  }
  else {
    con$gamma2.init <- log(gamma2.init)
    if (struct[2] == "CAR") 
      con$phi.init <- qlogis(phi.init)
    if (is.element(struct[2], c("SPEXP", "SPGAU", "SPLIN", 
                                "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD"))) 
      con$phi.init <- log(phi.init)
    if (!is.element(struct[2], c("CAR", "SPEXP", "SPGAU", 
                                 "SPLIN", "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD"))) 
      con$phi.init <- atanh(phi.init)
  }
  optimizer <- match.arg(con$optimizer, c("optim", "nlminb", 
                                          "uobyqa", "newuoa", "bobyqa", "nloptr", "nlm", "hjk", 
                                          "nmk", "mads", "ucminf", "lbfgsb3c", "subplex", "BBoptim", 
                                          "optimParallel", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", 
                                          "SANN", "Brent", "Rcgmin", "Rvmmin"))
  optmethod <- match.arg(con$optmethod, c("Nelder-Mead", "BFGS", 
                                          "CG", "L-BFGS-B", "SANN", "Brent"))
  if (optimizer %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", 
                       "SANN", "Brent")) {
    optmethod <- optimizer
    optimizer <- "optim"
  }
  nearpd <- con$nearpd
  cholesky <- con$cholesky
  parallel <- con$parallel
  cl <- con$cl
  ncpus <- con$ncpus
  optcontrol <- control[is.na(con.pos)]
  if (length(optcontrol) == 0L) 
    optcontrol <- list()
  if (ncpus > 1L) 
    optimizer <- "optimParallel"
  reml <- ifelse(method == "REML", TRUE, FALSE)
  con$hesspack <- match.arg(con$hesspack, c("numDeriv", "pracma", 
                                            "calculus"))
  if ((isTRUE(cvvc) || cvvc %in% c("varcor", "varcov", "transf") || 
       optbeta) && !requireNamespace(con$hesspack, quietly = TRUE)) 
    stop((paste0("Please install the '", con$hesspack, 
                            "' package to compute the Hessian.")))
  if (withS) {
    if (length(con$sigma2.init) != sigma2s) 
      stop((paste0("Length of the 'sigma2.init' argument (", 
                              length(con$sigma2.init), ") does not match the actual number of variance components (", 
                              sigma2s, ").")))
  }
  else {
    con$sigma2.init <- 0
  }
  if (withG) {
    if (length(con$tau2.init) != tau2s) 
      stop((paste0("Length of the 'tau2.init' argument (", 
                              length(con$tau2.init), ") does not match the actual number of variance components (", 
                              tau2s, ").")))
  }
  else {
    con$tau2.init <- 0
  }
  if (withG) {
    if (length(con$rho.init) != rhos) 
      stop((paste0("Length of the 'rho.init' argument (", 
                              length(con$rho.init), ") does not match the actual number of correlations (", 
                              rhos, ").")))
  }
  else {
    con$rho.init <- 0
  }
  if (withH) {
    if (length(con$gamma2.init) != gamma2s) 
      stop((paste0("Length of the 'gamma2.init' argument (", 
                              length(con$gamma2.init), ") does not match the actual number of variance components (", 
                              gamma2s, ").")))
  }
  else {
    con$gamma2.init <- 0
  }
  if (withH) {
    if (length(con$phi.init) != phis) 
      stop((paste0("Length of the 'phi.init' argument (", 
                              length(con$phi.init), ") does not match the actual number of correlations (", 
                              phis, ").")))
  }
  else {
    con$phi.init <- 0
  }
  if (withZslope) {
    if (length(con$tau2B.init) != tau2Bs) 
      stop((paste0("Length of the 'tau2B.init' argument (", 
                              length(con$tau2B.init), ") does not match the actual number of variance components (", 
                              tau2Bs, ").")))
  }
  else {
    con$tau2B.init <- 0
  }
  if (withS) {
    sigma2.fix <- !is.na(sigma2)
  }
  else {
    sigma2.fix <- NA
  }
  if (withG) {
    tau2.fix <- !is.na(tau2)
    rho.fix <- !is.na(rho)
  }
  else {
    tau2.fix <- NA
    rho.fix <- NA
  }
  if (withH) {
    gamma2.fix <- !is.na(gamma2)
    phi.fix <- !is.na(phi)
  }
  else {
    gamma2.fix <- NA
    phi.fix <- NA
  }
  if (withZslope) {
    tau2B.fix <- !is.na(tau2B)
  }
  else {
    tau2B.fix <- NA
  }
  vc.fix <- list(sigma2 = sigma2.fix, tau2 = tau2.fix, rho = rho.fix,tau2B = tau2B.fix, 
                 gamma2 = gamma2.fix, phi = phi.fix)
  tmp <- .chkopt(optimizer, optcontrol)
  optimizer <- tmp$optimizer
  optcontrol <- tmp$optcontrol
  par.arg <- tmp$par.arg
  ctrl.arg <- tmp$ctrl.arg
  if (optimizer == "optimParallel::optimParallel") {
    parallel$cl <- NULL
    if (is.null(cl)) {
      ncpus <- as.integer(ncpus)
      if (ncpus < 1L) 
        stop(("Control argument 'ncpus' must be >= 1."))
      cl <- parallel::makePSOCKcluster(ncpus)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }
    else {
      if (!inherits(cl, "SOCKcluster")) 
        stop(("Specified cluster is not of class 'SOCKcluster'."))
    }
    parallel$cl <- cl
    if (is.null(parallel$forward)) 
      parallel$forward <- FALSE
    if (is.null(parallel$loginfo)) {
        parallel$loginfo <- FALSE
    }
  }
  if (optbeta || (!is.element(method, c("FE", "EE", "CE")) && 
                  !is.null(random))) {
    if (optbeta) {
      par.val <- "c(rep(0,p), con$sigma2.init, con$tau2.init, con$rho.init, con$gamma2.init, con$phi.init)"
    }
    else {
      par.val <- "c(sigma2=con$sigma2.init, tau2=con$tau2.init, rho=con$rho.init, gamma2=con$gamma2.init, phi=con$phi.init, tau2B=con$tau2B.init)"
    }
    if (anyNA(c(sigma2, tau2, rho, gamma2, phi, tau2B)) || optbeta) {
      optcall <- paste0(optimizer, "(", par.arg, "=", 
                        par.val, ", .ll.rma.mv.exch, reml=reml, ", ifelse(optimizer == 
                                                                            "optim", "method=optmethod, ", ""), "tau2B.arg=tau2B,tau2Bs=tau2Bs,withZslope=withZslope,Z.slope=Z_slope,\n Y=Y, M=V, A=NULL, X=X, k=k, pX=p,\n            D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2, g.Dmat=g.Dmat, h.Dmat=h.Dmat,\n            sigma2.arg=sigma2, tau2.arg=tau2, rho.arg=rho, gamma2.arg=gamma2, phi.arg=phi, beta.arg=beta.arg,\n            sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,\n            withS=withS, withG=withG, withH=withH, struct=struct,\n            g.levels.r=g.levels.r, h.levels.r=h.levels.r, g.values=g.values, h.values=h.values,\n            sparse=sparse, cholesky=cholesky, nearpd=nearpd, vctransf=TRUE, vccov=FALSE, vccon=vccon,\n            verbose=verbose, digits=digits, REMLf=con$REMLf,\n            dofit=FALSE, hessian=FALSE, optbeta=", 
                        optbeta, ", lambda=", lambda, ", intercept=", 
                        intercept, ctrl.arg, ")\n")
      iteration <- 0
      try(assign("iteration", iteration, envir = .metafor), 
          silent = TRUE)
        opt.res <- try(suppressWarnings(eval(str2lang(optcall))), 
                       silent = !verbose)
      if (isTRUE(args$retopt)) 
        return(opt.res)
      opt.res$par <- .chkconv(optimizer = optimizer, opt.res = opt.res, 
                              optcontrol = optcontrol, fun = "rma.mv", verbose = FALSE)
      if (p == k) {
        sigma2[is.na(sigma2)] <- 0
        tau2[is.na(tau2)] <- 0
        rho[is.na(rho)] <- 0
        gamma2[is.na(gamma2)] <- 0
        phi[is.na(phi)] <- 0
        tau2B[is.na(tau2B)] <- 0
      }
    }
    else {
      opt.res <- list(par = c(sigma2, tau2, rho, gamma2, 
                              phi, tau2B))
    }
    sigma2.arg <- sigma2
    tau2.arg <- tau2
    rho.arg <- rho
    gamma2.arg <- gamma2
    tau2B.arg <- tau2B    
    phi.arg <- phi
  }
  else {
    opt.res <- list(par = c(0, 0, 0, 0, 0,0))
  }
  fitcall <- .ll.rma.mv.exch(opt.res$par, reml = reml, Y = Y, M = V, 
                             A = A, X = X, k = k, pX = p, D.S = D.S, Z.G1 = Z.G1,
                             tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                             Z.G2 = Z.G2, Z.H1 = Z.H1, Z.H2 = Z.H2, g.Dmat = g.Dmat, 
                             h.Dmat = h.Dmat, sigma2.arg = sigma2, tau2.arg = tau2, 
                             rho.arg = rho, gamma2.arg = gamma2, phi.arg = phi, beta.arg = beta.arg, 
                             sigma2s = sigma2s, tau2s = tau2s, rhos = rhos, gamma2s = gamma2s, 
                             phis = phis, withS = withS, withG = withG, withH = withH, 
                             struct = struct, g.levels.r = g.levels.r, h.levels.r = h.levels.r, 
                             g.values = g.values, h.values = h.values, sparse = sparse, 
                             cholesky = cholesky, nearpd = nearpd, vctransf = TRUE, 
                             vccov = FALSE, vccon = vccon, verbose = FALSE, digits = digits, 
                             REMLf = con$REMLf, dofit = TRUE, optbeta = optbeta, 
                             lambda = lambda, intercept = intercept)
  beta <- as.matrix(fitcall$beta)
  vb <- matrix(NA_real_, nrow = p, ncol = p)
  hessian <- NA_real_
  vvc <- NA_real_
  if (optbeta) {
    if (con$hesspack == "numDeriv") 
      hessian <- try(numDeriv::hessian(func = .ll.rma.mv.exch, 
                                       tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                       x = opt.res$par, method.args = con$hessianCtrl, 
                                       reml = reml, Y = Y, M = V, A = A, X = X, k = k, 
                                       pX = p, D.S = D.S, Z.G1 = Z.G1, Z.G2 = Z.G2, 
                                       Z.H1 = Z.H1, Z.H2 = Z.H2, g.Dmat = g.Dmat, h.Dmat = h.Dmat, 
                                       sigma2.arg = sigma2, tau2.arg = tau2, rho.arg = rho, 
                                       gamma2.arg = gamma2, phi.arg = phi, beta.arg = beta.arg, 
                                       sigma2s = sigma2s, tau2s = tau2s, rhos = rhos, 
                                       gamma2s = gamma2s, phis = phis, withS = withS, 
                                       withG = withG, withH = withH, struct = struct, 
                                       g.levels.r = g.levels.r, h.levels.r = h.levels.r, 
                                       g.values = g.values, h.values = h.values, sparse = sparse, 
                                       cholesky = cholesky, nearpd = nearpd, vctransf = TRUE, 
                                       vccov = FALSE, vccon = vccon, verbose = verbose, 
                                       digits = digits, REMLf = con$REMLf, dofit = FALSE, 
                                       hessian = TRUE, optbeta = optbeta, lambda = lambda, 
                                       intercept = intercept), silent = !verbose)
    if (con$hesspack == "pracma") 
      hessian <- try(pracma::hessian(f = .ll.rma.mv.exch, x0 = opt.res$par, 
                                     tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                     reml = reml, Y = Y, M = V, A = A, X = X, k = k, 
                                     pX = p, D.S = D.S, Z.G1 = Z.G1, Z.G2 = Z.G2, 
                                     Z.H1 = Z.H1, Z.H2 = Z.H2, g.Dmat = g.Dmat, h.Dmat = h.Dmat, 
                                     sigma2.arg = sigma2, tau2.arg = tau2, rho.arg = rho, 
                                     gamma2.arg = gamma2, phi.arg = phi, beta.arg = beta.arg, 
                                     sigma2s = sigma2s, tau2s = tau2s, rhos = rhos, 
                                     gamma2s = gamma2s, phis = phis, withS = withS, 
                                     withG = withG, withH = withH, struct = struct, 
                                     g.levels.r = g.levels.r, h.levels.r = h.levels.r, 
                                     g.values = g.values, h.values = h.values, sparse = sparse, 
                                     cholesky = cholesky, nearpd = nearpd, vctransf = TRUE, 
                                     vccov = FALSE, vccon = vccon, verbose = verbose, 
                                     digits = digits, REMLf = con$REMLf, dofit = FALSE, 
                                     hessian = TRUE, optbeta = optbeta, lambda = lambda, 
                                     intercept = intercept), silent = !verbose)
    if (con$hesspack == "calculus") 
      hessian <- try(calculus::hessian(f = .ll.rma.mv.exch, 
                                       var = opt.res$par, params = list(reml = reml,
                                                                        tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                                                        Y = Y, M = V, A = A, X = X, k = k, pX = p, 
                                                                        D.S = D.S, Z.G1 = Z.G1, Z.G2 = Z.G2, Z.H1 = Z.H1, 
                                                                        Z.H2 = Z.H2, g.Dmat = g.Dmat, h.Dmat = h.Dmat, 
                                                                        sigma2.arg = sigma2, tau2.arg = tau2, rho.arg = rho, 
                                                                        gamma2.arg = gamma2, phi.arg = phi, beta.arg = beta.arg, 
                                                                        sigma2s = sigma2s, tau2s = tau2s, rhos = rhos, 
                                                                        gamma2s = gamma2s, phis = phis, withS = withS, 
                                                                        withG = withG, withH = withH, struct = struct, 
                                                                        g.levels.r = g.levels.r, h.levels.r = h.levels.r, 
                                                                        g.values = g.values, h.values = h.values, 
                                                                        sparse = sparse, cholesky = cholesky, nearpd = nearpd, 
                                                                        vctransf = TRUE, vccov = FALSE, vccon = vccon, 
                                                                        verbose = verbose, digits = digits, REMLf = con$REMLf, 
                                                                        dofit = FALSE, hessian = TRUE, optbeta = optbeta, 
                                                                        lambda = lambda, intercept = intercept)), 
                     silent = !verbose)
    if (inherits(hessian, "try-error")) {
      warning(("Error when trying to compute the Hessian."), 
              call. = FALSE)
      hessian <- NA_real_
    }
    else {
      colnames(hessian) <- rep("", ncol(hessian))
      if (int.incl) {
        colnames(hessian)[1:p] <- paste0("beta", 0:(p - 
                                                      1))
      }
      else {
        colnames(hessian)[1:p] <- paste0("beta", 1:p)
      }
      rownames(hessian) <- colnames(hessian)
      hest <- !apply(hessian, 1, function(x) all(abs(x) <= 
                                                   con$hesstol))
      hessian <- hessian[hest, hest, drop = FALSE]
      if (any(hest)) {
        vvc <- try(suppressWarnings(chol2inv(chol(hessian))), 
                   silent = TRUE)
        if (inherits(vvc, "try-error") || anyNA(vvc) || 
            any(is.infinite(vvc))) 
          warning(("Error when trying to invert the Hessian."), 
                  call. = FALSE)
        sel <- grep("beta", colnames(hessian), fixed = TRUE)
        vb[hest[1:p], hest[1:p]] <- vvc[sel, sel, drop = FALSE]
      }
      else {
        vb <- matrix(NA_real_, nrow = p, ncol = p)
      }
    }
  }
  else {
    vb <- as.matrix(fitcall$vb)
    vb[!beta.est, ] <- NA_real_
    vb[, !beta.est] <- NA_real_
  }
  if (withS) 
    sigma2 <- fitcall$sigma2
  if (withZslope)
    tau2B <- fitcall$tau2B
  if (withG) {
    G <- as.matrix(fitcall$G)
    if (is.element(struct[1], c("SPEXP", "SPGAU", "SPLIN", 
                                "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD"))) 
      colnames(G) <- rownames(G) <- seq_len(nrow(G))
    if (is.element(struct[1], c("CS", "HCS", "UN", "UNR", 
                                "AR", "HAR", "CAR", "ID", "DIAG"))) 
      colnames(G) <- rownames(G) <- g.levels.f[[1]]
    if (is.element(struct[1], c("GEN", "GDIAG"))) 
      colnames(G) <- rownames(G) <- g.names[-length(g.names)]
    tau2 <- fitcall$tau2
    rho <- fitcall$rho
    cov1 <- G[lower.tri(G)]
  }
  else {
    cov1 <- 0
  }
  if (withH) {
    H <- as.matrix(fitcall$H)
    if (is.element(struct[2], c("SPEXP", "SPGAU", "SPLIN", 
                                "SPRAT", "SPSPH", "PHYBM", "PHYPL", "PHYPD"))) 
      colnames(H) <- rownames(H) <- seq_len(nrow(H))
    if (is.element(struct[2], c("CS", "HCS", "UN", "UNR", 
                                "AR", "HAR", "CAR", "ID", "DIAG"))) 
      colnames(H) <- rownames(H) <- h.levels.f[[1]]
    if (is.element(struct[2], c("GEN", "GDIAG"))) 
      colnames(H) <- rownames(H) <- h.names[-length(h.names)]
    gamma2 <- fitcall$gamma2
    phi <- fitcall$phi
    cov2 <- H[lower.tri(H)]
  }
  else {
    cov2 <- 0
  }
  M <- fitcall$M
  if (!is.null(dimnames(M))) 
    M <- unname(M)
  if (is.element(test, c("knha", "adhoc", "t"))) {
    ddf <- .ddf.calc(dfs, X = X, k = k, p = p, mf.s = mf.s, 
                     mf.g = mf.g, mf.h = mf.h)
  }
  else {
    ddf <- rep(NA_integer_, p)
  }
  s2w <- 1
  if (is.element(test, c("knha", "adhoc"))) {
    knha.rma.mv.warn <- .getfromenv("knha.rma.mv.warn", 
                                    default = TRUE)
    if (knha.rma.mv.warn) {
      warning(("Use of the Knapp and Hartung method for 'rma.mv()' models is experimental.\nNote: This warning is only issued once per session (ignore at your peril)."), 
              call. = FALSE)
      try(assign("knha.rma.mv.warn", FALSE, envir = .metafor), 
          silent = TRUE)
    }
    RSS <- try(as.vector(t(Y - X %*% beta) %*% chol2inv(chol(M)) %*% 
                           (Y - X %*% beta)), silent = TRUE)
    if (inherits(RSS, "try-error")) 
      stop((paste0("Failure when trying to compute adjustment factor for Knapp and Hartung method.")))
    if (RSS <= .Machine$double.eps) {
      s2w <- 0
    }
    else {
      s2w <- as.vector(RSS/(k - p))
    }
  }
  if (test == "adhoc") 
    s2w[s2w < 1] <- 1
  vb <- s2w * vb
  QM <- try(as.vector(t(beta)[btt] %*% chol2inv(chol(vb[btt, 
                                                        btt])) %*% beta[btt]), silent = TRUE)
  if (inherits(QM, "try-error")) 
    QM <- NA_real_
  if (isTRUE(args$abbrev)) {
    tmp <- colnames(X)
    tmp <- gsub("relevel(factor(", "", tmp, fixed = TRUE)
    tmp <- gsub("\\), ref = \"[[:alnum:]]*\")", "", tmp)
    tmp <- gsub("poly(", "", tmp, fixed = TRUE)
    tmp <- gsub(", degree = [[:digit:]], raw = TRUE)", "^", 
                tmp)
    tmp <- gsub(", degree = [[:digit:]], raw = T)", "^", 
                tmp)
    tmp <- gsub(", degree = [[:digit:]])", "^", tmp)
    tmp <- gsub("rcs\\([[:alnum:]]*, [[:digit:]]\\)", "", 
                tmp)
    tmp <- gsub("factor(", "", tmp, fixed = TRUE)
    tmp <- gsub("I(", "", tmp, fixed = TRUE)
    tmp <- gsub(")", "", tmp, fixed = TRUE)
    colnames(X) <- tmp
  }
  rownames(beta) <- rownames(vb) <- colnames(vb) <- colnames(X.f) <- colnames(X)
  se <- sqrt(diag(vb))
  names(se) <- NULL
  zval <- c(beta/se)
  if (is.element(test, c("knha", "adhoc", "t"))) {
    QM <- QM/m
    QMdf <- c(m, min(ddf[btt]))
    QMp <- if (QMdf[2] > 0) 
      pf(QM, df1 = QMdf[1], df2 = QMdf[2], lower.tail = FALSE)
    else NA_real_
    pval <- sapply(seq_along(ddf), function(j) if (ddf[j] > 
                                                   0) 
      2 * pt(abs(zval[j]), df = ddf[j], lower.tail = FALSE)
      else NA_real_)
    crit <- sapply(seq_along(ddf), function(j) if (ddf[j] > 
                                                   0) 
      qt(level/2, df = ddf[j], lower.tail = FALSE)
      else NA_real_)
  }
  else {
    QMdf <- c(m, NA_integer_)
    QMp <- pchisq(QM, df = QMdf[1], lower.tail = FALSE)
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    crit <- qnorm(level/2, lower.tail = FALSE)
  }
  ci.lb <- c(beta - crit * se)
  ci.ub <- c(beta + crit * se)
  QEdf <- k - p
  if (QEdf > 0L) {
    QEp <- pchisq(QE, df = QEdf, lower.tail = FALSE)
  }
  else {
    QE <- 0
    QEp <- 1
  }
  ll.QE <- -1/2 * (k) * log(2 * base::pi) - 1/2 * determinant(V, 
                                                              logarithm = TRUE)$modulus
  if (!optbeta && (!is.element(method, c("FE", "EE", "CE")) && 
                   !is.null(random)) && (isTRUE(cvvc) || cvvc %in% c("varcor", 
                                                                     "varcov", "transf"))) {
    if (cvvc == "varcov" && (any(sigma2.fix, na.rm = TRUE) || 
                             any(tau2.fix, na.rm = TRUE) || any(rho.fix, na.rm = TRUE) || 
                             any(gamma2.fix, na.rm = TRUE) || any(phi.fix, na.rm = TRUE))) {
      warning(("Cannot use cvvc='varcov' when one or more components are fixed. Setting cvvc='varcor'."), 
              call. = FALSE)
      cvvc <- "varcor"
    }
    if (cvvc == "varcov" && any(!is.element(struct, c("UN", 
                                                      "GEN")))) {
      warning(("Cannot use cvvc='varcov' for the specified structure(s). Setting cvvc='varcor'."), 
              call. = FALSE)
      cvvc <- "varcor"
    }
    if (cvvc == "varcov") {
      if (con$hesspack == "numDeriv") 
        hessian <- try(numDeriv::hessian(func = .ll.rma.mv.exch, 
                                         x = c(sigma2, tau2, cov1, gamma2, cov2), method.args = con$hessianCtrl,
                                         tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                         reml = reml, Y = Y, M = V, A = NULL, X = X, 
                                         k = k, pX = p, D.S = D.S, Z.G1 = Z.G1, Z.G2 = Z.G2, 
                                         Z.H1 = Z.H1, Z.H2 = Z.H2, g.Dmat = g.Dmat, 
                                         h.Dmat = h.Dmat, sigma2.arg = sigma2.arg, 
                                         tau2.arg = tau2.arg, rho.arg = rho.arg, gamma2.arg = gamma2.arg, 
                                         phi.arg = phi.arg, beta.arg = beta.arg, sigma2s = sigma2s, 
                                         tau2s = tau2s, rhos = rhos, gamma2s = gamma2s, 
                                         phis = phis, withS = withS, withG = withG, 
                                         withH = withH, struct = struct, g.levels.r = g.levels.r, 
                                         h.levels.r = h.levels.r, g.values = g.values, 
                                         h.values = h.values, sparse = sparse, cholesky = c(FALSE, 
                                                                                            FALSE), nearpd = nearpd, vctransf = FALSE, 
                                         vccov = TRUE, vccon = vccon, verbose = verbose, 
                                         digits = digits, REMLf = con$REMLf, hessian = TRUE), 
                       silent = TRUE)
      if (con$hesspack == "pracma") 
        hessian <- try(pracma::hessian(f = .ll.rma.mv.exch, 
                                       x0 = c(sigma2, tau2, cov1, gamma2, cov2),
                                       tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                       reml = reml, Y = Y, M = V, A = NULL, X = X, 
                                       k = k, pX = p, D.S = D.S, Z.G1 = Z.G1, Z.G2 = Z.G2, 
                                       Z.H1 = Z.H1, Z.H2 = Z.H2, g.Dmat = g.Dmat, 
                                       h.Dmat = h.Dmat, sigma2.arg = sigma2.arg, 
                                       tau2.arg = tau2.arg, rho.arg = rho.arg, gamma2.arg = gamma2.arg, 
                                       phi.arg = phi.arg, beta.arg = beta.arg, sigma2s = sigma2s, 
                                       tau2s = tau2s, rhos = rhos, gamma2s = gamma2s, 
                                       phis = phis, withS = withS, withG = withG, 
                                       withH = withH, struct = struct, g.levels.r = g.levels.r, 
                                       h.levels.r = h.levels.r, g.values = g.values, 
                                       h.values = h.values, sparse = sparse, cholesky = c(FALSE, 
                                                                                          FALSE), nearpd = nearpd, vctransf = FALSE, 
                                       vccov = TRUE, vccon = vccon, verbose = verbose, 
                                       digits = digits, REMLf = con$REMLf, hessian = TRUE), 
                       silent = TRUE)
      if (con$hesspack == "calculus") 
        hessian <- try(calculus::hessian(f = .ll.rma.mv.exch, 
                                         var = c(sigma2, tau2, cov1, gamma2, cov2), 
                                         params = list(reml = reml, Y = Y, M = V, A = NULL, 
                                                       tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                                       X = X, k = k, pX = p, D.S = D.S, Z.G1 = Z.G1, 
                                                       Z.G2 = Z.G2, Z.H1 = Z.H1, Z.H2 = Z.H2, g.Dmat = g.Dmat, 
                                                       h.Dmat = h.Dmat, sigma2.arg = sigma2.arg, 
                                                       tau2.arg = tau2.arg, rho.arg = rho.arg, 
                                                       gamma2.arg = gamma2.arg, phi.arg = phi.arg, 
                                                       beta.arg = beta.arg, sigma2s = sigma2s, 
                                                       tau2s = tau2s, rhos = rhos, gamma2s = gamma2s, 
                                                       phis = phis, withS = withS, withG = withG, 
                                                       withH = withH, struct = struct, g.levels.r = g.levels.r, 
                                                       h.levels.r = h.levels.r, g.values = g.values, 
                                                       h.values = h.values, sparse = sparse, cholesky = c(FALSE, 
                                                                                                          FALSE), nearpd = nearpd, vctransf = FALSE, 
                                                       vccov = TRUE, vccon = vccon, verbose = verbose, 
                                                       digits = digits, REMLf = con$REMLf, hessian = TRUE)), 
                       silent = TRUE)
    }
    else {
      if (con$hesspack == "numDeriv") 
        hessian <- try(numDeriv::hessian(func = .ll.rma.mv.exch, 
                                         x = if (cvvc == "transf") 
                                           opt.res$par
                                         else c(sigma2, tau2, rho, gamma2, phi), method.args = con$hessianCtrl, 
                                         tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                         reml = reml, Y = Y, M = V, A = NULL, X = X, 
                                         k = k, pX = p, D.S = D.S, Z.G1 = Z.G1, Z.G2 = Z.G2, 
                                         Z.H1 = Z.H1, Z.H2 = Z.H2, g.Dmat = g.Dmat, 
                                         h.Dmat = h.Dmat, sigma2.arg = sigma2.arg, 
                                         tau2.arg = tau2.arg, rho.arg = rho.arg, gamma2.arg = gamma2.arg, 
                                         phi.arg = phi.arg, beta.arg = beta.arg, sigma2s = sigma2s, 
                                         tau2s = tau2s, rhos = rhos, gamma2s = gamma2s, 
                                         phis = phis, withS = withS, withG = withG, 
                                         withH = withH, struct = struct, g.levels.r = g.levels.r, 
                                         h.levels.r = h.levels.r, g.values = g.values, 
                                         h.values = h.values, sparse = sparse, cholesky = ifelse(c(cvvc == 
                                                                                                     "transf", cvvc == "transf") & cholesky, 
                                                                                                 TRUE, FALSE), nearpd = nearpd, vctransf = cvvc == 
                                           "transf", vccov = FALSE, vccon = vccon, 
                                         verbose = verbose, digits = digits, REMLf = con$REMLf, 
                                         hessian = TRUE), silent = TRUE)
      if (con$hesspack == "pracma") 
        hessian <- try(pracma::hessian(f = .ll.rma.mv.exch, 
                                       x0 = if (cvvc == "transf") 
                                         opt.res$par
                                       else c(sigma2, tau2, rho, gamma2, phi), reml = reml, 
                                       tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                       Y = Y, M = V, A = NULL, X = X, k = k, pX = p, 
                                       D.S = D.S, Z.G1 = Z.G1, Z.G2 = Z.G2, Z.H1 = Z.H1, 
                                       Z.H2 = Z.H2, g.Dmat = g.Dmat, h.Dmat = h.Dmat, 
                                       sigma2.arg = sigma2.arg, tau2.arg = tau2.arg, 
                                       rho.arg = rho.arg, gamma2.arg = gamma2.arg, 
                                       phi.arg = phi.arg, beta.arg = beta.arg, sigma2s = sigma2s, 
                                       tau2s = tau2s, rhos = rhos, gamma2s = gamma2s, 
                                       phis = phis, withS = withS, withG = withG, 
                                       withH = withH, struct = struct, g.levels.r = g.levels.r, 
                                       h.levels.r = h.levels.r, g.values = g.values, 
                                       h.values = h.values, sparse = sparse, cholesky = ifelse(c(cvvc == 
                                                                                                   "transf", cvvc == "transf") & cholesky, 
                                                                                               TRUE, FALSE), nearpd = nearpd, vctransf = cvvc == 
                                         "transf", vccov = FALSE, vccon = vccon, 
                                       verbose = verbose, digits = digits, REMLf = con$REMLf, 
                                       hessian = TRUE), silent = TRUE)
      if (con$hesspack == "calculus") 
        hessian <- try(calculus::hessian(f = .ll.rma.mv.exch, 
                                         var = if (cvvc == "transf") 
                                           opt.res$par
                                         else c(sigma2, tau2, rho, gamma2, phi), params = list(reml = reml, 
                                                                                               tau2B.arg=tau2B, withZslope=withZslope,tau2Bs=tau2Bs,Z.slope=Z_slope,
                                                                                               Y = Y, M = V, A = NULL, X = X, k = k, pX = p, 
                                                                                               D.S = D.S, Z.G1 = Z.G1, Z.G2 = Z.G2, Z.H1 = Z.H1, 
                                                                                               Z.H2 = Z.H2, g.Dmat = g.Dmat, h.Dmat = h.Dmat, 
                                                                                               sigma2.arg = sigma2.arg, tau2.arg = tau2.arg, 
                                                                                               rho.arg = rho.arg, gamma2.arg = gamma2.arg, 
                                                                                               phi.arg = phi.arg, beta.arg = beta.arg, 
                                                                                               sigma2s = sigma2s, tau2s = tau2s, rhos = rhos, 
                                                                                               gamma2s = gamma2s, phis = phis, withS = withS, 
                                                                                               withG = withG, withH = withH, struct = struct, 
                                                                                               g.levels.r = g.levels.r, h.levels.r = h.levels.r, 
                                                                                               g.values = g.values, h.values = h.values, 
                                                                                               sparse = sparse, cholesky = ifelse(c(cvvc == 
                                                                                                                                      "transf", cvvc == "transf") & cholesky, 
                                                                                                                                  TRUE, FALSE), nearpd = nearpd, vctransf = cvvc == 
                                                                                                 "transf", vccov = FALSE, vccon = vccon, 
                                                                                               verbose = verbose, digits = digits, REMLf = con$REMLf, 
                                                                                               hessian = TRUE)), silent = TRUE)
    }
    if (inherits(hessian, "try-error")) {
      warning(("Error when trying to compute the Hessian."), 
              call. = FALSE)
      hessian <- NA_real_
    }
    else {
      colnames(hessian) <- seq_len(ncol(hessian))
      if (sigma2s == 1) {
        colnames(hessian)[1] <- "sigma^2"
      }
      else {
        colnames(hessian)[1:sigma2s] <- paste0("sigma^2.", 
                                               seq_len(sigma2s))
      }
      if (tau2s == 1) {
        colnames(hessian)[sigma2s + 1] <- "tau^2"
      }
      else {
        colnames(hessian)[(sigma2s + 1):(sigma2s + tau2s)] <- paste0("tau^2.", 
                                                                     seq_len(tau2s))
      }
      term <- ifelse(cvvc == "varcov", ifelse(withH, "cov1", 
                                              "cov"), "rho")
      if (rhos == 1) {
        colnames(hessian)[sigma2s + tau2s + 1] <- term
      }
      else {
        colnames(hessian)[(sigma2s + tau2s + 1):(sigma2s + 
                                                   tau2s + rhos)] <- paste0(term, ".", outer(seq_len(g.nlevels.f[1]), 
                                                                                             seq_len(g.nlevels.f[1]), paste, sep = ".")[lower.tri(matrix(NA, 
                                                                                                                                                         nrow = g.nlevels.f, ncol = g.nlevels.f))])
      }
      if (gamma2s == 1) {
        colnames(hessian)[sigma2s + tau2s + rhos + 1] <- "gamma^2"
      }
      else {
        colnames(hessian)[(sigma2s + tau2s + rhos + 
                             1):(sigma2s + tau2s + rhos + gamma2s)] <- paste0("gamma^2.", 
                                                                              seq_len(gamma2s))
      }
      term <- ifelse(cvvc == "varcov", "cov2", "phi")
      if (phis == 1) {
        colnames(hessian)[sigma2s + tau2s + rhos + gamma2s + 
                            1] <- term
      }
      else {
        colnames(hessian)[(sigma2s + tau2s + rhos + 
                             gamma2s + 1):(sigma2s + tau2s + rhos + gamma2s + 
                                             phis)] <- paste0(term, ".", outer(seq_len(h.nlevels.f[1]), 
                                                                               seq_len(h.nlevels.f[1]), paste, sep = ".")[lower.tri(matrix(NA, 
                                                                                                                                           nrow = h.nlevels.f, ncol = h.nlevels.f))])
      }
      rownames(hessian) <- colnames(hessian)
      if (withS && withG && !withH) 
        hessian <- hessian[1:(nrow(hessian) - 2), 1:(ncol(hessian) - 
                                                       2), drop = FALSE]
      if (withS && !withG && !withH) 
        hessian <- hessian[1:(nrow(hessian) - 4), 1:(ncol(hessian) - 
                                                       4), drop = FALSE]
      if (!withS && withG && withH) 
        hessian <- hessian[2:nrow(hessian), 2:ncol(hessian), 
                           drop = FALSE]
      if (!withS && withG && !withH) 
        hessian <- hessian[2:(nrow(hessian) - 2), 2:(ncol(hessian) - 
                                                       2), drop = FALSE]
      if (!withS && !withG && !withH) 
        hessian <- NA_real_
      if (cvvc == "varcov" && withG) {
        posG <- matrix(NA_real_, nrow = tau2s, ncol = tau2s)
        diag(posG) <- 1:tau2s
        posG[lower.tri(posG)] <- (tau2s + 1):(tau2s * 
                                                (tau2s + 1)/2)
        posG <- posG[lower.tri(posG, diag = TRUE)]
        if (withS) {
          pos <- c(1:sigma2s, sigma2s + posG)
        }
        else {
          pos <- posG
        }
        if (withH) {
          posH <- matrix(NA_real_, nrow = gamma2s, ncol = gamma2s)
          diag(posH) <- 1:gamma2s
          posH[lower.tri(posH)] <- (gamma2s + 1):(gamma2s * 
                                                    (gamma2s + 1)/2)
          posH <- posH[lower.tri(posH, diag = TRUE)]
          pos <- c(pos, max(pos) + posH)
        }
        hessian <- hessian[pos, pos]
      }
      hest <- !apply(hessian, 1, function(x) all(abs(x) <= 
                                                   con$hesstol))
      hessian <- hessian[hest, hest, drop = FALSE]
      vvc <- try(suppressWarnings(chol2inv(chol(hessian))), 
                 silent = TRUE)
      if (inherits(vvc, "try-error") || anyNA(vvc) || 
          any(is.infinite(vvc))) {
        warning(("Error when trying to invert the Hessian."), 
                call. = FALSE)
        vvc <- NA_real_
      }
      else {
        dimnames(vvc) <- dimnames(hessian)
      }
    }
  }
  p <- sum(beta.est)
  if (is.null(vccon)) {
    parms <- p + ifelse(withS, sum(ifelse(sigma2.fix, 0, 
                                          1)), 0) + ifelse(withG, sum(ifelse(tau2.fix, 0, 
                                                                             1)), 0) + ifelse(withG, sum(ifelse(rho.fix, 0, 1)), 
                                                                                              0) + ifelse(withH, sum(ifelse(gamma2.fix, 0, 1)), 
                                                                                                          0) + ifelse(withH, sum(ifelse(phi.fix, 0, 1)), 0)
  }
  else {
    parms <- p + ifelse(withS && !is.null(vccon$sigma2), 
                        length(unique(vccon$sigma2)) - sum(sigma2.fix), 
                        0) + ifelse(withG && !is.null(vccon$tau2), length(unique(vccon$tau2)) - 
                                      sum(tau2.fix), 0) + ifelse(withG && !is.null(vccon$rho), 
                                                                 length(unique(vccon$rho)) - sum(rho.fix), 0) + ifelse(withH && 
                                                                                                                         !is.null(vccon$gamma2), length(unique(vccon$gamma2)) - 
                                                                                                                         sum(gamma2.fix), 0) + ifelse(withH && !is.null(vccon$phi), 
                                                                                                                                                      length(unique(vccon$phi)) - sum(phi.fix), 0)
  }
  ll.ML <- fitcall$llvals[1]
  ll.REML <- fitcall$llvals[2]
  if (allvipos) {
    dev.ML <- -2 * (ll.ML - ll.QE)
  }
  else {
    dev.ML <- -2 * ll.ML
  }
  AIC.ML <- -2 * ll.ML + 2 * parms
  BIC.ML <- -2 * ll.ML + parms * log(k)
  AICc.ML <- -2 * ll.ML + 2 * parms * max(k, parms + 2)/(max(k, 
                                                             parms + 2) - parms - 1)
  dev.REML <- -2 * (ll.REML - 0)
  AIC.REML <- -2 * ll.REML + 2 * parms
  BIC.REML <- -2 * ll.REML + parms * log(k - p)
  AICc.REML <- -2 * ll.REML + 2 * parms * max(k - p, parms + 
                                                2)/(max(k - p, parms + 2) - parms - 1)
  fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, 
                        ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol = 2, 
                      byrow = FALSE)
  dimnames(fit.stats) <- list(c("ll", "dev", "AIC", "BIC", 
                                "AICc"), c("ML", "REML"))
  fit.stats <- data.frame(fit.stats)
  replfun <- function(x) {
    if (grepl("interaction(", x, fixed = TRUE) || grepl("paste(", 
                                                        x, fixed = TRUE) || grepl("paste0(", x, fixed = TRUE)) {
      x <- gsub("interaction\\((.*)\\)", "(\\1)", x)
      x <- gsub("paste[0]?\\((.*)\\)", "(\\1)", x)
      x <- gsub(",", ":", x, fixed = TRUE)
      x <- gsub(" ", "", x, fixed = TRUE)
      x <- gsub("^\\((.*)\\)$", "\\1", x)
    }
    return(x)
  }
  s.names <- sapply(s.names, replfun)
  g.names <- sapply(g.names, replfun)
  h.names <- sapply(h.names, replfun)
  p.eff <- p
  k.eff <- k
  weighted <- TRUE
  if (is.null(args$outlist) || args$outlist == "nodata") {
    res <- list(b = beta, beta = beta, se = se, zval = zval, 
                pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, vb = vb, 
                sigma2 = sigma2, tau2 = tau2, rho = rho, gamma2 = gamma2, 
                phi = phi, QE = QE, QEdf = QEdf, QEp = QEp, QM = QM, 
                QMdf = QMdf, QMp = QMp, k = k, k.f = k.f, k.eff = k.eff, 
                k.all = k.all, p = p, p.eff = p.eff, parms = parms, 
                int.only = int.only, int.incl = int.incl, intercept = intercept, 
                allvipos = allvipos, coef.na = coef.na, yi = yi, 
                vi = vi, V = V, W = A, X = X, yi.f = yi.f, vi.f = vi.f, 
                V.f = V.f, X.f = X.f, W.f = W.f, ni = ni, ni.f = ni.f, 
                M = M, G = G, H = H, hessian = hessian, vvc = vvc, 
                vccon = vccon, chksumyi = digest::digest(as.vector(yi)), 
                chksumV = digest::digest(as.matrix(V)), chksumX = digest::digest(X), 
                ids = ids, not.na = not.na, subset = subset, slab = slab, 
                slab.null = slab.null, measure = measure, method = method, 
                weighted = weighted, optbeta = optbeta, test = test, 
                dfs = dfs, ddf = ddf, s2w = s2w, btt = btt, m = m, 
                digits = digits, level = level, sparse = sparse, 
                dist = args$dist, control = control, verbose = verbose, 
                fit.stats = fit.stats, vc.fix = vc.fix, withS = withS, 
                withG = withG, withH = withH, withR = withR, formulas = formulas, 
                sigma2s = sigma2s, tau2s = tau2s, rhos = rhos, gamma2s = gamma2s, 
                phis = phis, s.names = s.names, g.names = g.names, 
                h.names = h.names, s.levels = s.levels, s.levels.f = s.levels.f, 
                s.nlevels = s.nlevels, s.nlevels.f = s.nlevels.f, 
                g.nlevels.f = g.nlevels.f, g.nlevels = g.nlevels, 
                h.nlevels.f = h.nlevels.f, h.nlevels = h.nlevels, 
                g.levels.f = g.levels.f, g.levels.k = g.levels.k, 
                g.levels.comb.k = g.levels.comb.k, h.levels.f = h.levels.f, 
                h.levels.k = h.levels.k, h.levels.comb.k = h.levels.comb.k, 
                struct = struct, Rfix = Rfix, R = R, Rscale = Rscale, 
                mf.r = mf.r, mf.s = mf.s, mf.g = mf.g, mf.g.f = mf.g.f, 
                mf.h = mf.h, mf.h.f = mf.h.f, Z.S = Z.S, Z.G1 = Z.G1, 
                Z.G2 = Z.G2, Z.H1 = Z.H1, Z.H2 = Z.H2, formula.yi = formula.yi, 
                Z_slope=Z_slope,withZslope=withZslope, tau2B=tau2B,
                formula.mods = formula.mods, random = random, version = packageVersion("metafor"), 
                call = mf)
    if (is.null(args$outlist)) 
      res <- append(res, list(data = data), which(names(res) == 
                                                    "fit.stats"))
  }
  else {
    if (args$outlist == "minimal") {
      res <- list(b = beta, beta = beta, se = se, zval = zval, 
                  pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, vb = vb, 
                  sigma2 = sigma2, tau2 = tau2, rho = rho, gamma2 = gamma2, 
                  phi = phi, QE = QE, QEdf = QEdf, QEp = QEp, 
                  QM = QM, QMdf = QMdf, QMp = QMp, k = k, k.f = k.f, 
                  k.eff = k.eff, k.all = k.all, p = p, p.eff = p.eff, 
                  parms = parms, int.only = int.only, int.incl = int.incl, 
                  intercept = intercept, chksumyi = digest::digest(as.vector(yi)), 
                  chksumV = digest::digest(as.matrix(V)), chksumX = digest::digest(X), 
                  measure = measure, method = method, weighted = weighted, 
                  optbeta = optbeta, test = test, dfs = dfs, ddf = ddf, 
                  btt = btt, m = m, digits = digits, level = level, 
                  fit.stats = fit.stats, vc.fix = vc.fix, withS = withS, 
                  withG = withG, withH = withH, withR = withR, 
                  s.names = s.names, g.names = g.names, h.names = h.names, 
                  s.nlevels = s.nlevels, g.nlevels.f = g.nlevels.f, 
                  g.nlevels = g.nlevels, h.nlevels.f = h.nlevels.f, 
                  h.nlevels = h.nlevels, g.levels.f = g.levels.f, 
                  g.levels.k = g.levels.k, g.levels.comb.k = g.levels.comb.k, 
                  h.levels.f = h.levels.f, h.levels.k = h.levels.k, 
                  h.levels.comb.k = h.levels.comb.k, struct = struct, 
                  Rfix = Rfix,
                  Z_slope=Z_slope,withZslope=withZslope, tau2B=tau2B)
    }
    else {
      res <- eval(str2lang(paste0("list(", args$outlist, 
                                  ")")))
    }
  }
  time.end <- proc.time()
  res$time <- unname(time.end - time.start)[3]
  if (isTRUE(args$time)) 
    .print.time(res$time)
  if (isTRUE(args$time)) 
    cat("\n")
  class(res) <- c("rma.mv", "rma")
  return(res)
}



.ll.rma.mv.exch<-function (par, reml, Y, M, A, X, k, pX, D.S, Z.G1, Z.G2, Z.H1, 
                           Z.slope = NULL, tau2B.arg,tau2Bs,withZslope,   # NEW arguments
                           Z.H2, g.Dmat, h.Dmat, sigma2.arg, tau2.arg, rho.arg, gamma2.arg, 
                           phi.arg, beta.arg, sigma2s, tau2s, rhos, gamma2s, phis, 
                           withS, withG, withH, struct, g.levels.r, h.levels.r, g.values, 
                           h.values, sparse, cholesky, nearpd, vctransf, vccov, vccon, 
                           verbose, digits, REMLf, dofit = FALSE, hessian = FALSE, 
                           optbeta = FALSE, lambda = 0, intercept = TRUE) 
{
  if (optbeta) {
    beta <- par[1:pX]
    par <- par[-c(1:pX)]
  }
  if (withS) {
    vars <- par[seq_len(sigma2s)]
    if (vctransf) {
      sigma2 <- ifelse(is.na(sigma2.arg), exp(vars), sigma2.arg)
    }
    else {
      sigma2 <- ifelse(is.na(sigma2.arg), vars, sigma2.arg)
      sigma2[sigma2 < 0] <- 0
    }
    sigma2 <- ifelse(sigma2 <= .Machine$double.eps * 10, 
                     0, sigma2)
    if (!is.null(vccon) && !is.null(vccon$sigma2)) {
      for (id in unique(vccon$sigma2)) sigma2[vccon$sigma2 == 
                                                id] <- mean(sigma2[vccon$sigma2 == id])
    }
    for (j in seq_len(sigma2s)) {
      M <- M + sigma2[j] * D.S[[j]]
    }
  }
  if (withG) {
    vars <- par[(sigma2s + 1):(sigma2s + tau2s)]
    cors <- par[(sigma2s + tau2s + 1):(sigma2s + tau2s + 
                                         rhos)]
    resG <- .con.E(v = vars, r = cors, v.arg = tau2.arg,
                   r.arg = rho.arg, Z1 = Z.G1, Z2 = Z.G2, levels.r = g.levels.r,
                   values = g.values, Dmat = g.Dmat, struct = struct[1],
                   cholesky = cholesky[1], vctransf = vctransf, vccov = vccov,
                   nearpd = nearpd, sparse = sparse)
    tau2 <- resG$v
    tau2 <- ifelse(tau2 <= .Machine$double.eps * 10, 0, tau2)
    rho <- resG$r#cors
    G <- resG$E
    if (!is.null(vccon)) {
      if (!is.null(vccon$tau2)) {
        for (id in unique(vccon$tau2)) tau2[vccon$tau2 == 
                                              id] <- mean(tau2[vccon$tau2 == id])
      }
      if (!is.null(vccon$rho)) {
        for (id in unique(vccon$rho)) {
          rho[vccon$rho == id] <- mean(rho[vccon$rho == 
                                             id])
        }
      }
      resG <- .con.E(v = tau2, r = rho, v.arg = tau2.arg, 
                     r.arg = rho.arg, Z1 = Z.G1, Z2 = Z.G2, levels.r = g.levels.r, 
                     values = g.values, Dmat = g.Dmat, struct = struct[1], 
                     cholesky = FALSE, vctransf = FALSE, vccov = vccov, 
                     nearpd = nearpd, sparse = sparse)
      tau2 <- resG$v
      rho <- resG$r
      G <- resG$E
    }
    M <- M + (Z.G1 %*% G %*% t(Z.G1)) * tcrossprod(Z.G2)
  }
  if (withH) {
    vars <- par[(sigma2s + tau2s + rhos + 1):(sigma2s + 
                                                tau2s + rhos + gamma2s)]
    cors <- par[(sigma2s + tau2s + rhos + gamma2s + 1):(sigma2s + 
                                                          tau2s + rhos + gamma2s + phis)]
    resH <- .con.E(v = vars, r = cors, v.arg = gamma2.arg, 
                   r.arg = phi.arg, Z1 = Z.H1, Z2 = Z.H2, levels.r = h.levels.r, 
                   values = h.values, Dmat = h.Dmat, struct = struct[2], 
                   cholesky = cholesky[2], vctransf = vctransf, vccov = vccov, 
                   nearpd = nearpd, sparse = sparse)
    gamma2 <- resH$v
    phi <- resH$r
    H <- resH$E
    if (!is.null(vccon)) {
      if (!is.null(vccon$gamma2)) {
        for (id in unique(vccon$gamma2)) {
          gamma2[vccon$gamma2 == id] <- mean(gamma2[vccon$gamma2 == 
                                                      id])
        }
      }
      if (!is.null(vccon$phi)) {
        for (id in unique(vccon$phi)) {
          phi[vccon$phi == id] <- mean(phi[vccon$phi == 
                                             id])
        }
      }
      resH <- .con.E(v = gamma2, r = phi, v.arg = gamma2.arg, 
                     r.arg = phi.arg, Z1 = Z.H1, Z2 = Z.H2, levels.r = h.levels.r, 
                     values = h.values, Dmat = h.Dmat, struct = struct[2], 
                     cholesky = FALSE, vctransf = FALSE, vccov = vccov, 
                     nearpd = nearpd, sparse = sparse)
      gamma2 <- resH$v
      phi <- resH$r
      H <- resH$E
    }
    M <- M + (Z.H1 %*% H %*% t(Z.H1)) * tcrossprod(Z.H2)
  }
  # -----------------------------------------------------------------
  # NEW BLOCK – handle optional Z.slope random‑effects matrix
  # -----------------------------------------------------------------
  if(withZslope){
    tau2B <- tau2B.arg
    vars <- par[(sigma2s + tau2s + rhos + gamma2s + phis + 1):(sigma2s + tau2s + rhos + gamma2s + phis + tau2Bs)]
    if (!is.null(Z.slope)) {
      n_slope <- ncol(Z.slope)
      if (!is.null(tau2B.arg)) {
        
        if (vctransf) tau2B <- log(tau2B)
      }
      if (vctransf) {
        tau2B <- ifelse(is.null(tau2B.arg), exp(vars), tau2B.arg)
      } else {
        tau2B <- ifelse(is.null(tau2B.arg), vars, tau2B.arg)
        tau2B[tau2B < 0] <- 0
      }
      tau2B <- ifelse(tau2B <= .Machine$double.eps * 10, 0, tau2B)
      Gb.slope <- diag(tau2B, n_slope)
      M <- M + Z.slope %*% Gb.slope %*% t(Z.slope)
    }
  }
  
  # -----------------------------------------------------------------
  # Continue with the original routine (unchanged from here on)
  # -----------------------------------------------------------------
  
  if (!hessian) {
    pars <- list(sigma2 = if (withS) sigma2 else NULL, tau2 = if (withG) tau2 else NULL, 
                 rho = if (withG) rho else NULL, gamma2 = if (withH) gamma2 else NULL, 
                 phi = if (withH) phi else NULL,
                 tau2B = if (withZslope) tau2B else NULL)
    try(assign("rma.mv", pars, envir = .metafor), silent = TRUE)
  }
  if (nearpd) 
    M <- as.matrix(nearPD(M)$mat)
  else {
    W <- try(suppressWarnings(chol2inv(chol(M))), silent = TRUE)
  }
  if (inherits(W, "try-error")) {
    if (dofit) {
      stop(("Final variance-covariance matrix is not positive definite."), 
           call. = FALSE)
    }
    else {
      llval <- -Inf
    }
  }else {
    if (!dofit || is.null(A)) {
      stXWX <- chol2inv(chol(as.matrix(t(X) %*% W %*% 
                                         X)))
      if (!optbeta) 
        beta <- matrix(stXWX %*% crossprod(X, W) %*% 
                         Y, ncol = 1)
      beta <- ifelse(is.na(beta.arg), beta, beta.arg)
      RSS <- as.vector(t(Y - X %*% beta) %*% W %*% (Y - 
                                                      X %*% beta))
      if (optbeta && lambda > 0) {
        if (intercept) {
          RSS <- RSS + c(lambda * crossprod(beta[-1]))
        }
        else {
          RSS <- RSS + c(lambda * crossprod(beta))
        }
      }
      vb <- stXWX
    }
    else {
      stXAX <- chol2inv(chol(as.matrix(t(X) %*% A %*% 
                                         X)))
      beta <- matrix(stXAX %*% crossprod(X, A) %*% Y, 
                     ncol = 1)
      beta <- ifelse(is.na(beta.arg), beta, beta.arg)
      RSS <- as.vector(t(Y - X %*% beta) %*% W %*% (Y - 
                                                      X %*% beta))
      if (optbeta && lambda > 0) {
        if (intercept) {
          RSS <- RSS + c(lambda * crossprod(beta[-1]))
        }
        else {
          RSS <- RSS + c(lambda * crossprod(beta))
        }
      }
      vb <- matrix(stXAX %*% t(X) %*% A %*% M %*% A %*% 
                     X %*% stXAX, nrow = pX, ncol = pX)
    }
    llvals <- c(NA_real_, NA_real_)
    if (dofit || !reml) 
      llvals[1] <- -1/2 * (k) * log(2 * base::pi) - 1/2 * 
      determinant(M, logarithm = TRUE)$modulus - 1/2 * 
      RSS
    if (dofit || reml) 
      llvals[2] <- -1/2 * (k - pX) * log(2 * base::pi) + 
      ifelse(REMLf, 1/2 * determinant(crossprod(X), 
                                      logarithm = TRUE)$modulus, 0) + -1/2 * determinant(M, 
                                                                                         logarithm = TRUE)$modulus - 1/2 * determinant(crossprod(X, 
                                                                                                                                                 W) %*% X, logarithm = TRUE)$modulus - 1/2 * 
      RSS
    if (dofit) {
      res <- list(beta = beta, vb = vb, M = M, llvals = llvals)
      if (withS) 
        res$sigma2 <- sigma2
      if (withG) {
        res$G <- G
        res$tau2 <- tau2
        res$rho <- rho
      }
      if (withH) {
        res$H <- H
        res$gamma2 <- gamma2
        res$phi <- phi
      }
      if (withZslope) { 
        res$Z.slope <- Z.slope
        res$tau2B <- tau2B 
      }
      return(res)
    }
    else {
      llval <- ifelse(reml, llvals[2], llvals[1])
    }
  }
  return(-1 * c(llval))
}
