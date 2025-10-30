# =========================
# 1. Build directed network
# =========================
build_directed_network_from_hat_row <- function(hat_matrix, row_label) {
  row_vals <- hat_matrix[row_label, , drop = FALSE]
  #
  edge_list <- c()
  for (col_label in colnames(row_vals)) {
    nodes <- unlist(strsplit(col_label, ":"))
    #
    if (length(nodes) == 2) {
      i <- nodes[1]; j <- nodes[2]; val <- row_vals[1, col_label]
      if (val > 0) {
        edge_list <- c(edge_list, i, j)
      }
      else if (val < 0) {
        edge_list <- c(edge_list, j, i)
      }
    }
  }
  #
  graph(edges = edge_list, directed = TRUE)
}

# ==========================================
# 2. Extract edge sets from igraph path list
# ==========================================
get_edges_from_path <- function(x) {
  nodes <- as.vector(x)
  if (length(nodes) < 2)
    return(character(0))
  #
  edges <- mapply(function(a, b) paste(a, b, sep = "-"),
                  nodes[-length(nodes)], nodes[-1])
  # Standardize to alphabetical order
  sapply(strsplit(edges, "-"), function(x) paste(sort(x), collapse = "-"))
}

# ======================================
# 3. Reduce matrix to full rank (if needed)
# ======================================
reduce_to_full_rank <- function(x, tol = 1e-8) {
  qr_decomp <- qr(x, tol = tol)
  independent <- qr_decomp$pivot[seq_len(qr_decomp$rank)]
  x[independent, independent, drop = FALSE]
}

# =======================================
# 4. Construct variance-covariance matrix V
# =======================================
construct_V_matrix <- function(x) {
  pairs <- t(combn(sort(x$trts), 2))
  pair_names <- apply(pairs, 1, paste0, collapse = ":")
  valid_pairs <- pair_names[pair_names %in% x$comparisons]
  #
  get_seTE_for_pair <- function(pair, seTEs) {
    trts <- unlist(strsplit(pair, ":"))
    seTEs[trts[1], trts[2]]
  }
  #
  V_diag <- sapply(valid_pairs, get_seTE_for_pair, 
                   seTEs = (x$seTE.direct.common)^2)
  #
  V <- diag(V_diag)
  rownames(V) <- colnames(V) <- valid_pairs
  #
  V
}

# ==========================================
# 5. Compute theta_p (path-specific effects)
# ==========================================
get_theta_for_path <- function(x, TEs) {
  nodes <- as.vector(x)
  #
  if (length(nodes) < 2)
    return(0)
  #
  edges <- mapply(function(a, b) TEs[a, b], nodes[-length(nodes)], nodes[-1])
  #
  sum(edges, na.rm = TRUE)
}

# =======================================
# 6. Core function: compute Q and p-value
# =======================================
run_path_inconsistency <- function(net, hat_common, node1, node2) {
  comp_intr <- paste0(node1, ":", node2)
  g <- build_directed_network_from_hat_row(hat_common, comp_intr)
  all_paths <- all_simple_paths(g, from = node1, to = node2)
  if (length(all_paths) == 0) {
    return(data.frame(Comparison = paste(node1, node2, sep = " - "),
                      Q = NA, df = NA, p_value = NA,
                      Note = "No path found"))
  }
  #
  path_list <- lapply(all_paths, function(x) V(g)[x]$name)

  # Build path overlap matrix
  edge_sets <- lapply(all_paths, get_edges_from_path)
  n <- length(edge_sets)
  path_matrix <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j)
        path_matrix[i, j] <- length(edge_sets[[i]])
      else
        path_matrix[i, j] <- length(intersect(edge_sets[[i]], edge_sets[[j]]))
    }
  }
  #
  rownames(path_matrix) <- colnames(path_matrix) <- paste0("Path", seq_len(n))
  
  # Rank reduction if singular
  det_val <- det(path_matrix)
  if (abs(det_val) < 1e-6) {
    reduced_matrix <- reduce_to_full_rank(path_matrix)
  }
  else {
    reduced_matrix <- path_matrix
  }
  #
  kept_paths <- rownames(reduced_matrix)
  kept_indices <- as.integer(sub("Path", "", kept_paths))
  kept_path_list <- path_list[kept_indices]
  
  # Construct C matrix
  comp_edges <- gsub(":", "-", net$comparisons)
  comp_edges_std <- sapply(strsplit(comp_edges, "-"), function(x) paste(sort(x), collapse = "-"))
  C <- matrix(0, nrow = length(kept_paths), ncol = length(comp_edges_std))
  rownames(C) <- kept_paths
  colnames(C) <- net$comparisons
  #
  edge_sets_std <- lapply(kept_path_list, get_edges_from_path)
  for (i in seq_along(kept_paths)) {
    C[i, comp_edges_std %in% edge_sets_std[[i]]] <- 1
  }

  # Construct V and align with C
  V <- construct_V_matrix(net)
  #
  common_edges <- intersect(colnames(C), colnames(V))
  #
  C <- C[, common_edges, drop = FALSE]
  V <- V[common_edges, common_edges, drop = FALSE]

  # Compute S
  S <- C %*% V %*% t(C)
  
  if (det(S) < 1e-6) {
    S_inv <- MASS::ginv(S)
    note <- "Matrix S singular, pseudo-inverse used"
  }
  else {
    S_inv <- solve(S)
    note <- ""
  }
  
  # Compute theta_p
  theta_p <- sapply(kept_path_list, get_theta_for_path, 
                    TEs = net$TE.direct.common)
  #write.csv(theta_p, file.path(dir_path, "theta_p.csv"), row.names = TRUE  )
  theta_diff <- theta_p - net$TE.common[node1, node2]

  # Compute Q
  Q <- t(theta_diff) %*% S_inv %*% theta_diff
  df <- nrow(S) - 1
  pval <- pchisq(Q, df = df, lower.tail = FALSE)
  
  res <-
    list(
      results =
        data.frame(comparison = paste(node1, node2, sep = " - "),
                   Q = as.numeric(Q), df = df,
                   pval = as.numeric(pval),
                   path_index = kept_indices, #I(list(kept_indices)),
                   note = note),
      path_matrix = path_matrix,
      S = S,
      theta_p = theta_p
  )
  #
  class(res) <- "netpath"
  #
  res
}

# ======================================
# 7. Wrapper for single or all comparisons
# ======================================
run_all_path_inconsistencies <- function(net, hat_common,
                                         node1 = NULL, node2 = NULL) {
  if (!is.null(node1) && !is.null(node2))
    return(run_path_inconsistency(net, hat_common, node1, node2))
  #
  stop("Error: A comparison of interest must be chosen. Please specify both 
          'node1' and 'node2'.")
}

# =============================================
# 8. Create standardized heatmap 
# ============================================
create_standardized_heatmap <- function(xhat, Sigma, indx) {
  # Input validation
  N <- length(xhat)
  if (!all(dim(Sigma) == c(N, N))) {
    stop("Covariance matrix dimensions do not match xhat length")
  }
  
  # Calculate pairwise standardized differences (z_ij)
  z <- matrix(0, nrow = N, ncol = N)
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      var <- Sigma[i,i] + Sigma[j,j] - 2 * Sigma[i,j]
      var <- max(var, 1e-12)  # regularize to avoid negatives
      z[i,j] <- (xhat[i] - xhat[j]) / sqrt(var)
    }
  }
  abs_z <- abs(z)
  
  # Create heatmap using ggplot2
  p_i <- p_j <- value <- NULL
  
  # Convert matrix to long format for ggplot
  melted_z <- melt(abs_z)
  colnames(melted_z) <- c("p_i", "p_j", "value")
  
  # Create the plot
  mylabs <- as.expression(lapply(indx, function(i) bquote(pi[.(i)])))
  
  p <- ggplot(melted_z, aes(x = p_j, y = p_i, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red",
                        name = expression(z[p*","*p[prime]]^{T[6]*T[9]})) +
    scale_x_continuous(
      breaks = seq_along(mylabs),
      labels = mylabs,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(mylabs),
      labels = mylabs,
      expand = c(0, 0)
    ) +
    theme_minimal(base_size = 15) +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 16), # increased legend title size
          legend.text = element_text(size = 16),  
          panel.grid = element_blank()) +
    coord_fixed()

  # Add text annotations
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      text_val <- if(i == j)
        sprintf("%.0f", abs_z[i, j])
      else
        sprintf("%.2f", abs_z[i, j])
      
      # Calculate text color based on background intensity
      text_color <- if(abs_z[i, j] > max(abs_z) / 2) "white" else "black"
      
      p <- p + annotate("text", x = j, y = i, 
                        label = text_val, color = text_color, size = 6)
    }
  }
  #
  p
}
