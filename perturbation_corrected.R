run_perturbation_corrected <- function(network, adj_m, expr_sub, maxSR,
                                       original_SR,
                                       perturbation_rates = c(0.05, 0.10, 0.20),
                                       n_iterations = 100,
                                       seed = 789,
                                       n_cores = parallel::detectCores() - 1) {
  
  n_samples <- ncol(expr_sub)
  mc_genes <- rownames(adj_m)
  n_genes <- length(mc_genes)
  
  # Convert to sparse once
  if (!inherits(adj_m, "sparseMatrix")) {
    adj_m <- Matrix::Matrix(adj_m, sparse = TRUE)
  }
  
  # Extract existing edges as (i,j) pairs from upper triangle
  adj_upper <- Matrix::triu(adj_m)
  edge_idx <- Matrix::which(adj_upper > 0, arr.ind = TRUE)
  original_edges <- nrow(edge_idx)
  
  # All possible non-edge pairs would be huge — instead we sample on the fly
  
  results_list <- list()
  
  for (rate in perturbation_rates) {
    
    cat(sprintf("\n=== Testing %.0f%% edge perturbation (corrected SR) ===\n", rate * 100))
    
    n_to_modify <- round(original_edges * rate)
    n_to_delete <- round(n_to_modify / 2)
    n_to_add <- round(n_to_modify / 2)
    
    set.seed(seed)
    iter_seeds <- sample.int(.Machine$integer.max, n_iterations)
    
    perturb_fn <- function(iter) {
      set.seed(iter_seeds[iter])
      
      # ---- Perturb edges directly on sparse matrix ----
      
      # Delete random edges
      keep_mask <- rep(TRUE, original_edges)
      if (n_to_delete > 0) {
        del_idx <- sample.int(original_edges, min(n_to_delete, original_edges))
        keep_mask[del_idx] <- FALSE
      }
      remaining <- edge_idx[keep_mask, , drop = FALSE]
      
      # Add random edges (avoiding self-loops and duplicates)
      new_i <- integer(0)
      new_j <- integer(0)
      if (n_to_add > 0) {
        # Oversample to account for self-loops and duplicates
        cand_i <- sample.int(n_genes, n_to_add * 3, replace = TRUE)
        cand_j <- sample.int(n_genes, n_to_add * 3, replace = TRUE)
        valid <- cand_i != cand_j
        cand_i <- cand_i[valid]
        cand_j <- cand_j[valid]
        # Canonicalize to upper triangle
        swap <- cand_i > cand_j
        tmp <- cand_i[swap]
        cand_i[swap] <- cand_j[swap]
        cand_j[swap] <- tmp
        # Deduplicate
        edge_keys <- paste0(cand_i, "_", cand_j)
        uniq <- !duplicated(edge_keys)
        cand_i <- cand_i[uniq]
        cand_j <- cand_j[uniq]
        n_take <- min(n_to_add, length(cand_i))
        new_i <- cand_i[seq_len(n_take)]
        new_j <- cand_j[seq_len(n_take)]
      }
      
      # Build perturbed sparse adjacency (symmetric)
      all_i <- c(remaining[, 1], new_i)
      all_j <- c(remaining[, 2], new_j)
      # Both directions for symmetric matrix
      adj_pert <- Matrix::sparseMatrix(
        i = c(all_i, all_j),
        j = c(all_j, all_i),
        x = 1,
        dims = c(n_genes, n_genes),
        dimnames = list(mc_genes, mc_genes)
      )
      # sparseMatrix sums duplicates, clamp to binary
      adj_pert@x[adj_pert@x > 1] <- 1
      
      edge_count <- Matrix::nnzero(Matrix::triu(adj_pert))
      
      # ---- maxSR for perturbed network ----
      maxSR_pert <- tryCatch({
        ev <- RSpectra::eigs_sym(adj_pert, k = 1, which = "LM")
        if (ev$values[1] <= 1) return(rep(NA, n_samples + 1))
        log(ev$values[1])
      }, error = function(e) {
        return(rep(NA, n_samples + 1))
      })
      
      if (length(maxSR_pert) > 1) return(maxSR_pert)  # NA fallback
      
      # ---- Vectorized SR computation (same pattern as permutation) ----
      sumexp_m <- as.matrix(adj_pert %*% expr_sub)
      
      invP_m <- expr_sub * sumexp_m
      nf <- colSums(invP_m)
      nf[nf == 0] <- 1
      invP_m <- t(t(invP_m) / nf)
      
      SR_vec <- numeric(n_samples)
      for (s in seq_len(n_samples)) {
        se <- sumexp_m[, s]
        se[se == 0] <- 1
        p_m <- t(t(adj_pert) * expr_sub[, s])
        p_m <- p_m / se
        
        p_vals <- p_m@x
        log_vals <- log(p_vals)
        p_log <- p_m
        p_log@x <- p_vals * log_vals
        S_v <- -Matrix::rowSums(p_log)
        
        SR_vec[s] <- sum(invP_m[, s] * S_v) / maxSR_pert
      }
      
      c(SR_vec, edge_count)
    }
    
    # ---- Parallel ----
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, varlist = c(
      "expr_sub", "adj_m", "edge_idx", "original_edges", "mc_genes",
      "n_genes", "n_samples", "n_to_delete", "n_to_add", "iter_seeds"
    ), envir = environment())
    parallel::clusterEvalQ(cl, { library(Matrix); library(RSpectra) })
    
    raw <- parallel::parLapply(cl, seq_len(n_iterations), perturb_fn)
    parallel::stopCluster(cl)
    
    raw_mat <- do.call(rbind, raw)
    perturbed_SR <- raw_mat[, seq_len(n_samples), drop = FALSE]
    edge_counts <- raw_mat[, n_samples + 1]
    colnames(perturbed_SR) <- colnames(expr_sub)
    
    # ---- Statistics ----
    cv_per_sample <- apply(perturbed_SR, 2, sd, na.rm = TRUE) /
                     colMeans(perturbed_SR, na.rm = TRUE)
    overall_cv <- mean(cv_per_sample, na.rm = TRUE)
    mean_shift <- mean(colMeans(perturbed_SR, na.rm = TRUE) - original_SR, na.rm = TRUE)
    
    cat(sprintf("  Overall CV: %.4f (%.2f%%)\n", overall_cv, overall_cv * 100))
    cat(sprintf("  Mean edges: %.0f (original: %d, change: %+.1f%%)\n",
                mean(edge_counts, na.rm = TRUE), original_edges,
                100 * (mean(edge_counts, na.rm = TRUE) - original_edges) / original_edges))
    cat(sprintf("  Mean SR shift from original: %+.4f\n", mean_shift))
    
    if (overall_cv < 0.01) {
      cat("  -> VERY ROBUST\n")
    } else if (overall_cv < 0.05) {
      cat("  -> ROBUST\n")
    } else if (overall_cv < 0.10) {
      cat("  -> MODERATE\n")
    } else {
      cat("  -> SENSITIVE\n")
    }
    
    results_list[[paste0("rate_", rate)]] <- list(
      perturbation_rate = rate,
      perturbed_means = perturbed_SR,
      original_means = original_SR,
      cv_per_sample = cv_per_sample,
      overall_cv = overall_cv,
      mean_shift = mean_shift,
      edge_counts = edge_counts
    )
  }
  
  cat("\n===== PERTURBATION SUMMARY =====\n")
  for (rate_name in names(results_list)) {
    r <- results_list[[rate_name]]
    cat(sprintf("  %s: CV = %.4f (%.2f%%), shift = %+.4f\n",
                rate_name, r$overall_cv, r$overall_cv * 100, r$mean_shift))
  }
  cat("================================\n\n")
  
  return(results_list)
}
perturb_results <- run_perturbation_corrected(
  network = ppi_net,
  adj_m = adj_m,
  expr_sub = expr_sub,
  maxSR = maxSR,
  original_SR = SR_per_sample,
  perturbation_rates = c(0.05, 0.10, 0.20),
  n_iterations = 100,
  seed = 789
)
