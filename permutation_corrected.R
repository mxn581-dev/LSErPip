run_permutation_corrected <- function(adj_m, expr_sub, maxSR,
                                      original_SR,
                                      n_permutations = 100,
                                      seed = 456,
                                      n_cores = parallel::detectCores() - 1) {
  
  # Convert to sparse once — ~99% of adj_m is zeros
  if (!inherits(adj_m, "sparseMatrix")) {
    adj_m <- Matrix::Matrix(adj_m, sparse = TRUE)
  }
  
  n_samples <- ncol(expr_sub)
  n_genes <- nrow(expr_sub)
  
  cat("Running", n_permutations, "permutations (corrected SR) on", n_cores, "cores...\n")
  
  # Reproducible parallel seeds
  set.seed(seed)
  perm_seeds <- sample.int(.Machine$integer.max, n_permutations)
  
  perm_fn <- function(i) {
    set.seed(perm_seeds[i])
    
    # Shuffle within each column via index matrix
    perm_idx <- vapply(seq_len(n_samples), function(j) sample.int(n_genes),
                       integer(n_genes))
    perm_exp <- expr_sub[cbind(as.vector(perm_idx),
                  rep(seq_len(n_samples), each = n_genes))]
    dim(perm_exp) <- dim(expr_sub)
    
    # Vectorized SR for ALL samples at once (matrix-matrix multiply)
    # Mass action: sumexp = A %*% X  (genes x samples, one shot)
    sumexp_m <- as.matrix(adj_m %*% perm_exp)
    
    # Stationary distribution: pi_i = x_i * sumexp_i / Z
    invP_m <- perm_exp * sumexp_m
    nf <- colSums(invP_m)
    nf[nf == 0] <- 1
    invP_m <- t(t(invP_m) / nf)
    
    # Stochastic matrix per sample and local entropy
    # Done sample-by-sample here because p_m is genes x genes per sample
    SR_vec <- numeric(n_samples)
    for (s in seq_len(n_samples)) {
      se <- sumexp_m[, s]
      se[se == 0] <- 1
      # WRONG:  adj_m * perm_exp[, s]          → row-wise (x_i)
# RIGHT:  t(t(adj_m) * perm_exp[, s])    → column-wise (x_j)
      p_m <- t(t(adj_m) * perm_exp[, s])
      p_m <- p_m / se
      
      # Local entropy via sparse operations
      p_vals <- p_m@x               # nonzero entries only
      log_vals <- log(p_vals)
      p_log <- p_m
      p_log@x <- p_vals * log_vals
      S_v <- -Matrix::rowSums(p_log)
      
      SR_vec[s] <- sum(invP_m[, s] * S_v) / maxSR
    }
    SR_vec
  }
  
  # Parallel across permutations
  # Replace the mclapply call with:
cl <- parallel::makeCluster(n_cores)
parallel::clusterExport(cl, varlist = c(
  "expr_sub", "adj_m", "maxSR", "n_samples", "n_genes", "perm_seeds"
), envir = environment())
parallel::clusterEvalQ(cl, library(Matrix))

permuted_SR <- do.call(rbind,
  parallel::parLapply(cl, seq_len(n_permutations), perm_fn)
)

parallel::stopCluster(cl)
  
  # Statistics (unchanged)
  perm_means <- colMeans(permuted_SR, na.rm = TRUE)
  perm_sds <- apply(permuted_SR, 2, sd, na.rm = TRUE)
  
  z_scores <- (original_SR - perm_means) / perm_sds
  
  p_values <- vapply(seq_len(n_samples), function(s) {
    perm_dist <- permuted_SR[, s]
    mu <- mean(perm_dist)
    mean(abs(perm_dist - mu) >= abs(original_SR[s] - mu))
  }, numeric(1))
  
  cat("=== PERMUTATION RESULTS (Corrected SR) ===\n")
  cat(sprintf("  Samples with p < 0.05: %d / %d (%.1f%%)\n",
              sum(p_values < 0.05), n_samples,
              100 * sum(p_values < 0.05) / n_samples))
  cat(sprintf("  Mean |Z-score|: %.3f\n", mean(abs(z_scores))))
  cat(sprintf("  Z-score range: %.3f to %.3f\n", min(z_scores), max(z_scores)))
  cat(sprintf("  Observed SR mean: %.4f\n", mean(original_SR)))
  cat(sprintf("  Permuted SR mean: %.4f\n\n", mean(perm_means)))
  
  frac_sig <- mean(p_values < 0.05)
  if (frac_sig > 0.8) {
    cat("  -> HIGHLY SIGNIFICANT: Network captures real biological signal\n\n")
  } else if (frac_sig > 0.5) {
    cat("  -> MODERATELY SIGNIFICANT\n\n")
  } else {
    cat("  -> NOT SIGNIFICANT\n\n")
  }
  
  list(
    permuted_means = permuted_SR,
    original_means = original_SR,
    p_values = p_values,
    z_scores = z_scores
  )
}
expr_sub <- expr_corrected[mc_genes, , drop = FALSE]

perm_results <- run_permutation_corrected(
  adj_m = adj_m,
  expr_sub = expr_sub,
  maxSR = maxSR,
  original_SR = SR_per_sample,
  n_permutations = 100,
  seed = 456
)

# Export
perm_summary <- data.frame(
  sample = colnames(expr_sub),
  original_SR = perm_results$original_means,
  permuted_mean = colMeans(perm_results$permuted_means, na.rm = TRUE),
  permuted_sd = apply(perm_results$permuted_means, 2, sd, na.rm = TRUE),
  z_score = perm_results$z_scores,
  p_value = perm_results$p_values
)
write.csv(perm_summary, "permutation_corrected_summary.csv", row.names = FALSE)
write.csv(perm_results$permuted_means, "permutation_corrected_means.csv", row.names = FALSE)
cat("Saved: permutation_corrected_summary.csv, permutation_corrected_means.csv\n\n")
