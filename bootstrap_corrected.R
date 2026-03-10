bootstrap_entropy <- function(entropy_mat, n_bootstrap = 100, seed = 123) {
  
  set.seed(seed)
  n_genes <- nrow(entropy_mat)
  n_samples <- ncol(entropy_mat)
  
  # Store bootstrap results
  bootstrap_means <- matrix(NA, nrow = n_bootstrap, ncol = n_samples)
  colnames(bootstrap_means) <- colnames(entropy_mat)
  
  bootstrap_sds <- matrix(NA, nrow = n_bootstrap, ncol = n_samples)
  colnames(bootstrap_sds) <- colnames(entropy_mat)
  
  cat("Running", n_bootstrap, "bootstrap iterations...\n")
  
  for (b in 1:n_bootstrap) {
    if (b %% 10 == 0) cat(sprintf("  Iteration %d/%d\r", b, n_bootstrap))
    
    # Sample genes WITH REPLACEMENT
    boot_indices <- sample(1:n_genes, size = n_genes, replace = TRUE)
    boot_matrix <- entropy_mat[boot_indices, , drop = FALSE]
    
    # Calculate mean and SD for each sample
    bootstrap_means[b, ] <- colMeans(boot_matrix, na.rm = TRUE)
    bootstrap_sds[b, ] <- apply(boot_matrix, 2, sd, na.rm = TRUE)
  }
  
  cat("\n\n")
  
  # Calculate coefficient of variation across bootstrap iterations
  cv_per_sample <- apply(bootstrap_means, 2, sd, na.rm = TRUE) / 
                   colMeans(bootstrap_means, na.rm = TRUE)
  
  results <- list(
    bootstrap_means = bootstrap_means,
    bootstrap_sds = bootstrap_sds,
    cv_per_sample = cv_per_sample,
    overall_cv = mean(cv_per_sample, na.rm = TRUE)
  )
  
  return(results)
}

# Run bootstrap
boot_results <- bootstrap_entropy(entropy_matrix_corrected, n_bootstrap = 100)

cat("Bootstrap Results:\n")
cat("Overall CV across samples:", round(boot_results$overall_cv, 4), "\n")
cat("Range of CVs:", round(range(boot_results$cv_per_sample), 4), "\n\n")

# Interpretation guide
if (boot_results$overall_cv < 0.05) {
  cat("Interpretation: ROBUST - Entropy estimates are very stable (CV < 5%)\n\n")
} else if (boot_results$overall_cv < 0.10) {
  cat("Interpretation: MODERATELY ROBUST - Entropy estimates are reasonably stable (CV < 10%)\n\n")
} else {
  cat("Interpretation: LOW ROBUSTNESS - Entropy estimates vary considerably (CV > 10%)\n\n")
}

# Save bootstrap results
write.csv(boot_results$bootstrap_means, "bootstrap_means.csv", row.names = FALSE)
write.csv(data.frame(sample = colnames(entropy_matrix_corrected), 
                     cv = boot_results$cv_per_sample),
          "bootstrap_cv_per_sample.csv", row.names = FALSE)
