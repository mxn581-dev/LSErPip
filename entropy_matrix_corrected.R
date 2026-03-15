library(tidyverse)
library(igraph)
library(entropy)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


# Memory management
gc(reset = TRUE)
options(expressions = 500000)

clean_mem <- function() {
  gc(verbose = FALSE)
  Sys.sleep(0.5)
}
# STEP 1: Load PPI Network with CONFIDENCE FILTERING

cat("STEP 1: Loading PPI network...\n")

# WINDOWS FILE PATH - Update this to your actual path
# Example: "C:/Users/YourName/Downloads/9606.protein.links.full.v12.0.txt"
ppi_file <- ""

ppi_raw <- read.table(ppi_file, header = TRUE, stringsAsFactors = FALSE)

cat("Original network size:", nrow(ppi_raw), "edges\n")
cat("Columns:", paste(colnames(ppi_raw), collapse = ", "), "\n\n")

# Extract protein pairs and score
edges <- data.frame(
  from = ppi_raw[, 1],
  to = ppi_raw[, 2],
  score = ppi_raw$combined_score,  # Use combined_score column
  stringsAsFactors = FALSE
)

cat("Score range:", min(edges$score, na.rm = TRUE), "-", 
    max(edges$score, na.rm = TRUE), "\n")

# CONFIDENCE FILTERING
confidence_threshold <- 400

cat("\n=== CONFIDENCE FILTERING ===\n")
cat("Using threshold:", confidence_threshold, "\n")

edges_filtered <- edges[edges$score >= confidence_threshold, ]

cat("Edges before:", nrow(edges), "\n")
cat("Edges after:", nrow(edges_filtered), "\n")
cat("Reduction:", round(100 * (1 - nrow(edges_filtered)/nrow(edges)), 1), "%\n\n")

# Increase threshold if network too large
if (nrow(edges_filtered) > 100000) {
  confidence_threshold <- 700
  cat("Network large. Increasing threshold to", confidence_threshold, "\n")
  edges_filtered <- edges[edges$score >= confidence_threshold, ]
  cat("Edges after filtering:", nrow(edges_filtered), "\n\n")
}

if (nrow(edges_filtered) > 50000) {
  confidence_threshold <- 900
  cat("Network still large. Using highest confidence:", confidence_threshold, "\n")
  edges_filtered <- edges[edges$score >= confidence_threshold, ]
  cat("Edges after filtering:", nrow(edges_filtered), "\n\n")
}

# Keep only protein columns
edges <- edges_filtered[, c("from", "to")]

# Remove 9606. prefix
edges$from <- gsub("^9606\\.", "", edges$from)
edges$to <- gsub("^9606\\.", "", edges$to)

# Remove duplicates and self-loops
edges <- edges[edges$from != edges$to, ]
edges <- unique(edges)

rm(ppi_raw, edges_filtered)
clean_mem()

cat("Final edges after cleaning:", nrow(edges), "\n\n")
# STEP 2: Load Expression Data (68 Samples)
cat("STEP 2: Loading expression data...\n")

expr_data <- read.csv("path/to/your/expression_matrix.csv",
                      header = TRUE, row.names = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)

cat("Expression (before filtering):", nrow(expr_data), "genes x", 
    ncol(expr_data), "samples\n")

# Filter to CR/PD samples only
all_samples <- colnames(expr_data)
cohort_samples <- all_samples

if (length(cohort_samples) > 0) {
  expr_data <- expr_data[, cohort_samples, drop = FALSE]
  cat("Expression (after filtering):", nrow(expr_data), "genes x", 
      ncol(expr_data), "samples\n")
} else {
  cat("WARNING: No CR/PD samples found. Using all samples.\n")
}

cat("Sample names:", paste(head(colnames(expr_data), 5), "..."), "\n\n")
clean_mem()
# STEP 3: ID Conversion (ENSP to Gene Symbols)

cat("STEP 3: Converting ENSP IDs to Gene Symbols...\n")

# Install biomaRt if needed
if (!require("biomaRt", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt", ask = FALSE)
  library(biomaRt)
}

# Get unique ENSP IDs
all_ensp <- unique(c(edges$from, edges$to))
cat("Converting", length(all_ensp), "ENSP IDs...\n")

# Conversion function
convert_ensp <- function(ensp_ids, batch_size = 200) {
  
  library(biomaRt)
  
  ensembl <- tryCatch({
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  }, error = function(e) {
    cat("Trying mirror...\n")
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
               mirror = "useast")
  })
  
  n_batches <- ceiling(length(ensp_ids) / batch_size)
  all_results <- data.frame()
  
  for (i in 1:n_batches) {
    cat(sprintf("  Batch %d/%d\r", i, n_batches))
    
    start <- (i - 1) * batch_size + 1
    end <- min(i * batch_size, length(ensp_ids))
    batch <- ensp_ids[start:end]
    
    result <- tryCatch({
      getBM(attributes = c('ensembl_peptide_id', 'hgnc_symbol'),
            filters = 'ensembl_peptide_id',
            values = batch,
            mart = ensembl)
    }, error = function(e) {
      data.frame()
    })
    
    if (nrow(result) > 0) {
      all_results <- rbind(all_results, result)
    }
    
    Sys.sleep(0.3)
  }
  
  cat("\n")
  return(all_results[all_results$hgnc_symbol != "", ])
}

id_map <- convert_ensp(all_ensp, batch_size = 200)

cat("Mapped:", nrow(id_map), "IDs\n")
write.csv(id_map, "id_mapping.csv", row.names = FALSE)

# Apply mapping to edges
map_vec <- setNames(id_map$hgnc_symbol, id_map$ensembl_peptide_id)
edges$from_symbol <- map_vec[edges$from]
edges$to_symbol <- map_vec[edges$to]

# Keep only mapped edges
edges_mapped <- edges[!is.na(edges$from_symbol) & !is.na(edges$to_symbol), ]
edges_mapped <- edges_mapped[, c("from_symbol", "to_symbol")]
colnames(edges_mapped) <- c("from", "to")

rm(edges, id_map, map_vec)
clean_mem()

cat("Mapped edges:", nrow(edges_mapped), "\n")

# Create network
ppi_net <- graph_from_data_frame(edges_mapped, directed = FALSE)
ppi_net <- simplify(ppi_net)

rm(edges_mapped)
clean_mem()

cat("Network:", vcount(ppi_net), "nodes,", ecount(ppi_net), "edges\n")

# Check overlap with expression data
overlap <- intersect(V(ppi_net)$name, rownames(expr_data))
cat("Genes in both network and expression:", length(overlap), "\n\n")

if (length(overlap) < 100) {
  stop("Too few overlapping genes! Check ID conversion.")
}

#CORRECTED FIXES

library(igraph)
library(tidyverse)
# FIX 1: Handle log-scale expression data

prepare_expression <- function(expr_matrix) {
  
  cat("Checking expression data scale...\n")
  cat(sprintf("  Range: %.3f to %.3f\n", min(expr_matrix), max(expr_matrix)))
  cat(sprintf("  Negative values: %d (%.1f%%)\n",
              sum(expr_matrix < 0),
              100 * sum(expr_matrix < 0) / length(expr_matrix)))
  
  # If there are negative values, data is log-transformed
  if (any(expr_matrix < 0)) {
    cat("  -> Detected log-scale data (negative values present)\n")
    cat("  -> Converting back to linear scale with 2^x\n")
    expr_matrix <- 2^expr_matrix
  }
  
  # Ensure all values are strictly positive (required for mass action)
  # SCENT uses an offset of 1.1 before log, so minimum is ~0.14
  # For safety, replace any zeros/near-zeros with a small pseudocount
  min_nonzero <- min(expr_matrix[expr_matrix > 0], na.rm = TRUE)
  pseudocount <- min_nonzero * 0.01
  expr_matrix[expr_matrix <= 0] <- pseudocount
  
  cat(sprintf("  Final range: %.4f to %.3f\n", min(expr_matrix), max(expr_matrix)))
  cat(sprintf("  All values positive: %s\n\n", all(expr_matrix > 0)))
  
  return(expr_matrix)
}

# STEP: Build adjacency matrix from igraph network

build_adjacency <- function(network, genes) {
  
  cat("Building adjacency matrix...\n")
  
  # Subset network to genes present in expression data
  net_sub <- induced_subgraph(network, genes)
  
  # Extract maximally connected component (as in SCENT DoIntegPPI)
  components <- components(net_sub)
  largest_comp <- which.max(components$csize)
  mc_nodes <- names(components$membership[components$membership == largest_comp])
  net_mc <- induced_subgraph(net_sub, mc_nodes)
  
  cat(sprintf("  Maximally connected component: %d nodes, %d edges\n",
              vcount(net_mc), ecount(net_mc)))
  
  # Get adjacency matrix (binary, undirected)
  adj_m <- as_adjacency_matrix(net_mc, sparse = FALSE)
  diag(adj_m) <- 0  # No self-loops
  
  return(list(
    adj = adj_m,
    network = net_mc,
    genes = mc_nodes
  ))
}

# FIX 3: Compute maximum entropy rate (maxSR)
# maxSR = log(lambda_max) where lambda_max is largest eigenvalue of A
# This is the theoretical maximum entropy rate on this network

compute_maxSR <- function(adj_m) {
  
  cat("Computing maximum entropy rate (maxSR)...\n")
  
  # Use ARPACK to find largest eigenvalue of adjacency matrix
  fa <- function(x, extra = NULL) {
    as.vector(adj_m %*% x)
  }
  
  ap <- arpack(fa, options = list(n = nrow(adj_m), nev = 1, which = "LM"), sym = TRUE)
  lambda_max <- ap$values[1]
  maxSR <- log(lambda_max)
  
  cat(sprintf("  Largest eigenvalue: %.4f\n", lambda_max))
  cat(sprintf("  maxSR = log(lambda): %.4f\n\n", maxSR))
  
  return(maxSR)
}


# ============================================================================
# LOCAL ENTROPY FUNCTION
# S_i = -sum_j p_ij * log(p_ij)  (natural log, following SCENT)
# ============================================================================

compute_local_entropy <- function(prob_vec) {
  # Remove zeros (log(0) undefined)
  p <- prob_vec[prob_vec > 0]
  return(-sum(p * log(p)))
}


# FIX 2 + 3: Compute entropy rate per sample
# Following SCENT CompSRanaPRL exactly:

#   For sample s with expression vector x:
#   1. sumexp_i = sum_j A_ij * x_j         (weighted degree)
#   2. p_ij = A_ij * x_j / sumexp_i        (stochastic/transition matrix)
#   3. pi_i = x_i * sumexp_i / sum(x * sumexp) (stationary distribution)
#   4. S_i = -sum_j p_ij * log(p_ij)       (local entropy)
#   5. SR = sum(pi_i * S_i)                (entropy rate)
#   6. SR_norm = SR / maxSR                (normalized, 0 to 1)

calc_entropy_rate_per_sample <- function(adj_m, expr_matrix, genes, maxSR) {
  
  n_genes <- length(genes)
  n_samples <- ncol(expr_matrix)
  sample_names <- colnames(expr_matrix)
  
  # Subset expression to network genes
  expr_sub <- expr_matrix[genes, , drop = FALSE]
  
  # Initialize output matrices
  SR_vec <- numeric(n_samples)                           # global SR per sample
  local_entropy_mat <- matrix(NA, n_genes, n_samples)    # local S_i per gene per sample
  norm_local_mat <- matrix(NA, n_genes, n_samples)       # normalized local entropy
  pi_mat <- matrix(NA, n_genes, n_samples)               # stationary distribution
  
  rownames(local_entropy_mat) <- genes
  colnames(local_entropy_mat) <- sample_names
  rownames(norm_local_mat) <- genes
  colnames(norm_local_mat) <- sample_names
  rownames(pi_mat) <- genes
  colnames(pi_mat) <- sample_names
  names(SR_vec) <- sample_names
  
  # Degree vector for normalization of local entropy
  degree_vec <- rowSums(adj_m)
  max_local_S <- ifelse(degree_vec > 0, log(degree_vec), 1)
  
  cat("Computing entropy rate for", n_samples, "samples...\n")
  
  for (s in 1:n_samples) {
    if (s %% 10 == 0) cat(sprintf("  Sample %d/%d\r", s, n_samples))
    
    # Expression vector for this sample
    exp_v <- as.numeric(expr_sub[, s])
    
    # ---- MASS ACTION PRINCIPLE (FIX 2) ----
    # sumexp_i = sum_j A_ij * x_j  (for each gene i, sum of neighbor expressions)
    sumexp_v <- as.vector(adj_m %*% matrix(exp_v, ncol = 1))
    
    # ---- STATIONARY DISTRIBUTION (FIX 3) ----
    # pi_i = x_i * sumexp_i / normalization
    # Under detailed balance, this is the invariant measure
    invP_v <- exp_v * sumexp_v
    nf <- sum(invP_v)
    
    if (nf == 0) {
      SR_vec[s] <- NA
      next
    }
    
    invP_v <- invP_v / nf  # Normalized: sum(invP_v) = 1
    pi_mat[, s] <- invP_v
    
    # ---- STOCHASTIC/TRANSITION MATRIX (FIX 2) ----
    # p_ij = A_ij * x_j / sumexp_i
    # This is the probability of signaling from gene i to gene j
    # Following SCENT: p.m <- t(t(adj.m)*exp.v)/sumexp.v
    
    # Handle genes with no neighbors (sumexp = 0)
    safe_sumexp <- sumexp_v
    safe_sumexp[safe_sumexp == 0] <- 1  # Avoid division by zero
    
    p_m <- t(t(adj_m) * exp_v) / safe_sumexp
    
    # ---- LOCAL ENTROPY (per gene) ----
    # S_i = -sum_j p_ij * log(p_ij)
    S_v <- apply(p_m, 1, compute_local_entropy)
    local_entropy_mat[, s] <- S_v
    
    # Normalized local entropy: S_i / log(k_i)
    norm_local_mat[, s] <- ifelse(degree_vec > 0, S_v / max_local_S, 0)
    
    # ---- ENTROPY RATE (FIX 3) ----
    # SR = sum(pi_i * S_i)
    SR <- sum(invP_v * S_v)
    
    # Normalize by maxSR
    SR_vec[s] <- SR / maxSR
  }
  
  cat("\n\n")
  
  cat("=== ENTROPY RATE RESULTS ===\n")
  cat(sprintf("  Samples: %d\n", n_samples))
  cat(sprintf("  Genes in network: %d\n", n_genes))
  cat(sprintf("  SR range: %.4f - %.4f\n", min(SR_vec, na.rm = TRUE), max(SR_vec, na.rm = TRUE)))
  cat(sprintf("  SR mean: %.4f\n", mean(SR_vec, na.rm = TRUE)))
  cat(sprintf("  SR sd: %.4f\n\n", sd(SR_vec, na.rm = TRUE)))
  
  return(list(
    SR = SR_vec,                          # Global entropy rate per sample (normalized)
    local_entropy = local_entropy_mat,    # Unnormalized local entropy (S_i)
    norm_local_entropy = norm_local_mat,  # Normalized local entropy (S_i / log(k_i))
    pi = pi_mat,                          # Stationary distribution per sample
    maxSR = maxSR,
    genes = genes
  ))
}

# MAIN PIPELINE
cat("CORRECTED SIGNALING ENTROPY RATE PIPELINE\n")
cat("Following Teschendorff et al. / SCENT framework\n")
# ---- 1. Prepare expression data ----
expr_corrected <- prepare_expression(expr_data)

# ---- 2. Build adjacency matrix from PPI network ----
overlap_genes <- intersect(V(ppi_net)$name, rownames(expr_corrected))
cat(sprintf("Genes overlapping network and expression: %d\n\n", length(overlap_genes)))

adj_result <- build_adjacency(ppi_net, overlap_genes)
adj_m <- adj_result$adj
mc_genes <- adj_result$genes

cat(sprintf("Genes in maximally connected component: %d\n\n", length(mc_genes)))

# ---- 3. Compute maxSR ----
maxSR <- compute_maxSR(adj_m)

# ---- 4. Compute entropy rate ----
sr_result <- calc_entropy_rate_per_sample(adj_m, expr_corrected, mc_genes, maxSR)

# The normalized local entropy matrix (genes x samples) replaces your old entropy_matrix
entropy_matrix_corrected <- sr_result$norm_local_entropy

# Global SR per sample
SR_per_sample <- sr_result$SR

cat("Output objects created:\n")
cat("  entropy_matrix_corrected: normalized local entropy (genes x samples)\n")
cat(sprintf("    Dimensions: %d genes x %d samples\n",
            nrow(entropy_matrix_corrected), ncol(entropy_matrix_corrected)))
cat("  SR_per_sample: global entropy rate per sample (single value per sample)\n\n")

write.csv(entropy_matrix_corrected, "entropy_matrix_corrected.csv", row.names = TRUE)
write.csv(data.frame(
  sample = names(SR_per_sample),
  SR = SR_per_sample
), "SR_per_sample.csv", row.names = FALSE)

# Gene-level summary
gene_summary <- data.frame(
  gene = mc_genes,
  degree = rowSums(adj_m),
  mean_local_entropy = rowMeans(sr_result$local_entropy, na.rm = TRUE),
  mean_norm_entropy = rowMeans(sr_result$norm_local_entropy, na.rm = TRUE),
  sd_norm_entropy = apply(sr_result$norm_local_entropy, 1, sd, na.rm = TRUE),
  mean_pi = rowMeans(sr_result$pi, na.rm = TRUE)
) %>% arrange(desc(mean_norm_entropy))

write.csv(gene_summary, "gene_entropy_summary_corrected.csv", row.names = FALSE)

cat("Saved:\n")
cat("  - entropy_matrix_corrected.csv\n")
cat("  - SR_per_sample.csv\n")
cat("  - gene_entropy_summary_corrected.csv\n\n")

cat("============================================================\n")
cat("DONE. You can now re-run your bootstrap/permutation/perturbation\n")
cat("analyses using 'entropy_matrix_corrected' in place of 'entropy_matrix'\n")
cat("and 'expr_corrected' in place of 'expr_data'\n")
cat("============================================================\n")
