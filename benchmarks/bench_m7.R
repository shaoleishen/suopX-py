#!/usr/bin/env Rscript
# Benchmark: R SoupX on M7 sample

library(SoupX)
library(Matrix)

DATA_DIR <- "/home/bioshen/CLI/snRNAseq/M7"
RAW_DIR <- file.path(DATA_DIR, "raw_feature_bc_matrix")
FILT_DIR <- file.path(DATA_DIR, "filtered_feature_bc_matrix")

cat(sprintf("[%s] Starting R SoupX benchmark on M7\n", Sys.time()))
t0 <- Sys.time()

# ---------------------------------------------------------------
# Helper: read 10X mtx directory
# ---------------------------------------------------------------
read_10x_dir <- function(dir_path) {
  cat(sprintf("  Reading mtx from %s...\n", dir_path))
  mtx <- readMM(gzfile(file.path(dir_path, "matrix.mtx.gz")))
  # 10X CellRanger v3+: genes × barcodes (no transpose needed)
  mtx <- as(mtx, "CsparseMatrix")
  genes <- read.table(gzfile(file.path(dir_path, "features.tsv.gz")),
                       header = FALSE, stringsAsFactors = FALSE)[, 1]
  barcodes <- read.table(gzfile(file.path(dir_path, "barcodes.tsv.gz")),
                          header = FALSE, stringsAsFactors = FALSE)[, 1]
  # Must set dimnames for SoupX setClusters to work (requires named clusters)
  rownames(mtx) <- genes
  colnames(mtx) <- barcodes
  list(mat = mtx, genes = genes, barcodes = barcodes)
}

# ---------------------------------------------------------------
# Load matrices
# ---------------------------------------------------------------
raw <- read_10x_dir(RAW_DIR)
tod <- raw$mat
raw_genes <- raw$genes
cat(sprintf("  Raw: %d genes × %d droplets\n", nrow(tod), ncol(tod)))

filt <- read_10x_dir(FILT_DIR)
toc <- filt$mat
filt_genes <- filt$genes
cat(sprintf("  Filt: %d genes × %d cells\n", nrow(toc), ncol(toc)))

# Align genes by position (both should be in same order for first N genes)
common_genes <- intersect(raw_genes, filt_genes)
common_idx_raw <- match(common_genes, raw_genes)
common_idx_filt <- match(common_genes, filt_genes)
cat(sprintf("  Common genes: %d\n", length(common_genes)))
tod <- tod[common_idx_raw, , drop = FALSE]
toc <- toc[common_idx_filt, , drop = FALSE]

# ---------------------------------------------------------------
# Create SoupChannel + estimate soup
# ---------------------------------------------------------------
t_load <- Sys.time()
# Manually compute soup profile (estimateSoup fails on dgCMatrix subclasses)
droplet_sums <- Matrix::colSums(tod)
empty_mask <- droplet_sums >= 0 & droplet_sums <= 100
soup_counts <- Matrix::rowSums(tod[, empty_mask, drop = FALSE])
sc <- SoupChannel(as(tod, "TsparseMatrix"), toc, channelName = "M7")
sc$soupProfile <- data.frame(
    est = soup_counts / max(sum(soup_counts), 1),
    counts = soup_counts,
    row.names = NULL
)
cat(sprintf("  SoupChannel ready: %d genes × %d cells, %d empty droplets (%.1fs)\n",
    nrow(sc$toc), ncol(sc$toc), sum(empty_mask), difftime(t_load, t0, units = "secs")))

# ---------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------
set.seed(42)
clusters <- as.character(sample(1:16, ncol(sc$toc), replace = TRUE))
names(clusters) <- colnames(sc$toc)
sc <- setClusters(sc, clusters)

# ---------------------------------------------------------------
# autoEstCont
# ---------------------------------------------------------------
t2 <- Sys.time()
suppressMessages(sc <- autoEstCont(sc, tfidfMin = 0.2, forceAccept = TRUE))
rho <- sc$fit$rhoEst
t_auto <- difftime(Sys.time(), t2, units = "secs")
cat(sprintf("  autoEstCont: ρ=%.4f (%.1fs)\n", rho, t_auto))

# ---------------------------------------------------------------
# adjustCounts - subtraction
# ---------------------------------------------------------------
t3 <- Sys.time()
suppressMessages(out_sub <- adjustCounts(sc, method = "subtraction"))
t_sub <- difftime(Sys.time(), t3, units = "secs")
cat(sprintf("  adjustCounts (subtraction): %.1fs\n", t_sub))

# ---------------------------------------------------------------
# adjustCounts - soupOnly
# ---------------------------------------------------------------
t4 <- Sys.time()
suppressMessages(out_soup <- adjustCounts(sc, method = "soupOnly"))
t_soup <- difftime(Sys.time(), t4, units = "secs")
cat(sprintf("  adjustCounts (soupOnly): %.1fs\n", t_soup))

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
t_total <- difftime(Sys.time(), t0, units = "secs")
cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat(sprintf("  R SOUPX BENCHMARK SUMMARY (M7)\n"))
cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))
cat(sprintf("  Data:         %d cells × %d genes\n", ncol(sc$toc), nrow(sc$toc)))
cat(sprintf("  Droplets:     %d\n", ncol(sc$tod)))
cat(sprintf("  ρ (auto):     %.4f\n", rho))
cat(sprintf("  Load+soup:    %.1fs\n", difftime(t_load, t0, units = "secs")))
cat(sprintf("  autoEstCont:  %.1fs\n", t_auto))
cat(sprintf("  subtraction:  %.1fs\n", t_sub))
cat(sprintf("  soupOnly:     %.1fs\n", t_soup))
cat(sprintf("  Total:        %.1fs\n", t_total))
cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))
