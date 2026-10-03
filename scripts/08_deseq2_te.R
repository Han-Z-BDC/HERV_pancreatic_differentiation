#!/usr/bin/env Rscript
# Pseudobulk DESeq2 for TE subfamilies (input from 07_pseudobulk_te.py), one modality per call.
#
# Units are sample x ATAC_WNN_annotation group. In the early stages each group comes from one sample
# (hPSC only from D0, DE only from D4), so blocking on sample would remove the stage/cell-type
# contrasts this project asks about. Design is therefore ~ group, with the sample x group pseudobulks
# as replicates; stage and differentiation batch stay confounded (only D7/D9 have two batches).
#
# Size factors come from each unit's library size (RNA: gene UMIs; ATAC: all cut sites), scaled to a
# geometric mean of 1 -- not from the TE counts, which are dominated by old intronic repeats.
#
#   groups: each group vs the mean of all other groups (marker-style contrast)
#   stages: units summed per sample, ~ stage, each stage vs the mean of the other stages
# Calls: padj < 0.05 and log2FC > 1; subfamilies with mean count < 10 dropped.
# Usage: Rscript 08_deseq2_te.R <deseq2 dir for one modality: counts.csv, coldata.csv, subfamily_info.csv>

suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
d <- args[1]
out <- file.path(d, "results"); dir.create(out, showWarnings = FALSE)
counts <- as.matrix(read.csv(file.path(d, "counts.csv"), row.names = 1, check.names = FALSE))
cd <- read.csv(file.path(d, "coldata.csv"), row.names = 1, stringsAsFactors = TRUE)[colnames(counts), ]
info <- read.csv(file.path(d, "subfamily_info.csv"), row.names = 1)
counts <- counts[rowMeans(counts) >= 10, ]
cat(sprintf("%s: %d subfamilies, %d units, %d groups\n", basename(d), nrow(counts), ncol(counts), nlevels(cd$group)))

vs_rest <- function(dds, factor_name, levels_, tag) {
  rn <- resultsNames(dds)
  res_all <- list()
  for (g in levels_) {
    w <- setNames(rep(-1 / (length(levels_) - 1), length(levels_)), paste0(factor_name, levels_))
    w[paste0(factor_name, g)] <- 1
    v <- setNames(rep(0, length(rn)), rn); v[names(w)] <- w
    r <- as.data.frame(results(dds, contrast = v))
    r$subfamily <- rownames(r); r$family <- info[r$subfamily, "family"]; r$class <- info[r$subfamily, "class"]
    r$level <- g
    res_all[[g]] <- r
  }
  all <- do.call(rbind, res_all)
  write.csv(all, file.path(out, sprintf("%s_vs_rest_all.csv", tag)), row.names = FALSE)
  up <- all[!is.na(all$padj) & all$padj < 0.05 & all$log2FoldChange > 1, ]
  up <- up[order(up$level, -up$log2FoldChange), ]
  write.csv(up, file.path(out, sprintf("%s_vs_rest_up.csv", tag)), row.names = FALSE)
  erv <- up[up$class == "LTR", ]
  cat(sprintf("\n== %s: subfamilies up vs rest (padj<0.05, log2FC>1): %d total, %d ERV/LTR\n", tag, nrow(up), nrow(erv)))
  for (g in levels_) {
    e <- erv[erv$level == g, ]
    e <- e[e$baseMean > 20, ]
    if (nrow(e)) cat(sprintf("  %s (%d): %s\n", g, nrow(e),
      paste(sprintf("%s %.1f", head(e$subfamily, 12), head(e$log2FoldChange, 12)), collapse = ", ")))
  }
  invisible(all)
}

sf <- function(lib) { s <- lib / exp(mean(log(lib))); s }

# groups
dds <- DESeqDataSetFromMatrix(counts, cd, ~ 0 + group)
sizeFactors(dds) <- sf(cd$lib_size)
dds <- DESeq(dds, quiet = TRUE)
vs_rest(dds, "group", levels(cd$group), "group")
saveRDS(dds, file.path(out, "dds_group.rds"))

# stages (sum units per sample)
smp <- factor(cd$sample)
sc <- sapply(levels(smp), function(s) rowSums(counts[, smp == s, drop = FALSE]))
scd <- data.frame(row.names = levels(smp), stage = factor(sapply(strsplit(levels(smp), "_"), `[`, 2),
                  levels = c("D0", "D4", "D7", "D9", "D12", "D15", "D22", "D45")),
                  lib_size = tapply(cd$lib_size, smp, sum)[levels(smp)])
scd$stage <- droplevels(scd$stage)
dds2 <- DESeqDataSetFromMatrix(sc, scd, ~ 0 + stage)
sizeFactors(dds2) <- sf(scd$lib_size)
dds2 <- DESeq(dds2, quiet = TRUE)
vs_rest(dds2, "stage", levels(scd$stage), "stage")
norm <- counts(dds2, normalized = TRUE)
write.csv(round(norm, 2), file.path(out, "stage_normalized_counts.csv"))
cat("\ndone:", out, "\n")
