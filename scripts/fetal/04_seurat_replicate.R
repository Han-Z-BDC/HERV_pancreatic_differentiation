#!/usr/bin/env Rscript
# Reproduce the Olaniru et al. 2023 (Cell Metab; bioRxiv 2022.02.04.478971) scRNA-seq analysis of human
# fetal pancreas from our own Cell Ranger output, and write cell tables for the TE analysis.
#
# Paper methods (preprint): Cell Ranger 6.0.0 / GRCh38; HTODemux (Seurat, CLR, positive.quantile 0.99;
# pool 3 by genotype, Vireo); QC nUMI > 500, nGene > 250, log10GenesPerUMI > 0.75, mitoRatio < 0.25;
# SCTransform; IntegrateData; PCA 20 dims; Louvain clustering; Wilcoxon markers; Seurat 4.0.2.
# Expected: 24,080 cells (12 PCW 4,099; 13 7,168; 14 3,530; 15 1,664; 18 2,171; 19 3,827; 20 1,621).
#
# Here (Seurat 5.3): per pool, gene + hashtag matrices from $FETAL_OUT/cellranger/<pool>/outs/
#   1. HTODemux as in the paper; their published per-cell PCW calls (GEO *_demultiplexing_*.csv) are the
#      final assignment (pool 3 was assigned by genotype); the HTO x PCW crosstab gives the hashtag map
#   2. their QC thresholds
#   3. SCTransform per pool (vst.flavor v1, as Seurat 4), CCA integration across the 4 pools
#      (IntegrateLayers = FindIntegrationAnchors + IntegrateData), PCA 20, Louvain, UMAP
#   4. annotation: the paper's marker panels (module scores per cluster) -> major types; endocrine
#      cells sub-clustered and labelled (NEUROG3 progenitors, beta, alpha, delta, epsilon, polyhormonal)
# Outputs ($FETAL_OUT/seurat/): fetal_seurat.rds, cells.csv (barcode, pool, pcw, donor, QC, cluster,
#   cell_type, endocrine subtype, UMAP), cluster_markers.csv, hto_vs_geo_<pool>.csv, qc_vs_paper.csv,
#   whitelists/<pool>.txt (kept barcodes, for IRescue), umap_*.png
# Usage: Rscript scripts/fetal/04_seurat_replicate.R  (via 04_seurat_replicate.sb)

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2)
})
options(future.globals.maxSize = 32 * 1024^3)
set.seed(0)
FETAL_OUT <- Sys.getenv("FETAL_OUT"); REPO <- Sys.getenv("REPO")
out <- file.path(FETAL_OUT, "seurat"); dir.create(file.path(out, "whitelists"), recursive = TRUE, showWarnings = FALSE)
samples <- read.delim(file.path(REPO, "configs/fetal_olaniru_samples.tsv"))

# ---------------------------------------------------------------- 1. per pool: read, HTODemux, GEO labels
objs <- list(); qc_rows <- list()
for (i in seq_len(nrow(samples))) {
  pool <- samples$pool[i]
  m <- Read10X(file.path(FETAL_OUT, "cellranger", pool, "outs/filtered_feature_bc_matrix"))
  # pool4 was run gene-expression only (no usable hashtag library -- see 03_cellranger_count.sb),
  # so Read10X returns a bare matrix there rather than a list of feature types.
  has_hto <- is.list(m) && "Antibody Capture" %in% names(m)
  gex <- if (is.list(m)) m[["Gene Expression"]] else m
  o <- CreateSeuratObject(gex, project = pool)
  if (has_hto) {
    o[["HTO"]] <- CreateAssayObject(m[["Antibody Capture"]][, colnames(o)])
    o <- NormalizeData(o, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
    o <- HTODemux(o, assay = "HTO", positive.quantile = 0.99, verbose = FALSE)
  } else {
    o$hash.ID <- NA_character_
    o$HTO_classification.global <- NA_character_
  }
  geo <- read.csv(file.path(FETAL_OUT, "geo_suppl", samples$geo_demux[i]), check.names = FALSE)
  colnames(geo) <- c("barcode", "geo_label")
  o$geo_label <- geo$geo_label[match(colnames(o), geo$barcode)]
  if (has_hto) {
    tab <- as.data.frame.matrix(table(o$hash.ID, ifelse(is.na(o$geo_label), "not_in_GEO", o$geo_label)))
    write.csv(tab, file.path(out, sprintf("hto_vs_geo_%s.csv", pool)))
    cat(sprintf("\n== %s: %d cells called by Cell Ranger, %d in the GEO table\nHTODemux classes x GEO labels:\n",
                pool, ncol(o), sum(!is.na(o$geo_label))))
    print(tab)
  } else {
    cat(sprintf("\n== %s: %d cells called by Cell Ranger, %d in the GEO table (no HTO; GEO labels only)\n",
                pool, ncol(o), sum(!is.na(o$geo_label))))
    print(table(ifelse(is.na(o$geo_label), "not_in_GEO", o$geo_label)))
  }
  o$pool <- pool
  o$pcw <- ifelse(grepl("PCW$", o$geo_label), sub("PCW", "", o$geo_label), NA)
  # donor: PCW alone is ambiguous in pool 1 (two 13 PCW embryos); add the HTO when it is a singlet.
  # pool4's three donors have distinct PCW, so the GEO label alone identifies them.
  o$donor <- ifelse(is.na(o$pcw), NA, paste0(pool, "_", o$pcw, "PCW",
                    ifelse(!is.na(o$HTO_classification.global) & o$HTO_classification.global == "Singlet",
                           paste0("_", o$hash.ID), "")))
  # ---- QC as in the paper
  o$mitoRatio <- PercentageFeatureSet(o, pattern = "^MT-") / 100
  o$log10GenesPerUMI <- log10(o$nFeature_RNA) / log10(o$nCount_RNA)
  o$pass_qc <- o$nCount_RNA > 500 & o$nFeature_RNA > 250 & o$log10GenesPerUMI > 0.75 & o$mitoRatio < 0.25
  keep <- !is.na(o$pcw) & o$pass_qc
  qc_rows[[pool]] <- data.frame(pool = pool, called = ncol(o), in_geo = sum(!is.na(o$geo_label)),
                                assigned = sum(!is.na(o$pcw)), pass_qc = sum(o$pass_qc), kept = sum(keep))
  objs[[pool]] <- subset(o, cells = colnames(o)[keep])
  objs[[pool]] <- RenameCells(objs[[pool]], new.names = paste0(pool, "#", colnames(objs[[pool]])))
}
qc <- bind_rows(qc_rows); print(qc)

# ---------------------------------------------------------------- 2. merge, SCT per pool, integrate
obj <- merge(objs[[1]], objs[-1])
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- split(JoinLayers(obj[["RNA"]]), f = obj$pool)
obj <- SCTransform(obj, vst.flavor = "v1", method = "glmGamPoi", verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- IntegrateLayers(obj, method = CCAIntegration, normalization.method = "SCT",
                       orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:20, verbose = FALSE)
for (r in c(0.4, 0.8, 1.2)) obj <- FindClusters(obj, resolution = r, algorithm = 1, verbose = FALSE)
obj$cluster <- obj$SCT_snn_res.0.8
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:20, verbose = FALSE)

pcw_counts <- table(obj$pcw)
paper <- c("12" = 4099, "13" = 7168, "14" = 3530, "15" = 1664, "18" = 2171, "19" = 3827, "20" = 1621)
cmp <- data.frame(pcw = names(paper), paper = as.integer(paper), ours = as.integer(pcw_counts[names(paper)]))
cmp <- rbind(cmp, data.frame(pcw = "total", paper = sum(paper), ours = ncol(obj)))
write.csv(cmp, file.path(out, "qc_vs_paper.csv"), row.names = FALSE); write.csv(qc, file.path(out, "qc_per_pool.csv"), row.names = FALSE)
cat("\ncells per PCW, ours vs paper:\n"); print(cmp)

# ---------------------------------------------------------------- 3. annotation with the paper's markers
obj <- PrepSCTFindMarkers(obj, verbose = FALSE)
MAJOR <- list(
  Endocrine = c("INS", "GCG", "SST", "PAX6", "CHGA", "NEUROD1"),
  Acinar = c("CPA1", "PRSS1", "CLPS", "CTRB2"),
  Ductal = c("SLC4A4", "CFTR", "ANXA4", "SOX9", "KRT19"),
  Mesenchymal = c("COL3A1", "DCN", "VIM", "COL1A1"),
  Immune = c("RAC2", "LYZ", "TRAC", "PTPRC"),
  Endothelial = c("VWF", "ADGRL4", "ANGPT2", "PECAM1"),
  Schwann = c("CRYAB", "CDH19", "SOX10"),
  Erythroblast = c("HBB", "HBG2", "HBA1"))
obj <- AddModuleScore(obj, features = MAJOR, name = "major_", assay = "SCT")
sc_cols <- paste0("major_", seq_along(MAJOR))
cl_score <- aggregate(obj@meta.data[, sc_cols], list(cluster = obj$cluster), mean)
colnames(cl_score)[-1] <- names(MAJOR)
cl_score$cell_type <- names(MAJOR)[apply(cl_score[, -1], 1, which.max)]
obj$cell_type <- cl_score$cell_type[match(obj$cluster, cl_score$cluster)]
write.csv(cl_score, file.path(out, "cluster_major_scores.csv"), row.names = FALSE)
cat("\ncluster -> major type:\n"); print(cl_score[, c("cluster", "cell_type")])
mk <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
write.csv(mk, file.path(out, "cluster_markers.csv"), row.names = FALSE)

# ---- endocrine sub-clusters (paper: 3 NEUROG3 progenitor types, 4 beta, 2 delta, alpha, epsilon)
ENDO <- list(EP_NEUROG3 = c("NEUROG3", "NEUROD1", "INSM1", "HES6", "SOX4"), Beta = c("INS", "IAPP", "NKX6-1", "PDX1"),
             Alpha = c("GCG", "ARX", "IRX2", "TTR"), Delta = c("SST", "HHEX", "LEPR"), Epsilon = c("GHRL", "ACSL1"),
             PP_cell = c("PPY"))
en <- subset(obj, subset = cell_type == "Endocrine")
en <- FindNeighbors(en, reduction = "integrated.cca", dims = 1:20, verbose = FALSE)
en <- FindClusters(en, resolution = 1.0, verbose = FALSE)
en <- AddModuleScore(en, features = ENDO, name = "endo_", assay = "SCT")
es <- aggregate(en@meta.data[, paste0("endo_", seq_along(ENDO))], list(sub = en$seurat_clusters), mean)
colnames(es)[-1] <- names(ENDO)
hormones <- c("INS", "GCG", "SST", "GHRL")
avg <- AverageExpression(en, features = c(hormones, "NEUROG3"), group.by = "seurat_clusters", assays = "SCT")$SCT
# Seurat 5's AverageExpression prefixes numeric group names with "g" ("g0", "g1", ...). Strip it so
# the columns can be indexed by the cluster labels in es$sub -- otherwise the polyhormonal lookup
# below silently returns NA and the cbind() for endocrine_subclusters.csv goes out of bounds.
colnames(avg) <- sub("^g", "", colnames(avg))
stopifnot(all(as.character(es$sub) %in% colnames(avg)))
lab <- names(ENDO)[apply(es[, -1], 1, which.max)]
# polyhormonal: >= 2 hormones each at >= 30% of the top hormone in that sub-cluster
poly <- apply(avg[hormones, , drop = FALSE], 2, function(v) sum(v >= 0.3 * max(v)) >= 2 & max(v) > 1)
lab[poly[as.character(es$sub)] & lab != "EP_NEUROG3"] <- "Polyhormonal"
es$endo_type <- lab
en$endo_type <- es$endo_type[match(en$seurat_clusters, es$sub)]
obj$endo_type <- NA; obj$endo_type[colnames(en)] <- en$endo_type
obj$endo_subcluster <- NA; obj$endo_subcluster[colnames(en)] <- as.character(en$seurat_clusters)
write.csv(cbind(es, t(avg)[as.character(es$sub), ]), file.path(out, "endocrine_subclusters.csv"), row.names = FALSE)
cat("\nendocrine sub-clusters:\n"); print(es[, c("sub", "endo_type")])
obj$cell_type_fine <- ifelse(is.na(obj$endo_type), obj$cell_type, obj$endo_type)

# ---------------------------------------------------------------- 4. outputs
um <- Embeddings(obj, "umap")
cells <- data.frame(cell = colnames(obj), pool = obj$pool, barcode = sub(".*#", "", colnames(obj)), pcw = obj$pcw,
                    donor = obj$donor, hash = obj$hash.ID, nCount_RNA = obj$nCount_RNA, nFeature_RNA = obj$nFeature_RNA,
                    mitoRatio = obj$mitoRatio, cluster = obj$cluster, cell_type = obj$cell_type, endo_type = obj$endo_type,
                    endo_subcluster = obj$endo_subcluster, cell_type_fine = obj$cell_type_fine, UMAP_1 = um[, 1], UMAP_2 = um[, 2])
write.csv(cells, file.path(out, "cells.csv"), row.names = FALSE)
for (p in unique(cells$pool)) writeLines(cells$barcode[cells$pool == p], file.path(out, "whitelists", paste0(p, ".txt")))
cat("\ncell types:\n"); print(table(obj$cell_type_fine, obj$pcw))
saveRDS(obj, file.path(out, "fetal_seurat.rds"))
for (g in c("cell_type", "cell_type_fine", "pcw", "pool")) {
  p <- DimPlot(obj, group.by = g, label = g %in% c("cell_type", "cell_type_fine"), repel = TRUE, raster = TRUE) + ggtitle(g)
  ggsave(file.path(out, sprintf("umap_%s.png", g)), p, width = 8, height = 6.5, dpi = 110)
}
p <- FeaturePlot(obj, c("NEUROG3", "INS", "GCG", "SST", "GHRL", "CHGA", "CPA1", "SOX9", "CDH6"), order = TRUE, raster = TRUE, ncol = 3)
ggsave(file.path(out, "umap_markers.png"), p, width = 13, height = 12, dpi = 100)
cat("done", out, "\n")
