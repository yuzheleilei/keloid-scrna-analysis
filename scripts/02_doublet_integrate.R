#!/usr/bin/env Rscript
## =============================================================================
## 02_doublet_integrate.R — doublet tagging + batch integration
## Main steps:
##   1. Restrict NS1-5 to CellRanger-filtered barcodes (processing-level parity)
##   2. Per-donor DoubletFinder on raw counts: normalize -> HVG2000 -> PCA15 ->
##      paramSweep(PCs 1:10) -> pK by max BCmetric -> pN=0.25, rate = 0.008*n/1000
##      clamped to [0.02, 0.08]. Cells are TAGGED (pANN / DF.classification), not removed.
##   3. LogNormalize(1e4) + per-donor HVG(2000) -> merge 15 donors (JoinLayers)
##   4. Integration: (a) uncorrected PCA, (b) Harmony (batch = dataset), (c) Seurat CCA
##   5. Metrics: mean iLISI (perplexity 30, label = dataset) + silhouette (cluster proxy)
##   6. Save integrated objects + metrics table
## =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(DoubletFinder); library(harmony)
  library(lisi); library(cluster)
})
set.seed(42)

## ---- CONFIG -----------------------------------------------------------------
DATA_DIR <- "rawdata"                       # for the NS filtered-barcode list
OUT_DIR  <- "results/rebuild"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DONORS <- c("AK01", "AK02", "AK03", "IAK01", "IAK02", "IAK03",
            "MS01", "MS02", "MS03", "MS04", "NS1", "NS2", "NS3", "NS4", "NS5")

obj_list <- readRDS(file.path("results/qc", "per_donor_qc.rds"))[DONORS]

## ---- 1. NS restriction to CellRanger-filtered barcodes -----------------------
## NOTE: supply the filtered barcode list from GEO (CellRanger filtered_feature_bc_matrix)
ref_bc <- scan(gzfile(file.path(DATA_DIR, "GSE130973_barcodes_filtered.tsv.gz")),
               what = "character", sep = "\t", quiet = TRUE)
for (d in c("NS1", "NS2", "NS3", "NS4", "NS5")) {
  keep <- obj_list[[d]]$orig.barcode %in% ref_bc
  obj_list[[d]] <- obj_list[[d]][, keep]
}

## ---- 2. Per-donor DoubletFinder (tag only, no removal) -----------------------
for (d in DONORS) {
  obj <- obj_list[[d]]
  n   <- ncol(obj)
  obj <- NormalizeData(obj, scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 15, verbose = FALSE)

  sweep  <- paramSweep(obj, PCs = 1:10, sct = FALSE, num.cores = 1)
  bcmvn  <- find.pK(summarizeSweep(sweep, GT = FALSE))
  pK     <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  rate   <- min(0.08, max(0.02, 0.008 * n / 1000))
  nExp   <- round(n * rate)
  obj    <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp, sct = FALSE)

  ## standardize metadata column names to pANN_0.25_{pK} / DF.classifications_0.25_{pK}
  old <- grep("^pANN_0\\.25_", colnames(obj@meta.data), value = TRUE)
  obj[[paste0("pANN_0.25_", pK)]] <- obj[[old[1]]]
  cls <- grep("^DF\\.classifications_0\\.25_", colnames(obj@meta.data), value = TRUE)
  obj[[paste0("DF.classifications_0.25_", pK)]] <- obj[[cls[1]]]

  obj <- DietSeurat(obj, layers = c("counts", "data"))   # keep memory low; pANN kept in meta
  obj_list[[d]] <- obj
  cat(sprintf("%s: pK=%.3f rate=%.3f flagged %d doublets\n", d, pK, rate,
              sum(obj[[paste0("DF.classifications_0.25_", pK)]][[1]] == "Doublet")))
}
df_tab <- do.call(rbind, lapply(DONORS, function(d) {
  o <- obj_list[[d]]
  pann_col <- grep("^pANN_0\\.25_", colnames(o@meta.data), value = TRUE)
  data.frame(donor = d, n_cells = ncol(o), pANN_mean = mean(o[[pann_col]][[1]]))
}))
write.csv(df_tab, file.path(OUT_DIR, "doubletfinder_summary.csv"), row.names = FALSE)

## ---- 3. Normalize + HVG + merge ----------------------------------------------
for (d in DONORS) {
  obj_list[[d]] <- NormalizeData(obj_list[[d]], scale.factor = 1e4, verbose = FALSE)
  obj_list[[d]] <- FindVariableFeatures(obj_list[[d]], nfeatures = 2000, verbose = FALSE)
}
merged <- merge(obj_list[[1]], obj_list[-1], add.cell.ids = DONORS)
merged <- JoinLayers(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
rm(obj_list); invisible(gc())

## ---- 4a. Uncorrected ----------------------------------------------------------
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.8, algorithm = 1, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, seed.use = 42, verbose = FALSE)

## ---- 4b. Harmony (main integration) -------------------------------------------
harm <- RunHarmony(merged, group.by.vars = "batch", dims = 1:30, verbose = FALSE)
harm <- FindNeighbors(harm, reduction = "harmony", dims = 1:30, verbose = FALSE)
harm <- FindClusters(harm, resolution = 0.8, algorithm = 1, verbose = FALSE)
harm <- RunUMAP(harm, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)

## ---- 4c. Seurat CCA (comparison; anchors restricted to genes common to all 4 datasets)
spl <- SplitObject(merged, split.by = "batch")
spl <- lapply(spl, function(x) FindVariableFeatures(DietSeurat(x, layers = c("counts", "data")),
                                                    nfeatures = 2000, verbose = FALSE))
common <- Reduce(intersect, lapply(spl, rownames))
sif    <- SelectIntegrationFeatures(spl, nfeatures = 2000)
spl    <- lapply(spl, function(x) ScaleData(x, features = intersect(sif, common), verbose = FALSE))
anchors <- FindIntegrationAnchors(spl, anchor.features = intersect(sif, common), dims = 1:30, verbose = FALSE)
cca   <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(cca) <- "integrated"
cca   <- ScaleData(cca, verbose = FALSE)
cca   <- RunPCA(cca, npcs = 30, verbose = FALSE)
cca   <- FindNeighbors(cca, dims = 1:30, verbose = FALSE)
cca   <- FindClusters(cca, resolution = 0.8, algorithm = 1, verbose = FALSE)
cca   <- RunUMAP(cca, dims = 1:30, seed.use = 42, verbose = FALSE)

## ---- 5. Integration metrics ----------------------------------------------------
set.seed(42)
metric <- function(emb, clusters) {
  ilisi <- mean(compute_lisi(emb, data.frame(dataset = merged$dataset[rownames(emb)]),
                             label_colnames = "dataset", perplexity = 30)[, 1])
  idx   <- if (nrow(emb) > 20000) sample.int(nrow(emb), 20000) else seq_len(nrow(emb))
  sil   <- mean(silhouette(as.integer(clusters)[idx], dist(emb[idx, , drop = FALSE]))[, "sil_width"])
  c(iLISI = ilisi, silhouette = sil)
}
metrics <- rbind(
  data.frame(run = "uncorrected", t(metric(Embeddings(merged, "pca")[, 1:30], merged$seurat_clusters))),
  data.frame(run = "harmony",     t(metric(Embeddings(harm, "harmony")[, 1:30], harm$seurat_clusters))),
  data.frame(run = "cca",         t(metric(Embeddings(cca, "pca")[, 1:30],   cca$seurat_clusters))))
write.csv(metrics, file.path(OUT_DIR, "ilisi_silhouette.csv"), row.names = FALSE)
print(metrics)

## ---- 6. Save -------------------------------------------------------------------
saveRDS(list(merged = merged, harmony = harm, cca = cca),
        file.path(OUT_DIR, "integrated_objects.rds"))
cat("saved", file.path(OUT_DIR, "integrated_objects.rds"), "\n")
