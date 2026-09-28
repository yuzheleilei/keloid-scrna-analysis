#!/usr/bin/env Rscript
## =============================================================================
## 03_annotate_atlas.R — major cell types + lineage subclustering
## Main steps:
##   1. Canonical marker panel -> per-cluster z-scored lineage score -> cell_type
##   2. Subcluster lineages on the integrated assay (resolution grid, algorithm=1):
##      VEC -> 8 subclusters, FB -> 13, immune compartment -> 7
##   3. Name subclusters by curated marker sets (tip / postnp / vwffb / cd45p /
##      arterial / venous / capillary / lymphatic; FB ks_high ... ; immune ...)
##   4. Write labels back -> atlas_annotated.rds + subcluster_mapping.csv
## =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(Seurat); library(ggplot2) })

## ---- CONFIG -----------------------------------------------------------------
IN_DIR  <- "results/rebuild"
OUT_DIR <- "results/atlas"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

objs <- readRDS(file.path(IN_DIR, "integrated_objects.rds"))
obj  <- objs$cca
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

## ---- 1. major cell types by canonical markers --------------------------------
markers <- list(
  VEC     = c("PECAM1", "VWF", "CLDN5", "CDH5"),
  FB      = c("COL1A1", "DCN", "LUM", "COL3A1"),
  EPI     = c("KRT5", "KRT14", "KRT1", "KRT10", "KRT15"),
  PERI    = c("ACTA2", "MYH11", "CSPG4", "RGS5", "PDGFRB"),
  MEL     = c("MLANA", "TYR", "PMEL", "DCT"),
  SCH     = c("S100B", "PLP1", "SOX10", "MPZ"),
  T       = c("CD3D", "CD3E", "CD2", "TRAC"),
  NK      = c("NKG7", "GNLY", "KLRD1", "PRF1"),
  B       = c("MS4A1", "CD79A", "CD74", "CD19"),
  MonoMac = c("LST1", "C1QA", "C1QB", "FCGR3A", "LYZ", "S100A8"),
  DC      = c("FCER1A", "CLEC10A", "LAMP3", "CD1C"),
  MAST    = c("TPSAB1", "TPSB2", "CPA3", "MS4A2"))
present <- intersect(unlist(markers), rownames(obj))

avg   <- AverageExpression(obj, features = present, group.by = "seurat_clusters", verbose = FALSE)[[1]]
z     <- t(scale(t(as.matrix(avg))))                       # z-score each gene across clusters
score <- t(sapply(markers, function(mk) colMeans(z[intersect(mk, rownames(z)), , drop = FALSE])))
assign <- apply(score, 2, function(s) { i <- which.max(s); if (s[i] <= 0) "Unassigned" else names(score)[i] })
obj$cell_type <- unname(assign[paste0("g", obj$seurat_clusters)])
obj$cell_type[is.na(obj$cell_type)] <- "Unassigned"
print(table(obj$cell_type))

## ---- 2. lineage subclustering (resolution grid hits target cluster count) -----
subcluster <- function(obj, lineage, target_n, prefix, res_grid = c(0.2, 0.3, 0.5, 0.8, 1.2)) {
  cells <- colnames(obj)[obj$cell_type == lineage]
  if (length(cells) < 200) return(NULL)
  sub <- subset(obj, cells = cells)
  DefaultAssay(sub) <- "integrated"
  sub <- RunPCA(sub, npcs = 20, verbose = FALSE)
  best <- NULL; best_d <- Inf
  for (r in res_grid) {                       # pick resolution closest to target k
    s <- FindClusters(sub, resolution = r, algorithm = 1, verbose = FALSE)
    d <- abs(length(unique(s$seurat_clusters)) - target_n)
    if (d < best_d) { best_d <- d; best <- s }
    if (d == 0) break
  }
  best
}

## curated marker sets for subtype naming
vec_sets <- list(tip = c("ESM1","ANGPT2","KDR","DLL4"), postnp = c("POSTN","MKI67","TOP2A","CENPF"),
                 vwffb = c("VWF","PLVAP","CLDN5"), cd45p = c("PTPRC","CD45","LAPTM5"),
                 arterial = c("GJA5","EFNB2","CXCL12"), venous = c("ACKR1","VWF","SELE"),
                 capillary = c("RGCC","CA4","CD36"), lymph = c("PROX1","LYVE1","PDPN"))
fb_sets  <- list(ks_high = c("POSTN","ASPN","THBS2","COMP"), papillary = c("APCDD1","SLC1A3","WIF1"),
                 reticular = c("DCN","LUM","COL1A1"), myofibro = c("ACTA2","MYH11","TAGLN"),
                 inflammatory = c("CCL2","IL6","CXCL8"), cycling = c("MKI67","TOP2A"),
                 secretory = c("MMP1","MMP3","PTGS2"), adipofibro = c("APOE","PLIN2","ADIPOQ"))

name_by_markers <- function(sub, marker_sets) {
  av <- AverageExpression(sub, features = unique(unlist(marker_sets)),
                          group.by = "subtype", verbose = FALSE)[[1]]
  z <- t(scale(t(as.matrix(av))))
  vapply(colnames(z), function(g) {
    sc <- sapply(names(marker_sets), function(tp) mean(z[intersect(marker_sets[[tp]], rownames(z)), g]))
    paste0(g, "_", tolower(names(which.max(sc))))
  }, character(1))
}

## ---- 3. run subclustering + naming --------------------------------------------
vec <- subcluster(obj, "VEC", 8, "VEC")
vec$subtype  <- paste0("VEC", as.numeric(vec$seurat_clusters))
vec$subtype_named <- unname(name_by_markers(vec, vec_sets))

fb  <- subcluster(obj, "FB", 13, "FB")
fb$subtype   <- paste0("FB", as.numeric(fb$seurat_clusters))
fb$subtype_named  <- unname(name_by_markers(fb, fb_sets))

imm_cells <- colnames(obj)[obj$cell_type %in% c("T","NK","B","MonoMac","DC","MAST")]
sub_i <- subset(obj, cells = imm_cells)
DefaultAssay(sub_i) <- "integrated"
sub_i <- RunPCA(sub_i, npcs = 20, verbose = FALSE)
sub_i <- FindClusters(sub_i, resolution = 0.5, algorithm = 1, verbose = FALSE)
sub_i$immune_subtype <- paste0("IMM", as.numeric(sub_i$seurat_clusters))

## ---- 4. write labels back to the atlas object ---------------------------------
sm <- data.frame(cell = colnames(obj), cell_type = obj$cell_type,
                 vec_subtype = NA, fb_subtype = NA, immune_subtype = NA)
sm$vec_subtype[match(colnames(vec), sm$cell)]   <- vec$subtype_named
sm$fb_subtype[match(colnames(fb), sm$cell)]     <- fb$subtype_named
sm$immune_subtype[match(colnames(sub_i), sm$cell)] <- sub_i$immune_subtype
write.csv(sm, file.path(OUT_DIR, "subcluster_mapping.csv"), row.names = FALSE)

obj$vec_subtype    <- sm$vec_subtype[match(colnames(obj), sm$cell)]
obj$fb_subtype     <- sm$fb_subtype[match(colnames(obj), sm$cell)]
obj$immune_subtype <- sm$immune_subtype[match(colnames(obj), sm$cell)]
saveRDS(obj, file.path(OUT_DIR, "atlas_annotated.rds"))
saveRDS(list(vec = vec, fb = fb, immune = sub_i), file.path(OUT_DIR, "lineage_objects.rds"))
cat("saved", file.path(OUT_DIR, "atlas_annotated.rds"), "\n")
