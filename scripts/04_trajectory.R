#!/usr/bin/env Rscript
## =============================================================================
## 04_trajectory_cellchat.R — developmental potency, pseudotime, full-data CellChat
## Main steps:
##   1. CytoTRACE2 on VEC and FB (developmental potency per subtype)
##   2. monocle3 on VEC: Seurat -> cds via new_cell_data_set() + reducedDims()
##      (Seurat UMAP imported; NOT as.cell_data_set()), learn_graph, order_cells
##      rooted at the highest-potency subtype; cross-lineage VEC+FB ordering
##   3. CellChat AK vs IAK (all 12 major cell types): triMean, population.size=TRUE
##   4. Pathway-level and key ligand-receptor-pair comparison AK vs IAK
## Key conventions (CellChat v2.1.2 / monocle3 1.3.1):
##   - pathway dimension of netP$prob is dimnames[[3]]
##   - computeCommunProbPathway() (not computePathway)
##   - learn_graph() has no reduction_method argument
## =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(monocle3); library(CellChat); library(CytoTRACE2)
})
set.seed(42)

## ---- CONFIG -----------------------------------------------------------------
IN_DIR  <- "results/atlas"
dir.create(file.path(IN_DIR, "monocle3"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(IN_DIR, "cellchat"), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(IN_DIR, "atlas_annotated.rds"))
lin <- readRDS(file.path(IN_DIR, "lineage_objects.rds"))
vec <- lin$vec; fb <- lin$fb

## ---- 1. CytoTRACE2 ------------------------------------------------------------
run_cyto <- function(sub, top_genes = 8000) {
  DefaultAssay(sub) <- "RNA"
  genes <- unique(c(head(VariableFeatures(sub), top_genes)))
  dat <- as(GetAssayData(sub, assay = "RNA", layer = "data")[genes, ], "dgCMatrix")
  cytotrace2(dat, species = "human")
}
cyto_vec <- run_cyto(vec)
cyto_fb  <- run_cyto(fb, top_genes = 6000)

attach_cyto <- function(sub, res) {
  df <- as.data.frame(res)
  sc <- df[[grep("CytoTRACE2_Score|score", colnames(df), ignore.case = TRUE, value = TRUE)[1]]]
  names(sc) <- rownames(df)
  sub$cytotrace2 <- as.numeric(sc[colnames(sub)])
  write.csv(data.frame(subtype = names(tapply(sub$cytotrace2, sub$subtype_named, mean)),
                       mean_potency = as.numeric(tapply(sub$cytotrace2, sub$subtype_named, mean))),
            file.path(IN_DIR, "cytotrace_summary.csv"))
  sub
}
vec <- attach_cyto(vec, cyto_vec)
fb  <- attach_cyto(fb,  cyto_fb)

## ---- 2. monocle3 (Seurat conversion via new_cell_data_set + reducedDims) -------
run_monocle <- function(sub) {
  DefaultAssay(sub) <- "RNA"
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  cds <- new_cell_data_set(counts,
                           cell_metadata = sub[[]],
                           gene_metadata = data.frame(gene_short_name = rownames(counts),
                                                      row.names = rownames(counts)))
  cds <- tryCatch(estimate_size_factors(cds), error = function(e) cds)
  reducedDims(cds)$UMAP <- Embeddings(sub, "umap")[colnames(cds), , drop = FALSE]
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, close_loop = FALSE)
  ## root at the highest-CytoTRACE2-potency subtype
  pot <- tapply(colData(cds)$cytotrace2, colData(cds)$subtype_named, mean, na.rm = TRUE)
  root <- head(colnames(cds)[colData(cds)$subtype_named == names(sort(pot, decreasing = TRUE))[1]], 100)
  order_cells(cds, root_cells = root)
}

cds_vec <- run_monocle(vec)
pt <- data.frame(cell = colnames(cds_vec), pseudotime = pseudotime(cds_vec),
                 subtype = colData(cds_vec)$subtype_named)
write.csv(pt, file.path(IN_DIR, "monocle3/vec_pseudotime.csv"), row.names = FALSE)
print(sort(tapply(pt$pseudotime, pt$subtype, median, na.rm = TRUE)))

## cross-lineage VEC + FB ordering (the putative VEC -> FB transition axis)
vf     <- subset(obj, cells = c(colnames(vec), colnames(fb)))
cds_vf <- run_monocle(vf)
lab <- ifelse(!is.na(colData(cds_vf)$vec_subtype), as.character(colData(cds_vf)$vec_subtype),
              as.character(colData(cds_vf)$fb_subtype))
write.csv(data.frame(cell = colnames(cds_vf), pseudotime = pseudotime(cds_vf), subtype = lab),
          file.path(IN_DIR, "monocle3/vecfb_pseudotime.csv"), row.names = FALSE)

## ---- 3. CellChat AK vs IAK (full data) -----------------------------------------
run_cellchat <- function(group_label) {
  sub <- subset(obj, cells = colnames(obj)[obj$group == group_label])
  DefaultAssay(sub) <- "RNA"
  cc <- createCellChat(GetAssayData(sub, assay = "RNA", layer = "data"),
                       meta = sub[[]][, c("cell_type", "group")], group.by = "cell_type")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}
cc_ak  <- run_cellchat("AK")
cc_iak <- run_cellchat("IAK")

## ---- 4. pathway + key pair comparison ------------------------------------------
extract_pw <- function(cc) {
  prob <- cc@netP$prob                          # [source, target, pathway] in CellChat 2.x
  data.frame(pathway = dimnames(prob)[[3]], prob = as.numeric(apply(prob, 3, sum)))
}
cmp <- merge(extract_pw(cc_ak), extract_pw(cc_iak), by = "pathway", suffixes = c("_AK", "_IAK"))
cmp$log2FC_AKvsIAK <- log2((cmp$prob_AK + 1e-6) / (cmp$prob_IAK + 1e-6))
write.csv(cmp[order(-abs(cmp$log2FC_AKvsIAK)), ],
          file.path(IN_DIR, "cellchat/pathway_comparison.csv"), row.names = FALSE)

key_pairs <- c("CXCL12_CXCR4", "COL1A1_ITGA2", "COL1A1_ITGB1", "IL1B_IL1R1", "CD40_CD40LG")
pairs <- do.call(rbind, Map(function(cc, g) {
  df <- subsetCommunication(cc)
  out <- df[df$interaction_name %in% key_pairs, c("source", "target", "interaction_name", "prob")]
  out$group <- g
  out
}, list(AK = cc_ak, IAK = cc_iak), c("AK", "IAK")))
write.csv(pairs, file.path(IN_DIR, "cellchat/key_interactions.csv"), row.names = FALSE)

saveRDS(list(AK = cc_ak, IAK = cc_iak), file.path(IN_DIR, "cellchat/cellchat_objects.rds"))
cat("saved trajectory + cellchat outputs under", IN_DIR, "\n")
