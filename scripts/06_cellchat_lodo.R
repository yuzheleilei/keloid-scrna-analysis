#!/usr/bin/env Rscript
## =============================================================================
## 06_cellchat_lodo.R — leave-one-donor-out CellChat (AK vs IAK per fold)
## Main steps:
##   1. Trim CellChatDB to the pathways of interest (IL1/CD40/CXCL/COLLAGEN/CD99/MIF)
##   2. Per fold x group: rebuild CellChat on the remaining donors
##      (triMean, population.size = TRUE, filterCommunication min.cells = 10)
##   3. Extract pathway-level communication probability (netP$prob, dim 3 = pathway)
##     + key ligand-receptor pair probabilities (net$prob, dim 3 = interaction)
##   4. Write per-fold tables (basis of the LODO robustness figures)
## =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(Seurat); library(CellChat) })
set.seed(42)

## ---- CONFIG -----------------------------------------------------------------
IN_DIR  <- "results/atlas"
OUT_DIR <- "results/robustness"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(IN_DIR, "atlas_annotated.rds"))
obj <- JoinLayers(obj)

KEL      <- c("AK01", "AK02", "AK03", "IAK01", "IAK02", "IAK03")
KEEP_PW  <- c("IL1", "CD40", "CXCL", "COLLAGEN", "CD99", "MIF")
KEY_PAIRS <- c("CXCL12_CXCR4", "COL1A1_ITGA2", "COL1A1_ITGB1", "COL1A1_ITGA2_ITGB1",
               "COL1A1_CD44", "IL1B_IL1R1", "CD40_CD40LG", "CD99_CD99", "MIF_CD74_CD44")

## trimmed database: only the pathways compared in the manuscript
db <- CellChatDB.human
db$interaction <- db$interaction[db$interaction$pathway_name %in% KEEP_PW, ]

run_fold <- function(group_label, holdout) {
  cells <- colnames(obj)[obj$group == group_label & obj$manuscript_id != holdout]
  sub   <- subset(obj, cells = cells)
  DefaultAssay(sub) <- "RNA"

  cc <- createCellChat(GetAssayData(sub, assay = "RNA", layer = "data"),
                       meta = sub[[]][, c("cell_type", "group")], group.by = "cell_type")
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", population.size = TRUE)

  ## key L-R pairs from net$prob BEFORE filtering (raw probabilities retained)
  np  <- cc@net$prob                                  # [source, target, interaction]
  pairs <- data.frame(
    holdout = holdout, group = group_label,
    pair = KEY_PAIRS,
    prob = vapply(KEY_PAIRS, function(p)
      if (p %in% dimnames(np)[[3]]) sum(np[, , p], na.rm = TRUE) else 0, numeric(1)))

  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  pw <- cc@netP$prob                                  # [source, target, pathway]
  pathways <- data.frame(
    holdout = holdout, group = group_label,
    pathway = KEEP_PW,
    prob  = vapply(KEEP_PW, function(p) {
      sel <- grepl(p, dimnames(pw)[[3]], ignore.case = TRUE)
      if (any(sel)) sum(pw[, , sel, drop = FALSE], na.rm = TRUE) else 0
    }, numeric(1)))
  list(pathways = pathways, pairs = pairs)
}

## run all 6 folds x 2 groups, collect pathway and pair tables
pw_all <- do.call(rbind, lapply(KEL, function(holdout)
  do.call(rbind, lapply(c("AK", "IAK"), function(g)
    tryCatch(run_fold(g, holdout)$pathways, error = function(e) NULL)))))
pr_all <- do.call(rbind, lapply(KEL, function(holdout)
  do.call(rbind, lapply(c("AK", "IAK"), function(g)
    tryCatch(run_fold(g, holdout)$pairs, error = function(e) NULL)))))

write.csv(pw_all, file.path(OUT_DIR, "cellchat_pathway_perfold.csv"), row.names = FALSE)
write.csv(pr_all, file.path(OUT_DIR, "cellchat_pair_perfold.csv"), row.names = FALSE)
cat("LODO CellChat tables written to", OUT_DIR, "\n")
