#!/usr/bin/env Rscript
## =============================================================================
## 05_robustness_lodo.R — donor-level robustness analyses
## Main steps:
##   1. Donor-level proportions of key subtypes within their compartment
##      (POSTN+P / VWF+FB / CD45+P / Tip(arterial proxy) in VEC; FB_KS_high in FB)
##   2. Leave-one-donor-out (LODO, 6 folds): direction (AK vs IAK) per fold
##   3. Donor-level bootstrap (1000x, 3 donors resampled with replacement/group) -> 95% CI
##   4. Jackknife (leave-one-out bias / SE)
##   5. Pseudobulk edgeR (donor x subtype sums, IAK vs AK)
##   6. LODO monocle3 state preference (rooted at the CD45+P subset)
## Replicate unit = donor (cells are nested within donors, not independent).
## =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(Seurat); library(edgeR) })
set.seed(42)

## ---- CONFIG -----------------------------------------------------------------
IN_DIR  <- "results/atlas"
OUT_DIR <- "results/robustness"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(IN_DIR, "atlas_annotated.rds"))
obj <- JoinLayers(obj)

KEL <- c("AK01", "AK02", "AK03", "IAK01", "IAK02", "IAK03")
KEY <- list(
  postnp = list(col = "vec_subtype", pat = "postnp",   label = "POSTN+P",             comp = "VEC"),
  vwffb  = list(col = "vec_subtype", pat = "vwffb",    label = "VWF+FB",              comp = "VEC"),
  cd45p  = list(col = "vec_subtype", pat = "cd45p",    label = "CD45+P",              comp = "VEC"),
  tip    = list(col = "vec_subtype", pat = "arterial", label = "Tip(proxy=arterial)", comp = "VEC"),
  kshigh = list(col = "fb_subtype",  pat = "ks_high",  label = "FB_KS_high",          comp = "FB"))

## ---- 1. donor-level proportions -----------------------------------------------
donor_props <- function(o, key) {
  comp <- colnames(o)[!is.na(o$cell_type) & o$cell_type == key$comp]
  sub  <- o[[]][comp, , drop = FALSE]
  hit  <- !is.na(sub[[key$col]]) & grepl(key$pat, sub[[key$col]], ignore.case = TRUE)
  d <- data.frame(donor = names(table(sub$manuscript_id)),
                  prop  = as.numeric(tapply(hit, sub$manuscript_id, sum)) / as.numeric(table(sub$manuscript_id)))
  d$group <- ifelse(grepl("^AK", d$donor), "AK", "IAK")
  d
}

## full-cohort direction (donor-level Wilcoxon, 3 v 3)
full <- do.call(rbind, lapply(names(KEY), function(k) {
  dp <- donor_props(obj, KEY[[k]])
  w  <- suppressWarnings(wilcox.test(prop ~ group, dp, exact = FALSE))
  data.frame(key = KEY[[k]]$label, direction = ifelse(mean(dp$prop[dp$group=="AK"]) > mean(dp$prop[dp$group=="IAK"]),
                                                      "AK>IAK", "IAK>AK"), p = w$p.value)
}))
print(full)
write.csv(full, file.path(OUT_DIR, "full_cohort_direction.csv"), row.names = FALSE)

## ---- 2. LODO 6-fold: is the direction preserved without each donor? ------------
lodo <- do.call(rbind, lapply(KEL, function(holdout) {
  keep <- KEL[KEL != holdout]
  do.call(rbind, lapply(names(KEY), function(k) {
    dp <- donor_props(obj, KEY[[k]]); dp <- dp[dp$donor %in% keep, ]
    ak <- dp$prop[dp$group == "AK"]; iak <- dp$prop[dp$group == "IAK"]
    data.frame(holdout = holdout, key = KEY[[k]]$label,
               direction = ifelse(mean(ak) > mean(iak), "AK>IAK", "IAK>AK"),
               p = suppressWarnings(wilcox.test(ak, iak, exact = FALSE)$p.value))
  }))
}))
write.csv(lodo, file.path(OUT_DIR, "lodo_direction_perfold.csv"), row.names = FALSE)

## ---- 3. donor bootstrap (1000x) -------------------------------------------------
boot <- do.call(rbind, lapply(names(KEY), function(k) {
  dp  <- donor_props(obj, KEY[[k]])
  ak  <- dp$prop[dp$group == "AK"]; iak <- dp$prop[dp$group == "IAK"]
  d   <- replicate(1000, mean(sample(iak, replace = TRUE)) - mean(sample(ak, replace = TRUE)))
  data.frame(key = KEY[[k]]$label, delta_mean = mean(iak) - mean(ak),
             ci_low = quantile(d, 0.025), ci_high = quantile(d, 0.975))
}))
write.csv(boot, file.path(OUT_DIR, "bootstrap_delta_ci.csv"), row.names = FALSE)

## ---- 4. jackknife ----------------------------------------------------------------
jack <- do.call(rbind, lapply(names(KEY), function(k) {
  dp <- donor_props(obj, KEY[[k]]); p <- setNames(dp$prop, dp$donor)
  d  <- vapply(KEL, function(h) mean(p[grepl("^IAK", setdiff(KEL, h))]) - mean(p[grepl("^AK", setdiff(KEL, h))]), numeric(1))
  data.frame(key = KEY[[k]]$label, jackknife_mean = mean(d), jackknife_se = sd(d))
}))
write.csv(jack, file.path(OUT_DIR, "jackknife_delta.csv"), row.names = FALSE)

## ---- 5. pseudobulk edgeR (donor-level aggregate counts, IAK vs AK) --------------
pseudobulk_de <- function(key) {
  comp    <- colnames(obj)[obj$cell_type == key$comp]
  subtype <- obj[[key$col]][comp]
  hit     <- !is.na(subtype) & grepl(key$pat, subtype, ignore.case = TRUE)
  counts  <- GetAssayData(obj, assay = "RNA", layer = "counts")[, comp]
  agg     <- t(rowsum(t(counts[, hit, drop = FALSE]), group = obj$manuscript_id[comp][hit]))
  meta    <- data.frame(group = ifelse(grepl("^AK", colnames(agg)), "AK", "IAK"),
                        row.names = colnames(agg))               # donors as replicates
  design <- model.matrix(~ group, meta)
  y <- DGEList(agg)
  y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  topTags(glmQLFTest(fit, coef = "groupIAK"), n = Inf)$table
}
de_kshigh <- pseudobulk_de(KEY$kshigh)          # e.g. FB_KS_high: ECM/collagen genes up in IAK
write.csv(de_kshigh, file.path(OUT_DIR, "pseudobulk_edger_FKS_high_IAK_vs_AK.csv"))

## ---- 6. LODO monocle3 state preference (rooted at CD45+P subset) -----------------
## per fold: rebuild cds on VEC, root at cd45p cells, compare median pseudotime of
## Tip-like and POSTN+P cells between remaining AK and IAK donors (Wilcoxon on cells).
library(monocle3)
lodo_pt <- do.call(rbind, lapply(KEL, function(holdout) {
  cells <- colnames(obj)[obj$group %in% c("AK", "IAK") & obj$manuscript_id != holdout &
                         obj$cell_type == "VEC"]
  sub <- subset(obj, cells = cells)
  DefaultAssay(sub) <- "RNA"
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  cds <- new_cell_data_set(counts, cell_metadata = sub[[]],
                           gene_metadata = data.frame(gene_short_name = rownames(counts),
                                                      row.names = rownames(counts)))
  cds <- tryCatch(estimate_size_factors(cds), error = function(e) cds)
  reducedDims(cds)$UMAP <- Embeddings(sub, "umap")[colnames(cds), , drop = FALSE]
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, close_loop = FALSE)
  root <- head(colnames(cds)[grepl("cd45p", colData(cds)$subtype_named, ignore.case = TRUE)], 100)
  if (length(root) == 0) return(NULL)
  cds <- order_cells(cds, root_cells = root)
  pt  <- pseudotime(cds)
  do.call(rbind, lapply(c("arterial", "postnp"), function(pat) {
    sel <- grepl(pat, colData(cds)$subtype_named, ignore.case = TRUE)
    df <- data.frame(pt = pt[sel], subtype = pat,
                     group = colData(cds)$group[sel], donor = colData(cds)$manuscript_id[sel])
    ak <- df$pt[df$group == "AK"]; iak <- df$pt[df$group == "IAK"]
    if (length(ak) < 10 || length(iak) < 10) return(NULL)
    data.frame(holdout = holdout, subtype = pat,
               median_pt_AK = median(ak), median_pt_IAK = median(iak),
               direction = ifelse(median(ak) > median(iak), "AK_higher_pt", "IAK_higher_pt"),
               p = suppressWarnings(wilcox.test(ak, iak, exact = FALSE)$p.value),
               n_ak = length(ak), n_iak = length(iak))
  }))
}))
write.csv(lodo_pt, file.path(OUT_DIR, "lodo_monocle_state_preference.csv"), row.names = FALSE)
cat("robustness outputs written to", OUT_DIR, "\n")
