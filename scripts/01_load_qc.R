#!/usr/bin/env Rscript
## =============================================================================
## 01_load_qc.R — load 15 donors, unify gene IDs to symbols, per-donor QC
## Main steps:
##   1. Read 10x triplets (plain or .gz) for 10 donors (AK/IAK/MS)
##   2. Split GSE130973 raw combined matrix into NS1-NS5 by barcode suffix -1..-5
##   3. Per-donor QC: 500 <= nFeature <= 4000, nCount <= 8000, percent.mt < 10
##   4. Save per-donor raw-count Seurat objects (unintegrated) + QC table
## =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(Seurat); library(data.table) })

## ---- CONFIG -----------------------------------------------------------------
DATA_DIR <- "rawdata"                       # see data/README.md
OUT_DIR  <- "results/qc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## QC thresholds (identical for all 15 donors)
QC <- list(min_feat = 500, max_feat = 4000, max_count = 8000, max_mt = 10)

## ---- 10x reader with gene-symbol unification --------------------------------
## Duplicated symbols: keep the row with the highest total count.
load10x <- function(mtx, bc, feat) {
  f <- read.table(gzfile(feat), header = FALSE, sep = "\t",
                  comment.char = "", fill = TRUE, colClasses = "character")
  sym <- if (!is.null(f[[2]])) f[[2]] else f[[1]]
  bad <- is.na(sym) | sym %in% c("", ".")
  sym[bad] <- f[[1]][bad]
  bcv <- scan(gzfile(bc), what = "character", sep = "\t", quiet = TRUE)
  m   <- as(Matrix::readMM(mtx), "CsparseMatrix")
  if (ncol(m) == length(sym) && nrow(m) == length(bcv)) m <- Matrix::t(m)  # safety transpose
  dimnames(m) <- list(sym, bcv)
  rs   <- Matrix::rowSums(m)
  keep <- !duplicated(sym[order(sym, -rs)])
  list(m = m[keep, , drop = FALSE])
}

## ---- QC + Seurat object (raw counts only) -----------------------------------
qc_build <- function(m, id, group, dataset, sample) {
  genes  <- rownames(m)
  mt_idx <- grepl("^MT-", genes)
  nCount   <- Matrix::colSums(m)
  nFeature <- diff(m@p)
  pct_mt <- if (any(mt_idx)) 100 * Matrix::colSums(m[mt_idx, , drop = FALSE]) / nCount else rep(0, ncol(m))
  keep <- nFeature >= QC$min_feat & nFeature <= QC$max_feat &
          nCount  <= QC$max_count   & pct_mt < QC$max_mt
  obj <- CreateSeuratObject(m[, keep, drop = FALSE], project = id)
  obj$donor <- obj$manuscript_id <- id
  obj$group <- group; obj$batch <- obj$dataset <- dataset; obj$sample <- sample
  obj$percent.mt <- pct_mt[keep]
  obj
}

## ---- donor definitions (paths under DATA_DIR) --------------------------------
donors <- list(
  list(id = "AK01",  group = "AK",  dataset = "inhouse",
       mtx = file.path(DATA_DIR, "inhouse/AK/matrix.mtx"),
       bc  = file.path(DATA_DIR, "inhouse/AK/barcodes.tsv"),
       feat= file.path(DATA_DIR, "inhouse/AK/genes.tsv")),
  list(id = "AK02",  group = "AK",  dataset = "GSE220300",
       mtx = file.path(DATA_DIR, "extracted/GSE220300/GSM6797956_AC01_matrix.mtx.gz"),
       bc  = file.path(DATA_DIR, "extracted/GSE220300/GSM6797956_AC01_barcodes.tsv.gz"),
       feat= file.path(DATA_DIR, "extracted/GSE220300/GSM6797956_AC01_features.tsv.gz")),
  list(id = "AK03",  group = "AK",  dataset = "GSE220300",
       mtx = file.path(DATA_DIR, "extracted/GSE220300/GSM6797958_AC02_matrix.mtx.gz"),
       bc  = file.path(DATA_DIR, "extracted/GSE220300/GSM6797958_AC02_barcodes.tsv.gz"),
       feat= file.path(DATA_DIR, "extracted/GSE220300/GSM6797958_AC02_features.tsv.gz")),
  list(id = "IAK01", group = "IAK", dataset = "inhouse",
       mtx = file.path(DATA_DIR, "inhouse/IAK/matrix.mtx"),
       bc  = file.path(DATA_DIR, "inhouse/IAK/barcodes.tsv"),
       feat= file.path(DATA_DIR, "inhouse/IAK/features.tsv")),
  list(id = "IAK02", group = "IAK", dataset = "GSE220300",
       mtx = file.path(DATA_DIR, "extracted/GSE220300/GSM6797963_IC02_matrix.mtx.gz"),
       bc  = file.path(DATA_DIR, "extracted/GSE220300/GSM6797963_IC02_barcodes.tsv.gz"),
       feat= file.path(DATA_DIR, "extracted/GSE220300/GSM6797963_IC02_features.tsv.gz")),
  list(id = "IAK03", group = "IAK", dataset = "GSE220300",
       mtx = file.path(DATA_DIR, "extracted/GSE220300/GSM6797960_IC01_matrix.mtx.gz"),
       bc  = file.path(DATA_DIR, "extracted/GSE220300/GSM6797960_IC01_barcodes.tsv.gz"),
       feat= file.path(DATA_DIR, "extracted/GSE220300/GSM6797960_IC01_features.tsv.gz")),
  list(id = "MS01",  group = "MS",  dataset = "GSE163973",
       mtx = file.path(DATA_DIR, "extracted/GSE163973/GSM4994382_NS1_matrix/NF1_matrix/matrix.mtx.gz"),
       bc  = file.path(DATA_DIR, "extracted/GSE163973/GSM4994382_NS1_matrix/NF1_matrix/barcodes.tsv.gz"),
       feat= file.path(DATA_DIR, "extracted/GSE163973/GSM4994382_NS1_matrix/NF1_matrix/features.tsv.gz")),
  list(id = "MS02",  group = "MS",  dataset = "GSE163973",
       mtx = file.path(DATA_DIR, "extracted/GSE163973/GSM4994383_NS2_matrix/NF2_matrix/matrix.mtx.gz"),
       bc  = file.path(DATA_DIR, "extracted/GSE163973/GSM4994383_NS2_matrix/NF2_matrix/barcodes.tsv.gz"),
       feat= file.path(DATA_DIR, "extracted/GSE163973/GSM4994383_NS2_matrix/NF2_matrix/features.tsv.gz")),
  list(id = "MS03",  group = "MS",  dataset = "GSE163973",
       mtx = file.path(DATA_DIR, "extracted/GSE163973/GSM4994384_NS3_matrix/NF3_matrix/matrix.mtx.gz"),
       bc  = file.path(DATA_DIR, "extracted/GSE163973/GSM4994384_NS3_matrix/NF3_matrix/barcodes.tsv.gz"),
       feat= file.path(DATA_DIR, "extracted/GSE163973/GSM4994384_NS3_matrix/NF3_matrix/features.tsv.gz")),
  list(id = "MS04",  group = "MS",  dataset = "GSE220300",
       mtx = file.path(DATA_DIR, "extracted/GSE220300/GSM6797965_MS02_matrix.mtx.gz"),
       bc  = file.path(DATA_DIR, "extracted/GSE220300/GSM6797965_MS02_barcodes.tsv.gz"),
       feat= file.path(DATA_DIR, "extracted/GSE220300/GSM6797965_MS02_features.tsv.gz"))
)

## ---- load the 10 non-NS donors ------------------------------------------------
obj_list <- list()
for (d in donors) {
  m <- load10x(d$mtx, d$bc, d$feat)$m
  obj_list[[d$id]] <- qc_build(m, d$id, d$group, d$dataset, d$id)
  cat(sprintf("%s: %d cells pass QC\n", d$id, ncol(obj_list[[d$id]])))
  rm(m); invisible(gc())
}

## ---- GSE130973: raw combined matrix -> 5 donors by barcode suffix -------------
raw <- load10x(file.path(DATA_DIR, "inhouse/GSE130973_matrix_raw.mtx.gz"),
               file.path(DATA_DIR, "inhouse/GSE130973_barcodes_raw.tsv.gz"),
               file.path(DATA_DIR, "inhouse/GSE130973_genes_raw.tsv.gz"))$m
suf <- sub(".*(-[0-9]+)$", "\\1", colnames(raw))
ns_samples <- c("y1", "y2", "o1", "o2", "o3")
for (i in 1:5) {
  sub_m <- raw[, suf == paste0("-", i), drop = FALSE]
  id <- paste0("NS", i)
  obj_list[[id]] <- qc_build(sub_m, id, "NS", "GSE130973", ns_samples[i])
  cat(sprintf("%s: %d cells pass QC\n", id, ncol(obj_list[[id]])))
  rm(sub_m); invisible(gc())
}
rm(raw); invisible(gc())

## ---- save ----------------------------------------------------------------------
saveRDS(obj_list, file.path(OUT_DIR, "per_donor_qc.rds"))
qc_table <- data.frame(
  donor    = names(obj_list),
  group    = vapply(obj_list, function(o) o$group[1], ""),
  n_pass   = vapply(obj_list, ncol, integer(1)),
  med_feat = vapply(obj_list, function(o) median(o$nFeature_RNA), numeric(1)),
  med_umi  = vapply(obj_list, function(o) median(o$nCount_RNA), numeric(1)),
  med_mt   = vapply(obj_list, function(o) median(o$percent.mt), numeric(1)))
fwrite(qc_table, file.path(OUT_DIR, "qc_table.csv"))
cat("saved", file.path(OUT_DIR, "per_donor_qc.rds"), "and qc_table.csv\n")
