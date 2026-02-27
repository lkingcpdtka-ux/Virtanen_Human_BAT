## =========================================================
## Part 1 — Main Analysis: Virtanen Human BAT RNA-seq
## =========================================================
## Goal : Determine whether CLDN1 is increased in human BAT
##        relative to paired subcutaneous WAT
## Data : GSE113764 (Virtanen et al. Cell Metabolism 2018)
## Design: Paired — 14 subjects × 2 tissues (BAT / WAT)
## Model : ~ subject + tissue  (paired DESeq2)
## =========================================================

cat("=== Part 1: Virtanen Human BAT — CLDN1 validation ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


is_gzip_file <- function(path) {
  if (!file.exists(path)) return(FALSE)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  sig <- readBin(con, what = "raw", n = 2)
  length(sig) == 2 && identical(as.integer(sig), c(31L, 139L))
}

is_zip_file <- function(path) {
  if (!file.exists(path)) return(FALSE)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  sig <- readBin(con, what = "raw", n = 4)
  length(sig) == 4 && identical(as.integer(sig), c(80L, 75L, 3L, 4L))
}

prepare_text_input_file <- function(path, force_gzip = FALSE) {
  raw_size <- file.info(path)$size
  if (is.na(raw_size) || raw_size <= 0) {
    stop("Count file is empty or unreadable: ", path)
  }

  if (is_zip_file(path)) {
    zfiles <- unzip(path, list = TRUE)
    if (nrow(zfiles) == 0) stop("ZIP file contains no entries: ", path)
    entry <- zfiles$Name[1]
    tmp_zip_dir <- tempfile("unzipped_counts_")
    dir.create(tmp_zip_dir, recursive = TRUE, showWarnings = FALSE)
    unzip(path, files = entry, exdir = tmp_zip_dir)
    return(file.path(tmp_zip_dir, entry))
  }

  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  raw_bytes <- readBin(con, what = "raw", n = raw_size)

  use_gzip <- isTRUE(force_gzip) || is_gzip_file(path)
  if (use_gzip) {
    raw_bytes <- tryCatch(memDecompress(raw_bytes, type = "gzip"),
                          error = function(e) {
                            stop("Could not decompress gzip count file: ", path,
                                 " | ", conditionMessage(e))
                          })
  }

  tmp <- tempfile(fileext = ".txt")
  writeBin(raw_bytes, tmp)
  return(tmp)
}

read_geo_counts_with_fallback <- function(path, is_gz = FALSE) {
  txt_path <- prepare_text_input_file(path, force_gzip = is_gz)
  on.exit(unlink(txt_path, recursive = TRUE, force = TRUE), add = TRUE)

  encodings <- c("UTF-8", "latin1", "CP1252", "UTF-16LE", "UTF-16BE", "unknown")
  seps <- c("	", ",", ";")
  skip_options <- c(0, 1, 2, 3, 4, 5, 10, 20)

  for (enc in encodings) {
    for (sep in seps) {
      for (skip_n in skip_options) {
        dat <- tryCatch(
          suppressWarnings(
            read.table(
              txt_path,
              header = TRUE,
              sep = sep,
              quote = "\"",
              comment.char = "",
              fill = TRUE,
              check.names = FALSE,
              stringsAsFactors = FALSE,
              fileEncoding = enc,
              skip = skip_n
            )
          ),
          error = function(e) e
        )

        if (!inherits(dat, "error") && ncol(dat) >= 2 && nrow(dat) > 0) {
          return(as.data.frame(dat, stringsAsFactors = FALSE))
        }
      }
    }
  }

  stop("Unable to read GEO counts file with delimiter/encoding fallbacks: ", basename(path),
       ". Please provide a plain text count table via LOCAL_COUNTS_PATH.")
}

verify_output_file <- function(path, label = "output", must_exist = TRUE) {
  ok <- file.exists(path)
  if (ok) {
    size_bytes <- file.size(path)
    cat(sprintf("[SAVE-CHECK] %-24s OK  %s (%.1f KB)\n",
                label, basename(path), as.numeric(size_bytes) / 1024))
    return(invisible(TRUE))
  }

  msg <- paste0("[SAVE-CHECK] ", label, " MISSING: ", path)
  if (must_exist) stop(msg) else warning(msg)
  invisible(FALSE)
}

print_run_structure_sanity <- function(outdir) {
  cat("\n[SANITY] savepoint structure check\n")
  expected <- c("tables", "plots", "logs", "cache")
  for (d in expected) {
    dpath <- file.path(outdir, d)
    cat(sprintf("  - %-6s : %s\n", d, ifelse(dir.exists(dpath), "present", "MISSING")))
  }
}

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

coerce_numeric_like <- function(x) {
  vals <- as.character(x)
  vals <- gsub(",", "", vals, fixed = TRUE)
  vals <- gsub("[^0-9eE+.-]", "", vals)
  suppressWarnings(as.numeric(vals))
}

select_numeric_count_columns <- function(df, min_numeric_fraction = 0.5) {
  if (ncol(df) == 0) return(df)

  cn <- colnames(df)
  cn[is.na(cn)] <- ""

  ## Prefer explicit sample-like columns if present (Virtanen BAT/WAT columns)
  sample_like <- grepl("(BAT|WAT)$", cn, ignore.case = TRUE)
  if (sum(sample_like) >= 6) {
    out <- df[, sample_like, drop = FALSE]
    for (j in seq_len(ncol(out))) out[[j]] <- coerce_numeric_like(out[[j]])

    non_empty <- colSums(!is.na(out)) > 0
    out <- out[, non_empty, drop = FALSE]

    if (ncol(out) > 0) {
      dropped <- cn[!sample_like]
      if (length(dropped) > 0) {
        cat("[SANITY] Dropping non-sample columns:", paste(dropped, collapse = ", "), "\n")
      }
      return(out)
    }
  }

  numeric_fraction <- vapply(df, function(col) {
    vals <- coerce_numeric_like(col)
    mean(!is.na(vals), na.rm = TRUE)
  }, numeric(1))

  numeric_fraction[is.na(numeric_fraction)] <- 0
  keep <- numeric_fraction >= min_numeric_fraction

  if (!any(keep)) {
    ranked <- sort(numeric_fraction, decreasing = TRUE)
    top_names <- names(head(ranked, 10))
    msg <- paste0("No numeric-like columns detected in count table after removing gene column. ",
                  "Top candidate columns by numeric fraction: ",
                  paste(sprintf("%s=%.2f", top_names, head(ranked, 10)), collapse = ", "))
    stop(msg)
  }

  dropped <- names(df)[!keep]
  if (length(dropped) > 0) {
    cat("[SANITY] Dropping non-count columns:", paste(dropped, collapse = ", "), "\n")
  }

  out <- df[, keep, drop = FALSE]
  for (j in seq_len(ncol(out))) out[[j]] <- coerce_numeric_like(out[[j]])

  non_empty <- colSums(!is.na(out)) > 0
  out <- out[, non_empty, drop = FALSE]
  if (ncol(out) == 0) {
    stop("Count columns were detected but all converted values are NA after numeric coercion")
  }

  out
}

find_local_count_file <- function(geo_accession, workdir) {
  if (exists("LOCAL_COUNTS_PATH") &&
      !is.null(LOCAL_COUNTS_PATH) &&
      !is.na(LOCAL_COUNTS_PATH) &&
      nzchar(LOCAL_COUNTS_PATH)) {
    explicit <- if (file.exists(LOCAL_COUNTS_PATH)) LOCAL_COUNTS_PATH else file.path(workdir, LOCAL_COUNTS_PATH)
    if (file.exists(explicit)) {
      cat("[INFO] Using LOCAL_COUNTS_PATH:", normalizePath(explicit, mustWork = FALSE), "\n")
      return(explicit)
    }
    warning("[WARN] LOCAL_COUNTS_PATH was set but not found: ", LOCAL_COUNTS_PATH)
  }

  candidates <- list.files(
    workdir,
    pattern = paste0("^", geo_accession, ".*(humanBATWAT|count|raw).*(txt|csv)(\\.gz)?$"),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(candidates) > 0) {
    cat("[INFO] Found local count file candidate:", basename(candidates[1]), "\n")
    return(candidates[1])
  }

  return(NULL)
}

## ---- 1) Packages -----------------------------------------
required_pkgs <- c(
  "DESeq2", "edgeR", "GEOquery", "SummarizedExperiment",
  "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
  "EnhancedVolcano", "org.Hs.eg.db", "AnnotationDbi",
  "dplyr", "tidyr", "tibble", "stringr"
)

## BiocManager install if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
cat("[OK] All packages loaded\n")

## ---- 1.5) Source parameters ------------------------------
source(file.path(getwd(), "parameters.R"))

## ---- 2) save_core utilities ------------------------------
utils_dir <- file.path(getwd(), "save_core")
if (file.exists(file.path(utils_dir, "save.R"))) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))
  source(file.path(utils_dir, "purge.R"))
  source(file.path(utils_dir, "dedupe.R"))
  source(file.path(utils_dir, "cache_utils.R"))
  use_save_core <- TRUE
  cat("[OK] save_core utilities loaded\n")
} else {
  use_save_core <- FALSE
  cat("[WARN] save_core not found at: ", utils_dir, "\n", sep = "")
  cat("[WARN] Creating minimal output directory structure manually...\n")

  init_run <- function(script_name, species, data_type, keywords, notes, ...) {
    run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
    outdir  <- file.path(getwd(), "savepoints", paste0("RUN_", run_tag))
    dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "plots"),  recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "logs"),   recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "cache"),  recursive = TRUE, showWarnings = FALSE)
    return(list(outdir = outdir, run_tag = run_tag,
                cache_dir = file.path(outdir, "cache")))
  }
}

## ---- 2.5) Run metadata -----------------------------------
run_info <- list(
  script_name = "Part1_main_analysis.R",
  species     = SPECIES,
  data_type   = DATA_TYPE,
  keywords    = c("CLDN1", "Claudin-1", "BAT", "WAT", "human",
                   "brown adipose", "Virtanen", "GSE113764", "DESeq2"),
  notes       = "Primary goal: test whether CLDN1 is elevated in human BAT vs WAT",
  message     = "DESeq2 paired design + VST QC + CLDN1 focused validation"
)

## ---- 3) Initialise run -----------------------------------
run_ctx <- init_run(
  script_name   = run_info$script_name,
  species       = run_info$species,
  data_type     = run_info$data_type,
  keywords      = run_info$keywords,
  notes         = run_info$notes,
  message       = run_info$message,
  mode          = SAVECORE_MODE,
  run_tag       = SAVECORE_RUN_TAG,
  savepoint_dir = file.path(getwd(), "savepoints")
)

outdir    <- run_ctx$outdir
run_tag   <- run_ctx$run_tag
cache_dir <- run_ctx$cache_dir
cat("Run directory:", normalizePath(outdir, mustWork = FALSE), "\n")
cat("[SANITY] save_core mode:", SAVECORE_MODE, "| run_tag override:", ifelse(is.null(SAVECORE_RUN_TAG), "NULL", SAVECORE_RUN_TAG), "\n")
print_run_structure_sanity(outdir)
cat("\n")

## ---- 3.5) Analysis flags ---------------------------------
set.seed(SEED)

## Sanity-check tracker
sanity <- list()

## ==========================================================
## 4) MAIN ANALYSIS — wrapped in tryCatch
## ==========================================================
tryCatch({

  ## ---- 4.1) Download / load GEO data ---------------------
  cat("== 4.1  Fetching GEO data:", GEO_ACCESSION, "==\n")

  geo_cache <- file.path(cache_dir, paste0(GEO_ACCESSION, "_SE.rds"))

  if (file.exists(geo_cache)) {
    cat("[CACHE] Loading cached GEO data\n")
    se <- readRDS(geo_cache)
  } else {
    ## Download the supplementary count matrix from GEO
    gse <- getGEO(GEO_ACCESSION, GSEMatrix = TRUE, getGPL = FALSE)
    gse <- gse[[1]]

    ## Try local counts first; otherwise download GEO supplementary file
    supp_dir <- file.path(cache_dir, "geo_supp")
    dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

    local_count_file <- find_local_count_file(GEO_ACCESSION, getwd())
    count_file <- NULL

    if (!is.null(local_count_file)) {
      count_file <- local_count_file
    } else {
      supp_files <- getGEOSuppFiles(GEO_ACCESSION, baseDir = supp_dir,
                                     makeDirectory = FALSE)

      count_candidates <- rownames(supp_files)[grepl("count|raw|humanBATWAT|txt|csv",
                                                     rownames(supp_files),
                                                     ignore.case = TRUE)]
      if (length(count_candidates) > 0) {
        count_file <- count_candidates[1]
      }
    }

    if (!is.null(count_file)) {
      cat("[INFO] Reading count file:", basename(count_file), "\n")
      counts_raw <- read_geo_counts_with_fallback(
        count_file,
        is_gz = isTRUE(grepl("\\.gz$", count_file, ignore.case = TRUE))
      )

      if (ncol(counts_raw) < 2) {
        stop("Count file appears malformed (fewer than 2 columns): ", count_file,
             ". If this is a gzipped file without .gz extension, recompress/rename or set LOCAL_COUNTS_PATH correctly.")
      }

      ## First column is gene identifier; remaining columns should be counts
      gene_col <- counts_raw[[1]]
      counts_raw <- counts_raw[, -1, drop = FALSE]
      counts_raw <- select_numeric_count_columns(as.data.frame(counts_raw))

      ## Handle duplicate gene IDs by summing counts
      if (anyDuplicated(gene_col)) {
        n_dup <- sum(duplicated(gene_col))
        cat("[INFO] Found", n_dup, "duplicate gene names — aggregating by sum\n")
        counts_raw$..gene.. <- gene_col
        counts_raw <- aggregate(. ~ ..gene.., data = counts_raw, FUN = sum)
        gene_col   <- counts_raw$..gene..
        counts_raw <- counts_raw[, colnames(counts_raw) != "..gene..", drop = FALSE]
      }

      rownames(counts_raw) <- gene_col
      cat("[SANITY] Count matrix dimensions after cleanup:", nrow(counts_raw), "x", ncol(counts_raw), "\n")

      ## ---- Sanity block A: Characterise raw count matrix ----
      cat("\n--- SANITY BLOCK A: Raw count matrix diagnostics ---\n")

      ## A1. Library sizes
      lib_sizes_raw <- colSums(counts_raw, na.rm = TRUE)
      cat("[A1] Library sizes (total counts per sample):\n")
      cat("     min =", format(min(lib_sizes_raw), big.mark = ","),
          " median =", format(median(lib_sizes_raw), big.mark = ","),
          " max =", format(max(lib_sizes_raw), big.mark = ","), "\n")
      if (any(lib_sizes_raw < 1e5)) {
        cat("[A1] WARNING: Some samples have very low library sizes (<100K).\n")
        cat("     Low samples:", paste(names(which(lib_sizes_raw < 1e5)), collapse = ", "), "\n")
      }

      ## A2. Gene ID classification
      gene_ids <- rownames(counts_raw)
      n_total_genes <- length(gene_ids)
      n_hgnc_like   <- sum(grepl("^[A-Z][A-Z0-9-]*$", gene_ids) &
                            nchar(gene_ids) <= 15 &
                            !grepl("^(AK|BC|AJ|AF|NM_|NR_|XM_|XR_)", gene_ids))
      n_genbank     <- sum(grepl("^(AK|BC|AJ|AF|AB|AL|CR|BX)[0-9]{5,}", gene_ids))
      n_refseq      <- sum(grepl("^(NM_|NR_|XM_|XR_)", gene_ids))
      n_ensg        <- sum(grepl("^ENSG", gene_ids))
      n_numeric     <- sum(grepl("^[0-9]+$", gene_ids))
      n_other       <- n_total_genes - n_hgnc_like - n_genbank - n_refseq - n_ensg - n_numeric

      cat("[A2] Gene ID classification (", n_total_genes, " total):\n", sep = "")
      cat("     HGNC-like symbols : ", n_hgnc_like, " (", round(100 * n_hgnc_like / n_total_genes, 1), "%)\n", sep = "")
      cat("     GenBank accessions: ", n_genbank, " (", round(100 * n_genbank / n_total_genes, 1), "%)\n", sep = "")
      cat("     RefSeq IDs        : ", n_refseq, " (", round(100 * n_refseq / n_total_genes, 1), "%)\n", sep = "")
      cat("     ENSEMBL IDs       : ", n_ensg, " (", round(100 * n_ensg / n_total_genes, 1), "%)\n", sep = "")
      cat("     Numeric only      : ", n_numeric, " (", round(100 * n_numeric / n_total_genes, 1), "%)\n", sep = "")
      cat("     Other/ambiguous   : ", n_other, " (", round(100 * n_other / n_total_genes, 1), "%)\n", sep = "")
      cat("     Sample gene IDs (first 10):", paste(head(gene_ids, 10), collapse = ", "), "\n")
      cat("     Sample gene IDs (random 10):", paste(sample(gene_ids, min(10, length(gene_ids))), collapse = ", "), "\n")

      ## A3. Zero-count genes
      row_sums <- rowSums(counts_raw, na.rm = TRUE)
      n_zero   <- sum(row_sums == 0)
      n_tiny   <- sum(row_sums > 0 & row_sums < 10)
      cat("[A3] Count distribution by gene:\n")
      cat("     Zero across all samples: ", n_zero, " (", round(100 * n_zero / n_total_genes, 1), "%)\n", sep = "")
      cat("     Total counts 1-9       : ", n_tiny, " (", round(100 * n_tiny / n_total_genes, 1), "%)\n", sep = "")
      cat("     Total counts >= 10     : ", sum(row_sums >= 10), "\n", sep = "")
      cat("     Row-sum quantiles: ", paste(names(quantile(row_sums, c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1))),
          "=", format(quantile(row_sums, c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1)), big.mark = ","),
          collapse = "  "), "\n")

      ## A4. Flag non-standard IDs for potential removal
      is_standard <- grepl("^[A-Z][A-Z0-9-]*$", gene_ids) &
                     nchar(gene_ids) <= 15 &
                     !grepl("^(AK|BC|AJ|AF|AB|AL|CR|BX)[0-9]{5,}", gene_ids)
      n_standard <- sum(is_standard)
      cat("[A4] Standard gene symbols (potential HGNC): ", n_standard, " of ", n_total_genes, "\n", sep = "")
      if (n_standard > 0 && n_standard < n_total_genes) {
        ## Show what non-standard IDs contribute in terms of counts
        nonstandard_total <- sum(row_sums[!is_standard])
        standard_total    <- sum(row_sums[is_standard])
        cat("     Standard symbols carry ",
            round(100 * standard_total / (standard_total + nonstandard_total), 1),
            "% of all counts\n", sep = "")
        cat("     Non-standard IDs carry ",
            round(100 * nonstandard_total / (standard_total + nonstandard_total), 1),
            "% of all counts\n", sep = "")
      }

      ## A5. Key gene presence check (before any filtering)
      key_genes <- c("UCP1", "CLDN1", "ACTB", "GAPDH", "PPARG", "LEP")
      for (kg in key_genes) {
        idx <- which(gene_ids == kg)
        if (length(idx) > 0) {
          cat("[A5] ", kg, " present — row sum = ",
              format(row_sums[idx[1]], big.mark = ","),
              ", mean = ", round(mean(as.numeric(counts_raw[idx[1], ])), 1), "\n", sep = "")
        } else {
          cat("[A5] ", kg, " NOT FOUND in gene IDs\n", sep = "")
        }
      }
      cat("--- END SANITY BLOCK A ---\n\n")

    } else {
      ## Fallback: use exprs() from the GEO Series Matrix
      cat("[INFO] No count file found; using Series Matrix expression data\n")
      counts_raw <- as.data.frame(exprs(gse))
    }

    ## Build sample metadata from GEO phenoData
    pdata <- pData(gse)

    ## Parse tissue type and subject ID from characteristics columns
    ## Virtanen GEO samples have characteristics like:
    ##   "tissue: BAT" / "tissue: WAT" and "subject: S1", etc.
    char_cols <- grep("characteristics_ch", colnames(pdata), value = TRUE)

    ## Extract tissue
    tissue_col <- NULL
    for (cc in char_cols) {
      vals <- as.character(pdata[[cc]])
      if (any(grepl("tissue|BAT|WAT", vals, ignore.case = TRUE))) {
        tissue_col <- cc
        break
      }
    }

    if (!is.null(tissue_col)) {
      pdata$tissue <- gsub("^.*:\\s*", "", as.character(pdata[[tissue_col]]))
    } else {
      ## Try to infer from title
      pdata$tissue <- ifelse(grepl("BAT|brown", pdata$title, ignore.case = TRUE),
                              "BAT", "WAT")
    }

    ## Extract subject
    subject_col <- NULL
    for (cc in char_cols) {
      vals <- as.character(pdata[[cc]])
      if (any(grepl("subject|patient|individual|donor", vals, ignore.case = TRUE))) {
        subject_col <- cc
        break
      }
    }

    if (!is.null(subject_col)) {
      pdata$subject <- gsub("^.*:\\s*", "", as.character(pdata[[subject_col]]))
    } else {
      ## Infer paired subject from sample ordering
      n_bat <- sum(pdata$tissue == "BAT")
      pdata$subject <- paste0("S", rep(seq_len(n_bat), 2))
    }

    ## Ensure column alignment between counts and metadata
    pdata$geo_accession <- as.character(pdata$geo_accession)
    pdata_id <- ifelse(is.na(pdata$geo_accession) | pdata$geo_accession == "",
                       rownames(pdata), pdata$geo_accession)
    pdata_norm <- normalize_sample_id(pdata_id)
    counts_norm <- normalize_sample_id(colnames(counts_raw))

    map_idx <- match(counts_norm, pdata_norm)
    matched <- !is.na(map_idx)

    cat("[SANITY] Sample-name matching:", sum(matched), "of", length(counts_norm), "count columns matched to metadata\n")

    if (sum(matched) < 2) {
      ## Last resort: strict order match only if dimensions agree
      cat("[WARN] Weak sample-name matching; attempting order-based fallback\n")
      if (ncol(counts_raw) != nrow(pdata)) {
        stop("Unable to align count columns to metadata: matched ", sum(matched),
             ", ncol(counts_raw)=", ncol(counts_raw),
             ", nrow(pdata)=", nrow(pdata),
             ". Consider setting LOCAL_COUNTS_PATH to the correct file.")
      }
      colnames(counts_raw) <- pdata_id
      map_idx <- seq_len(nrow(pdata))
      matched <- rep(TRUE, ncol(counts_raw))
    }

    counts_raw <- counts_raw[, matched, drop = FALSE]
    map_idx <- map_idx[matched]
    pdata <- pdata[map_idx, , drop = FALSE]

    sample_ids <- pdata_id[map_idx]
    colnames(counts_raw) <- sample_ids
    rownames(pdata) <- sample_ids

    ## Diagnostic: show sample alignment so user can verify tissue assignment
    cat("[DIAG] Sample alignment (first 6):\n")
    diag_n <- min(6, ncol(counts_raw))
    diag_df <- data.frame(
      geo_id  = head(sample_ids, diag_n),
      tissue  = head(as.character(pdata$tissue), diag_n),
      subject = head(as.character(pdata$subject), diag_n),
      stringsAsFactors = FALSE
    )
    print(diag_df, row.names = FALSE)

    ## Keep only BAT/WAT samples used in the contrast
    keep_samples <- !is.na(pdata$tissue) & pdata$tissue %in% c(CONDITION_REF, CONDITION_TEST)
    counts_raw <- counts_raw[, keep_samples, drop = FALSE]
    pdata <- pdata[keep_samples, , drop = FALSE]

    ## Build SummarizedExperiment
    col_data <- DataFrame(
      tissue  = factor(pdata$tissue, levels = c(CONDITION_REF, CONDITION_TEST)),
      subject = factor(pdata$subject),
      row.names = rownames(pdata)
    )

    se <- SummarizedExperiment(
      assays  = list(counts = as.matrix(counts_raw)),
      colData = col_data
    )

    saveRDS(se, geo_cache)
    cat("[CACHE] Saved SummarizedExperiment to cache\n")
  }

  cat("[OK] SummarizedExperiment:",
      nrow(se), "genes x", ncol(se), "samples\n")
  cat("     Tissues:", paste(levels(se$tissue), collapse = " vs "), "\n")
  cat("     Subjects:", nlevels(se$subject), "\n")
  cat("[SANITY] First 6 sample IDs:", paste(head(colnames(se), 6), collapse = ", "), "\n")
  cat("[SANITY] Tissue counts:\n")
  print(table(se$tissue, useNA = "ifany"))
  cat("[SANITY] Subjects with paired samples (top 10):\n")
  subj_tab <- sort(table(se$subject), decreasing = TRUE)
  print(head(subj_tab, 10))

  ## Exclude subjects if specified (to match paper's 14-subject cohort)
  if (exists("SUBJECTS_TO_EXCLUDE") && !is.null(SUBJECTS_TO_EXCLUDE) &&
      length(SUBJECTS_TO_EXCLUDE) > 0) {
    keep_subj <- !(se$subject %in% SUBJECTS_TO_EXCLUDE)
    cat("[FILTER] Excluding subjects:", paste(SUBJECTS_TO_EXCLUDE, collapse = ", "),
        "->", sum(keep_subj), "of", ncol(se), "samples remain\n")
    se <- se[, keep_subj]
    se$subject <- droplevels(se$subject)
  }

  sanity$n_genes_raw   <- nrow(se)
  sanity$n_samples     <- ncol(se)
  sanity$n_subjects    <- nlevels(se$subject)
  sanity$tissues       <- levels(se$tissue)

  ## ---- 4.2) Pre-filtering --------------------------------
  ## Standard approach (Nature/Cell Metabolism convention):
  ##   Keep genes with CPM >= 1 in at least the smallest group size.
  ## This is purely count-driven — no gene-ID pattern matching needed.
  ## Non-coding junk (GenBank accessions etc.) naturally drops out
  ## because those features have zero/near-zero counts.
  ##
  ## References:
  ##   - edgeR::filterByExpr() (Chen et al. 2016, F1000Research)
  ##   - DESeq2 vignette: rowSums(counts >= 10) >= smallest_group
  ##   - Robinson et al. 2010, Genome Biology (edgeR)
  cat("\n== 4.2  Pre-filtering low-count genes ==\n")

  counts_mat <- assay(se, "counts")

  ## Ensure finite non-negative integer-like counts for DESeq2
  counts_mat <- as.matrix(counts_mat)
  storage.mode(counts_mat) <- "numeric"

  n_na_before <- sum(is.na(counts_mat) | !is.finite(counts_mat))
  if (n_na_before > 0) {
    cat("[SANITY] Replacing", n_na_before, "NA/Inf count values with 0 before filtering\n")
    counts_mat[is.na(counts_mat) | !is.finite(counts_mat)] <- 0
  }

  counts_mat[counts_mat < 0] <- 0
  counts_mat <- round(counts_mat)

  int_max <- .Machine$integer.max
  n_overflow <- sum(counts_mat > int_max)
  if (n_overflow > 0) {
    overflow_idx <- which(counts_mat > int_max, arr.ind = TRUE)
    cat("[SANITY] Capping", n_overflow, "count values above integer max at", int_max, "\n")
    cat("[SANITY] Overflow in gene(s):",
        paste(head(rownames(counts_mat)[overflow_idx[, 1]], 3), collapse = ", "), "\n")
    cat("[SANITY] Overflow in sample(s):",
        paste(head(colnames(counts_mat)[overflow_idx[, 2]], 3), collapse = ", "), "\n")
    counts_mat[counts_mat > int_max] <- int_max
  }

  storage.mode(counts_mat) <- "integer"
  assay(se, "counts") <- counts_mat

  if (LOW_COUNT_FILTER) {
    design_group <- se$tissue
    min_group    <- min(table(design_group))
    lib_sizes    <- colSums(counts_mat)
    med_lib      <- median(lib_sizes)

    cat("[FILTER] Smallest group size:", min_group, "samples\n")
    cat("[FILTER] Median library size:", format(med_lib, big.mark = ","), "\n")

    ## Step 1: fast pre-screen — remove genes with zero counts in ALL samples.
    ##         This is instant and eliminates the bulk of ~217K features.
    nonzero <- rowSums(counts_mat > 0) > 0
    n_allzero <- sum(!nonzero)
    cat("[FILTER] Removing", n_allzero, "genes with zero counts across all samples\n")
    counts_mat_nz <- counts_mat[nonzero, , drop = FALSE]
    cat("[FILTER] After zero removal:", nrow(counts_mat_nz), "genes remain\n")

    ## Step 2: CPM >= 1 in at least min_group samples (standard threshold).
    ##         This is the workhorse filter used in most Nature/Cell papers.
    cpm_mat  <- t(t(counts_mat_nz) / lib_sizes * 1e6)
    keep_cpm <- rowSums(cpm_mat >= 1) >= min_group
    n_pass_cpm <- sum(keep_cpm)
    cat("[FILTER] Genes with CPM >= 1 in >=", min_group, "samples:", n_pass_cpm, "\n")

    ## Step 3 (optional): if edgeR available, refine with filterByExpr
    ##         which also considers lib sizes and a min total count.
    if (requireNamespace("edgeR", quietly = TRUE)) {
      min_count_cpm1 <- max(MIN_COUNTS_PER_GENE, ceiling(med_lib / 1e6))
      keep_expr <- edgeR::filterByExpr(counts_mat_nz, group = design_group,
                                       min.count = min_count_cpm1)
      ## Use the intersection of both filters (most conservative)
      keep <- keep_cpm & keep_expr
      cat("[FILTER] edgeR::filterByExpr (min.count=", min_count_cpm1,
          ") keeps:", sum(keep_expr), "\n", sep = "")
      cat("[FILTER] Intersection of CPM + filterByExpr:", sum(keep), "\n")
    } else {
      keep <- keep_cpm
      cat("[FILTER] edgeR not available; using CPM filter only\n")
    }

    ## Apply filter back to full SE (must index into original row space)
    keep_names <- rownames(counts_mat_nz)[keep]
    se_filt <- se[keep_names, ]
    cat("[FILTER] Final: kept", nrow(se_filt), "of", nrow(se), "genes",
        "(removed", nrow(se) - nrow(se_filt), ")\n")
  } else {
    se_filt <- se
  }

  sanity$n_genes_after_id_filter <- nrow(se)  # no separate ID filter; same as raw

  ## ---- Sanity block B: Post-filter diagnostics ----
  cat("\n--- SANITY BLOCK B: Post-filter diagnostics ---\n")

  sanity$n_genes_filtered <- nrow(se_filt)
  filt_counts <- assay(se_filt, "counts")

  ## B1. Dimension check
  cat("[B1] Filtered matrix: ", nrow(se_filt), " genes x ", ncol(se_filt), " samples\n", sep = "")

  ## B2. Hard cap — if still >30K genes, something is wrong; do a second-pass
  MAX_GENES_FOR_DESEQ <- 30000
  if (nrow(se_filt) > MAX_GENES_FOR_DESEQ) {
    cat("[B2] WARNING: ", nrow(se_filt), " genes exceeds hard cap of ",
        MAX_GENES_FOR_DESEQ, "\n", sep = "")
    cat("[B2] Applying second-pass filter: keeping top genes by mean expression\n")
    row_means <- rowMeans(filt_counts)
    keep_top  <- order(row_means, decreasing = TRUE)[seq_len(MAX_GENES_FOR_DESEQ)]
    cutoff_mean <- row_means[keep_top[MAX_GENES_FOR_DESEQ]]
    se_filt <- se_filt[keep_top, ]
    filt_counts <- assay(se_filt, "counts")
    cat("[B2] Second-pass filter: kept ", nrow(se_filt), " genes (mean count cutoff >= ",
        round(cutoff_mean, 1), ")\n", sep = "")
    sanity$n_genes_filtered <- nrow(se_filt)
    sanity$second_pass_filter <- TRUE
  } else {
    cat("[B2] Gene count (", nrow(se_filt), ") is within safe range for DESeq2\n", sep = "")
    sanity$second_pass_filter <- FALSE
  }

  ## B3. Library sizes after filtering
  filt_lib <- colSums(filt_counts, na.rm = TRUE)
  cat("[B3] Post-filter library sizes:\n")
  cat("     min =", format(min(filt_lib), big.mark = ","),
      " median =", format(median(filt_lib), big.mark = ","),
      " max =", format(max(filt_lib), big.mark = ","), "\n")

  ## B4. Key gene survival check
  cat("[B4] Key gene survival after filtering:\n")
  filt_gene_ids <- rownames(se_filt)
  for (kg in c("UCP1", "CLDN1", "ACTB", "GAPDH", "PPARG", "LEP",
               BAT_MARKER_GENES[1:3], WAT_MARKER_GENES[1:3])) {
    status <- if (kg %in% filt_gene_ids) "PRESENT" else "REMOVED"
    cat("     ", kg, ": ", status, "\n", sep = "")
  }

  ## B5. Count distribution summary
  filt_row_sums <- rowSums(filt_counts, na.rm = TRUE)
  cat("[B5] Filtered gene count distribution:\n")
  cat("     Quantiles: ", paste(
    names(quantile(filt_row_sums, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
    "=", format(quantile(filt_row_sums, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)), big.mark = ","),
    collapse = "  "), "\n")

  cat("--- END SANITY BLOCK B ---\n\n")

  ## ---- 4.3) DESeq2 — paired design ----------------------
  cat("\n== 4.3  DESeq2 paired analysis ==\n")

  ## Sanity block C: Pre-DESeq2 diagnostics
  cat("--- SANITY BLOCK C: Pre-DESeq2 diagnostics ---\n")

  n_genes_for_deseq <- nrow(se_filt)
  n_samples_for_deseq <- ncol(se_filt)
  n_subjects_for_deseq <- nlevels(se_filt$subject)

  ## C1. Design matrix dimensions
  ## Design ~ subject + tissue: ncol = n_subjects + 1 (tissue)
  ## (intercept absorbed into subject coding)
  design_cols <- n_subjects_for_deseq  # subject levels - 1 + intercept + tissue = n_subjects + 1
  cat("[C1] Design: ~ subject + tissue\n")
  cat("     Subjects: ", n_subjects_for_deseq, " | Samples: ", n_samples_for_deseq, "\n", sep = "")
  cat("     Approximate design matrix: ", n_samples_for_deseq, " x ", design_cols + 1, "\n", sep = "")

  ## C2. Memory estimate (rough: DESeq2 needs ~8 bytes * genes * samples * 4 matrices)
  mem_est_mb <- (8 * n_genes_for_deseq * n_samples_for_deseq * 6) / (1024^2)
  cat("[C2] Estimated DESeq2 memory: ~", round(mem_est_mb, 0), " MB for core matrices\n", sep = "")
  avail_mem <- tryCatch({
    meminfo <- readLines("/proc/meminfo", n = 3)
    avail_line <- grep("MemAvailable", meminfo, value = TRUE)
    if (length(avail_line) > 0) {
      avail_kb <- as.numeric(gsub("[^0-9]", "", avail_line))
      round(avail_kb / 1024)
    } else NA
  }, error = function(e) NA)
  if (!is.na(avail_mem)) {
    cat("[C2] Available system memory: ~", avail_mem, " MB\n", sep = "")
    if (mem_est_mb > avail_mem * 0.5) {
      cat("[C2] WARNING: DESeq2 may consume >50% of available memory!\n")
    }
  }

  ## C3. Time estimate (empirical: ~0.5-2 sec per 1000 genes for paired design)
  time_est_min <- n_genes_for_deseq / 1000 * 1.5 / 60
  cat("[C3] Estimated DESeq2 runtime: ~", round(time_est_min, 1), " minutes\n", sep = "")
  if (time_est_min > 30) {
    cat("[C3] WARNING: This will take a LONG time. Consider more aggressive filtering.\n")
  }

  cat("--- END SANITY BLOCK C ---\n\n")

  dds <- DESeqDataSet(se_filt, design = ~ subject + tissue)
  dds$tissue <- relevel(dds$tissue, ref = CONDITION_REF)

  cat("[INFO] Starting DESeq2 at:", format(Sys.time(), "%H:%M:%S"), "\n")
  deseq_start <- Sys.time()

  ## Run DESeq2 (or load from cache)
  ## C4. Cache validation: invalidate if gene count changed (e.g. new filters)
  if (use_save_core) {
    cached_dds_path <- file.path(cache_dir, "dds_BAT_vs_WAT.rds")
    if (file.exists(cached_dds_path)) {
      cached_dds <- readRDS(cached_dds_path)
      cached_ngenes <- nrow(cached_dds)
      if (cached_ngenes != nrow(dds)) {
        cat("[C4] CACHE STALE: cached DESeq2 has", cached_ngenes,
            "genes but current filter yields", nrow(dds), "\n")
        cat("[C4] Removing stale cache to force recompute\n")
        file.remove(cached_dds_path)
        ## Also remove stale VST cache
        vst_cache_path <- file.path(cache_dir, "vst_counts.rds")
        if (file.exists(vst_cache_path)) file.remove(vst_cache_path)
      } else {
        cat("[C4] Cache gene count matches (", cached_ngenes, ") — using cache\n", sep = "")
      }
      rm(cached_dds)
    }
    dds <- cache_deseq(dds, cache_dir, contrast_name = "BAT_vs_WAT")
  } else {
    dds <- DESeq(dds)
  }

  deseq_elapsed <- difftime(Sys.time(), deseq_start, units = "mins")
  cat("[INFO] DESeq2 finished at:", format(Sys.time(), "%H:%M:%S"),
      "(", round(as.numeric(deseq_elapsed), 1), "min )\n")
  sanity$deseq_runtime_min <- round(as.numeric(deseq_elapsed), 2)

  ## Extract results
  res <- results(dds,
                  contrast = c("tissue", CONDITION_TEST, CONDITION_REF),
                  alpha = PADJ_CUTOFF)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

  cat("[OK] DESeq2 complete:", nrow(res_df), "genes tested\n")
  cat("     Up in BAT (padj <", PADJ_CUTOFF, ", LFC >=", LFC_CUTOFF, "):",
      sum(res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange >= LFC_CUTOFF,
          na.rm = TRUE), "\n")
  cat("     Up in WAT (padj <", PADJ_CUTOFF, ", LFC <=", -LFC_CUTOFF, "):",
      sum(res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange <= -LFC_CUTOFF,
          na.rm = TRUE), "\n")

  sanity$n_up_bat <- sum(res_df$padj < PADJ_CUTOFF &
                          res_df$log2FoldChange >= LFC_CUTOFF, na.rm = TRUE)
  sanity$n_up_wat <- sum(res_df$padj < PADJ_CUTOFF &
                          res_df$log2FoldChange <= -LFC_CUTOFF, na.rm = TRUE)

  ## Compare to paper
  cat("[PAPER] Virtanen et al. report", EXPECTED_UP_BAT,
      "BAT-enriched genes (14 subjects)\n")
  cat("[PAPER] This run:", sanity$n_up_bat, "BAT-enriched genes (",
      nlevels(se_filt$subject), "subjects)\n")
  if (abs(sanity$n_up_bat - EXPECTED_UP_BAT) / EXPECTED_UP_BAT <= DEG_COUNT_TOLERANCE) {
    cat("[PAPER] Within", DEG_COUNT_TOLERANCE * 100, "% tolerance\n")
  } else {
    cat("[PAPER] Outside tolerance — may reflect different subject count or filtering\n")
  }

  ## ---- 4.4) Gene symbol mapping --------------------------
  cat("\n== 4.4  Mapping gene IDs to symbols ==\n")

  ## Check if gene IDs are already symbols or ENSEMBL
  sample_ids <- head(res_df$gene_id, 20)
  ids_are_ensembl <- any(grepl("^ENSG", sample_ids))

  if (ids_are_ensembl) {
    ## Strip version suffix (ENSG00000123456.7 -> ENSG00000123456)
    res_df$gene_id_clean <- gsub("\\.\\d+$", "", res_df$gene_id)

    gene_symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys    = res_df$gene_id_clean,
      keytype = "ENSEMBL",
      column  = "SYMBOL",
      multiVals = "first"
    )
    res_df$symbol <- gene_symbols[res_df$gene_id_clean]
    n_mapped <- sum(!is.na(res_df$symbol))
  } else {
    ## Gene IDs are already symbols (e.g. GEO humanBATWAT.txt)
    res_df$symbol <- res_df$gene_id
    n_mapped <- sum(res_df$symbol != "")
  }

  ## Fill NA symbols with gene_id
  res_df$symbol[is.na(res_df$symbol)] <- res_df$gene_id[is.na(res_df$symbol)]

  cat("[OK] Mapped", n_mapped, "of", nrow(res_df), "genes to HGNC symbols\n")

  ## ---- 4.5) CLDN1 — the primary question -----------------
  cat("\n")
  cat("##########################################################\n")
  cat("##  PRIMARY QUERY: Is CLDN1 increased in human BAT?     ##\n")
  cat("##########################################################\n\n")

  cldn1_row <- res_df[res_df$symbol == "CLDN1", ]

  if (nrow(cldn1_row) == 0) {
    cat("[!!] CLDN1 NOT FOUND in the expression dataset\n")
    cat("     Possible reasons: filtered out (low counts),\n")
    cat("     not annotated, or not detected in this tissue panel.\n")
    sanity$cldn1_found <- FALSE
    sanity$cldn1_verdict <- "NOT_DETECTED"
  } else {
    cldn1 <- cldn1_row[1, ]
    cat("  Gene ID      :", cldn1$gene_id, "\n")
    cat("  Symbol       : CLDN1\n")
    cat("  log2FC (BAT/WAT):", round(cldn1$log2FoldChange, 3), "\n")
    cat("  p-value      :", signif(cldn1$pvalue, 4), "\n")
    cat("  padj (BH)    :", signif(cldn1$padj, 4), "\n")
    cat("  baseMean     :", round(cldn1$baseMean, 1), "\n\n")

    sanity$cldn1_found    <- TRUE
    sanity$cldn1_lfc      <- cldn1$log2FoldChange
    sanity$cldn1_padj     <- cldn1$padj
    sanity$cldn1_baseMean <- cldn1$baseMean

    if (!is.na(cldn1$padj) && cldn1$padj < PADJ_CUTOFF &&
        cldn1$log2FoldChange > 0) {
      cat("  >>> VERDICT: YES — CLDN1 is SIGNIFICANTLY INCREASED in BAT\n")
      if (cldn1$log2FoldChange >= LFC_CUTOFF) {
        cat("      (meets both padj and log2FC thresholds)\n")
      } else {
        cat("      (significant but below LFC cutoff of", LFC_CUTOFF, ")\n")
      }
      sanity$cldn1_verdict <- "UP_IN_BAT_SIGNIFICANT"
    } else if (!is.na(cldn1$log2FoldChange) && cldn1$log2FoldChange > 0) {
      cat("  >>> VERDICT: TREND — CLDN1 shows a positive fold-change in BAT\n")
      cat("      but does NOT reach statistical significance (padj =",
          signif(cldn1$padj, 3), ")\n")
      sanity$cldn1_verdict <- "UP_TREND_NOT_SIGNIFICANT"
    } else if (!is.na(cldn1$log2FoldChange) && cldn1$log2FoldChange < 0) {
      cat("  >>> VERDICT: NO — CLDN1 is LOWER in BAT than WAT\n")
      cat("      (log2FC =", round(cldn1$log2FoldChange, 3), ")\n")
      sanity$cldn1_verdict <- "DOWN_IN_BAT"
    } else {
      cat("  >>> VERDICT: INCONCLUSIVE (NA values in results)\n")
      sanity$cldn1_verdict <- "INCONCLUSIVE"
    }
  }

  cat("\n##########################################################\n\n")

  ## ---- 4.5b) UCP1 — positive control ----------------------
  cat("\n")
  cat("##########################################################\n")
  cat("##  CONTROL: UCP1 expression (canonical BAT marker)     ##\n")
  cat("##########################################################\n\n")

  ucp1_row <- res_df[res_df$symbol == "UCP1", ]

  if (nrow(ucp1_row) == 0) {
    cat("[!!] UCP1 NOT FOUND — this is unexpected for BAT/WAT data\n")
    sanity$ucp1_found   <- FALSE
    sanity$ucp1_verdict <- "NOT_DETECTED"
  } else {
    ucp1 <- ucp1_row[1, ]
    cat("  Gene ID      :", ucp1$gene_id, "\n")
    cat("  Symbol       : UCP1\n")
    cat("  log2FC (BAT/WAT):", round(ucp1$log2FoldChange, 3), "\n")
    cat("  p-value      :", signif(ucp1$pvalue, 4), "\n")
    cat("  padj (BH)    :", signif(ucp1$padj, 4), "\n")
    cat("  baseMean     :", round(ucp1$baseMean, 1), "\n\n")

    sanity$ucp1_found    <- TRUE
    sanity$ucp1_lfc      <- ucp1$log2FoldChange
    sanity$ucp1_padj     <- ucp1$padj
    sanity$ucp1_baseMean <- ucp1$baseMean

    if (!is.na(ucp1$padj) && ucp1$padj < PADJ_CUTOFF &&
        ucp1$log2FoldChange > 0) {
      cat("  >>> CONTROL PASS: UCP1 is SIGNIFICANTLY UP in BAT\n")
      if (ucp1$log2FoldChange >= LFC_CUTOFF) {
        cat("      (meets both padj and log2FC thresholds)\n")
      } else {
        cat("      (significant but below LFC cutoff of", LFC_CUTOFF, ")\n")
      }
      sanity$ucp1_verdict <- "UP_IN_BAT_SIGNIFICANT"
    } else if (!is.na(ucp1$log2FoldChange) && ucp1$log2FoldChange > 0) {
      cat("  >>> CONTROL WARNING: UCP1 trends up but is not significant\n")
      sanity$ucp1_verdict <- "UP_TREND_NOT_SIGNIFICANT"
    } else {
      cat("  >>> CONTROL WARNING: UCP1 not up in BAT — check data integrity\n")
      sanity$ucp1_verdict <- "UNEXPECTED"
    }
  }

  cat("\n##########################################################\n\n")

  ## ---- 4.6) Validate BAT markers -------------------------
  cat("== 4.6  BAT marker gene validation ==\n")

  marker_check <- function(gene_list, label) {
    hits <- res_df[res_df$symbol %in% gene_list, ]
    cat("\n--- ", label, " (", nrow(hits), "/", length(gene_list),
        " found) ---\n", sep = "")
    if (nrow(hits) > 0) {
      for (i in seq_len(nrow(hits))) {
        sig_flag <- ""
        if (!is.na(hits$padj[i]) && hits$padj[i] < PADJ_CUTOFF)
          sig_flag <- " ***"
        cat(sprintf("  %-12s  LFC = %+7.3f   padj = %-10s  baseMean = %8.1f%s\n",
                    hits$symbol[i],
                    hits$log2FoldChange[i],
                    ifelse(is.na(hits$padj[i]), "NA",
                           formatC(hits$padj[i], format = "e", digits = 2)),
                    hits$baseMean[i],
                    sig_flag))
      }
    }
    return(hits)
  }

  bat_hits <- marker_check(BAT_MARKER_GENES, "BAT markers")
  wat_hits <- marker_check(WAT_MARKER_GENES, "WAT markers")

  ## Also show CLDN1 among claudin family
  claudins <- grep("^CLDN", res_df$symbol, value = TRUE)
  if (length(claudins) > 0) {
    cat("\n--- Claudin family members detected ---\n")
    cldn_df <- res_df[res_df$symbol %in% claudins, ] %>% arrange(padj)
    for (i in seq_len(nrow(cldn_df))) {
      sig_flag <- ""
      if (!is.na(cldn_df$padj[i]) && cldn_df$padj[i] < PADJ_CUTOFF)
        sig_flag <- " ***"
      cat(sprintf("  %-12s  LFC = %+7.3f   padj = %-10s  baseMean = %8.1f%s\n",
                  cldn_df$symbol[i],
                  cldn_df$log2FoldChange[i],
                  ifelse(is.na(cldn_df$padj[i]), "NA",
                         formatC(cldn_df$padj[i], format = "e", digits = 2)),
                  cldn_df$baseMean[i],
                  sig_flag))
    }
  }

  sanity$n_bat_markers_found <- nrow(bat_hits)
  sanity$n_bat_markers_sig   <- sum(bat_hits$padj < PADJ_CUTOFF &
                                     bat_hits$log2FoldChange > 0,
                                     na.rm = TRUE)

  ## ---- 4.7) VST transformation for QC --------------------
  cat("\n== 4.7  VST transformation ==\n")

  if (use_save_core) {
    vsd <- cache_vst(dds, cache_dir)
  } else {
    vsd <- vst(dds, blind = FALSE)
  }

  vst_mat <- assay(vsd)
  cat("[OK] VST matrix:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n")

  ## ---- 4.8) Save full DE table ---------------------------
  cat("\n== 4.8  Saving results tables ==\n")

  de_file <- paste0("DE_BAT_vs_WAT_", run_tag, ".csv")
  write.csv(res_df, file = file.path(outdir, "tables", de_file),
            row.names = FALSE)
  cat("[SAVED]", de_file, "\n")
  verify_output_file(file.path(outdir, "tables", de_file), "DE full CSV")

  ## DEG-only table
  deg_df <- res_df %>%
    filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) >= LFC_CUTOFF)

  deg_file <- paste0("DEGs_BAT_vs_WAT_padj", PADJ_CUTOFF,
                      "_lfc", LFC_CUTOFF, "_", run_tag, ".csv")
  write.csv(deg_df, file = file.path(outdir, "tables", deg_file),
            row.names = FALSE)
  cat("[SAVED]", deg_file, "(", nrow(deg_df), "DEGs )\n")
  verify_output_file(file.path(outdir, "tables", deg_file), "DEG CSV")

  ## CLDN1-specific table
  cldn1_file <- paste0("CLDN1_result_", run_tag, ".csv")
  write.csv(cldn1_row, file = file.path(outdir, "tables", cldn1_file),
            row.names = FALSE)
  cat("[SAVED]", cldn1_file, "\n")
  verify_output_file(file.path(outdir, "tables", cldn1_file), "CLDN1 CSV")

  ## UCP1-specific table
  ucp1_file <- paste0("UCP1_result_", run_tag, ".csv")
  write.csv(ucp1_row, file = file.path(outdir, "tables", ucp1_file),
            row.names = FALSE)
  cat("[SAVED]", ucp1_file, "\n")
  verify_output_file(file.path(outdir, "tables", ucp1_file), "UCP1 CSV")

  ## Authoritative DEG lists for downstream scripts
  deg_lists <- list(
    up_in_BAT = deg_df %>% filter(log2FoldChange > 0) %>% pull(symbol),
    up_in_WAT = deg_df %>% filter(log2FoldChange < 0) %>% pull(symbol),
    all_degs  = deg_df$symbol
  )
  deg_lists_file <- paste0("DEG_lists_authoritative_", run_tag, ".rds")
  saveRDS(deg_lists, file = file.path(outdir, "tables", deg_lists_file))
  cat("[SAVED]", deg_lists_file, "\n")
  verify_output_file(file.path(outdir, "tables", deg_lists_file), "DEG list RDS")

  ## Save VST matrix for downstream parts
  vst_file <- paste0("VST_matrix_", run_tag, ".rds")
  saveRDS(list(vst_mat = vst_mat, col_data = colData(dds)),
          file = file.path(outdir, "tables", vst_file))
  cat("[SAVED]", vst_file, "\n")
  verify_output_file(file.path(outdir, "tables", vst_file), "VST RDS")

  ## ---- 4.9) Plots ----------------------------------------
  if (GENERATE_PLOTS) {
    cat("\n== 4.9  Generating plots ==\n")

    ## -- 4.9a) Volcano plot highlighting CLDN1 + BAT markers --
    cat("[PLOT] Volcano plot\n")

    highlight_genes <- c("CLDN1", BAT_MARKER_GENES[1:5], WAT_MARKER_GENES[1:3])
    res_df$highlight <- ifelse(res_df$symbol %in% highlight_genes,
                                res_df$symbol, NA)

    volcano <- EnhancedVolcano(
      res_df,
      lab             = res_df$symbol,
      selectLab       = highlight_genes,
      x               = "log2FoldChange",
      y               = "padj",
      pCutoff         = PADJ_CUTOFF,
      FCcutoff        = LFC_CUTOFF,
      title           = "BAT vs WAT (Virtanen et al. 2018)",
      subtitle        = paste("CLDN1 highlighted | n =", N_SUBJECTS, "paired subjects"),
      caption          = paste("Total genes:", nrow(res_df)),
      drawConnectors  = TRUE,
      widthConnectors = 0.5,
      colConnectors   = "grey30",
      max.overlaps    = 30,
      pointSize       = 1.5,
      labSize         = 4.0,
      col             = c("grey80", "forestgreen", "royalblue", "red2")
    )

    volcano_file <- paste0("Volcano_BAT_vs_WAT_CLDN1_", run_tag, ".png")
    ggsave(file.path(outdir, "plots", volcano_file), plot = volcano,
           width = VOLCANO_W, height = VOLCANO_H, dpi = PLOT_DPI)
    cat("[SAVED]", volcano_file, "\n")
    verify_output_file(file.path(outdir, "plots", volcano_file), "Volcano plot")

    ## -- 4.9b) PCA (with subject labels for outlier identification) --
    cat("[PLOT] PCA biplot (subject-labelled)\n")

    pca_data <- plotPCA(vsd, intgroup = c("tissue", "subject"), returnData = TRUE)
    pct_var  <- round(100 * attr(pca_data, "percentVar"))

    p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = tissue,
                                   label = subject)) +
      geom_point(size = 4, alpha = 0.8) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20,
                                show.legend = FALSE) +
      scale_colour_manual(values = c("WAT" = "#2166AC", "BAT" = "#B2182B")) +
      labs(
        title = "PCA: BAT vs WAT (Virtanen 2018)",
        x = paste0("PC1 (", pct_var[1], "% variance)"),
        y = paste0("PC2 (", pct_var[2], "% variance)")
      ) +
      theme_bw(base_size = 14) +
      theme(legend.position = "bottom")

    pca_file <- paste0("PCA_BAT_vs_WAT_", run_tag, ".png")
    ggsave(file.path(outdir, "plots", pca_file), plot = p_pca,
           width = PCA_W, height = PCA_H, dpi = PLOT_DPI)
    cat("[SAVED]", pca_file, "\n")
    verify_output_file(file.path(outdir, "plots", pca_file), "PCA plot")

    ## -- 4.9c) CLDN1 expression box plot (per-sample) --
    cat("[PLOT] CLDN1 per-sample boxplot\n")

    cldn1_id <- res_df$gene_id[res_df$symbol == "CLDN1"]
    if (length(cldn1_id) > 0 && cldn1_id[1] %in% rownames(vst_mat)) {
      cldn1_vst <- data.frame(
        expression = vst_mat[cldn1_id[1], ],
        tissue     = colData(dds)$tissue,
        subject    = colData(dds)$subject
      )

      p_cldn1 <- ggplot(cldn1_vst, aes(x = tissue, y = expression, fill = tissue)) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
        geom_line(aes(group = subject), colour = "grey50", linewidth = 0.4) +
        geom_point(aes(colour = tissue), size = 3) +
        scale_fill_manual(values = c("WAT" = "#2166AC", "BAT" = "#B2182B")) +
        scale_colour_manual(values = c("WAT" = "#2166AC", "BAT" = "#B2182B")) +
        labs(
          title    = "CLDN1 Expression: BAT vs WAT",
          subtitle = paste0("log2FC = ",
                            round(cldn1_row$log2FoldChange[1], 2),
                            ", padj = ",
                            signif(cldn1_row$padj[1], 3)),
          y = "VST-normalised expression",
          x = NULL
        ) +
        theme_bw(base_size = 14) +
        theme(legend.position = "none")

      cldn1_plot_file <- paste0("CLDN1_boxplot_BAT_vs_WAT_", run_tag, ".png")
      ggsave(file.path(outdir, "plots", cldn1_plot_file), plot = p_cldn1,
             width = 5, height = 6, dpi = PLOT_DPI)
      cat("[SAVED]", cldn1_plot_file, "\n")
      verify_output_file(file.path(outdir, "plots", cldn1_plot_file), "CLDN1 boxplot")
    } else {
      cat("[SKIP] CLDN1 not in VST matrix — cannot plot\n")
    }

    ## -- 4.9c2) UCP1 expression box plot (positive control) --
    cat("[PLOT] UCP1 per-sample boxplot\n")

    ucp1_id <- res_df$gene_id[res_df$symbol == "UCP1"]
    if (length(ucp1_id) > 0 && ucp1_id[1] %in% rownames(vst_mat)) {
      ucp1_vst <- data.frame(
        expression = vst_mat[ucp1_id[1], ],
        tissue     = colData(dds)$tissue,
        subject    = colData(dds)$subject
      )

      p_ucp1 <- ggplot(ucp1_vst, aes(x = tissue, y = expression, fill = tissue)) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
        geom_line(aes(group = subject), colour = "grey50", linewidth = 0.4) +
        geom_point(aes(colour = tissue), size = 3) +
        scale_fill_manual(values = c("WAT" = "#2166AC", "BAT" = "#B2182B")) +
        scale_colour_manual(values = c("WAT" = "#2166AC", "BAT" = "#B2182B")) +
        labs(
          title    = "UCP1 Expression: BAT vs WAT (Positive Control)",
          subtitle = paste0("log2FC = ", round(ucp1_row$log2FoldChange[1], 2),
                            ", padj = ", signif(ucp1_row$padj[1], 3)),
          y = "VST-normalised expression", x = NULL
        ) +
        theme_bw(base_size = 14) + theme(legend.position = "none")

      ucp1_plot_file <- paste0("UCP1_boxplot_BAT_vs_WAT_", run_tag, ".png")
      ggsave(file.path(outdir, "plots", ucp1_plot_file), plot = p_ucp1,
             width = 5, height = 6, dpi = PLOT_DPI)
      cat("[SAVED]", ucp1_plot_file, "\n")
      verify_output_file(file.path(outdir, "plots", ucp1_plot_file), "UCP1 boxplot")
    } else {
      cat("[SKIP] UCP1 not in VST matrix — cannot plot\n")
    }

    ## -- 4.9c3) Correlation: CLDN1 vs UCP1 in BAT samples --
    cat("[PLOT] CLDN1 vs UCP1 Spearman correlation (BAT only)\n")

    cor_cldn1_id <- res_df$gene_id[res_df$symbol == "CLDN1"][1]
    cor_ucp1_id  <- res_df$gene_id[res_df$symbol == "UCP1"][1]

    if (!is.na(cor_cldn1_id) && !is.na(cor_ucp1_id) &&
        cor_cldn1_id %in% rownames(vst_mat) &&
        cor_ucp1_id %in% rownames(vst_mat)) {

      bat_idx <- which(colData(dds)$tissue == "BAT")
      cor_df <- data.frame(
        CLDN1   = vst_mat[cor_cldn1_id, bat_idx],
        UCP1    = vst_mat[cor_ucp1_id,  bat_idx],
        subject = colData(dds)$subject[bat_idx]
      )

      cor_test <- cor.test(cor_df$CLDN1, cor_df$UCP1, method = "spearman",
                           exact = FALSE)
      cat("[CORR] Spearman rho =", round(cor_test$estimate, 3),
          ", p =", signif(cor_test$p.value, 3), "\n")

      sanity$cldn1_ucp1_rho  <- cor_test$estimate
      sanity$cldn1_ucp1_pval <- cor_test$p.value

      p_cor <- ggplot(cor_df, aes(x = UCP1, y = CLDN1)) +
        geom_point(colour = "#B2182B", size = 3) +
        geom_smooth(method = "lm", se = TRUE, colour = "grey30",
                    linewidth = 0.8) +
        ggrepel::geom_text_repel(aes(label = subject), size = 2.5,
                                  max.overlaps = 15) +
        labs(
          title = "CLDN1 vs UCP1 in BAT (Spearman)",
          subtitle = paste0("rho = ", round(cor_test$estimate, 3),
                            ", p = ", signif(cor_test$p.value, 3)),
          x = "UCP1 (VST-normalised)", y = "CLDN1 (VST-normalised)"
        ) +
        theme_bw(base_size = 14)

      cor_file <- paste0("CLDN1_vs_UCP1_BAT_correlation_", run_tag, ".png")
      ggsave(file.path(outdir, "plots", cor_file), plot = p_cor,
             width = 6, height = 6, dpi = PLOT_DPI)
      cat("[SAVED]", cor_file, "\n")
      verify_output_file(file.path(outdir, "plots", cor_file), "Correlation plot")
    } else {
      cat("[SKIP] CLDN1 or UCP1 not found — cannot compute correlation\n")
    }

    ## -- 4.9d) Heatmap of top DEGs with CLDN1 marked ------
    cat("[PLOT] Top DEG heatmap\n")

    top_degs <- res_df %>%
      filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) >= LFC_CUTOFF) %>%
      slice_min(padj, n = HEATMAP_TOP_N)

    ## Force-include CLDN1 if it's a DEG but not in top N
    if (nrow(cldn1_row) > 0 && !("CLDN1" %in% top_degs$symbol)) {
      if (!is.na(cldn1_row$padj[1]) && cldn1_row$padj[1] < PADJ_CUTOFF) {
        top_degs <- bind_rows(top_degs, cldn1_row[1, ])
      }
    }

    hm_ids <- top_degs$gene_id[top_degs$gene_id %in% rownames(vst_mat)]
    if (length(hm_ids) >= 5) {
      hm_mat <- vst_mat[hm_ids, ]

      ## Use symbols as row labels
      rownames(hm_mat) <- top_degs$symbol[match(rownames(hm_mat), top_degs$gene_id)]

      ## Annotation
      annot_col <- data.frame(
        Tissue = colData(dds)$tissue,
        row.names = colnames(hm_mat)
      )
      annot_colors <- list(
        Tissue = c("BAT" = "#B2182B", "WAT" = "#2166AC")
      )

      heat_file <- paste0("Heatmap_topDEGs_CLDN1_", run_tag, ".png")
      png(file.path(outdir, "plots", heat_file),
          width = HEATMAP_W, height = HEATMAP_H, res = HEATMAP_RES)
      pheatmap(
        hm_mat,
        scale            = HEATMAP_SCALE,
        clustering_method = "ward.D2",
        annotation_col   = annot_col,
        annotation_colors = annot_colors,
        show_colnames    = FALSE,
        fontsize_row     = 8,
        main             = paste0("Top ", HEATMAP_TOP_N,
                                   " DEGs (BAT vs WAT) + CLDN1"),
        color            = colorRampPalette(
                             rev(brewer.pal(11, "RdBu")))(100)
      )
      dev.off()
      cat("[SAVED]", heat_file, "\n")
      verify_output_file(file.path(outdir, "plots", heat_file), "Top DEG heatmap")
    } else {
      cat("[SKIP] Too few DEGs for heatmap (", length(hm_ids), ")\n")
    }

    ## -- 4.9e) BAT marker + CLDN1 bar chart ----------------
    cat("[PLOT] BAT markers LFC bar chart\n")

    focus_genes <- c("CLDN1", BAT_MARKER_GENES)
    focus_df <- res_df %>%
      filter(symbol %in% focus_genes) %>%
      mutate(
        sig = ifelse(!is.na(padj) & padj < PADJ_CUTOFF, "Significant", "NS"),
        symbol = factor(symbol, levels = rev(
          c("CLDN1", BAT_MARKER_GENES[BAT_MARKER_GENES %in% symbol])
        ))
      )

    if (nrow(focus_df) > 0) {
      p_bar <- ggplot(focus_df, aes(x = symbol, y = log2FoldChange, fill = sig)) +
        geom_col(width = 0.7) +
        geom_hline(yintercept = 0, linewidth = 0.5) +
        geom_hline(yintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
                   linetype = "dashed", colour = "grey50") +
        scale_fill_manual(values = c("Significant" = "#B2182B", "NS" = "grey60")) +
        coord_flip() +
        labs(
          title = "CLDN1 among BAT marker genes",
          subtitle = "log2 fold-change (BAT / WAT)",
          x = NULL, y = "log2 Fold-Change",
          fill = NULL
        ) +
        theme_bw(base_size = 13)

      bar_file <- paste0("BAT_markers_CLDN1_barplot_", run_tag, ".png")
      ggsave(file.path(outdir, "plots", bar_file), plot = p_bar,
             width = 7, height = 5, dpi = PLOT_DPI)
      cat("[SAVED]", bar_file, "\n")
      verify_output_file(file.path(outdir, "plots", bar_file), "Marker barplot")
    }

    ## Register hero plot
    if (use_save_core && exists("save_run_file")) {
      hero_path <- file.path(outdir, "plots", volcano_file)
      if (file.exists(hero_path)) {
        tryCatch(save_run_file(hero_path, tag = "hero_volcano"),
                 error = function(e)
                   cat("[WARN] save_run_file failed:", conditionMessage(e), "\n"))
      }
    }
  }

  ## ---- 4.10) DEG summary table ---------------------------
  cat("\n== 4.10  DEG summary ==\n")

  deg_summary <- data.frame(
    contrast      = "BAT_vs_WAT",
    n_tested      = nrow(res_df),
    n_up_bat      = sanity$n_up_bat,
    n_up_wat      = sanity$n_up_wat,
    n_total_deg   = sanity$n_up_bat + sanity$n_up_wat,
    paper_up_bat  = EXPECTED_UP_BAT,
    ucp1_lfc      = ifelse(isTRUE(sanity$ucp1_found), round(sanity$ucp1_lfc, 4), NA),
    ucp1_padj     = ifelse(isTRUE(sanity$ucp1_found), signif(sanity$ucp1_padj, 4), NA),
    ucp1_verdict  = ifelse(is.null(sanity$ucp1_verdict), NA, sanity$ucp1_verdict),
    cldn1_lfc     = ifelse(sanity$cldn1_found, round(sanity$cldn1_lfc, 4), NA),
    cldn1_padj    = ifelse(sanity$cldn1_found, signif(sanity$cldn1_padj, 4), NA),
    cldn1_verdict = sanity$cldn1_verdict,
    cldn1_ucp1_rho = ifelse(!is.null(sanity$cldn1_ucp1_rho),
                             round(sanity$cldn1_ucp1_rho, 4), NA),
    stringsAsFactors = FALSE
  )

  deg_summary_file <- paste0("DEG_summary_", run_tag, ".csv")
  write.csv(deg_summary, file = file.path(outdir, "tables", deg_summary_file),
            row.names = FALSE)
  cat("[SAVED]", deg_summary_file, "\n")
  verify_output_file(file.path(outdir, "tables", deg_summary_file), "DE summary CSV")
  print(deg_summary)

  ## ---- 4.11) Sanity check table --------------------------
  cat("\n== 4.11  Validation summary ==\n")

  sanity_df <- data.frame(
    check = c(
      "Samples loaded",
      "Subjects detected",
      "Genes raw (before any filter)",
      "Genes after count filter (CPM>=1)",
      "Second-pass filter applied",
      "DESeq2 runtime (min)",
      "BAT-enriched genes (paper: 463)",
      "UCP1 verdict",
      "UCP1 log2FC",
      "UCP1 padj",
      "CLDN1 found",
      "CLDN1 verdict",
      "CLDN1 vs UCP1 rho (BAT)"
    ),
    expected = c(
      N_SAMPLES_EXPECTED,
      N_SUBJECTS,
      "~217000 (GEO file)",
      "15000-25000",
      "FALSE (ideally)",
      "< 10",
      paste0(EXPECTED_UP_BAT, " +/- ", DEG_COUNT_TOLERANCE * 100, "%"),
      "UP_IN_BAT_SIGNIFICANT",
      "> 0",
      "< 0.05",
      "TRUE",
      "user-defined",
      "positive"
    ),
    observed = c(
      sanity$n_samples,
      sanity$n_subjects,
      sanity$n_genes_raw,
      sanity$n_genes_filtered,
      ifelse(!is.null(sanity$second_pass_filter), sanity$second_pass_filter, "N/A"),
      ifelse(!is.null(sanity$deseq_runtime_min), sanity$deseq_runtime_min, "N/A"),
      sanity$n_up_bat,
      ifelse(is.null(sanity$ucp1_verdict), "NA", sanity$ucp1_verdict),
      ifelse(isTRUE(sanity$ucp1_found), round(sanity$ucp1_lfc, 3), "NOT_FOUND"),
      ifelse(isTRUE(sanity$ucp1_found), signif(sanity$ucp1_padj, 4), "NOT_FOUND"),
      sanity$cldn1_found,
      sanity$cldn1_verdict,
      ifelse(!is.null(sanity$cldn1_ucp1_rho),
             round(sanity$cldn1_ucp1_rho, 3), "NA")
    ),
    stringsAsFactors = FALSE
  )

  sanity_file <- paste0("sanity_checks_", run_tag, ".csv")
  write.csv(sanity_df, file = file.path(outdir, "tables", sanity_file),
            row.names = FALSE)
  cat("[SAVED]", sanity_file, "\n")
  verify_output_file(file.path(outdir, "tables", sanity_file), "Sanity CSV")
  print(sanity_df)

  ## ---- 4.12) Session info --------------------------------
  session_file <- file.path(outdir, "logs", paste0("sessionInfo_", run_tag, ".txt"))
  sink(session_file)
  cat("=== Session Info ===\n")
  cat("Run tag:", run_tag, "\n")
  cat("CLDN1 verdict:", sanity$cldn1_verdict, "\n")
  cat("UCP1  verdict:", ifelse(is.null(sanity$ucp1_verdict), "NA",
                                sanity$ucp1_verdict), "\n\n")
  print(sessionInfo())
  sink()
  verify_output_file(session_file, "Session info")

  ## ---- 4.13) Parameter log -------------------------------
  params_lines <- c(
    "=== Parameters ===",
    paste0("GEO_ACCESSION     = ", GEO_ACCESSION),
    paste0("SPECIES           = ", SPECIES),
    paste0("PADJ_CUTOFF       = ", PADJ_CUTOFF),
    paste0("LFC_CUTOFF        = ", LFC_CUTOFF),
    paste0("MIN_COUNTS_PER_GENE = ", MIN_COUNTS_PER_GENE),
    paste0("MIN_SAMPLES_DETECTED = ", MIN_SAMPLES_DETECTED),
    paste0("SEED              = ", SEED),
    paste0("Design            = ~ subject + tissue (paired)")
  )
  params_file <- file.path(outdir, "logs", paste0("params_", run_tag, ".txt"))
  writeLines(params_lines, con = params_file)
  verify_output_file(params_file, "Params log")

  if (use_save_core && exists("finalize_run")) {
    finalize_run(
      outdir,
      status = "success",
      summary_lines = c(
        paste0("run_tag: ", run_tag),
        paste0("cldn1_verdict: ", sanity$cldn1_verdict),
        paste0("n_samples: ", sanity$n_samples),
        paste0("n_total_deg: ", sanity$n_up_bat + sanity$n_up_wat)
      )
    )
  }

}, error = function(e) {
  ## ---- Error handler -------------------------------------
  cat("\n[ERROR] ", conditionMessage(e), "\n\n")
  err_file <- paste0("ERROR_", run_tag, ".txt")
  writeLines(
    c("Script error:", conditionMessage(e), "",
      "Traceback:", capture.output(traceback())),
    con = file.path(outdir, "logs", err_file)
  )
  cat("[SAVED] Error log:", err_file, "\n")

  if (use_save_core && exists("finalize_run")) {
    finalize_run(
      outdir,
      status = "error",
      summary_lines = c(
        paste0("run_tag: ", run_tag),
        paste0("error: ", conditionMessage(e))
      )
    )
  }
})

## ---- 5) Cleanup ------------------------------------------
if (use_save_core && exists("dedupe")) {
  try(dedupe(outdir), silent = TRUE)
}

cat("\n[SANITY] Output inventory\n")
tables_written <- list.files(file.path(outdir, "tables"), full.names = FALSE)
plots_written  <- list.files(file.path(outdir, "plots"), full.names = FALSE)
cat("  Tables:", length(tables_written), "files\n")
if (length(tables_written) > 0) {
  cat("   -", paste(head(tables_written, 10), collapse = "\n   - "), "\n")
}
cat("  Plots :", length(plots_written), "files\n")
if (length(plots_written) > 0) {
  cat("   -", paste(head(plots_written, 10), collapse = "\n   - "), "\n")
}

cat("\n")
cat("##########################################################\n")
cat("##  DONE — Part 1 complete                              ##\n")
cat("##  Run directory: ", basename(outdir), "\n", sep = "")
cat("##  UCP1  verdict: ", ifelse(is.null(sanity$ucp1_verdict), "NA",
                                   sanity$ucp1_verdict), "\n", sep = "")
cat("##  CLDN1 verdict: ", sanity$cldn1_verdict, "\n", sep = "")
cat("##########################################################\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
