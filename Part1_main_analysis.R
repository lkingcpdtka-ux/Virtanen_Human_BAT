## =========================================================
## Part 1  Main Analysis: Virtanen Human BAT RNA-seq
## =========================================================
## Goal : Determine whether CLDN1 is increased in human BAT
##        relative to paired subcutaneous WAT
## Data : GSE113764 (Virtanen et al. Cell Metabolism 2018)
## Design: Paired  14 subjects  2 tissues (BAT / WAT)
## Model : ~ subject + tissue  (paired DESeq2)
## =========================================================

cat("=== Part 1: Virtanen Human BAT  CLDN1 validation ===\n")
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
  ## Transliterate accented characters to ASCII (a, o, etc.)
  ## before stripping, so Finnish names like "Poi" match correctly.
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
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
## 4) MAIN ANALYSIS  wrapped in tryCatch
## ==========================================================
tryCatch({

  ## ---- 4.1) Download / load GEO data ---------------------
  cat("== 4.1  Fetching GEO data:", GEO_ACCESSION, "==\n")

  ## Clean slate  prevent stale objects from a previous run in the same
  ## R session from silently skipping the rebuild.
  if (exists("se", inherits = FALSE))  rm(se)
  if (exists("dds", inherits = FALSE)) rm(dds)
  if (exists("vsd", inherits = FALSE)) rm(vsd)

  geo_cache <- file.path(cache_dir, paste0(GEO_ACCESSION, "_SE.rds"))

  ## Also check for stale SE caches in ALL previous savepoint directories
  ## (when re-running after a parsing fix, old caches cause silent failures)
  savepoints_root <- file.path(getwd(), "savepoints")
  if (dir.exists(savepoints_root)) {
    all_se_caches <- list.files(savepoints_root,
                                 pattern = paste0(GEO_ACCESSION, "_SE\\.rds$"),
                                 recursive = TRUE, full.names = TRUE)
    for (old_cache in all_se_caches) {
      old_se <- tryCatch(readRDS(old_cache), error = function(e) NULL)
      if (!is.null(old_se) && (nrow(old_se) > 50000 || !("UCP1" %in% rownames(old_se)))) {
        cat("[CACHE] Purging stale SE cache:", old_cache, "\n")
        file.remove(old_cache)
        ## Also remove associated DESeq2 and VST caches in the same directory
        old_cache_dir <- dirname(old_cache)
        stale_dds <- list.files(old_cache_dir, pattern = "^dds_.*\\.rds$", full.names = TRUE)
        stale_vst <- list.files(old_cache_dir, pattern = "^vst_.*\\.rds$", full.names = TRUE)
        for (f in c(stale_dds, stale_vst)) {
          cat("[CACHE] Purging stale:", basename(f), "\n")
          file.remove(f)
        }
      }
      if (!is.null(old_se)) rm(old_se)
    }
  }

  if (file.exists(geo_cache)) {
    cat("[CACHE] Loading cached GEO data\n")
    se <- readRDS(geo_cache)
    ## Validate: must have recognisable gene symbols and reasonable dimensions
    if (nrow(se) > 50000 || !("UCP1" %in% rownames(se))) {
      cat("[CACHE] STALE: cached SE has", nrow(se), "features",
          "(UCP1 present:", "UCP1" %in% rownames(se), ")\n")
      cat("[CACHE] Removing stale cache and rebuilding...\n")
      file.remove(geo_cache)
      ## Also clear downstream caches
      for (f in list.files(cache_dir, pattern = "^(dds_|vst_).*\\.rds$", full.names = TRUE)) {
        file.remove(f)
      }
      rm(se)
    } else {
      cat("[CACHE] Valid: ", nrow(se), " genes, UCP1 present\n", sep = "")
    }
  }

  if (!exists("se")) {
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

      ## ---- Smart gene-ID extraction ----
      ## GEO count files vary in structure.  Common layouts:
      ##   (a) col1 = ENSG/symbol, rest = counts
      ##   (b) col1 = transcript accession, "Symbol" col has gene names, rest = counts
      ## We prefer a "Symbol" column if it contains recognisable gene names.
      cn <- colnames(counts_raw)

      ## Look for a dedicated Symbol / GeneName column
      symbol_col_idx <- grep("^(symbol|gene[._]?name|gene[._]?symbol|hgnc)$",
                              cn, ignore.case = TRUE)

      if (length(symbol_col_idx) > 0) {
        ## Found a Symbol column  use it as the gene identifier
        sym_idx <- symbol_col_idx[1]
        gene_col <- as.character(counts_raw[[sym_idx]])
        cat("[INFO] Using column '", cn[sym_idx], "' as gene identifier\n", sep = "")

        ## Drop the first column (accession) and the Symbol column.
        ## All remaining column filtering (numeric detection, BAT/WAT
        ## pattern matching) is handled by select_numeric_count_columns().
        cols_to_drop <- sort(unique(c(1, sym_idx)))
        dropped_names <- cn[cols_to_drop]
        cat("[INFO] Dropping identifier columns: ", paste(dropped_names, collapse = ", "), "\n", sep = "")
        counts_raw <- counts_raw[, -cols_to_drop, drop = FALSE]
      } else {
        ## No Symbol column  fall back to col1 as gene identifier
        gene_col <- as.character(counts_raw[[1]])
        counts_raw <- counts_raw[, -1, drop = FALSE]
        cat("[INFO] No 'Symbol' column found; using column 1 as gene identifier\n")
      }

      counts_raw <- select_numeric_count_columns(as.data.frame(counts_raw))

      ## ---- Gene-ID quality check ----
      ## Warn if gene_col looks like transcript accessions, not gene symbols
      sample_ids_check <- head(gene_col[gene_col != "" & !is.na(gene_col)], 50)
      n_look_like_accessions <- sum(grepl("^(AK|BC|AJ|AF|AB|AL|ENST|NM_|NR_|XM_|XR_)[0-9]",
                                           sample_ids_check))
      if (n_look_like_accessions > length(sample_ids_check) * 0.3) {
        cat("[WARNING] Gene IDs look like transcript accessions, not gene symbols!\n")
        cat("          Sample IDs: ", paste(head(sample_ids_check, 5), collapse = ", "), "\n", sep = "")
        cat("          This will prevent matching to gene names like UCP1/CLDN1.\n")
        cat("          Check that the correct column was used as gene identifier.\n")
      }

      ## Remove rows with empty/NA gene symbols
      valid_gene <- !is.na(gene_col) & nzchar(gene_col) & gene_col != "---" & gene_col != "NA"
      if (sum(!valid_gene) > 0) {
        cat("[INFO] Removing", sum(!valid_gene), "rows with empty/NA gene symbols\n")
        gene_col   <- gene_col[valid_gene]
        counts_raw <- counts_raw[valid_gene, , drop = FALSE]
      }

      ## Handle duplicate gene IDs by summing counts (transcript -> gene aggregation)
      if (anyDuplicated(gene_col)) {
        n_dup <- sum(duplicated(gene_col))
        n_unique <- length(unique(gene_col))
        cat("[INFO] Found", n_dup, "duplicate gene names (",
            length(gene_col), "transcripts -> ", n_unique, "unique genes)\n")
        cat("[INFO] Aggregating to gene level by summing transcript counts\n")
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
      ## Standard RNA-seq: 10M-200M reads per sample -> lib sizes ~10M-200M.
      ## If lib sizes are >1 billion, data may be transcript-level or pre-scaled.
      if (median(lib_sizes_raw) > 1e9) {
        cat("[A1] WARNING: Median library size is ", format(median(lib_sizes_raw), big.mark = ","),
            "  this is FAR above typical RNA-seq (10M-200M).\n", sep = "")
        cat("     Possible causes:\n")
        cat("     - Data is at transcript level (needs gene-level aggregation)\n")
        cat("     - Genomatix weighted counts (not raw read counts)\n")
        cat("     - Pre-scaled or multiplied values\n")
        cat("     DESeq2 may still work, but interpret size factors carefully.\n")
      } else if (median(lib_sizes_raw) >= 1e6 && median(lib_sizes_raw) <= 5e8) {
        cat("[A1] OK  library sizes in typical RNA-seq range\n")
      }

      ## A2. Gene ID classification
      gene_ids <- rownames(counts_raw)
      n_total_genes <- length(gene_ids)
      n_hgnc_like   <- sum(grepl("^[A-Z][A-Z0-9-]*$", gene_ids) &
                            nchar(gene_ids) <= 15 &
                            !grepl("^(AK|BC|AJ|AF|NM_|NR_|XM_|XR_|ENST|ENSG)", gene_ids))
      n_genbank     <- sum(grepl("^(AK|BC|AJ|AF|AB|AL|CR|BX)[0-9]{5,}", gene_ids))
      n_refseq      <- sum(grepl("^(NM_|NR_|XM_|XR_)", gene_ids))
      n_enst        <- sum(grepl("^ENST", gene_ids))
      n_ensg        <- sum(grepl("^ENSG", gene_ids))
      n_numeric     <- sum(grepl("^[0-9]+$", gene_ids))
      n_other       <- n_total_genes - n_hgnc_like - n_genbank - n_refseq - n_enst - n_ensg - n_numeric

      cat("[A2] Gene ID classification (", n_total_genes, " total):\n", sep = "")
      cat("     HGNC-like symbols : ", n_hgnc_like, " (", round(100 * n_hgnc_like / n_total_genes, 1), "%)\n", sep = "")
      cat("     GenBank accessions: ", n_genbank, " (", round(100 * n_genbank / n_total_genes, 1), "%)\n", sep = "")
      cat("     RefSeq IDs        : ", n_refseq, " (", round(100 * n_refseq / n_total_genes, 1), "%)\n", sep = "")
      cat("     ENSEMBL transcr.  : ", n_enst, " (", round(100 * n_enst / n_total_genes, 1), "%)\n", sep = "")
      cat("     ENSEMBL genes     : ", n_ensg, " (", round(100 * n_ensg / n_total_genes, 1), "%)\n", sep = "")
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
                     !grepl("^(AK|BC|AJ|AF|AB|AL|CR|BX)[0-9]{5,}", gene_ids) &
                     !grepl("^(ENST|ENSG)[0-9]", gene_ids)
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
          cat("[A5] ", kg, " present  row sum = ",
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

    ## Step 1: match count column names -> GEO accession IDs
    map_idx <- match(counts_norm, pdata_norm)
    matched <- !is.na(map_idx)
    cat("[SANITY] Sample-name matching (GEO accession):", sum(matched), "of",
        length(counts_norm), "count columns matched\n")

    ## Step 2: if that failed, try matching against GEO sample titles
    ##   GEO pdata$title often contains the original sample names used
    ##   as column headers in the count file (e.g. "JojuBAT", "JojuWAT").
    if (sum(matched) < 2 && "title" %in% colnames(pdata)) {
      title_norm <- normalize_sample_id(pdata$title)
      map_idx <- match(counts_norm, title_norm)
      matched <- !is.na(map_idx)
      cat("[SANITY] Sample-name matching (GEO title):", sum(matched), "of",
          length(counts_norm), "count columns matched\n")
    }

    ## Step 3: order-based fallback (last resort)  only if dimensions agree
    if (sum(matched) < 2) {
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

    original_count_colnames <- colnames(counts_raw)
    counts_raw <- counts_raw[, matched, drop = FALSE]
    map_idx <- map_idx[matched]
    pdata <- pdata[map_idx, , drop = FALSE]

    original_count_colnames <- original_count_colnames[matched]
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
    ## Use safe subject IDs for model-matrix column names (avoid DESeq2 note
    ## about non-standard characters, e.g. diacritics in subject labels).
    subject_safe <- normalize_sample_id(pdata$subject)
    col_data <- DataFrame(
      tissue  = factor(pdata$tissue, levels = c(CONDITION_REF, CONDITION_TEST)),
      subject = factor(subject_safe),
      subject_label = as.character(pdata$subject),
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

  ## Metadata mapping export (GSM -> subject/tissue + original count column)
  sample_map_df <- data.frame(
    sample_id = colnames(se),
    gsm = colnames(se),
    subject = as.character(colData(se)$subject_label),
    tissue = as.character(colData(se)$tissue),
    original_count_column = if (exists("original_count_colnames")) original_count_colnames else colnames(se),
    stringsAsFactors = FALSE
  )
  sample_map_file <- file.path(outdir, "tables", "sample_metadata_mapping.csv")
  write.csv(sample_map_df, sample_map_file, row.names = FALSE)
  cat("[SAVED] sample_metadata_mapping.csv\n")
  verify_output_file(sample_map_file, "Sample metadata map")
  cat("[META] Subject count from mapping:", length(unique(sample_map_df$subject)), "\n")
  cat("[META] Tissue counts from mapping:\n")
  print(table(sample_map_df$tissue, useNA = "ifany"))
  meta_cols <- grep("characteristics|description|title|source", colnames(pdata),
                    ignore.case = TRUE, value = TRUE)
  if (length(meta_cols) > 0) {
    meta_txt <- unlist(lapply(meta_cols, function(cc) as.character(pdata[[cc]])))
    qc_hits <- grep("exclude|fail|rin|quality", meta_txt, ignore.case = TRUE, value = TRUE)
    if (length(qc_hits) > 0) {
      cat("[META] Potential QC/exclusion metadata entries detected (review only):\n")
      print(unique(head(qc_hits, 10)))
    } else {
      cat("[META] No explicit QC-fail/exclude metadata flag detected in GEO annotations.\n")
    }
  }

  ## ---- Sanity: are these raw counts or pre-normalised? ----
  ## Raw counts are integers >= 0 with most values in hundreds-thousands range.
  ## RPKM/FPKM/TPM are floats, typically 0-100 with many decimals.
  ## If the data looks pre-normalised, DESeq2 will give wrong results.
  cat("\n--- SANITY: Raw count validation ---\n")
  count_sample <- assay(se, "counts")
  count_vals   <- as.numeric(count_sample[sample(nrow(count_sample), min(5000, nrow(count_sample))),
                                           sample(ncol(count_sample), min(5, ncol(count_sample)))])
  count_vals   <- count_vals[!is.na(count_vals)]

  frac_integer <- mean(count_vals == round(count_vals), na.rm = TRUE)
  frac_zero    <- mean(count_vals == 0, na.rm = TRUE)
  count_max    <- max(count_vals, na.rm = TRUE)
  count_median <- median(count_vals, na.rm = TRUE)

  cat("[RAW-CHECK] Fraction integer-valued: ", round(frac_integer, 3), "\n", sep = "")
  cat("[RAW-CHECK] Fraction zero          : ", round(frac_zero, 3), "\n", sep = "")
  cat("[RAW-CHECK] Median value           : ", round(count_median, 2), "\n", sep = "")
  cat("[RAW-CHECK] Max value              : ", format(count_max, big.mark = ","), "\n", sep = "")

  if (frac_integer < 0.9) {
    cat("[RAW-CHECK] WARNING: Data appears to be NORMALISED (RPKM/FPKM/TPM),\n")
    cat("            not raw counts. DESeq2 requires raw integer counts.\n")
    cat("            Fraction integer:", round(frac_integer, 3), "\n")
    sanity$data_type_warning <- "LIKELY_NORMALIZED"
  } else if (count_max < 500 && count_median < 10) {
    cat("[RAW-CHECK] WARNING: Values look unusually low for raw counts.\n")
    cat("            Possible CPM or log-transformed data.\n")
    sanity$data_type_warning <- "UNUSUALLY_LOW"
  } else {
    cat("[RAW-CHECK] OK  data looks like raw integer counts\n")
    sanity$data_type_warning <- "NONE"
  }

  ## ---- Sanity: paired design validation ----
  ## For a paired design (~ subject + tissue), every subject must have
  ## both BAT and WAT.  Unpaired subjects break the model.
  cat("\n--- SANITY: Paired design validation ---\n")
  subj_tissue <- table(se$subject, se$tissue)
  has_both    <- rowSums(subj_tissue > 0) == 2
  unpaired    <- names(has_both)[!has_both]

  if (length(unpaired) > 0) {
    cat("[PAIRED] WARNING: ", length(unpaired), " subject(s) lack one tissue:\n", sep = "")
    for (u in unpaired) {
      present <- colnames(subj_tissue)[subj_tissue[u, ] > 0]
      cat("         ", u, "  only has: ", paste(present, collapse = ", "), "\n", sep = "")
    }
    cat("[PAIRED] Removing unpaired subjects to avoid DESeq2 model errors\n")
    keep_paired <- !(se$subject %in% unpaired)
    se <- se[, keep_paired]
    se$subject <- droplevels(se$subject)
    cat("[PAIRED] Remaining: ", ncol(se), " samples from ",
        nlevels(se$subject), " complete pairs\n", sep = "")
  } else {
    cat("[PAIRED] OK  all ", nlevels(se$subject),
        " subjects have both ", CONDITION_REF, " and ", CONDITION_TEST, "\n", sep = "")
  }

  ## Check for duplicate tissue per subject (would also break paired model)
  max_per_cell <- max(subj_tissue[has_both, ])
  if (max_per_cell > 1) {
    cat("[PAIRED] WARNING: Some subject-tissue combos have >1 sample.\n")
    cat("         Max samples per cell:", max_per_cell, "\n")
    cat("         DESeq2 paired design expects exactly 1 per cell.\n")
  }

  ## Optional explicit sample exclusions (disabled by default)
  excluded_samples_df <- data.frame(sample_id = character(0), reason = character(0), stringsAsFactors = FALSE)
  if (exists("EXCLUDE_SAMPLES") && length(EXCLUDE_SAMPLES) > 0) {
    drop_idx <- colnames(se) %in% EXCLUDE_SAMPLES
    if (any(drop_idx)) {
      excluded_samples_df <- data.frame(
        sample_id = colnames(se)[drop_idx],
        reason = "manual_EXCLUDE_SAMPLES",
        stringsAsFactors = FALSE
      )
      se <- se[, !drop_idx]
      se$subject <- droplevels(se$subject)
      cat("[EXCLUDE] Applied EXCLUDE_SAMPLES:", paste(excluded_samples_df$sample_id, collapse = ", "), "\n")
    }
  }
  if (nrow(excluded_samples_df) == 0) {
    cat("[EXCLUDE] No exclusions applied.\n")
  }

  sanity$n_genes_raw   <- nrow(se)
  sanity$n_samples     <- ncol(se)
  sanity$n_subjects    <- nlevels(se$subject)
  sanity$tissues       <- levels(se$tissue)

  ## ---- 4.2) Pre-filtering --------------------------------
  ## Standard approach (Nature/Cell Metabolism convention):
  ##   Keep genes with CPM >= 1 in at least the smallest group size.
  ## This is purely count-driven  no gene-ID pattern matching needed.
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

  ## ---- Scaling gate (off by default unless final matrix is absurdly large) ----
  lib_sizes_unscaled <- colSums(counts_mat)
  median_lib_unscaled <- median(lib_sizes_unscaled)
  scale_gate <- if (exists("LIBSIZE_SCALE_GATE")) LIBSIZE_SCALE_GATE else 2e8
  target_lib_size <- if (exists("TARGET_LIB_SIZE")) TARGET_LIB_SIZE else 2e7

  auto_scale <- median_lib_unscaled > scale_gate
  if (exists("DO_SCALE_COUNTS") && !is.null(DO_SCALE_COUNTS)) {
    do_scale_counts <- isTRUE(DO_SCALE_COUNTS)
    cat("[SCALE] DO_SCALE_COUNTS explicitly set to", do_scale_counts, "\n")
  } else {
    do_scale_counts <- auto_scale
    cat("[SCALE] Auto gate: median library size =",
        format(median_lib_unscaled, big.mark = ","),
        "| threshold =", format(scale_gate, big.mark = ","),
        "-> DO_SCALE_COUNTS =", do_scale_counts, "\n")
  }

  sanity$libsize_unscaled_min <- min(lib_sizes_unscaled)
  sanity$libsize_unscaled_median <- median_lib_unscaled
  sanity$libsize_unscaled_max <- max(lib_sizes_unscaled)

  if (do_scale_counts) {
    scale_factor <- median_lib_unscaled / target_lib_size
    cat("[SCALE] WARNING: Extremely large final gene-level libraries detected.\n")
    cat("[SCALE] Scaling counts by factor", round(scale_factor, 4),
        "to target median", format(target_lib_size, big.mark = ","), "\n")
    counts_mat <- round(counts_mat / scale_factor)
    sanity$count_scale_factor <- round(scale_factor, 4)
  } else {
    sanity$count_scale_factor <- 1
    cat("[SCALE] No scaling applied (preferred for integer count-like data).\n")
  }

  lib_sizes_scaled <- colSums(counts_mat)
  sanity$libsize_scaled_min <- min(lib_sizes_scaled)
  sanity$libsize_scaled_median <- median(lib_sizes_scaled)
  sanity$libsize_scaled_max <- max(lib_sizes_scaled)

  ## Final integer overflow check (should be rare/none after scaling)
  int_max <- .Machine$integer.max
  n_overflow <- sum(counts_mat > int_max)
  if (n_overflow > 0) {
    cat("[SANITY] WARNING: ", n_overflow,
        " values still exceed integer max after scaling  capping\n", sep = "")
    counts_mat[counts_mat > int_max] <- int_max
  }

  storage.mode(counts_mat) <- "integer"
  assay(se, "counts") <- counts_mat

  filter_mode <- if (exists("FILTER_MODE")) FILTER_MODE else "filterByExpr"
  if (LOW_COUNT_FILTER && filter_mode != "deseq2_default") {
    design_group <- se$tissue
    min_group    <- min(table(design_group))
    lib_sizes    <- colSums(counts_mat)
    med_lib      <- median(lib_sizes)

    cat("[FILTER] Smallest group size:", min_group, "samples\n")
    cat("[FILTER] Median library size:", format(med_lib, big.mark = ","), "\n")

    ## Step 1: fast pre-screen  remove genes with zero counts in ALL samples.
    ##         This is instant and eliminates the bulk of ~217K features.
    nonzero <- rowSums(counts_mat > 0) > 0
    n_allzero <- sum(!nonzero)
    cat("[FILTER] Removing", n_allzero, "genes with zero counts across all samples\n")
    counts_mat_nz <- counts_mat[nonzero, , drop = FALSE]
    cat("[FILTER] After zero removal:", nrow(counts_mat_nz), "genes remain\n")

    cpm_mat  <- t(t(counts_mat_nz) / lib_sizes * 1e6)
    keep_cpm <- rowSums(cpm_mat >= 1) >= min_group
    cat("[FILTER] Genes with CPM >= 1 in >=", min_group, "samples:", sum(keep_cpm), "\n")

    keep <- keep_cpm
    if (filter_mode %in% c("filterByExpr", "intersection") && requireNamespace("edgeR", quietly = TRUE)) {
      min_count_cpm1 <- max(MIN_COUNTS_PER_GENE, ceiling(med_lib / 1e6))
      keep_expr <- edgeR::filterByExpr(counts_mat_nz, group = design_group,
                                       min.count = min_count_cpm1)
      cat("[FILTER] edgeR::filterByExpr (min.count=", min_count_cpm1,
          ") keeps:", sum(keep_expr), "\n", sep = "")
      if (filter_mode == "filterByExpr") {
        keep <- keep_expr
      } else if (filter_mode == "intersection") {
        keep <- keep_cpm & keep_expr
      }
    } else if (filter_mode == "filterByExpr") {
      cat("[FILTER] edgeR not available; falling back to CPM filter\n")
    }

    cat("[FILTER] Mode:", filter_mode, "| kept", sum(keep), "genes after filtering\n")

    ## Apply filter back to full SE (must index into original row space)
    keep_names <- rownames(counts_mat_nz)[keep]
    se_filt <- se[keep_names, ]
    cat("[FILTER] Final: kept", nrow(se_filt), "of", nrow(se), "genes",
        "(removed", nrow(se) - nrow(se_filt), ")\n")
  } else {
    se_filt <- se
    cat("[FILTER] Mode:", ifelse(LOW_COUNT_FILTER, filter_mode, "disabled"),
        "-> using unfiltered gene set\n")
  }

  sanity$n_genes_after_id_filter <- nrow(se)  # no separate ID filter; same as raw

  ## ---- Sanity block B: Post-filter diagnostics ----
  cat("\n--- SANITY BLOCK B: Post-filter diagnostics ---\n")

  sanity$n_genes_filtered <- nrow(se_filt)
  filt_counts <- assay(se_filt, "counts")

  ## B1. Dimension check
  cat("[B1] Filtered matrix: ", nrow(se_filt), " genes x ", ncol(se_filt), " samples\n", sep = "")

  ## B2. Hard cap  if still >30K genes, something is wrong; do a second-pass
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

  ## ---- 4.3) DESeq2  paired design ----------------------
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
    if (!file.exists("/proc/meminfo")) {
      NA_real_
    } else {
      meminfo <- suppressWarnings(readLines("/proc/meminfo", n = 3))
      avail_line <- grep("MemAvailable", meminfo, value = TRUE)
      if (length(avail_line) > 0) {
        avail_kb <- as.numeric(gsub("[^0-9]", "", avail_line))
        round(avail_kb / 1024)
      } else {
        NA_real_
      }
    }
  }, error = function(e) NA_real_)
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
        cat("[C4] Cache gene count matches (", cached_ngenes, ")  using cache\n", sep = "")
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

  ## ---- Sanity block D: Post-DESeq2 model diagnostics ----
  cat("\n--- SANITY BLOCK D: DESeq2 model diagnostics ---\n")

  ## D1. Size factors  should be ~0.5-2.0 for most samples.
  ##     Extreme values indicate library-size or composition problems.
  sf <- sizeFactors(dds)
  cat("[D1] Size factors:\n")
  cat("     min =", round(min(sf), 3), " median =", round(median(sf), 3),
      " max =", round(max(sf), 3), "\n")
  sf_outliers <- names(sf)[sf < 0.1 | sf > 10]
  if (length(sf_outliers) > 0) {
    cat("[D1] WARNING: Extreme size factors in:", paste(sf_outliers, collapse = ", "), "\n")
    cat("     This may indicate failed library prep or contamination.\n")
  } else {
    cat("[D1] OK  all size factors within normal range (0.1-10)\n")
  }
  sanity$size_factor_range <- paste0(round(min(sf), 3), "-", round(max(sf), 3))

  ## D2. Dispersion estimates  check for unusual patterns.
  ##     Median dispersion for human RNA-seq is typically 0.01-0.3.
  dispersions <- mcols(dds)$dispGeneEst
  dispersions <- dispersions[!is.na(dispersions)]
  if (length(dispersions) > 0) {
    cat("[D2] Gene-level dispersions:\n")
    cat("     median =", round(median(dispersions), 4),
        " mean =", round(mean(dispersions), 4), "\n")
    cat("     Quantiles: 5%=", round(quantile(dispersions, 0.05), 4),
        " 25%=", round(quantile(dispersions, 0.25), 4),
        " 75%=", round(quantile(dispersions, 0.75), 4),
        " 95%=", round(quantile(dispersions, 0.95), 4), "\n")
    if (median(dispersions) > 1) {
      cat("[D2] WARNING: High median dispersion suggests noisy data or poor model fit.\n")
    } else {
      cat("[D2] OK  dispersion range typical for human tissue RNA-seq\n")
    }
  }

  ## D3. Convergence  check for genes where DESeq2 maxed out iterations.
  n_not_converged <- sum(mcols(dds)$betaConv == FALSE, na.rm = TRUE)
  cat("[D3] Genes that did not converge:", n_not_converged, "of", nrow(dds), "\n")
  if (n_not_converged > nrow(dds) * 0.05) {
    cat("[D3] WARNING: >5% non-convergence. Model may be mis-specified.\n")
  } else {
    cat("[D3] OK  convergence rate normal\n")
  }

  ## D4. Cook's distance  flag outlier samples that dominate results.
  ##     High Cook's for many genes in one sample = probable outlier.
  cooks_mat <- assays(dds)[["cooks"]]
  if (!is.null(cooks_mat)) {
    cooks_per_sample <- apply(as.matrix(cooks_mat), 2, median, na.rm = TRUE)
    names(cooks_per_sample) <- colnames(dds)
    cooks_threshold <- qf(0.99, ncol(model.matrix(design(dds), colData(dds))),
                          ncol(dds) - ncol(model.matrix(design(dds), colData(dds))))
    cat("[D4] Cook's distance (median per sample):\n")
    cooks_df <- data.frame(
      sample = names(cooks_per_sample),
      tissue = as.character(dds$tissue),
      median_cooks = round(cooks_per_sample, 4),
      stringsAsFactors = FALSE
    )
    cooks_df <- cooks_df[order(-cooks_df$median_cooks), ]
    print(head(cooks_df, 6), row.names = FALSE)
    high_cooks <- cooks_df$sample[cooks_df$median_cooks > 1]
    if (length(high_cooks) > 0) {
      cat("[D4] WARNING: Samples with high Cook's distance: ",
          paste(high_cooks, collapse = ", "), "\n", sep = "")
      cat("     Consider removing these outliers and re-running.\n")
    } else {
      cat("[D4] OK  no sample-level outliers by Cook's distance\n")
    }
  }

  cat("--- END SANITY BLOCK D ---\n\n")

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

  ## ---- Sanity block E: P-value distribution ----
  ## A healthy DE analysis has a uniform distribution of p-values
  ## (null genes) with a spike near 0 (true DE genes).
  ## Anti-conservative (U-shape) = batch effects or model problems.
  ## All p ~ 1 = no signal or wrong comparison.
  cat("\n--- SANITY BLOCK E: P-value distribution ---\n")
  pvals <- res_df$pvalue[!is.na(res_df$pvalue)]
  n_total_pvals <- length(pvals)
  n_na_pvals    <- sum(is.na(res_df$pvalue))

  cat("[E1] Total genes with p-values:", n_total_pvals,
      " (NA:", n_na_pvals, ")\n")

  ## Bin p-values into deciles
  pval_bins <- cut(pvals, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  pval_tab  <- table(pval_bins)
  cat("[E2] P-value histogram (deciles):\n")
  for (b in names(pval_tab)) {
    bar_len <- round(pval_tab[b] / n_total_pvals * 50)
    cat(sprintf("     %-12s %5d (%4.1f%%) %s\n",
                b, pval_tab[b],
                100 * pval_tab[b] / n_total_pvals,
                paste(rep("#", bar_len), collapse = "")))
  }

  ## Diagnostic: expect bin[0,0.1] to be the largest (DE signal).
  ## If bin[0.9,1] is larger, something is wrong.
  first_bin <- pval_tab[1]
  last_bin  <- pval_tab[length(pval_tab)]
  mid_bins_mean <- mean(pval_tab[2:(length(pval_tab) - 1)])

  if (first_bin < mid_bins_mean) {
    cat("[E3] WARNING: No enrichment of small p-values. Possible issues:\n")
    cat("     - Wrong contrast or reference level\n")
    cat("     - Pre-normalised data passed to DESeq2\n")
    cat("     - Severe batch effects dominating the signal\n")
    sanity$pval_distribution <- "NO_SIGNAL"
  } else if (last_bin > 2 * mid_bins_mean) {
    cat("[E3] WARNING: Anti-conservative p-value distribution (spike near 1).\n")
    cat("     Possible batch effects, misspecified model, or violated assumptions.\n")
    sanity$pval_distribution <- "ANTI_CONSERVATIVE"
  } else {
    cat("[E3] OK  p-value distribution looks healthy (spike near 0, flat elsewhere)\n")
    sanity$pval_distribution <- "HEALTHY"
  }

  ## E4. NA/filtered gene accounting
  n_padj_na      <- sum(is.na(res_df$padj))
  n_outlier_na   <- sum(is.na(res_df$pvalue) & !is.na(res_df$baseMean) & res_df$baseMean > 0)
  n_low_count_na <- sum(is.na(res_df$pvalue) & (is.na(res_df$baseMean) | res_df$baseMean == 0))
  cat("[E4] padj NA (independent filtering):", n_padj_na, "of", nrow(res_df), "\n")
  cat("     pvalue NA (outlier/low count)   :", sum(is.na(res_df$pvalue)), "\n")
  cat("     Effectively tested              :", n_total_pvals, "\n")

  cat("--- END SANITY BLOCK E ---\n\n")

  sanity$n_up_bat <- sum(res_df$padj < PADJ_CUTOFF &
                          res_df$log2FoldChange >= LFC_CUTOFF, na.rm = TRUE)
  sanity$n_up_wat <- sum(res_df$padj < PADJ_CUTOFF &
                          res_df$log2FoldChange <= -LFC_CUTOFF, na.rm = TRUE)

  ## ---- 4.3b) Sensitivity model (unpaired) -----------------
  cat("\n== 4.3b  DESeq2 sensitivity model (unpaired) ==\n")
  dds_unpaired <- DESeqDataSet(se_filt, design = ~ tissue)
  dds_unpaired$tissue <- relevel(dds_unpaired$tissue, ref = CONDITION_REF)
  dds_unpaired <- DESeq(dds_unpaired, quiet = TRUE)
  res_unpaired <- results(dds_unpaired,
                          contrast = c("tissue", CONDITION_TEST, CONDITION_REF),
                          alpha = PADJ_CUTOFF)
  res_unpaired_df <- as.data.frame(res_unpaired) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

  sanity$n_up_bat_unpaired <- sum(res_unpaired_df$padj < PADJ_CUTOFF &
                                   res_unpaired_df$log2FoldChange >= LFC_CUTOFF, na.rm = TRUE)
  sanity$n_up_wat_unpaired <- sum(res_unpaired_df$padj < PADJ_CUTOFF &
                                   res_unpaired_df$log2FoldChange <= -LFC_CUTOFF, na.rm = TRUE)
  cat("[OK] Unpaired model complete:", nrow(res_unpaired_df), "genes tested\n")
  cat("     Up in BAT (padj <", PADJ_CUTOFF, ", LFC >=", LFC_CUTOFF, "):",
      sanity$n_up_bat_unpaired, "\n")
  cat("     Up in WAT (padj <", PADJ_CUTOFF, ", LFC <=", -LFC_CUTOFF, "):",
      sanity$n_up_wat_unpaired, "\n")

  ## Report paper mismatch without data-driven subject exclusion
  cat("[PAPER] Virtanen et al. report", EXPECTED_UP_BAT,
      "BAT-enriched genes from 14 subjects; this script keeps all paired subjects by default.\n")
  cat("[PAPER] Current cohort:", nlevels(se_filt$subject), "subjects |",
      ncol(se_filt), "samples\n")

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
    res_unpaired_df$gene_id_clean <- gsub("\\.\\d+$", "", res_unpaired_df$gene_id)
    res_unpaired_df$symbol <- gene_symbols[res_unpaired_df$gene_id_clean]
    n_mapped <- sum(!is.na(res_df$symbol))
  } else {
    ## Gene IDs are already symbols (e.g. GEO humanBATWAT.txt)
    res_df$symbol <- res_df$gene_id
    res_unpaired_df$symbol <- res_unpaired_df$gene_id
    n_mapped <- sum(res_df$symbol != "")
  }

  ## Fill NA symbols with gene_id
  res_df$symbol[is.na(res_df$symbol)] <- res_df$gene_id[is.na(res_df$symbol)]
  res_unpaired_df$symbol[is.na(res_unpaired_df$symbol)] <- res_unpaired_df$gene_id[is.na(res_unpaired_df$symbol)]

  cat("[OK] Mapped", n_mapped, "of", nrow(res_df), "genes to HGNC symbols\n")

  ## ---- 4.5) CLDN1  the primary question -----------------
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
      cat("  >>> VERDICT: YES  CLDN1 is SIGNIFICANTLY INCREASED in BAT\n")
      if (cldn1$log2FoldChange >= LFC_CUTOFF) {
        cat("      (meets both padj and log2FC thresholds)\n")
      } else {
        cat("      (significant but below LFC cutoff of", LFC_CUTOFF, ")\n")
      }
      sanity$cldn1_verdict <- "UP_IN_BAT_SIGNIFICANT"
    } else if (!is.na(cldn1$log2FoldChange) && cldn1$log2FoldChange > 0) {
      cat("  >>> VERDICT: TREND  CLDN1 shows a positive fold-change in BAT\n")
      cat("      but does NOT reach statistical significance (padj =",
          signif(cldn1$padj, 3), ")\n")
      sanity$cldn1_verdict <- "UP_TREND_NOT_SIGNIFICANT"
    } else if (!is.na(cldn1$log2FoldChange) && cldn1$log2FoldChange < 0) {
      cat("  >>> VERDICT: NO  CLDN1 is LOWER in BAT than WAT\n")
      cat("      (log2FC =", round(cldn1$log2FoldChange, 3), ")\n")
      sanity$cldn1_verdict <- "DOWN_IN_BAT"
    } else {
      cat("  >>> VERDICT: INCONCLUSIVE (NA values in results)\n")
      sanity$cldn1_verdict <- "INCONCLUSIVE"
    }
  }

  cat("\n##########################################################\n\n")

  ## Unpaired sensitivity key-gene lines
  cat("[SENS] Unpaired model key genes:\n")
  for (kg in c("CLDN1", "UCP1")) {
    row_u <- res_unpaired_df[res_unpaired_df$symbol == kg, ]
    if (nrow(row_u) > 0) {
      row_u <- row_u[1, ]
      cat(sprintf("  %-6s | log2FC=%+0.3f | p=%s | padj=%s | baseMean=%0.1f\n",
                  kg,
                  row_u$log2FoldChange,
                  formatC(row_u$pvalue, format = "e", digits = 3),
                  formatC(row_u$padj, format = "e", digits = 3),
                  row_u$baseMean))
    } else {
      cat(" ", kg, "| not found in unpaired model output\n")
    }
  }

  ## ---- 4.5b) UCP1  positive control ----------------------
  cat("\n")
  cat("##########################################################\n")
  cat("##  CONTROL: UCP1 expression (canonical BAT marker)     ##\n")
  cat("##########################################################\n\n")

  ucp1_row <- res_df[res_df$symbol == "UCP1", ]

  if (nrow(ucp1_row) == 0) {
    cat("[!!] UCP1 NOT FOUND  this is unexpected for BAT/WAT data\n")
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
      cat("  >>> CONTROL WARNING: UCP1 not up in BAT  check data integrity\n")
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

  ## ---- Sanity block F: Sample-level QC (on VST data) ----
  cat("\n--- SANITY BLOCK F: Sample-level QC ---\n")

  ## F1. Inter-sample Spearman correlation  should be high (>0.8)
  ##     within tissue, lower between tissues.
  if (ncol(vst_mat) <= 100) {  # skip for very large designs
    cor_mat <- cor(vst_mat, method = "spearman")
    diag(cor_mat) <- NA  # remove self-correlations
    tissue_vec <- as.character(colData(dds)$tissue)

    ## Within-group correlations
    for (tis in c(CONDITION_REF, CONDITION_TEST)) {
      idx <- which(tissue_vec == tis)
      if (length(idx) >= 2) {
        within_cors <- cor_mat[idx, idx][upper.tri(cor_mat[idx, idx])]
        cat("[F1] ", tis, " within-group Spearman: min=", round(min(within_cors, na.rm = TRUE), 3),
            " median=", round(median(within_cors, na.rm = TRUE), 3),
            " max=", round(max(within_cors, na.rm = TRUE), 3), "\n", sep = "")
        if (min(within_cors, na.rm = TRUE) < 0.8) {
          low_pairs <- which(cor_mat[idx, idx] < 0.8 & upper.tri(cor_mat[idx, idx]), arr.ind = TRUE)
          cat("[F1] WARNING: Low within-group correlation in ", tis,
              "  possible outlier sample(s)\n", sep = "")
        }
      }
    }

    ## Between-group correlations
    bat_idx <- which(tissue_vec == CONDITION_TEST)
    wat_idx <- which(tissue_vec == CONDITION_REF)
    if (length(bat_idx) > 0 && length(wat_idx) > 0) {
      between_cors <- as.vector(cor_mat[bat_idx, wat_idx])
      between_cors <- between_cors[!is.na(between_cors)]
      cat("[F1] Between-group Spearman: min=", round(min(between_cors), 3),
          " median=", round(median(between_cors), 3), "\n", sep = "")
    }
  }

  ## F2. PCA variance  PC1 should capture tissue effect.
  ##     If PC1 variance is <10%, the tissue signal may be weak.
  pca_obj <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
  pct_var <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)), 1)
  cat("[F2] PCA variance explained: PC1=", pct_var[1], "% PC2=", pct_var[2],
      "% PC3=", pct_var[3], "%\n", sep = "")

  ## Check if PC1 separates tissues
  pc1_by_tissue <- split(pca_obj$x[, 1], colData(dds)$tissue)
  if (length(pc1_by_tissue) == 2) {
    pc1_ttest <- tryCatch(
      t.test(pc1_by_tissue[[1]], pc1_by_tissue[[2]]),
      error = function(e) NULL
    )
    if (!is.null(pc1_ttest)) {
      cat("[F2] PC1 separates tissues: t-test p =", signif(pc1_ttest$p.value, 3), "\n")
      if (pc1_ttest$p.value > 0.05) {
        cat("[F2] WARNING: PC1 does NOT significantly separate tissues.\n")
        cat("     Tissue effect may be weak or confounded by batch/subject effects.\n")
      } else {
        cat("[F2] OK  PC1 clearly separates BAT from WAT\n")
      }
    }
  }
  sanity$pca_pc1_var <- pct_var[1]

  ## F3. Sample distance outlier detection
  ##     Flag any sample >3 median absolute deviations from group centroid.
  for (tis in c(CONDITION_REF, CONDITION_TEST)) {
    idx <- which(tissue_vec == tis)
    if (length(idx) >= 3) {
      group_mat <- t(vst_mat[, idx])
      centroid  <- colMeans(group_mat)
      dists     <- sqrt(rowSums((group_mat - rep(centroid, each = nrow(group_mat)))^2))
      med_dist  <- median(dists)
      mad_dist  <- mad(dists)
      outliers  <- names(dists)[dists > med_dist + 3 * mad_dist]
      if (length(outliers) > 0) {
        cat("[F3] ", tis, " outlier(s) by Euclidean distance: ",
            paste(outliers, collapse = ", "), "\n", sep = "")
      } else {
        cat("[F3] ", tis, ": no distance outliers detected\n", sep = "")
      }
    }
  }

  cat("--- END SANITY BLOCK F ---\n\n")

  ## ---- 4.8) Save full DE table ---------------------------
  cat("\n== 4.8  Saving results tables ==\n")

  de_file_paired <- "DE_BAT_vs_WAT_paired.csv"
  de_file_unpaired <- "DE_BAT_vs_WAT_unpaired.csv"
  write.csv(res_df, file = file.path(outdir, "tables", de_file_paired), row.names = FALSE)
  write.csv(res_unpaired_df, file = file.path(outdir, "tables", de_file_unpaired), row.names = FALSE)
  cat("[SAVED]", de_file_paired, "\n")
  cat("[SAVED]", de_file_unpaired, "\n")
  verify_output_file(file.path(outdir, "tables", de_file_paired), "DE paired CSV")
  verify_output_file(file.path(outdir, "tables", de_file_unpaired), "DE unpaired CSV")

  deg_df <- res_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) >= LFC_CUTOFF)
  deg_df_unpaired <- res_unpaired_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) >= LFC_CUTOFF)
  deg_df_padj_only <- res_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF)
  deg_df_unpaired_padj_only <- res_unpaired_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF)

  write.csv(deg_df, file.path(outdir, "tables", "DEGs_paired_padj0.05_lfc1.csv"), row.names = FALSE)
  write.csv(deg_df_unpaired, file.path(outdir, "tables", "DEGs_unpaired_padj0.05_lfc1.csv"), row.names = FALSE)
  write.csv(deg_df_padj_only, file.path(outdir, "tables", "DEGs_paired_padj0.05.csv"), row.names = FALSE)
  write.csv(deg_df_unpaired_padj_only, file.path(outdir, "tables", "DEGs_unpaired_padj0.05.csv"), row.names = FALSE)
  verify_output_file(file.path(outdir, "tables", "DEGs_paired_padj0.05_lfc1.csv"), "DEGs paired")
  verify_output_file(file.path(outdir, "tables", "DEGs_unpaired_padj0.05_lfc1.csv"), "DEGs unpaired")
  verify_output_file(file.path(outdir, "tables", "DEGs_paired_padj0.05.csv"), "DEGs paired padj")
  verify_output_file(file.path(outdir, "tables", "DEGs_unpaired_padj0.05.csv"), "DEGs unpaired padj")

  if (nrow(excluded_samples_df) > 0) {
    write.csv(excluded_samples_df, file.path(outdir, "tables", "excluded_samples.csv"), row.names = FALSE)
    verify_output_file(file.path(outdir, "tables", "excluded_samples.csv"), "Excluded samples")
  }

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

  ## Key genes summary across paired and unpaired models
  key_genes <- unique(c("CLDN1", "UCP1", "PPARG", "LEP", "CLDN11", BAT_MARKER_GENES, WAT_MARKER_GENES))
  key_paired <- res_df %>% filter(symbol %in% key_genes) %>% mutate(model = "paired")
  key_unpaired <- res_unpaired_df %>% filter(symbol %in% key_genes) %>% mutate(model = "unpaired")
  key_summary <- bind_rows(key_paired, key_unpaired) %>%
    select(model, gene_id, symbol, baseMean, log2FoldChange, pvalue, padj)
  write.csv(key_summary, file.path(outdir, "tables", "key_genes_summary.csv"), row.names = FALSE)
  verify_output_file(file.path(outdir, "tables", "key_genes_summary.csv"), "Key genes summary")

  ## QC summary table
  qc_summary <- data.frame(
    sample_id = colnames(dds),
    subject = as.character(colData(dds)$subject_label),
    tissue = as.character(colData(dds)$tissue),
    lib_size_unscaled = as.numeric(lib_sizes_unscaled[colnames(dds)]),
    lib_size_final = as.numeric(lib_sizes_scaled[colnames(dds)]),
    size_factor = as.numeric(sizeFactors(dds)),
    median_cooks = if (exists("cooks_df")) cooks_df$median_cooks[match(colnames(dds), cooks_df$sample)] else NA_real_,
    mean_spearman = if (exists("cor_mat")) round(rowMeans(cor_mat, na.rm = TRUE), 4)[colnames(dds)] else NA_real_,
    stringsAsFactors = FALSE
  )
  write.csv(qc_summary, file.path(outdir, "tables", "qc_summary.csv"), row.names = FALSE)
  verify_output_file(file.path(outdir, "tables", "qc_summary.csv"), "QC summary")

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

    volcano_file <- "Volcano_BAT_vs_WAT_paired.png"
    ggsave(file.path(outdir, "plots", volcano_file), plot = volcano,
           width = VOLCANO_W, height = VOLCANO_H, dpi = PLOT_DPI)
    cat("[SAVED]", volcano_file, "\n")
    verify_output_file(file.path(outdir, "plots", volcano_file), "Volcano plot")

    volcano_unpaired <- EnhancedVolcano(
      res_unpaired_df,
      lab             = res_unpaired_df$symbol,
      selectLab       = highlight_genes,
      x               = "log2FoldChange",
      y               = "padj",
      pCutoff         = PADJ_CUTOFF,
      FCcutoff        = LFC_CUTOFF,
      title           = "BAT vs WAT (unpaired sensitivity)",
      subtitle        = paste("CLDN1 highlighted | n =", ncol(se_filt), "samples"),
      caption         = paste("Total genes:", nrow(res_unpaired_df)),
      drawConnectors  = TRUE,
      widthConnectors = 0.5,
      colConnectors   = "grey30",
      max.overlaps    = 30,
      pointSize       = 1.5,
      labSize         = 4.0,
      col             = c("grey80", "forestgreen", "royalblue", "red2")
    )
    volcano_unpaired_file <- "Volcano_BAT_vs_WAT_unpaired.png"
    ggsave(file.path(outdir, "plots", volcano_unpaired_file), plot = volcano_unpaired,
           width = VOLCANO_W, height = VOLCANO_H, dpi = PLOT_DPI)
    verify_output_file(file.path(outdir, "plots", volcano_unpaired_file), "Volcano unpaired")

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
      cat("[SKIP] CLDN1 not in VST matrix  cannot plot\n")
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
      cat("[SKIP] UCP1 not in VST matrix  cannot plot\n")
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
      cat("[SKIP] CLDN1 or UCP1 not found  cannot compute correlation\n")
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
    n_tested_unpaired = nrow(res_unpaired_df),
    n_up_bat_unpaired = sanity$n_up_bat_unpaired,
    n_up_wat_unpaired = sanity$n_up_wat_unpaired,
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
      "Count scale factor applied",
      "Genes raw (before any filter)",
      "Genes after count filter (CPM>=1)",
      "Second-pass filter applied",
      "DESeq2 runtime (min)",
      "BAT-enriched genes paired (paper: 463)",
      "BAT-enriched genes unpaired",
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
      "1 (raw counts) or >1 (scaled)",
      "~217000 (GEO file)",
      "15000-25000",
      "FALSE (ideally)",
      "< 10",
      paste0(EXPECTED_UP_BAT, " +/- ", DEG_COUNT_TOLERANCE * 100, "%"),
      "context-only sensitivity",
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
      ifelse(!is.null(sanity$count_scale_factor), sanity$count_scale_factor, 1),
      sanity$n_genes_raw,
      sanity$n_genes_filtered,
      ifelse(!is.null(sanity$second_pass_filter), sanity$second_pass_filter, "N/A"),
      ifelse(!is.null(sanity$deseq_runtime_min), sanity$deseq_runtime_min, "N/A"),
      sanity$n_up_bat,
      sanity$n_up_bat_unpaired,
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
cat("##  DONE  Part 1 complete                              ##\n")
cat("##  Run directory: ", basename(outdir), "\n", sep = "")
cat("##  UCP1  verdict: ", ifelse(is.null(sanity$ucp1_verdict), "NA",
                                   sanity$ucp1_verdict), "\n", sep = "")
cat("##  CLDN1 verdict: ", sanity$cldn1_verdict, "\n", sep = "")
cat("##########################################################\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
