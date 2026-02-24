## =========================================================
## CACHE UTILITIES FOR KAT8 PIPELINE
## =========================================================
## Provides caching of expensive computations (DESeq2, GSEA, etc.)
## Integrates with save_core and savepoints directory structure
##
## Usage:
##   source("save_core/cache_utils.R")
##
##   # Check if cached result exists, otherwise compute and cache
##   dds <- cache_load_or_compute(
##     cache_key = "dds_fitted",
##     compute_fn = function() DESeq(dds),
##     cache_dir = file.path(outdir, "cache"),
##     force_recompute = FALSE
##   )
## =========================================================

## Initialize cache directory within a run
init_cache <- function(outdir) {
  cache_dir <- file.path(outdir, "cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cat("[CACHE] Initialized cache directory: ", cache_dir, "\n", sep = "")
  }
  return(cache_dir)
}

## Generate cache file path
get_cache_path <- function(cache_key, cache_dir) {
  file.path(cache_dir, paste0(cache_key, ".rds"))
}

## Check if cache exists and is valid
cache_exists <- function(cache_key, cache_dir) {
  cache_path <- get_cache_path(cache_key, cache_dir)
  file.exists(cache_path)
}

## Load from cache
cache_load <- function(cache_key, cache_dir, verbose = TRUE) {
  cache_path <- get_cache_path(cache_key, cache_dir)
  if (!file.exists(cache_path)) {
    if (verbose) cat("[CACHE] Not found: ", cache_key, "\n", sep = "")
    return(NULL)
  }

  if (verbose) cat("[CACHE] Loading: ", cache_key, " ... ", sep = "")
  result <- readRDS(cache_path)
  if (verbose) cat("done\n")
  return(result)
}

## Save to cache
cache_save <- function(obj, cache_key, cache_dir, verbose = TRUE) {
  cache_path <- get_cache_path(cache_key, cache_dir)
  if (verbose) cat("[CACHE] Saving: ", cache_key, " ... ", sep = "")
  saveRDS(obj, cache_path)
  if (verbose) cat("done (", round(file.size(cache_path) / 1024 / 1024, 1), " MB)\n", sep = "")
  return(invisible(cache_path))
}

## Main function: load from cache or compute and save
## This is the primary function you'll use
cache_load_or_compute <- function(cache_key,
                                   compute_fn,
                                   cache_dir,
                                   force_recompute = FALSE,
                                   verbose = TRUE) {

  cache_path <- get_cache_path(cache_key, cache_dir)

  ## Check if we should use cache
  if (!force_recompute && file.exists(cache_path)) {
    if (verbose) cat("[CACHE] HIT - Loading: ", cache_key, "\n", sep = "")
    result <- readRDS(cache_path)
    if (verbose) cat("[CACHE] Loaded successfully\n")
    return(result)
  }

  ## Compute the result
  if (verbose) {
    if (force_recompute) {
      cat("[CACHE] FORCED RECOMPUTE: ", cache_key, "\n", sep = "")
    } else {
      cat("[CACHE] MISS - Computing: ", cache_key, "\n", sep = "")
    }
  }

  start_time <- Sys.time()
  result <- compute_fn()
  elapsed <- difftime(Sys.time(), start_time, units = "secs")

  if (verbose) cat("[CACHE] Computed in ", round(elapsed, 1), " seconds\n", sep = "")

  ## Save to cache
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  saveRDS(result, cache_path)
  if (verbose) {
    cat("[CACHE] Saved: ", basename(cache_path),
        " (", round(file.size(cache_path) / 1024 / 1024, 1), " MB)\n", sep = "")
  }

  return(result)
}

## List all cached objects in a directory
cache_list <- function(cache_dir) {
  if (!dir.exists(cache_dir)) {
    cat("[CACHE] No cache directory found\n")
    return(character(0))
  }

  cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(cache_files) == 0) {
    cat("[CACHE] No cached objects found\n")
    return(character(0))
  }

  ## Get info about each cached file
  cache_info <- data.frame(
    key = gsub("\\.rds$", "", basename(cache_files)),
    size_mb = round(file.size(cache_files) / 1024 / 1024, 1),
    modified = file.mtime(cache_files),
    stringsAsFactors = FALSE
  )

  cache_info <- cache_info[order(cache_info$modified, decreasing = TRUE), ]

  cat("[CACHE] Found ", nrow(cache_info), " cached objects:\n", sep = "")
  for (i in seq_len(nrow(cache_info))) {
    cat("  - ", cache_info$key[i], " (", cache_info$size_mb[i], " MB, ",
        format(cache_info$modified[i], "%Y-%m-%d %H:%M"), ")\n", sep = "")
  }

  return(cache_info)
}

## Clear specific cache or all caches
cache_clear <- function(cache_dir, cache_key = NULL) {
  if (!dir.exists(cache_dir)) {
    cat("[CACHE] No cache directory found\n")
    return(invisible(FALSE))
  }

  if (is.null(cache_key)) {
    ## Clear all
    cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(cache_files) > 0) {
      file.remove(cache_files)
      cat("[CACHE] Cleared ", length(cache_files), " cached objects\n", sep = "")
    }
  } else {
    ## Clear specific key
    cache_path <- get_cache_path(cache_key, cache_dir)
    if (file.exists(cache_path)) {
      file.remove(cache_path)
      cat("[CACHE] Cleared: ", cache_key, "\n", sep = "")
    } else {
      cat("[CACHE] Not found: ", cache_key, "\n", sep = "")
    }
  }

  return(invisible(TRUE))
}

## =========================================================
## CONVENIENCE WRAPPERS FOR COMMON EXPENSIVE OPERATIONS
## =========================================================

## Cache DESeq2 object (Part 1)
cache_deseq <- function(dds_unfitted, cache_dir, contrast_name = "main", force = FALSE) {
  cache_key <- paste0("dds_", contrast_name)

  cache_load_or_compute(
    cache_key = cache_key,
    compute_fn = function() {
      cat("[INFO] Running DESeq2 (this may take a few minutes)...\n")
      DESeq(dds_unfitted)
    },
    cache_dir = cache_dir,
    force_recompute = force
  )
}

## Cache VST-transformed counts (Part 1)
cache_vst <- function(dds, cache_dir, force = FALSE) {
  cache_load_or_compute(
    cache_key = "vst_counts",
    compute_fn = function() {
      cat("[INFO] Computing VST transformation...\n")
      vst(dds, blind = FALSE)
    },
    cache_dir = cache_dir,
    force_recompute = force
  )
}

## Cache gseGO result (Part 3)
cache_gsego <- function(genelist, orgdb, ont, contrast_name, cache_dir,
                        min_gs = 5, max_gs = 500, force = FALSE) {
  cache_key <- paste0("gsego_", ont, "_", contrast_name)

  cache_load_or_compute(
    cache_key = cache_key,
    compute_fn = function() {
      cat("[INFO] Running gseGO (", ont, ") - this may take 1-3 minutes...\n", sep = "")
      gseGO(
        geneList     = genelist,
        OrgDb        = orgdb,
        ont          = ont,
        keyType      = "ENTREZID",
        minGSSize    = min_gs,
        maxGSSize    = max_gs,
        pvalueCutoff = 1.0,
        pAdjustMethod = "BH",
        verbose      = TRUE,
        seed         = TRUE,
        nPermSimple  = 1000
      )
    },
    cache_dir = cache_dir,
    force_recompute = force
  )
}

## Cache gseKEGG result (Part 3)
cache_gsekegg <- function(genelist, organism, contrast_name, cache_dir,
                          min_gs = 5, max_gs = 500, force = FALSE) {
  cache_key <- paste0("gsekegg_", contrast_name)

  cache_load_or_compute(
    cache_key = cache_key,
    compute_fn = function() {
      cat("[INFO] Running gseKEGG - this may take 30 seconds to 1 minute...\n")
      gseKEGG(
        geneList     = genelist,
        organism     = organism,
        minGSSize    = min_gs,
        maxGSSize    = max_gs,
        pvalueCutoff = 0.1,
        pAdjustMethod = "BH",
        verbose      = TRUE,
        seed         = TRUE,
        nPermSimple  = 1000
      )
    },
    cache_dir = cache_dir,
    force_recompute = force
  )
}

## Cache GO simplification result (Parts 2/3)
cache_simplify_go <- function(enrich_result, cutoff, contrast_name, direction, cache_dir, force = FALSE) {
  cache_key <- paste0("simplify_go_", contrast_name, "_", direction, "_cut", cutoff)

  cache_load_or_compute(
    cache_key = cache_key,
    compute_fn = function() {
      cat("[INFO] Simplifying GO terms (cutoff=", cutoff, ") - this may take 1-3 minutes...\n", sep = "")
      clusterProfiler::simplify(
        enrich_result,
        cutoff = cutoff,
        by = "p.adjust",
        select_fun = min
      )
    },
    cache_dir = cache_dir,
    force_recompute = force
  )
}

## =========================================================
## GLOBAL CACHE CONTROL
## =========================================================

## Set this to TRUE to force all caches to recompute
FORCE_RECOMPUTE_ALL <- FALSE

## Set this to skip caching entirely (for debugging)
DISABLE_CACHE <- FALSE

## Helper to check if caching is enabled
use_cache <- function() {
  !exists("DISABLE_CACHE") || !DISABLE_CACHE
}

## Helper to check if force recompute is set
force_recompute <- function() {
  exists("FORCE_RECOMPUTE_ALL") && FORCE_RECOMPUTE_ALL
}

cat("[INFO] Loaded cache_utils.R - use cache_load_or_compute() for expensive operations\n")
cat("[INFO] Set FORCE_RECOMPUTE_ALL <- TRUE to regenerate all cached results\n")
