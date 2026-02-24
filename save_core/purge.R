## =========================================================
## SAVE_CORE: Purge Old Run Directories
## =========================================================
## Provides:
##   - purge(): Remove old run directories based on retention policy
##   - purge_cache(): Clear cache from specific runs
## =========================================================

## Purge old run directories
## keep_recent: Number of recent runs to keep (default 5)
## age_days: Remove runs older than this many days (default NULL = no age limit)
## dry_run: If TRUE, only show what would be deleted without deleting
purge <- function(keep_recent = 5,
                  age_days = NULL,
                  dry_run = TRUE,
                  savepoint_dir = NULL) {

  if (is.null(savepoint_dir)) {
    savepoint_dir <- file.path(getwd(), "savepoints")
  }

  if (!dir.exists(savepoint_dir)) {
    cat("[purge] No savepoints directory found\n")
    return(invisible(NULL))
  }

  ## Find all run directories
  run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]

  if (length(run_dirs) == 0) {
    cat("[purge] No run directories found\n")
    return(invisible(NULL))
  }

  ## Get modification times
  mtimes <- file.mtime(run_dirs)
  run_order <- order(mtimes, decreasing = TRUE)
  run_dirs <- run_dirs[run_order]
  mtimes <- mtimes[run_order]

  ## Identify runs to keep
  keep_mask <- rep(FALSE, length(run_dirs))

  ## Keep the N most recent
  if (!is.null(keep_recent) && keep_recent > 0) {
    keep_mask[1:min(keep_recent, length(run_dirs))] <- TRUE
  }

  ## Apply age filter
  if (!is.null(age_days)) {
    age_cutoff <- Sys.time() - as.difftime(age_days, units = "days")
    age_keep <- mtimes >= age_cutoff
    keep_mask <- keep_mask | age_keep
  }

  ## Identify runs to purge
  purge_dirs <- run_dirs[!keep_mask]
  keep_dirs  <- run_dirs[keep_mask]

  if (length(purge_dirs) == 0) {
    cat("[purge] No runs to purge (keeping ", length(keep_dirs), " runs)\n", sep = "")
    return(invisible(NULL))
  }

  ## Calculate disk space
  purge_sizes <- sapply(purge_dirs, function(d) {
    files <- list.files(d, recursive = TRUE, full.names = TRUE)
    sum(file.size(files), na.rm = TRUE)
  })
  total_size_mb <- sum(purge_sizes) / 1024 / 1024

  cat("\n=== PURGE SUMMARY ===\n")
  cat("Keeping: ", length(keep_dirs), " runs\n", sep = "")
  cat("Purging: ", length(purge_dirs), " runs (", round(total_size_mb, 1), " MB)\n", sep = "")
  cat("\nRuns to purge:\n")
  for (i in seq_along(purge_dirs)) {
    cat("  - ", basename(purge_dirs[i]),
        " (", round(purge_sizes[i] / 1024 / 1024, 1), " MB)\n", sep = "")
  }

  if (dry_run) {
    cat("\n[DRY RUN] No files deleted. Set dry_run=FALSE to actually purge.\n")
    cat("=====================\n\n")
    return(invisible(data.frame(
      run_dir = purge_dirs,
      size_mb = round(purge_sizes / 1024 / 1024, 1),
      stringsAsFactors = FALSE
    )))
  }

  ## Actually delete
  cat("\nDeleting...\n")
  for (d in purge_dirs) {
    tryCatch({
      unlink(d, recursive = TRUE)
      cat("  Deleted: ", basename(d), "\n", sep = "")
    }, error = function(e) {
      cat("  FAILED: ", basename(d), " - ", conditionMessage(e), "\n", sep = "")
    })
  }

  cat("\n[purge] Freed ", round(total_size_mb, 1), " MB\n", sep = "")
  cat("=====================\n\n")

  return(invisible(NULL))
}

## Purge cache directories from all runs (keeps tables/plots/logs)
## Useful to free up disk space while keeping final outputs
purge_cache <- function(keep_recent = 3,
                        dry_run = TRUE,
                        savepoint_dir = NULL) {

  if (is.null(savepoint_dir)) {
    savepoint_dir <- file.path(getwd(), "savepoints")
  }

  ## Find all run directories
  run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]

  if (length(run_dirs) == 0) {
    cat("[purge_cache] No run directories found\n")
    return(invisible(NULL))
  }

  ## Sort by modification time
  mtimes <- file.mtime(run_dirs)
  run_order <- order(mtimes, decreasing = TRUE)
  run_dirs <- run_dirs[run_order]

  ## Skip the most recent N runs
  if (length(run_dirs) <= keep_recent) {
    cat("[purge_cache] All runs are recent (keeping cache for ", length(run_dirs), " runs)\n", sep = "")
    return(invisible(NULL))
  }

  old_runs <- run_dirs[(keep_recent + 1):length(run_dirs)]

  ## Find cache directories
  cache_dirs <- file.path(old_runs, "cache")
  cache_dirs <- cache_dirs[dir.exists(cache_dirs)]

  if (length(cache_dirs) == 0) {
    cat("[purge_cache] No cache directories to purge\n")
    return(invisible(NULL))
  }

  ## Calculate sizes
  cache_sizes <- sapply(cache_dirs, function(d) {
    files <- list.files(d, recursive = TRUE, full.names = TRUE)
    sum(file.size(files), na.rm = TRUE)
  })
  total_size_mb <- sum(cache_sizes) / 1024 / 1024

  cat("\n=== CACHE PURGE SUMMARY ===\n")
  cat("Keeping cache for: ", keep_recent, " most recent runs\n", sep = "")
  cat("Purging cache from: ", length(cache_dirs), " runs (", round(total_size_mb, 1), " MB)\n", sep = "")

  if (dry_run) {
    cat("\n[DRY RUN] No files deleted. Set dry_run=FALSE to actually purge.\n")
    cat("===========================\n\n")
    return(invisible(NULL))
  }

  ## Delete cache directories
  for (d in cache_dirs) {
    tryCatch({
      unlink(d, recursive = TRUE)
      cat("  Deleted: ", d, "\n", sep = "")
    }, error = function(e) {
      cat("  FAILED: ", d, " - ", conditionMessage(e), "\n", sep = "")
    })
  }

  cat("\n[purge_cache] Freed ", round(total_size_mb, 1), " MB\n", sep = "")
  cat("===========================\n\n")

  return(invisible(NULL))
}

cat("[INFO] Loaded purge.R - use purge() to remove old runs (dry_run=TRUE by default)\n")
