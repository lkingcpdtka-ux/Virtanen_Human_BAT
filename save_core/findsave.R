## =========================================================
## SAVE_CORE: Find Previous Run Directories
## =========================================================
## Provides:
##   - findsave(): Find previous run directories by various criteria
##   - find_latest_run(): Get most recent run directory
##   - find_run_by_tag(): Find run by specific tag
## =========================================================

## Find run directories matching criteria
## Returns: data.frame with run_dir, run_tag, mtime, and metadata
findsave <- function(keyword = NULL,
                     script_name = NULL,
                     species = NULL,
                     recent = NULL,
                     savepoint_dir = NULL) {

  ## Default savepoint directory
  if (is.null(savepoint_dir)) {
    savepoint_dir <- file.path(getwd(), "savepoints")
  }

  if (!dir.exists(savepoint_dir)) {
    cat("[findsave] No savepoints directory found\n")
    return(NULL)
  }

  ## Find all run directories
  run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]

  if (length(run_dirs) == 0) {
    cat("[findsave] No run directories found\n")
    return(NULL)
  }

  ## Build info data frame
  run_info <- data.frame(
    run_dir  = run_dirs,
    run_tag  = gsub("^RUN_", "", basename(run_dirs)),
    mtime    = file.mtime(run_dirs),
    stringsAsFactors = FALSE
  )

  ## Load metadata for filtering (if available)
  run_info$script_name <- NA_character_
  run_info$species     <- NA_character_
  run_info$keywords    <- NA_character_

  for (i in seq_len(nrow(run_info))) {
    meta_files <- list.files(
      file.path(run_info$run_dir[i], "logs"),
      pattern = "^run_metadata_.*\\.rds$",
      full.names = TRUE
    )
    if (length(meta_files) > 0) {
      tryCatch({
        meta <- readRDS(meta_files[1])
        run_info$script_name[i] <- meta$script_name
        run_info$species[i]     <- meta$species
        run_info$keywords[i]    <- paste(meta$keywords, collapse = ",")
      }, error = function(e) NULL)
    }
  }

  ## Apply filters
  if (!is.null(keyword)) {
    match_idx <- grep(keyword, run_info$keywords, ignore.case = TRUE)
    if (length(match_idx) > 0) {
      run_info <- run_info[match_idx, , drop = FALSE]
    } else {
      cat("[findsave] No runs matching keyword: ", keyword, "\n", sep = "")
      return(NULL)
    }
  }

  if (!is.null(script_name)) {
    match_idx <- grep(script_name, run_info$script_name, ignore.case = TRUE)
    if (length(match_idx) > 0) {
      run_info <- run_info[match_idx, , drop = FALSE]
    } else {
      cat("[findsave] No runs matching script: ", script_name, "\n", sep = "")
      return(NULL)
    }
  }

  if (!is.null(species)) {
    match_idx <- which(tolower(run_info$species) == tolower(species))
    if (length(match_idx) > 0) {
      run_info <- run_info[match_idx, , drop = FALSE]
    } else {
      cat("[findsave] No runs matching species: ", species, "\n", sep = "")
      return(NULL)
    }
  }

  ## Sort by modification time (most recent first)
  run_info <- run_info[order(run_info$mtime, decreasing = TRUE), ]

  ## Limit to recent N
  if (!is.null(recent) && recent > 0 && nrow(run_info) > recent) {
    run_info <- run_info[1:recent, , drop = FALSE]
  }

  return(run_info)
}

## Get the most recent run directory
## This is the most common use case
find_latest_run <- function(savepoint_dir = NULL, script_name = NULL) {
  runs <- findsave(savepoint_dir = savepoint_dir, script_name = script_name, recent = 1)

  if (is.null(runs) || nrow(runs) == 0) {
    return(NULL)
  }

  run_dir <- runs$run_dir[1]
  run_tag <- runs$run_tag[1]

  cat("[findsave] Found latest run: ", basename(run_dir), "\n", sep = "")
  cat("[findsave] Run tag: ", run_tag, "\n", sep = "")

  return(list(
    outdir    = run_dir,
    run_tag   = run_tag,
    cache_dir = file.path(run_dir, "cache")
  ))
}

## Find run by specific tag
find_run_by_tag <- function(run_tag, savepoint_dir = NULL) {
  if (is.null(savepoint_dir)) {
    savepoint_dir <- file.path(getwd(), "savepoints")
  }

  run_dir <- file.path(savepoint_dir, paste0("RUN_", run_tag))

  if (!dir.exists(run_dir)) {
    cat("[findsave] Run not found: ", run_tag, "\n", sep = "")
    return(NULL)
  }

  cat("[findsave] Found run: ", run_tag, "\n", sep = "")

  return(list(
    outdir    = run_dir,
    run_tag   = run_tag,
    cache_dir = file.path(run_dir, "cache")
  ))
}

## List all runs with summary
list_runs <- function(savepoint_dir = NULL, n = 10) {
  runs <- findsave(savepoint_dir = savepoint_dir, recent = n)

  if (is.null(runs) || nrow(runs) == 0) {
    cat("[findsave] No runs found\n")
    return(invisible(NULL))
  }

  cat("\n=== Recent Runs ===\n")
  for (i in seq_len(nrow(runs))) {
    cat(sprintf("%2d. %s (%s)\n", i, runs$run_tag[i],
                format(runs$mtime[i], "%Y-%m-%d %H:%M")))
    if (!is.na(runs$script_name[i])) {
      cat(sprintf("    Script: %s\n", runs$script_name[i]))
    }
  }
  cat("===================\n\n")

  return(invisible(runs))
}

cat("[INFO] Loaded findsave.R - use find_latest_run() or findsave() to locate previous runs\n")
