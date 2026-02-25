## =========================================================
## SAVE_CORE: Run Initialization and File Registration
## =========================================================
## Provides:
##   - init_run(): Create or resume timestamped run directories
##   - save_run_file(): Register important files for quick access
##   - get_run_file(): Retrieve indexed files by tag
##   - finalize_run(): Write run completion metadata
## =========================================================

savecore_default_savepoint_dir <- function(savepoint_dir = NULL) {
  if (!is.null(savepoint_dir)) return(savepoint_dir)
  file.path(getwd(), "savepoints")
}

savecore_write_metadata <- function(outdir, run_metadata, run_tag) {
  log_dir <- file.path(outdir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  metadata_file_txt <- file.path(log_dir, paste0("run_metadata_", run_tag, ".txt"))
  metadata_lines <- c(
    "=== RUN METADATA ===",
    paste0("Script: ", run_metadata$script_name),
    paste0("Species: ", run_metadata$species),
    paste0("Data Type: ", run_metadata$data_type),
    paste0("Keywords: ", paste(run_metadata$keywords, collapse = ", ")),
    paste0("Notes: ", run_metadata$notes),
    paste0("Message: ", run_metadata$message),
    paste0("Run Tag: ", run_metadata$run_tag),
    paste0("Timestamp: ", format(run_metadata$timestamp, "%Y-%m-%d %H:%M:%S")),
    paste0("R Version: ", run_metadata$r_version),
    paste0("Platform: ", run_metadata$platform),
    "===================="
  )

  writeLines(metadata_lines, metadata_file_txt)
  saveRDS(run_metadata, file.path(log_dir, paste0("run_metadata_", run_tag, ".rds")))
}

savecore_ensure_layout <- function(outdir) {
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "cache"), recursive = TRUE, showWarnings = FALSE)
}

## Initialize or resume an analysis run
## mode:
##   - "new"           : create RUN_<timestamp>
##   - "resume_latest" : reuse most recently modified RUN_ directory
##   - "resume_tag"    : reuse RUN_<run_tag>
init_run <- function(script_name,
                     species = "mouse",
                     data_type = "bulkRNAseq",
                     keywords = character(0),
                     notes = "",
                     message = "",
                     mode = c("new", "resume_latest", "resume_tag"),
                     run_tag = NULL,
                     savepoint_dir = NULL) {

  mode <- match.arg(mode)
  savepoint_dir <- savecore_default_savepoint_dir(savepoint_dir)
  dir.create(savepoint_dir, recursive = TRUE, showWarnings = FALSE)

  if (mode == "new") {
    run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
    outdir <- file.path(savepoint_dir, paste0("RUN_", run_tag))
    savecore_ensure_layout(outdir)
    cat("[INIT] Created run directory: ", outdir, "\n", sep = "")
  } else {
    if (mode == "resume_tag") {
      if (is.null(run_tag) || !nzchar(run_tag)) {
        stop("[init_run] run_tag is required when mode='resume_tag'")
      }
      outdir <- file.path(savepoint_dir, paste0("RUN_", run_tag))
      if (!dir.exists(outdir)) {
        stop("[init_run] Run directory not found for tag: ", run_tag)
      }
    } else {
      run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
      run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
      if (length(run_dirs) == 0) {
        stop("[init_run] No previous runs found to resume")
      }
      outdir <- run_dirs[which.max(file.mtime(run_dirs))]
      run_tag <- sub("^RUN_", "", basename(outdir))
    }

    savecore_ensure_layout(outdir)
    cat("[INIT] Resuming run directory: ", outdir, "\n", sep = "")
  }

  if (is.null(run_tag) || !nzchar(run_tag)) {
    run_tag <- sub("^RUN_", "", basename(outdir))
  }

  run_metadata <- list(
    script_name = script_name,
    species     = species,
    data_type   = data_type,
    keywords    = keywords,
    notes       = notes,
    message     = message,
    run_tag     = run_tag,
    timestamp   = Sys.time(),
    mode        = mode,
    r_version   = R.version.string,
    platform    = .Platform$OS.type
  )

  savecore_write_metadata(outdir, run_metadata, run_tag)

  cat("[INIT] Run tag: ", run_tag, "\n", sep = "")

  return(list(
    outdir    = outdir,
    run_tag   = run_tag,
    cache_dir = file.path(outdir, "cache"),
    metadata  = run_metadata
  ))
}

## Register an important file for quick access
save_run_file <- function(file_path, tag = NULL, description = "", run_dir = NULL) {

  if (!file.exists(file_path)) {
    warning("[save_run_file] File does not exist: ", file_path)
    return(invisible(FALSE))
  }

  file_path <- normalizePath(file_path, winslash = "/", mustWork = TRUE)

  if (is.null(run_dir)) {
    path_parts <- strsplit(file_path, "/")[[1]]
    run_dir_idx <- grep("^RUN_", path_parts)
    if (length(run_dir_idx) == 0) {
      warning("[save_run_file] Cannot determine run directory from path; pass run_dir explicitly")
      return(invisible(FALSE))
    }
    run_dir <- do.call(file.path, as.list(path_parts[1:max(run_dir_idx)]))
  }

  index_file <- file.path(run_dir, "logs", "file_index.tsv")

  if (is.null(tag)) {
    tag <- tools::file_path_sans_ext(basename(file_path))
  }

  entry <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    tag = tag,
    path = file_path,
    description = description,
    stringsAsFactors = FALSE
  )

  if (!file.exists(index_file)) {
    write.table(entry, index_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE)
  } else {
    write.table(entry, index_file, sep = "\t", row.names = FALSE, col.names = FALSE,
                quote = TRUE, append = TRUE)
  }

  cat("[INDEX] Registered: ", tag, " -> ", basename(file_path), "\n", sep = "")
  return(invisible(TRUE))
}

## Get a registered file by tag
get_run_file <- function(run_dir, tag) {
  index_file <- file.path(run_dir, "logs", "file_index.tsv")

  if (!file.exists(index_file)) {
    warning("[get_run_file] No file index found in run directory")
    return(NULL)
  }

  index <- read.delim(index_file, stringsAsFactors = FALSE)
  match_row <- index[index$tag == tag, , drop = FALSE]

  if (nrow(match_row) == 0) {
    warning("[get_run_file] Tag not found: ", tag)
    return(NULL)
  }

  match_row <- match_row[order(match_row$timestamp, decreasing = TRUE), , drop = FALSE]
  return(match_row$path[1])
}

## Write run completion marker
finalize_run <- function(outdir, status = c("success", "error"), summary_lines = character(0)) {
  status <- match.arg(status)
  status_file <- file.path(outdir, "logs", "run_status.txt")

  lines <- c(
    paste0("status: ", status),
    paste0("completed_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    summary_lines
  )

  writeLines(lines, status_file)
  cat("[FINALIZE] Wrote run status: ", status, "\n", sep = "")
  invisible(status_file)
}

cat("[INFO] Loaded save.R - use init_run() to create or resume a run directory\n")
