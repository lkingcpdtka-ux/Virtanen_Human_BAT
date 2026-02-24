## =========================================================
## SAVE_CORE: Run Initialization and File Registration
## =========================================================
## Provides:
##   - init_run(): Create timestamped run directory with metadata
##   - save_run_file(): Register important files for quick access
## =========================================================

## Initialize a new analysis run
## Returns: list(outdir, run_tag, cache_dir)
init_run <- function(script_name,
                     species = "mouse",
                     data_type = "bulkRNAseq",
                     keywords = character(0),
                     notes = "",
                     message = "") {

  ## Generate timestamp-based run tag

  run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outdir <- file.path(getwd(), "savepoints", paste0("RUN_", run_tag))

  ## Create directory structure
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "cache"), recursive = TRUE, showWarnings = FALSE)

  ## Save run metadata
  run_metadata <- list(
    script_name = script_name,
    species     = species,
    data_type   = data_type,
    keywords    = keywords,
    notes       = notes,
    message     = message,
    run_tag     = run_tag,
    timestamp   = Sys.time(),
    r_version   = R.version.string,
    platform    = .Platform$OS.type
  )

  ## Save metadata to JSON-like text file
  metadata_file <- file.path(outdir, "logs", paste0("run_metadata_", run_tag, ".txt"))
  metadata_lines <- c(
    "=== RUN METADATA ===",
    paste0("Script: ", script_name),
    paste0("Species: ", species),
    paste0("Data Type: ", data_type),
    paste0("Keywords: ", paste(keywords, collapse = ", ")),
    paste0("Notes: ", notes),
    paste0("Message: ", message),
    paste0("Run Tag: ", run_tag),
    paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("R Version: ", R.version.string),
    paste0("Platform: ", .Platform$OS.type),
    "===================="
  )
  writeLines(metadata_lines, metadata_file)

  ## Also save as RDS for programmatic access
  saveRDS(run_metadata, file.path(outdir, "logs", paste0("run_metadata_", run_tag, ".rds")))

  cat("[INIT] Created run directory: ", outdir, "\n", sep = "")
  cat("[INIT] Run tag: ", run_tag, "\n", sep = "")

  return(list(
    outdir    = outdir,
    run_tag   = run_tag,
    cache_dir = file.path(outdir, "cache"),
    metadata  = run_metadata
  ))
}

## Register an important file for quick access
## Creates an index of "hero" files that can be quickly retrieved
save_run_file <- function(file_path, tag = NULL, description = "") {

  if (!file.exists(file_path)) {
    warning("[save_run_file] File does not exist: ", file_path)
    return(invisible(FALSE))
  }

  ## Determine the run directory from file path
  ## Expects path like: savepoints/RUN_YYYYMMDD_HHMMSS/subdir/file.ext
  path_parts <- strsplit(file_path, .Platform$file.sep)[[1]]
  run_dir_idx <- grep("^RUN_", path_parts)

  if (length(run_dir_idx) == 0) {
    warning("[save_run_file] Cannot determine run directory from path")
    return(invisible(FALSE))
  }

  run_dir <- do.call(file.path, as.list(path_parts[1:run_dir_idx]))
  index_file <- file.path(run_dir, "logs", "file_index.txt")

  ## Use filename as tag if not provided
  if (is.null(tag)) {
    tag <- tools::file_path_sans_ext(basename(file_path))
  }

  ## Append to index
  entry <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\t",
    tag, "\t",
    file_path, "\t",
    description
  )

  ## Create header if file doesn't exist
  if (!file.exists(index_file)) {
    writeLines("timestamp\ttag\tpath\tdescription", index_file)
  }

  cat(entry, "\n", file = index_file, append = TRUE)
  cat("[INDEX] Registered: ", tag, " -> ", basename(file_path), "\n", sep = "")

  return(invisible(TRUE))
}

## Get a registered file by tag
get_run_file <- function(run_dir, tag) {
  index_file <- file.path(run_dir, "logs", "file_index.txt")

  if (!file.exists(index_file)) {
    warning("[get_run_file] No file index found in run directory")
    return(NULL)
  }

  index <- read.delim(index_file, stringsAsFactors = FALSE)
  match_row <- index[index$tag == tag, ]

  if (nrow(match_row) == 0) {
    warning("[get_run_file] Tag not found: ", tag)
    return(NULL)
  }

  ## Return most recent if multiple matches
  return(match_row$path[nrow(match_row)])
}

cat("[INFO] Loaded save.R - use init_run() to create a new run directory\n")
