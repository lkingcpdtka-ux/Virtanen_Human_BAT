## =========================================================
## SAVE_CORE: Deduplicate Files in Run Directory
## =========================================================
## Provides:
##   - dedupe(): Remove duplicate/versioned files, keeping only latest
## =========================================================

## Remove duplicate files from a run directory
## Looks for patterns like: file_v1.png, file_v2.png and keeps latest
## Also removes common temp file patterns
dedupe <- function(run_dir, dry_run = FALSE) {

  if (!dir.exists(run_dir)) {
    cat("[dedupe] Directory not found: ", run_dir, "\n", sep = "")
    return(invisible(NULL))
  }

  cat("[dedupe] Scanning: ", run_dir, "\n", sep = "")

  ## Common duplicate patterns to look for
  duplicate_patterns <- list(
    ## Version suffixes: file_v1.png, file_v2.png, etc.
    version = "_v[0-9]+\\.(png|pdf|csv|rds|txt)$",
    ## Numbered suffixes: file_1.png, file_2.png, etc.
    numbered = "_[0-9]+\\.(png|pdf|csv|rds|txt)$",
    ## Temp files
    temp = "^~\\$|^\\.~|^#.*#$|\\.(tmp|temp|bak)$"
  )

  ## Get all files
  all_files <- list.files(run_dir, recursive = TRUE, full.names = TRUE)

  files_to_remove <- character(0)

  ## Check for versioned files
  for (pattern_name in names(duplicate_patterns)) {
    pattern <- duplicate_patterns[[pattern_name]]
    matches <- grep(pattern, all_files, value = TRUE)

    if (length(matches) > 0 && pattern_name %in% c("version", "numbered")) {
      ## Group by base name (without version suffix)
      base_names <- gsub(pattern, "", matches)
      base_names <- gsub("\\.[^.]+$", "", base_names)  ## Remove extension too

      for (base in unique(base_names)) {
        group_files <- matches[startsWith(gsub("\\.[^.]+$", "", matches), base)]
        if (length(group_files) > 1) {
          ## Keep the most recently modified
          mtimes <- file.mtime(group_files)
          keep_idx <- which.max(mtimes)
          remove_files <- group_files[-keep_idx]
          files_to_remove <- c(files_to_remove, remove_files)
        }
      }
    } else if (pattern_name == "temp") {
      ## Remove all temp files
      files_to_remove <- c(files_to_remove, matches)
    }
  }

  ## Remove duplicates from list
  files_to_remove <- unique(files_to_remove)

  if (length(files_to_remove) == 0) {
    cat("[dedupe] No duplicate files found\n")
    return(invisible(NULL))
  }

  ## Calculate size
  total_size <- sum(file.size(files_to_remove), na.rm = TRUE)
  total_size_mb <- total_size / 1024 / 1024

  cat("[dedupe] Found ", length(files_to_remove), " duplicate/temp files (",
      round(total_size_mb, 2), " MB)\n", sep = "")

  if (dry_run) {
    cat("[dedupe] DRY RUN - files that would be removed:\n")
    for (f in files_to_remove) {
      cat("  - ", basename(f), "\n", sep = "")
    }
    return(invisible(files_to_remove))
  }

  ## Actually remove
  removed_count <- 0
  for (f in files_to_remove) {
    tryCatch({
      file.remove(f)
      removed_count <- removed_count + 1
    }, error = function(e) {
      cat("[dedupe] Failed to remove: ", basename(f), "\n", sep = "")
    })
  }

  cat("[dedupe] Removed ", removed_count, " files (", round(total_size_mb, 2), " MB)\n", sep = "")

  return(invisible(files_to_remove))
}

cat("[INFO] Loaded dedupe.R - use dedupe(run_dir) to remove duplicate files\n")
