# Fast species lookup with persisted reference-sequence caching.
# Annotation symbol matches are NOT assumed to be verified orthologs.
# This file is sourced after R/ensembl.R.

fast_cache_path <- function(species, symbol, window_bp, cache_dir = "cache/ensembl-v2") {
  key <- sprintf("%s_%s_%s.rds", species,
                 gsub("[^A-Za-z0-9_.-]", "_", toupper(symbol)), as.integer(window_bp))
  file.path(cache_dir, key)
}

fast_cache_read <- function(species, symbol, window_bp, cache_dir = "cache/ensembl-v2",
                            ttl_days = 14) {
  path <- fast_cache_path(species, symbol, window_bp, cache_dir)
  if (!file.exists(path)) return(NULL)
  info <- file.info(path)
  if (is.na(info$mtime) || difftime(Sys.time(), info$mtime, units = "days") > ttl_days) return(NULL)
  rec <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(rec) || !identical(rec$species, species) ||
      !is.character(rec$sequence) || length(rec$sequence) != 1L ||
      nchar(rec$sequence) != window_bp ||
      (length(rec$source) != 1L || !rec$source %in% c("ensembl", "ucsc"))) return(NULL)
  if (!is.character(rec$assembly) || length(rec$assembly) != 1L ||
      !nzchar(rec$assembly) || is.null(rec$gene_id)) return(NULL)
  # Reject invalid records cached by older versions where the assembly was
  # parsed from a response ID and could become GRCh38 for Mouse or Rat.
  if (!identical(species, "homo_sapiens") &&
      rec$assembly %in% c("GRCh37", "GRCh38")) return(NULL)
  rec
}

fast_cache_write <- function(record, symbol, window_bp, cache_dir = "cache/ensembl-v2") {
  if (!is.list(record) || (length(record$source) != 1L || !record$source %in% c("ensembl", "ucsc"))) return(invisible(FALSE))
  tryCatch({
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    path <- fast_cache_path(record$species, symbol, window_bp, cache_dir)
    temp <- tempfile(pattern = "cache-", tmpdir = cache_dir)
    on.exit(unlink(temp), add = TRUE)
    saveRDS(record, temp)
    if (file.exists(path)) unlink(path)
    file.rename(temp, path)
  }, error = function(e) FALSE)
}

# Prefer an organism's conventional gene-symbol casing before querying Ensembl.
# For unknown genes these are candidate identifiers, not orthology evidence.
preferred_symbol <- function(symbol, species) {
  s <- toupper(trimws(symbol))
  known <- LOCAL_SYMBOL_FALLBACKS[[s]]
  if (!is.null(known) && species %in% names(known)) return(unname(known[[species]]))
  if (species == "danio_rerio") return(tolower(symbol))
  if (species %in% c("mus_musculus", "rattus_norvegicus")) {
    return(paste0(toupper(substr(symbol, 1L, 1L)), tolower(substr(symbol, 2L, nchar(symbol)))))
  }
  toupper(symbol)
}

fast_gene_valid <- function(gene, species) {
  if (is.null(gene) || !identical(gene$object_type, "Gene") || is.null(gene$id)) return(FALSE)
  if (!is.null(gene$species) && !identical(as.character(gene$species), species)) return(FALSE)
  TRUE
}

fast_endpoints <- function(species, symbol) {
  vapply(seq_along(species), function(i) {
    sprintf("/lookup/symbol/%s/%s?expand=0", species[[i]],
            utils::URLencode(symbol[[i]], reserved = TRUE))
  }, character(1))
}

fetch_comparative_regions_fast <- function(symbol, selected, window_bp = 10000L,
                                          cache = NULL, api_batch = ensembl_json_batch,
                                          cache_dir = "cache/ensembl-v2", disk_cache = NULL, fallback_fetch = NULL) {
  if (!valid_gene_symbol(symbol)) stop("Invalid human gene symbol")
  if (!length(selected) || anyNA(selected) || any(!selected %in% unname(SPECIES))) {
    stop("Unsupported or missing species")
  }
  if (length(window_bp) != 1L || !is.finite(window_bp) ||
      window_bp < 1000L || window_bp > 1000000L) stop("Invalid window size")
  window_bp <- as.integer(window_bp)
  symbol <- trimws(symbol)
  selected <- unique(as.character(selected))
  # Mocked requests never read/write real-data caches unless explicitly enabled.
  if (is.null(disk_cache)) disk_cache <- identical(api_batch, ensembl_json_batch)
  # Offline/mock callers must never hit the network implicitly.
  if (is.null(fallback_fetch) && identical(api_batch, ensembl_json_batch) &&
      exists("fetch_ucsc_fallback", mode = "function")) {
    fallback_fetch <- fetch_ucsc_fallback
  }

  records <- stats::setNames(vector("list", length(selected)), selected)
  problems <- character()
  for (s in selected) {
    key <- paste(s, toupper(symbol), window_bp, sep = "|")
    if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
      records[[s]] <- get(key, envir = cache, inherits = FALSE)
    } else if (isTRUE(disk_cache)) {
      records[[s]] <- fast_cache_read(s, symbol, window_bp, cache_dir)
      if (!is.null(records[[s]]) && !is.null(cache)) assign(key, records[[s]], envir = cache)
    }
  }
  missing <- selected[vapply(records[selected], is.null, logical(1))]
  if (!length(missing)) return(list(records = unname(records), errors = problems))

  # One parallel round of organism-specific symbol lookup, including the human
  # anchor needed only if direct matches are incomplete.
  candidates <- unique(c("homo_sapiens", missing))
  preferred <- vapply(candidates, function(s) preferred_symbol(symbol, s), character(1))
  first <- api_batch(fast_endpoints(candidates, preferred))
  metadata <- stats::setNames(vector("list", length(candidates)), candidates)
  methods <- stats::setNames(rep("symbol_match_unverified", length(candidates)), candidates)
  for (i in seq_along(candidates)) {
    if (fast_gene_valid(first$results[[i]], candidates[[i]])) {
      metadata[[candidates[[i]]]] <- first$results[[i]]
    }
  }
  if (!is.null(metadata[["homo_sapiens"]])) methods[["homo_sapiens"]] <- "reference"

  # One additional parallel round for alternative spelling, only after 404.
  unresolved <- missing[vapply(metadata[missing], is.null, logical(1))]
  retry_species <- character()
  retry_symbols <- character()
  for (s in unresolved) {
    ix <- match(s, candidates)
    if (!identical(first$errors[[ix]], "Not found (HTTP 404)")) next
    alternatives <- setdiff(symbol_candidates(symbol, s), preferred[[s]])
    if (!length(alternatives)) next
    retry_species <- c(retry_species, s)
    retry_symbols <- c(retry_symbols, alternatives[[1L]])
  }
  if (length(retry_species)) {
    second <- api_batch(fast_endpoints(retry_species, retry_symbols))
    for (i in seq_along(retry_species)) {
      s <- retry_species[[i]]
      if (fast_gene_valid(second$results[[i]], s)) metadata[[s]] <- second$results[[i]]
    }
  }

  # Only unresolved species enter the slower, orthology-validated route.
  unresolved <- missing[vapply(metadata[missing], is.null, logical(1))]
  if (length(unresolved) && !is.null(metadata[["homo_sapiens"]])) {
    ref_id <- as.character(metadata[["homo_sapiens"]]$id)
    homology <- api_batch(sprintf("/homology/id/homo_sapiens/%s?type=orthologues&sequence=none&aligned=0",
                                  utils::URLencode(ref_id, reserved = TRUE)))
    if (!is.null(homology$results[[1L]])) {
      mapped <- resolve_orthologs(homology$results[[1L]], unresolved)
      problems <- c(problems, mapped$warnings)
      available <- unresolved[!is.na(mapped$ids[unresolved])]
      if (length(available)) {
        lookups <- api_batch(sprintf("/lookup/id/%s?expand=0", unname(mapped$ids[available])))
        for (i in seq_along(available)) {
          s <- available[[i]]
          if (fast_gene_valid(lookups$results[[i]], s)) {
            metadata[[s]] <- lookups$results[[i]]
            methods[[s]] <- as.character(mapped$types[[s]])
          }
        }
      }
    } else {
      problems <- c(problems, sprintf("Ortholog query: %s", homology$errors[[1L]]))
    }
  }

  ready <- missing[!vapply(metadata[missing], is.null, logical(1))]
  failed <- setdiff(missing, ready)
  for (s in failed) {
    index <- match(s, candidates)
    detail <- if (!is.na(index)) first$errors[[index]] else ""
    if (is.null(detail) || !length(detail) || !nzchar(detail)) detail <- "no gene or ortholog match"
    problems <- c(problems, sprintf("%s: gene lookup unavailable (%s)", species_label(s), detail))
  }

  intervals <- list()
  for (s in ready) {
    intervals[[s]] <- tryCatch(gene_interval(metadata[[s]], window_bp), error = function(e) {
      problems <<- c(problems, sprintf("%s: %s", species_label(s), conditionMessage(e)))
      NULL
    })
  }
  ready <- ready[!vapply(intervals[ready], is.null, logical(1))]
  if (length(ready)) {
    endpoints <- vapply(ready, function(s) {
      sprintf("/sequence/region/%s/%s", s,
              utils::URLencode(intervals[[s]]$region, reserved = FALSE))
    }, character(1))
    sequence_batch <- api_batch(endpoints)
    for (i in seq_along(ready)) {
      s <- ready[[i]]
      if (is.null(sequence_batch$results[[i]])) {
        problems <- c(problems, sprintf("%s: sequence request: %s", species_label(s), sequence_batch$errors[[i]]))
        next
      }
      rec <- tryCatch(record_from_sequence(s, metadata[[s]], intervals[[s]],
                                            sequence_batch$results[[i]], window_bp,
                                            orthology_type = methods[[s]]), error = function(e) {
        problems <<- c(problems, sprintf("%s: %s", species_label(s), conditionMessage(e)))
        NULL
      })
      if (is.null(rec)) next
      records[[s]] <- rec
      key <- paste(s, toupper(symbol), window_bp, sep = "|")
      if (!is.null(cache)) assign(key, rec, envir = cache)
      if (isTRUE(disk_cache)) fast_cache_write(rec, symbol, window_bp, cache_dir)
    }
  }
  # Independent public provider only for species that Ensembl could not load.
  # Never claim an annotation match proves orthology; record assembly/source.
  no_seq <- selected[vapply(records[selected], is.null, logical(1))]
  fallback_errors <- character()
  if (length(no_seq) && !is.null(fallback_fetch)) {
    alternative <- tryCatch(fallback_fetch(symbol, no_seq, window_bp),
                            error = function(e) list(records = list(), errors = conditionMessage(e)))
    for (s in intersect(names(alternative$records), no_seq)) {
      rec <- alternative$records[[s]]
      if (is.null(rec)) next
      records[[s]] <- rec
      key <- paste(s, toupper(symbol), window_bp, sep = "|")
      if (!is.null(cache)) assign(key, rec, envir = cache)
      if (isTRUE(disk_cache)) fast_cache_write(rec, symbol, window_bp, cache_dir)
    }
    fallback_errors <- alternative$errors
  }
  still_missing <- selected[vapply(records[selected], is.null, logical(1))]
  # A successful independent fetch supersedes the original Ensembl failure.
  # Keep only a short provider summary unless a species remains unavailable.
  recovered <- setdiff(no_seq, still_missing)
  if (length(recovered)) {
    # Do not leave warnings about an Ensembl request that was recovered by UCSC.
    problems <- Filter(function(msg) {
      any(vapply(still_missing, function(s) {
        startsWith(msg, paste0(species_label(s), ":"))
      }, logical(1)))
    }, problems)
  }
  if (length(still_missing)) {
    problems <- c(problems,
                  sprintf("Missing after both providers: %s",
                          paste(vapply(still_missing, species_label, character(1)), collapse = ", ")),
                  fallback_errors)
  }
  list(records = unname(Filter(Negate(is.null), records[selected])),
       errors = unique(problems[nzchar(problems)]))
}
