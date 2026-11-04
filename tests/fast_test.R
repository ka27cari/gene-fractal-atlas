# Offline contract for quick one-round symbol lookup + one-round sequence fetch.
source("R/fractal.R")
source("R/ensembl.R")
source("R/fast.R")

species <- unname(SPECIES)
ids <- stats::setNames(paste0("TEST_", seq_along(species)), species)
labels <- stats::setNames(c("TP53", "TP53", "Trp53", "Tp53", "TP53", "TP53", "tp53"), species)
mock_gene <- function(s) list(id = ids[[s]], species = s, display_name = labels[[s]],
                              object_type = "Gene", seq_region_name = "17", strand = 1L,
                              start = 1001L, end = 1700L, assembly_name = "TEST")

calls <- list()
mock_api <- function(endpoints) {
  calls[[length(calls) + 1L]] <<- endpoints
  results <- vector("list", length(endpoints))
  errors <- rep("", length(endpoints))
  for (i in seq_along(endpoints)) {
    endpoint <- endpoints[[i]]
    if (grepl("^/lookup/symbol/", endpoint)) {
      s <- sub("^/lookup/symbol/([^/]+)/.*", "\\1", endpoint)
      sym <- sub("^/lookup/symbol/[^/]+/([^?]+)\\?.*", "\\1", endpoint)
      if (identical(sym, labels[[s]])) results[i] <- list(mock_gene(s))
      else errors[[i]] <- "Not found (HTTP 404)"
    } else if (grepl("^/sequence/region/", endpoint)) {
      s <- sub("^/sequence/region/([^/]+)/.*", "\\1", endpoint)
      results[i] <- list(list(seq = paste0(rep("ACGT", 250L), collapse = ""),
                              id = "chromosome:TEST:17:501:1500:1"))
    } else stop(sprintf("Unexpected request: %s", endpoint))
  }
  list(results = results, errors = errors)
}

cache <- new.env(parent = emptyenv())
r <- fetch_comparative_regions_fast("TP53", species, 1000L,
                                    cache = cache, api_batch = mock_api, disk_cache = FALSE)
stopifnot(length(r$records) == length(species), length(calls) == 2L,
          all(vapply(r$records, function(z) nchar(z$sequence) == 1000L, logical(1))),
          identical(vapply(r$records, `[[`, character(1), "species"), species),
          identical(r$records[[1L]]$orthology, "reference"),
          all(vapply(r$records[-1L], function(z) identical(z$orthology, "symbol_match_unverified"), logical(1))))

again <- fetch_comparative_regions_fast("TP53", species, 1000L,
                                        cache = cache, api_batch = mock_api, disk_cache = FALSE)
stopifnot(length(again$records) == length(species), length(calls) == 2L)

# Persistent cache is independent of the in-process environment and expires.
root <- tempfile("gene-fractal-cache-")
dir.create(root)
for (rec in r$records) fast_cache_write(rec, "TP53", 1000L, cache_dir = root)
stopifnot(!is.null(fast_cache_read(species[[2L]], "TP53", 1000L, cache_dir = root)))
calls <- list()
on_disk <- fetch_comparative_regions_fast("TP53", species, 1000L,
                                          cache = new.env(parent = emptyenv()),
                                          api_batch = mock_api, cache_dir = root, disk_cache = TRUE)
stopifnot(length(on_disk$records) == length(species), length(calls) == 0L)
old_path <- fast_cache_path(species[[1L]], "TP53", 1000L, cache_dir = root)
Sys.setFileTime(old_path, Sys.time() - 16 * 86400)
stopifnot(is.null(fast_cache_read(species[[1L]], "TP53", 1000L, cache_dir = root)))
unlink(root, recursive = TRUE)

# A missing symbol should enter an ortholog fallback, not be shown as a verified
# direct symbol hit. Avoid checking a remote endpoint here.
mock_with_failure <- function(endpoints) {
  if (any(grepl("/homology/id/", endpoints, fixed = TRUE))) {
    return(list(results = list(NULL), errors = "Not found (HTTP 404)"))
  }
  out <- mock_api(endpoints)
  for (i in seq_along(endpoints)) {
    if (grepl("^/lookup/symbol/rattus_norvegicus/", endpoints[[i]])) {
      out$results[i] <- list(NULL)
      out$errors[[i]] <- "Not found (HTTP 404)"
    }
  }
  out
}
failure <- fetch_comparative_regions_fast("TP53", c("homo_sapiens", "rattus_norvegicus"), 1000L,
                                          api_batch = mock_with_failure, disk_cache = FALSE)
stopifnot(length(failure$records) == 1L, failure$records[[1L]]$species == "homo_sapiens",
          any(grepl("Rat: gene lookup unavailable", failure$errors, fixed = TRUE)))
cat("Fast fetch tests passed: two-round common path, memory/disk reuse, expiration, failure.\n")
