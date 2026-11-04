# Offline Ensembl API contract and species mapping tests.
source("R/fractal.R")
source("R/ensembl.R")

species <- unname(SPECIES)
stopifnot(length(species) == 7L)
ids <- stats::setNames(c("ENSG_TEST", "ENSPTR_TEST", "ENSMUS_TEST", "ENSRNO_TEST",
                         "ENSCAF_TEST", "ENSGAL_TEST", "ENSDAN_TEST"), species)
labels <- stats::setNames(c("TP53", "TP53", "Trp53", "Tp53", "TP53", "TP53", "tp53"), species)
assemblies <- stats::setNames(c("GRCh38", "Pan_tro_3.0", "GRCm39", "mRatBN7.2",
                                "ROS_Cfam_1.0", "GRCg7b", "GRCz11"), species)
mock_gene <- function(s) list(id = ids[[s]], object_type = "Gene", display_name = labels[[s]],
                              strand = if (s == "mus_musculus") -1L else 1L,
                              start = 1000L, end = 1700L, seq_region_name = "17")

calls <- list()
mock_api <- function(endpoints) {
  calls[[length(calls) + 1L]] <<- endpoints
  result <- vector("list", length(endpoints))
  errors <- rep("", length(endpoints))
  for (i in seq_along(endpoints)) {
    e <- endpoints[[i]]
    if (grepl("^/lookup/symbol/homo_sapiens/TP53", e)) {
      result[i] <- list(mock_gene("homo_sapiens"))
    } else if (grepl("^/homology/id/homo_sapiens/ENSG_TEST", e)) {
      homologies <- lapply(species[-1L], function(s) list(
        type = "ortholog_one2one", confidence = 1L,
        target = list(id = ids[[s]], species = s, perc_id = 92)))
      result[i] <- list(list(data = list(list(id = "ENSG_TEST", homologies = homologies))))
    } else if (grepl("^/lookup/id/", e)) {
      id <- sub("^/lookup/id/([^?]+).*", "\\1", e)
      s <- names(ids)[match(id, ids)]
      stopifnot(length(s) == 1L, !is.na(s))
      result[i] <- list(mock_gene(s))
    } else if (grepl("^/sequence/region/", e)) {
      s <- sub("^/sequence/region/([^/]+)/.*", "\\1", e)
      result[i] <- list(list(seq = paste(rep("ACGT", 250L), collapse = ""),
                             id = paste("chromosome", assemblies[[s]], "17", "500", "1499", "1", sep = ":")))
    } else stop(paste("Unexpected endpoint", e))
  }
  list(results = result, errors = errors)
}

cache <- new.env(parent = emptyenv())
r <- fetch_comparative_regions("TP53", species, 1000L, cache = cache, api_batch = mock_api)
stopifnot(length(r$records) == 7L, length(r$errors) == 0L,
          identical(vapply(r$records, `[[`, character(1), "species"), species),
          identical(vapply(r$records, `[[`, character(1), "symbol"), unname(labels)),
          length(calls) == 4L)

again <- fetch_comparative_regions("TP53", species, 1000L, cache = cache, api_batch = mock_api)
stopifnot(length(again$records) == 7L, length(calls) == 4L)

# Verify the corrected homology URL and stable test order.
stopifnot(grepl("/homology/id/homo_sapiens/ENSG_TEST", calls[[2L]][[1L]], fixed = TRUE))

# A valid direct-symbol fallback must be surfaced as fallback, not verified orthology.
fallback_api <- function(endpoints) {
  if (grepl("^/lookup/symbol", endpoints[[1L]])) {
    if (grepl("^/lookup/symbol/homo_sapiens", endpoints[[1L]])) {
      return(list(results = list(mock_gene("homo_sapiens")), errors = ""))
    }
    results <- lapply(endpoints, function(z) {
      if (grepl("/Trp53\\?", z)) mock_gene("mus_musculus") else NULL
    })
    return(list(results = results, errors = rep("Not found", length(endpoints))))
  }
  if (grepl("^/homology/id/", endpoints[[1L]])) {
    return(list(results = list(NULL), errors = "Not found"))
  }
  if (grepl("^/sequence/region/", endpoints[[1L]])) {
    return(list(results = list(list(seq = paste(rep("ACGT", 250), collapse = ""),
                                    id = "chromosome:GRCm39:17:500:1499:1")), errors = ""))
  }
  stop("Unexpected fallback request")
}
fall <- fetch_comparative_regions("TP53", "mus_musculus", 1000L,
                                  cache = new.env(parent = emptyenv()), api_batch = fallback_api)
stopifnot(length(fall$records) == 1L, grepl("symbol_fallback", fall$records[[1L]]$orthology))
cat("Ortholog lookup and fallback tests passed.\n")
