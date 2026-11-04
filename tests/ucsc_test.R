source("R/fractal.R")
source("R/ensembl.R")
source("R/fast.R")
source("R/ucsc.R")

# Independent provider mock: deterministic data, no network calls.
queries <- list()
mock_api <- function(urls) {
  queries[[length(queries) + 1L]] <<- urls
  r <- lapply(urls, function(url) {
    if (grepl("mygene.info", url, fixed = TRUE)) {
      taxid <- as.integer(sub(".*species=([0-9]+).*", "\\1", url))
      symbol <- if (taxid == 10090L) "Trp53" else if (taxid == 10116L) "Tp53" else "TP53"
      return(list(hits = list(list(taxid = taxid, symbol = symbol, entrezgene = 7157,
                                   genomic_pos = list(chr = "17", start = 1001L,
                                                      end = 2000L, strand = -1L),
                                   genomic_pos_hg19 = list(chr = "17", start = 1001L,
                                                           end = 2000L, strand = -1L)))))
    }
    if (grepl("api.genome.ucsc.edu", url, fixed = TRUE)) {
      return(list(dna = paste0(rep("ACGT", 250), collapse = ""),
                  genome = if (grepl("genome=hg19", url, fixed = TRUE)) "hg19" else
                    if (grepl("genome=mm10", url, fixed = TRUE)) "mm10" else
                      if (grepl("genome=rn6", url, fixed = TRUE)) "rn6" else "hg38"))
    }
    stop("unexpected URL")
  })
  list(results = r, errors = rep("", length(urls)))
}

selected <- c("homo_sapiens", "mus_musculus", "rattus_norvegicus")
result <- fetch_ucsc_fallback("TP53", selected, 1000L, api = mock_api)
stopifnot(length(result$records) == 3L,
          identical(names(result$records), selected),
          identical(result$records[[1L]]$assembly, "GRCh38"),
          identical(result$records[[2L]]$assembly, "GRCm38"),
          identical(result$records[[3L]]$assembly, "Rnor_6.0"),
          all(vapply(result$records, function(r) nchar(r$sequence) == 1000L &&
                       identical(r$source, "ucsc") && r$strand == -1L, logical(1))),
          length(queries) == 2L)
# A complete Ensembl timeout should enter the independent provider path.
all_timeout <- function(endpoints) {
  list(results = vector("list", length(endpoints)),
       errors = rep("Timeout was reached", length(endpoints)))
}
local_cache <- new.env(parent = emptyenv())
recovery <- fetch_comparative_regions_fast(
  "TP53", selected, 1000L, cache = local_cache,
  api_batch = all_timeout, disk_cache = FALSE,
  fallback_fetch = function(symbol, selected, window_bp) {
    fetch_ucsc_fallback(symbol, selected, window_bp, api = mock_api)
  }
)
stopifnot(length(recovery$records) == length(selected),
          !length(recovery$errors),
          all(vapply(recovery$records, function(r) identical(r$source, "ucsc"), logical(1))))
# Subsequent load uses memory cache even though both providers were mocked.
calls_before <- length(queries)
again <- fetch_comparative_regions_fast(
  "TP53", selected, 1000L, cache = local_cache,
  api_batch = all_timeout, disk_cache = FALSE,
  fallback_fetch = function(...) stop("Cache was not used")
)
stopifnot(length(again$records) == length(selected), length(queries) == calls_before)

# One-based inclusive to zero-based half-open conversion.
stopifnot(grepl("start=1499&end=2499", queries[[2]][[1L]], fixed = TRUE))
# No assembly guessing for chimpanzee.
stopifnot(!length(fetch_ucsc_fallback("TP53", "pan_troglodytes", 1000L, api = mock_api)$records))
# GRCh37 is separate from GRCh38 and matches Xena's hg19 coordinates.
grch37 <- fetch_ucsc_grch37_human("TP53", 1000L, api = mock_api)
stopifnot(identical(grch37$assembly, "GRCh37"),
          any(grepl("genome=hg19", queries[[length(queries)]], fixed = TRUE)))
# Refuse an annotation with the wrong target taxon.
bad <- list(hits = list(list(symbol = "TP53", taxid = 9606L,
                            genomic_pos = list(chr = "17", start = 100L, end = 200L, strand = 1L))))
stopifnot(is.null(mygene_pick_gene(bad, "TP53", "mus_musculus")))
# Fail closed when UCSC data length is wrong or the response assembly is wrong.
stopifnot(inherits(try(ucsc_record(list(dna = "ACGT"),
                                  list(chr = "17", start = 1000L, end = 2000L,
                                       strand = 1L, gene_id = "G", symbol = "TP53"),
                                  "homo_sapiens", 1000L), silent = TRUE), "try-error"))
cat("UCSC/MyGene provider contracts passed.\n")
