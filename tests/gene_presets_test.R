# Offline contract for popular symbols and expanded windows.
# Network-free fixture tests do NOT establish live Ensembl availability.
source("R/fractal.R")
source("R/ensembl.R")
source("R/fast.R")

presets <- c("TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "BRAF",
             "APC", "MYC", "TET2", "DNMT3A", "ASXL1", "PIK3CA")
stopifnot(identical(toupper("brca"), "BRCA"))
stopifnot(grepl("BRCA1 or BRCA2", gene_input_issue("BRCA"), fixed = TRUE))
stopifnot(all(vapply(presets, function(s) is.null(gene_input_issue(s)), logical(1))))
stopifnot(!is.null(gene_input_issue("TP53/BRCA1")))

species <- c("homo_sapiens", "mus_musculus", "rattus_norvegicus", "danio_rerio")
for (symbol in presets) {
  mock <- function(endpoints) {
    results <- vector("list", length(endpoints))
    errors <- rep("", length(endpoints))
    for (i in seq_along(endpoints)) {
      ep <- endpoints[[i]]
      if (grepl("^/lookup/symbol/", ep)) {
        s <- sub("^/lookup/symbol/([^/]+)/.*", "\\1", ep)
        results[i] <- list(list(id = paste0("TEST_", s), species = s,
                                object_type = "Gene", display_name = preferred_symbol(symbol, s),
                                seq_region_name = "17", strand = 1L, start = 1001L,
                                end = 4000000L, assembly_name = "TEST"))
      } else if (grepl("^/sequence/region/", ep)) {
        results[i] <- list(list(seq = paste0(rep("ACGT", 250L), collapse = ""),
                                id = "chromosome:TEST:17:501:1500:1"))
      } else stop(paste("Unexpected mock endpoint:", ep))
    }
    list(results = results, errors = errors)
  }
  loaded <- fetch_comparative_regions_fast(symbol, species, 1000L,
                                            api_batch = mock, disk_cache = FALSE)
  stopifnot(length(loaded$records) == length(species))
  stopifnot(all(vapply(loaded$records, function(x) identical(x$assembly, "TEST"), logical(1))))
}

# Reject a cross-species assembly conflict; do not silently label mouse GRCh38.
mouse <- list(id = "ENSMUSG_TEST", object_type = "Gene", species = "mus_musculus",
              seq_region_name = "11", strand = 1L, start = 1001L, end = 2001L,
              assembly_name = "GRCm39", display_name = "Trp53")
interval <- gene_interval(mouse, 1000L)
stopifnot(inherits(try(record_from_sequence("mus_musculus", mouse, interval,
            list(seq = paste(rep("A", 1000L), collapse = ""),
                 id = "chromosome:GRCh38:11:501:1500:1"), 1000L), silent = TRUE), "try-error"))

# Expanded interval, up to 1Mb (still below Ensembl region API's 10Mb limit).
for (window_bp in c(5000L, 50000L, 100000L, 250000L, 500000L, 1000000L)) {
  iv <- gene_interval(mouse, window_bp)
  stopifnot(iv$end - iv$start + 1L == window_bp)
}

# 100kb midpoint recurrence with finite positions and no color/geometry change.
rec <- list(sequence = paste0(rep("ACGT", 25000L), collapse = ""),
            region_start = 1L, region_end = 100000L, strand = 1L,
            species = "homo_sapiens", assembly = "GRCh38", chromosome = "1",
            source = "ensembl")
p <- make_cgr(rec, "square")
stopifnot(nrow(p) == 100000L, all(is.finite(p$x)), all(is.finite(p$y)),
          all(p$x >= 0 & p$x <= 1), all(p$y >= 0 & p$y <= 1))
cat("Gene fixture checks passed: 12 symbols, ambiguous BRCA, windows and assembly conflict.\n")
