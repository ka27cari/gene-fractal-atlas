# Cohort delta cutoffs are applied to real matching-probe summary values,
# not pixel positions or different CpGs in the same spatial CGR bin.
source("R/fractal.R")
source("R/svg.R")
source("R/xena.R")
source("R/atlas_cohorts.R")
source("R/data_catalog.R")

stopifnot("ASXL1" %in% CANDIDATE_GENES,
          "ASXL1" %in% unname(catalog_gene_choices()))
stopifnot(identical(catalog_external_cohorts("BRCA")$cohort, character(0)))
stopifnot(all(c("GSE131013", "GSE48684", "GSE53051_COAD") %in%
              unname(catalog_cohorts("COAD"))))
stopifnot("GSE42752" %in% catalog_external_cohorts("COAD")$cohort)
stopifnot("GSE76938" %in% unname(catalog_cohorts("PRAD")))
stopifnot("GSE69914" %in% unname(catalog_cohorts("BRCA")))

summary <- data.frame(
  probe = c("cg00000001", "cg00000001", "cg00000002", "cg00000002",
            "cg00000003", "cg00000003"),
  position = c(103, 103, 106, 106, 109, 109),
  group = c("Tumor", "Normal", "Tumor", "Normal", "Tumor", "Normal"),
  mean_beta = c(.90, .20, .45, .40, .95, .20),
  prevalence = c(1, 0, .5, .5, 1, 0),
  n = c(6L, 5L, 4L, 4L, 1L, 5L), stringsAsFactors = FALSE)

eligible <- atlas_filter_delta_probes(summary, min_samples = 2L, cutoff = 0)
stopifnot(nrow(eligible) == 2L,
          base::all(sort(round(eligible$delta, 2)) == c(.05, .70)))
strong <- atlas_filter_delta_probes(summary, 2L, 0.15)
stopifnot(nrow(strong) == 1L, strong$probe[[1L]] == "cg00000001",
          abs(strong$delta[[1L]] - 0.70) < 1e-12)
nohits <- atlas_filter_delta_probes(summary, 2L, 0.9)
stopifnot(nrow(nohits) == 0L, all(c("probe", "position", "delta") %in% names(nohits)))
stopifnot(inherits(try(atlas_filter_delta_probes(summary, 2L, 1.2), silent = TRUE), "try-error"))

rec <- list(sequence = "AACGATCGGCGTACG", strand = 1L,
            region_start = 101L, region_end = 115L, chromosome = "1")
coords <- make_cgr(rec, "triangle")
# A no-hit delta SVG still preserves reference triangle instead of inventing dots.
svg <- atlas_svg(coords, nohits, "triangle", delta = TRUE, min_samples = 1L)
stopifnot(grepl("<svg", svg, fixed = TRUE), grepl("#000000", svg, fixed = TRUE))
cat("Delta cutoff, ASXL1 preset, and multicancer GEO catalog checks passed\n")
