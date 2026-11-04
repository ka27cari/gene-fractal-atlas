# The configured TCGA cohort list contains 18 unique entries.
source("R/ensembl.R")
source("R/fractal.R")
source("R/methylation.R")
source("R/xena.R")
stopifnot(length(TCGA_COHORTS) == 18L, length(unique(TCGA_COHORTS)) == 18L,
          identical(xena_dataset("LUAD"), "TCGA.LUAD.sampleMap/HumanMethylation450"),
          identical(xena_dataset("THCA"), "TCGA.THCA.sampleMap/HumanMethylation450"),
          inherits(try(xena_dataset("FAKE"), silent = TRUE), "try-error"))

# Entire workflow with a small synthetic API fixture; no network required.
samples <- c("TCGA-AB-1234-01", "TCGA-AB-1234-11")
rec <- list(species = "homo_sapiens", species_label = "Human", symbol = "TP53",
            gene_id = "ENSG_TEST", assembly = "GRCh37", chromosome = "17",
            region_start = 100L, region_end = 1099L, strand = 1L,
            sequence = paste(rep("ACGT", 250), collapse = ""), source = "ensembl")
mock_post <- function(query) {
  stopifnot(grepl('TCGA.BRCA.sampleMap/HumanMethylation450', query, fixed = TRUE))
  if (grepl('sampleID', query, fixed = TRUE)) return(as.list(samples))
  list(
    list(name = list("cgEXAMPLE1"), position = list(list("chr17", 100, 101, "+"))),
    list(list(.21, .82))
  )
}
mock_ref <- function(symbol, window_bp) rec
out <- fetch_tcga("TP53", "BRCA", window_bp = 1000L, post = mock_post,
                  reference_fetch = mock_ref)
stopifnot(identical(out$cohort, "BRCA"), nrow(out$samples) == 2L,
          nrow(out$methylation) == 2L,
          identical(out$dataset, "TCGA.BRCA.sampleMap/HumanMethylation450"))
cat("18-cohort routing and sample parsing checks passed.\n")
