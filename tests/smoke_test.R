# Run from repository root: Rscript tests/smoke_test.R
source("R/fractal.R")
source("R/methylation.R")

# Exact geometry and strand-aware CpG anchoring on a controlled DNA fragment.
forward <- list(sequence = "ACGTACG", region_start = 100L, region_end = 106L,
                strand = 1L, source = "ensembl", species = "homo_sapiens",
                assembly = "GRCh38", chromosome = "1")
p <- make_cgr(forward)
stopifnot(nrow(p) == 7L, identical(which(p$is_cpg), c(2L, 6L)),
          identical(as.numeric(p$cpg_position[p$is_cpg]), c(101, 105)))
stopifnot(identical(reverse_complement("ACGTACG"), "CGTACGT"))

reverse <- forward
reverse$sequence <- reverse_complement(forward$sequence)
reverse$strand <- -1L
q <- make_cgr(reverse)
stopifnot(identical(as.numeric(q$cpg_position[q$is_cpg]), c(105, 101)))

stopifnot(all(is.finite(p$x)), all(is.finite(p$y)))
stopifnot(max(abs(make_cgr(forward, "square")$x)) <= 1)

csv <- tempfile(fileext = ".csv")
writeLines(c("species,assembly,chrom,position,beta,sample",
             "human,GRCh38,chr1,101,0.8,Blood",
             "homo_sapiens,GRCh38,1,105,0.2,Blood"), csv)
m <- read_methylation(csv)
stopifnot(identical(m$species, c("homo_sapiens", "homo_sapiens")))
matched <- match_methylation(forward, p, m, "Blood")
stopifnot(matched$matched == 2L, identical(matched$points$beta, c(0.8, 0.2)))
matched_rev <- match_methylation(reverse, q, m, "Blood")
stopifnot(matched_rev$matched == 2L)

writeLines(c("species,assembly,chrom,position,beta",
             "human,GRCh38,1,101,1.5"), csv)
stopifnot(inherits(try(read_methylation(csv), silent = TRUE), "try-error"))

cat("Smoke test passed: geometry, orientation, CpG coordinates and beta matching.\n")
