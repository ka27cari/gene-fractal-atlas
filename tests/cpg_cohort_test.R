# Run from the repository root: Rscript tests/cpg_cohort_test.R
source("R/xena.R")
source("R/cpg_cohorts.R")
source("R/svg.R")
source("R/fractal.R")

samples <- c("TCGA-AA-0001-01", "TCGA-AA-0001-11", "TCGA-AA-0002-01",
             "TCGA-AA-0002-11", "TCGA-AA-0003-01", "TCGA-AA-0004-01")
values <- c("0.89", "0.23", "0.73", "0.22", "NaN", "0.61")

# Two independent cohort-matrix calls per 3-sample chunk. Never subsample.
queries <- character()
mock_post <- function(query) {
  queries <<- c(queries, query)
  if (grepl(":field.name \"sampleID\"", query, fixed = TRUE)) return(as.list(samples))
  selected <- regmatches(query, gregexpr("TCGA-[A-Z]{2}-[0-9]{4}-[0-9]{2}", query))[[1L]]
  list(as.list(values[match(selected, samples)]))
}

out <- fetch_tcga_probe_cohort("cg01234567", "LUAD", post = mock_post, chunk_size = 3L)
stopifnot(nrow(out) == length(samples), length(queries) == 3L,
          identical(out$sample, samples),
          sum(is.na(out$beta)) == 1L, all(out$beta[is.finite(out$beta)] <= 1))
stopifnot(identical(tcga_sample_type(out$sample), out$group))
summary <- cpg_cohort_summary(out)
paired <- paired_tcga_probe(out)
stopifnot(summary$tumor_samples == 3L, summary$normal_samples == 2L,
          summary$paired_patients == 2L, nrow(paired) == 2L,
          abs(summary$paired_median_delta - mean(c(.89 - .23, .73 - .22))) < 1e-10)
stopifnot(grepl("cg01234567", xena_probe_query(xena_dataset("BRCA"), samples, "cg01234567")))

# No accidental silent resampling: query every sample even with 1253 records.
large <- sprintf("TCGA-AA-%04d-01", seq_len(1253))
mock_large <- function(query) {
  if (grepl(":field.name \"sampleID\"", query, fixed = TRUE)) return(as.list(large))
  selected <- regmatches(query, gregexpr("TCGA-AA-[0-9]{4}-01", query))[[1L]]
  list(as.list(rep("0.50", length(selected))))
}
large_out <- fetch_tcga_probe_cohort("cg01234567", "LUAD", post = mock_large, chunk_size = 350L)
stopifnot(nrow(large_out) == 1253L, all(large_out$beta == .5))

mini <- list(sequence = "ACGTCG", region_start = 101L, region_end = 106L,
             strand = 1L, chromosome = "7", species = "homo_sapiens",
             assembly = "GRCh37", source = "ensembl")
coords <- make_cgr(mini, "triangle")
marked <- coords[coords$is_cpg, , drop = FALSE]
svg <- svg_click_targets(marked, "triangle", "TCGA-AA-0001-01")
stopifnot(grepl('class="cpg-hit"', svg, fixed = TRUE),
          grepl('data-cpg-position="102"', svg, fixed = TRUE))

meth <- data.frame(position = c(102, 102, 105),
                   probe = c("cg01234567", "cg01234567", "cg76543210"))
stopifnot(identical(cpg_probe_for_position(meth, 102), "cg01234567"),
          is.null(cpg_probe_for_position(meth, 999)))
cat("Cohort-wide CpG and SVG interaction regression checks passed.\n")
