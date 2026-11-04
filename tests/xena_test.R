# Offline TCGA sample classification and Xena subset parsing contract.
source("R/fractal.R")
source("R/ensembl.R")
source("R/methylation.R")
source("R/xena.R")

samples <- c("TCGA-AB-1234-01", "TCGA-AB-1234-11", "TCGA-AB-5678-11", "TCGA-AB-5678-01",
             "TCGA-AB-9999-01")
choice <- choose_tcga_samples(samples)
stopifnot(nrow(choice) == 4L, all(choice$paired),
          identical(choice$group, c("Tumor", "Normal", "Tumor", "Normal")),
          identical(as.character(tcga_sample_type(c("TCGA-AB-5678-01", "TCGA-AB-5678-11"))),
                    c("Tumor", "Normal")))

rec <- list(species = "homo_sapiens", species_label = "Human", symbol = "TP53",
            gene_id = "ENSG_TEST", assembly = "GRCh37", chromosome = "17",
            region_start = 100L, region_end = 1099L, strand = 1L,
            sequence = paste(rep("ACGT", 250), collapse = ""), source = "ensembl")
stopifnot(101 %in% make_cgr(rec)$cpg_position)

mock_post <- function(query) {
  if (grepl("sampleID", query, fixed = TRUE)) return(as.list(samples))
  if (grepl("HumanMethylation450", query, fixed = TRUE)) {
    return(list(
      list(name = list("cgEXAMPLE1", "cgEXAMPLE2"),
           position = list(list("chr17", 100, 101, "+"), list("chr17", 104, 105, "+"))),
      list(list(.20, .90, .74, .12), list(.81, .08, .90, .36))
    ))
  }
  stop("Unexpected Xena query")
}
mock_reference <- function(symbol, window_bp) {
  stopifnot(identical(symbol, "TP53"), identical(as.integer(window_bp), 1000L))
  rec
}

out <- fetch_tcga_luad("TP53", window_bp = 1000L, post = mock_post,
                       reference_fetch = mock_reference)
stopifnot(nrow(out$samples) == 4L, nrow(out$methylation) == 8L,
          identical(out$record$assembly, "GRCh37"),
          all(out$methylation$position %in% c(101, 105)),
          all(out$methylation$beta >= 0 & out$methylation$beta <= 1))

# The C/mC split should change CGR geometry only at observed methylated CpGs.
sites <- methylation_sites_for_record(rec, out$methylation, out$samples$sample[[2L]])
coords <- make_cgr(rec, "mc_triangle", methylation_sites = sites,
                   methylation_threshold = 0.6)
stopifnot(any(coords$is_methyl_c), all(coords$draw_base[coords$is_methyl_c] == "mC"))

stopifnot(identical(parse_probe_position(list("chr17", 100, 101, "+"))$start, 100))
cat("TCGA barcode, Xena beta, genomic alignment and mC geometry tests passed.\n")
