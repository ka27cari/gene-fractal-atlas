# Offline tests for disease-specific GEO cohorts, sample-label guards,
# visible reference nucleotide traces and same-probe spatial concordance.
source("R/fractal.R")
source("R/svg.R")
source("R/atlas_cohorts.R")
source("R/data_catalog.R")

stopifnot("ASXL1" %in% CANDIDATE_GENES)
stopifnot(nrow(GEO_COHORT_CATALOG) == length(unique(GEO_COHORT_CATALOG$cohort)))
for (cancer in unique(GEO_COHORT_CATALOG$cancer)) {
  ids <- unname(catalog_cohorts(cancer))
  stopifnot(paste0("TCGA_", cancer) %in% ids, length(ids) > 1L)
  actual <- GEO_COHORT_CATALOG[GEO_COHORT_CATALOG$cohort %in% ids &
                                !grepl("^TCGA_", GEO_COHORT_CATALOG$cohort), , drop = FALSE]
  stopifnot(nrow(actual) > 0L, all(actual$cancer == cancer))
}
stopifnot(!"GSE32861" %in% unname(catalog_cohorts("BRCA")))
stopifnot(identical(catalog_geo_record("GSE53051_COAD")$study[[1L]], "GSE53051"))
stopifnot(identical(geo_matrix_url("GSE53051_COAD"),
                    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53051/matrix/GSE53051_series_matrix.txt.gz"))

check_groups <- function(titles, cohort, expected) {
  d <- data.frame(sample = paste0("GSM", seq_along(titles)), title = titles,
                  group = rep(NA_character_, length(titles)), stringsAsFactors = FALSE)
  z <- geo_restrict_study_groups(d, cohort)
  stopifnot(identical(z$group, expected))
}
check_groups(c("NK1", "RT163", "NB9", "TB198", "CP1", "P241"),
             "GSE52955_KIRC", c("Normal", "Tumor", NA, NA, NA, NA))
check_groups(c("NB9", "TB198", "NK1", "RT163"),
             "GSE52955_BLCA", c("Normal", "Tumor", NA, NA))
check_groups(c("CP1", "P241", "NB9", "TB198"),
             "GSE52955_PRAD", c("Normal", "Tumor", NA, NA))
check_groups(c("Genomic DNA from 2014ccRCC01T", "Genomic DNA from 2014ccRCC01N"),
             "GSE61441", c("Tumor", "Normal"))
check_groups(c("SAMPLE 1: genomic DNA from SCC_13", "SAMPLE 21: genomic DNA from Lung IPF1"),
             "GSE121849", c("Tumor", "Normal"))
check_groups(c("Breast_Cancer_1", "Breast_Normal_1", "Colon_Cancer_1", "Colon_Normal_1"),
             "GSE53051_BRCA", c("Tumor", "Normal", NA, NA))
check_groups(c("DNA methylation of normal esophagus tissues [C199N]", "DNA methylation of ESCC tissues [C199T]"),
             "GSE121930", c("Normal", "Tumor"))
stopifnot(is.na(geo_group_from_text("normal breast and breast cancer tissues")),
          identical(geo_group_from_text("lung, nontumor adjacent_12839"), "Normal"))

rec <- list(sequence = paste(rep("AACGATCGGCGTACG", 15L), collapse = ""),
            strand = 1L, region_start = 101L, region_end = 325L,
            chromosome = "1", assembly = "GRCh37")
coords <- make_cgr(rec, "triangle")
cp <- which(coords$is_cpg)[seq_len(3L)]
map <- data.frame(probe = sprintf("cg%08d", 1:3), position = coords$cpg_position[cp])

fixture <- c(
  '!Sample_title\t"Lung adenocarcinoma 1"\t"Adjacent non-tumor 1"\t"Lung adenocarcinoma 2"\t"Adjacent non-tumor 2"',
  '!series_matrix_table_begin',
  '"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"',
  sprintf('"%s"\t"0.91"\t"0.12"\t"0.85"\t"0.19"', map$probe[[1L]]),
  sprintf('"%s"\t"0.18"\t"0.72"\t"0.25"\t"0.65"', map$probe[[2L]]),
  '!series_matrix_table_end'
)
gz_file <- tempfile(fileext = ".gz")
w <- gzfile(gz_file, "wt")
writeLines(fixture, w)
close(w)
r <- gzfile(gz_file, "rt")
loaded <- geo_parse_series_connection(r, "GSE32861", map$probe)
close(r)
unlink(gz_file)
stopifnot(nrow(loaded$measurements) == 8L,
          all(c("Tumor", "Normal") %in% loaded$measurements$group))

make_summaries <- function(tumor, normal) {
  data.frame(probe = rep(map$probe, each = 2L),
             position = rep(map$position, each = 2L),
             group = rep(c("Tumor", "Normal"), 3L),
             mean_beta = c(rbind(tumor, normal)),
             prevalence = 0.5, n = 4L, stringsAsFactors = FALSE)
}
a <- make_summaries(c(.90, .15, .60), c(.10, .80, .50))
b <- make_summaries(c(.85, .85, .40), c(.15, .20, .55))
# Probes 1 and 2 qualify for both cohorts; probe 3 falls below 0.15 in cohort A.
agree <- atlas_spatial_concordance(list(`GSE-A` = a, `GSE-B` = b),
                                   coords, bins = 48L, min_samples = 2L,
                                   delta_cutoff = .15)
stopifnot(nrow(agree$probes) == 4L,
          all(agree$probes$probe %in% map$probe[1:2]),
          nrow(agree$cells) > 0L,
          all(agree$cells$fraction >= 0.5 & agree$cells$fraction <= 1))
svg <- atlas_concordance_svg(coords, "triangle", agree, bins = 48L)
stopifnot(grepl("#baff48", svg, fixed = TRUE),
          grepl("#e1ff98", svg, fixed = TRUE),
          grepl("atlas-cell-hit", svg, fixed = TRUE),
          grepl('data-cohort="Agreement"', svg, fixed = TRUE))
stopifnot(grepl("#baff48", atlas_svg(coords, data.frame(), "triangle"), fixed = TRUE))
cat("Multicancer catalog, GEO grouping, gzip reader and cohort-tick atlas passed\n")
