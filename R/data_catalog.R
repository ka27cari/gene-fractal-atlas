# Data catalog: registered connections are not live evidence of assay coverage.
# A dataset becomes "loaded" only after actual measured data are received.
# No synthetic patient values and no claims about WGS/WGBS integration.

CANDIDATE_GENES <- c("TP53", "EGFR", "KRAS", "BRAF", "PIK3CA", "STK11", "KEAP1",
                     "ALK", "ROS1", "MET", "ERBB2", "BRCA1", "BRCA2", "APC",
                     "MYC", "TET2", "DNMT3A", "ASXL1", "PTEN", "RB1", "CDKN2A")

catalog_cancer_codes <- function() unname(TCGA_COHORTS)

# The repository entry confirms a public study with tumor and normal tissue.
# It DOES NOT confirm that the GEO series matrix is usable by this live adapter.
# Each load must check the actual sample labels, beta scale and CpG IDs.
# Cohort IDs with a suffix select only the corresponding disease from a
# MULTI-CANCER GEO series (no cross-cancer contamination).
GEO_COHORT_CATALOG <- data.frame(
  cancer = c(
    "LUAD", "LUAD", "LUAD", "LUAD",
    "LUSC", "LUSC", "BRCA", "BRCA",
    "COAD", "COAD", "COAD", "KIRC", "KIRC",
    "LIHC", "HNSC", "PRAD", "PRAD", "THCA",
    "BLCA", "CHOL", "PAAD", "PAAD", "CESC", "STAD", "ESCA"
  ),
  cohort = c(
    "GSE32861", "GSE62948", "GSE63384", "GSE32866",
    "GSE121849", "GSE56044_LUSC", "GSE69914", "GSE53051_BRCA",
    "GSE131013", "GSE48684", "GSE53051_COAD", "GSE61441", "GSE52955_KIRC",
    "GSE54503", "GSE87053", "GSE76938", "GSE52955_PRAD", "GSE53051_THCA",
    "GSE52955_BLCA", "GSE49656", "GSE181740", "GSE53051_PAAD", "GSE211668", "GSE164988", "GSE121930"
  ),
  study = c(
    "GSE32861", "GSE62948", "GSE63384", "GSE32866",
    "GSE121849", "GSE56044", "GSE69914", "GSE53051",
    "GSE131013", "GSE48684", "GSE53051", "GSE61441", "GSE52955",
    "GSE54503", "GSE87053", "GSE76938", "GSE52955", "GSE53051",
    "GSE52955", "GSE49656", "GSE181740", "GSE53051", "GSE211668", "GSE164988", "GSE121930"
  ),
  platform = c(
    rep("27K", 4), "450K", "450K", "450K", "450K",
    "450K", "450K", "450K", "450K", "450K",
    "450K", "450K", "450K", "450K", "450K", "450K",
    "450K", "450K", "450K", "450K", "450K", "450K"
  ),
  note = c(
    "Selamat LUAD", "Mansfield LUAD", "Robles LUAD", "same parent series as GSE32861",
    "LUSC + surrounding lung", "mixed lung study; squamous tumors only", "breast tumor and adjacent/healthy normal", "breast samples ONLY from mixed study",
    "paired colon tumor and normals + healthy", "colon tumor/normal/adenoma", "colon tumor/normal ONLY from mixed study",
    "46 KIRC tumor-normal pairs", "kidney subset of mixed urological study",
    "66 liver tumor-normal pairs", "oral HNSC tumor and adjacent normal", "prostate tumor vs benign",
    "prostate subset of mixed urological study", "thyroid PTC vs normal subset ONLY",
    "bladder subset of mixed urological study", "cholangiocarcinoma with 4 adjacent normals",
    "7 matched pancreatic ductal adenocarcinoma pairs", "pancreas tumor/normal only from mixed GSE53051 study", "cervical cancer and histologically normal tissue",
    "12 paired stomach carcinoma and normal tissues", "15 matched esophageal squamous carcinoma / normal pairs"
  ),
  stringsAsFactors = FALSE
)
GEO_COHORT_CATALOG$url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GEO_COHORT_CATALOG$study)

catalog_geo_record <- function(cohort) {
  z <- GEO_COHORT_CATALOG[GEO_COHORT_CATALOG$cohort == cohort, , drop = FALSE]
  if (nrow(z) != 1L) stop("Unregistered GEO cohort: ", cohort)
  z
}

catalog_cohorts <- function(cancer) {
  if (length(cancer) != 1L || is.na(cancer) || !cancer %in% catalog_cancer_codes())
    stop("Unknown cancer type")
  tcga <- stats::setNames(paste0("TCGA_", cancer), paste0("TCGA-", cancer, " · 450K"))
  z <- GEO_COHORT_CATALOG[GEO_COHORT_CATALOG$cancer == cancer, , drop = FALSE]
  if (!nrow(z)) return(tcga)
  choices <- stats::setNames(z$cohort, paste0(z$cohort, " · GEO ", z$platform))
  c(tcga, choices)
}

# These additional studies are public but their matrices need a separate
# normalization/parsing workflow; never advertise one-click connection.
EXTERNAL_COHORTS <- data.frame(
  cancer = c("COAD", "COAD"),
  cohort = c("GSE42752", "GSE101764"),
  assay = rep("Methylation array · 450K", 2L),
  note = c("Tumor/normal series; manual beta import", "Large tumor/normal series; manual beta import"),
  url = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", c("GSE42752", "GSE101764")),
  stringsAsFactors = FALSE
)

catalog_external_cohorts <- function(cancer) {
  EXTERNAL_COHORTS[EXTERNAL_COHORTS$cancer == cancer, , drop = FALSE]
}

catalog_rows <- function(cancer, loaded = character(), loaded_probes = integer()) {
  cohorts <- catalog_cohorts(cancer)
  ids <- unname(cohorts)
  z <- GEO_COHORT_CATALOG[match(ids, GEO_COHORT_CATALOG$cohort), , drop = FALSE]
  assay <- ifelse(grepl("^TCGA_", ids), "Methylation array · 450K",
                  paste0("Methylation array · ", z$platform))
  assay[is.na(assay)] <- "Methylation array · 450K"
  provider <- ifelse(grepl("^TCGA_", ids), "UCSC Xena / TCGA", "NCBI GEO")
  status <- ifelse(ids %in% loaded, "Loaded · beta measurements present",
                   "Public study · adapter configured; online load unverified")
  for (i in seq_along(ids)) {
    if (ids[[i]] %in% names(loaded_probes)) {
      count <- loaded_probes[[ids[[i]]]]
      if (length(count) == 1L && is.finite(count))
        status[[i]] <- sprintf("Loaded · %s measured probes", format(as.integer(count), big.mark = ","))
    }
  }
  out <- data.frame(Dataset = names(cohorts), Assay = assay, Provider = provider,
                    Status = status, stringsAsFactors = FALSE)
  extra <- catalog_external_cohorts(cancer)
  if (nrow(extra)) out <- rbind(out, data.frame(
    Dataset = paste(extra$cohort, "· manual import"), Assay = extra$assay,
    Provider = "NCBI GEO", Status = "Public; CSV preparation required", stringsAsFactors = FALSE))
  out
}

catalog_gene_choices <- function() {
  stats::setNames(CANDIDATE_GENES, CANDIDATE_GENES)
}

catalog_supported_dataset <- function(cancer, id) {
  is.character(id) && length(id) == 1L && !is.na(id) &&
    id %in% unname(catalog_cohorts(cancer))
}
