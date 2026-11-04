# User-provided methylation beta values are never inferred from the DNA sequence.
# CSV requirements: species, assembly, chrom, position, beta; optional sample.
# Coordinates: 1-based, forward-reference C of a CpG, in the declared assembly.

normalize_species <- function(x) {
  x <- tolower(gsub("[[:space:]-]+", "_", trimws(as.character(x))))
  aliases <- c(human = "homo_sapiens", chimp = "pan_troglodytes",
               chimpanzee = "pan_troglodytes", mouse = "mus_musculus",
               rat = "rattus_norvegicus", dog = "canis_lupus_familiaris",
               chicken = "gallus_gallus", zebrafish = "danio_rerio")
  matched <- x %in% names(aliases)
  x[matched] <- unname(aliases[x[matched]])
  x
}

read_methylation <- function(path) {
  if (file.info(path)$size > 15 * 1024^2) stop("CSV upload limit: 15 MB")
  data <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                          strip.white = TRUE, fileEncoding = "UTF-8")
  names(data) <- tolower(trimws(names(data)))
  required <- c("species", "assembly", "chrom", "position", "beta")
  if (!all(required %in% names(data))) {
    stop(sprintf("Missing CSV columns: %s", paste(setdiff(required, names(data)), collapse = ", ")))
  }
  if (!"sample" %in% names(data)) data$sample <- "sample_1"
  if (nrow(data) == 0L) stop("Methylation CSV contains no rows")

  data$species <- normalize_species(data$species)
  data$assembly <- trimws(as.character(data$assembly))
  data$chrom <- sub("^chr", "", trimws(as.character(data$chrom)), ignore.case = TRUE)
  data$sample <- trimws(as.character(data$sample))
  data$position <- suppressWarnings(as.numeric(data$position))
  data$beta <- suppressWarnings(as.numeric(data$beta))
  if (anyNA(data[, c("species", "assembly", "chrom", "sample", "position", "beta")])) {
    stop("Missing or non-numeric methylation values found")
  }
  if (any(!nzchar(data$species) | !nzchar(data$assembly) |
          !nzchar(data$chrom) | !nzchar(data$sample))) stop("Empty metadata field")
  if (any(!is.finite(data$position) | data$position < 1 | data$position != floor(data$position))) {
    stop("Positions must be positive 1-based integers")
  }
  if (any(!is.finite(data$beta) | data$beta < 0 | data$beta > 1)) {
    stop("Beta values must be between 0 and 1")
  }
  key <- paste(data$sample, data$species, data$assembly, data$chrom, data$position, sep = "\r")
  if (anyDuplicated(key)) stop("Duplicate CpG positions within one sample: resolve duplicates before uploading")
  data
}

methylation_sites_for_record <- function(record, methylation = NULL, sample = NULL) {
  if (identical(record$source, "demo")) {
    sites <- record$embedded_methylation
  } else if (is.null(methylation) || is.null(sample)) {
    sites <- data.frame(position = numeric(), beta = numeric())
  } else {
    sites <- methylation[
      methylation$species == record$species &
        tolower(methylation$assembly) == tolower(record$assembly) &
        tolower(methylation$chrom) == tolower(sub("^chr", "", record$chromosome, ignore.case = TRUE)) &
        methylation$sample == sample,
      c("position", "beta"), drop = FALSE
    ]
  }
  if (!nrow(sites)) sites <- data.frame(position = numeric(), beta = numeric())
  sites
}

match_methylation <- function(record, geometry_data, methylation = NULL, sample = NULL) {
  sites <- methylation_sites_for_record(record, methylation, sample)
  cpg_rows <- which(geometry_data$is_cpg)
  if (nrow(sites) == 0L || length(cpg_rows) == 0L) {
    return(list(points = data.frame(), supplied = nrow(sites), matched = 0L))
  }
  matched <- match(sites$position, geometry_data$cpg_position[cpg_rows])
  valid <- !is.na(matched)
  rows <- cpg_rows[matched[valid]]
  points <- geometry_data[rows, , drop = FALSE]
  points$beta <- sites$beta[valid]
  list(points = points, supplied = nrow(sites), matched = nrow(points))
}
