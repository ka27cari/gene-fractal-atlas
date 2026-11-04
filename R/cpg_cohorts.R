# Exploratory single-probe cohort browser.
# The source is a single HM450 probe in legacy TCGA matrices (hg19/GRCh37).
# No patient-level values are imputed; failures are reported per cohort.

xena_probe_query <- function(dataset, samples, probe) {
  if (length(probe) != 1L || is.na(probe) || !grepl("^cg[0-9]{8}$", probe))
    stop("A valid HM450 cg probe identifier is required")
  if (!length(samples)) stop("No samples to query")
  sprintf('(fetch [{:table %s :columns [%s] :samples %s}])',
          xena_literal(dataset), xena_literal(probe),
          paste0('[', paste(vapply(samples, xena_literal, character(1)), collapse = ' '), ']'))
}

parse_xena_probe_values <- function(response, n) {
  if (!is.list(response) || length(response) != 1L) stop("Unexpected Xena probe response")
  raw_list <- response[[1L]]
  if (length(raw_list) != n) stop(sprintf("Xena returned %s values for %s samples", length(raw_list), n))
  raw <- vapply(raw_list, function(x) {
    if (is.null(x) || !length(x)) return(NA_character_)
    as.character(x[[1L]])
  }, character(1))
  out <- suppressWarnings(as.numeric(raw))
  out[!is.finite(out) | out < 0 | out > 1] <- NA_real_
  out
}

all_tcga_sample_info <- function(all_samples) {
  samples <- sort(unique(as.character(unlist(all_samples, use.names = FALSE))))
  group <- tcga_sample_type(samples)
  keep <- !is.na(group)
  data.frame(sample = samples[keep], group = group[keep],
             patient = substr(samples[keep], 1L, 12L), stringsAsFactors = FALSE)
}

fetch_tcga_probe_cohort <- function(probe, cancer, post = NULL, chunk_size = 350L) {
  dataset <- xena_dataset(cancer)
  sample_query <- sprintf("(%s %s)", XENA_SAMPLE_QUERY, xena_literal(dataset))
  metadata <- all_tcga_sample_info(xena_post(sample_query, post = post))
  if (!nrow(metadata)) stop("No primary-tumor or solid-tissue-normal samples in dataset")
  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size < 1L || chunk_size > 1000L)
    stop("Invalid chunk size")
  chunks <- split(seq_len(nrow(metadata)), ceiling(seq_len(nrow(metadata)) / chunk_size))
  values <- rep(NA_real_, nrow(metadata))
  # Each chunk is queried in the same sample order. No sample-size cap.
  for (ids in chunks) {
    res <- xena_post(xena_probe_query(dataset, metadata$sample[ids], probe), post = post)
    values[ids] <- parse_xena_probe_values(res, length(ids))
  }
  metadata$beta <- values
  metadata$cohort <- cancer
  metadata$probe <- probe
  # Return the full matrix-eligible sample inventory, including missing measurements.
  metadata
}

cpg_cohort_summary <- function(frame) {
  if (!nrow(frame)) stop("No cohort rows")
  t <- frame$beta[frame$group == "Tumor" & is.finite(frame$beta)]
  n <- frame$beta[frame$group == "Normal" & is.finite(frame$beta)]
  tf <- frame[frame$group == "Tumor" & is.finite(frame$beta), , drop = FALSE]
  nf <- frame[frame$group == "Normal" & is.finite(frame$beta), , drop = FALSE]
  matched <- intersect(unique(tf$patient), unique(nf$patient))
  data.frame(
    cohort = frame$cohort[[1L]],
    tumor_samples = length(t), tumor_patients = length(unique(tf$patient)),
    normal_samples = length(n), normal_patients = length(unique(nf$patient)),
    paired_patients = length(matched),
    paired_median_delta = if (length(matched)) {
      paired <- paired_tcga_probe(frame)
      if (nrow(paired)) stats::median(paired$delta) else NA_real_
    } else NA_real_,
    tumor_median = if (length(t)) stats::median(t) else NA_real_,
    normal_median = if (length(n)) stats::median(n) else NA_real_,
    mean_delta = if (length(t) && length(n)) mean(t) - mean(n) else NA_real_,
    missing_beta = sum(!is.finite(frame$beta)),
    stringsAsFactors = FALSE
  )
}

paired_tcga_probe <- function(frame) {
  measured <- frame[is.finite(frame$beta), , drop = FALSE]
  if (!nrow(measured)) return(data.frame())
  # Average technical/sample replicates within patient and tissue type.
  z <- stats::aggregate(beta ~ patient + group, data = measured, FUN = mean)
  t <- z[z$group == "Tumor", c("patient", "beta"), drop = FALSE]
  n <- z[z$group == "Normal", c("patient", "beta"), drop = FALSE]
  if (!nrow(t) || !nrow(n)) return(data.frame())
  out <- merge(t, n, by = "patient", suffixes = c("_tumor", "_normal"))
  out$delta <- out$beta_tumor - out$beta_normal
  out
}

# The browser exposes only the CpG probes already verified to map to the current
# reference CpG by the preview loader, never inferred from image proximity.
cpg_probe_for_position <- function(methylation, position) {
  pos <- suppressWarnings(as.numeric(position))
  if (length(pos) != 1L || !is.finite(pos)) return(NULL)
  hits <- methylation[methylation$position == pos &
                        grepl("^cg[0-9]{8}$", methylation$probe), , drop = FALSE]
  if (!nrow(hits)) return(NULL)
  hits$probe[[1L]]
}

# Compact row-wise dot-strip: every assayed, nonmissing sample is shown.
# Horizontal position is beta; small deterministic vertical jitter separates dots.
cpg_cohort_strip <- function(frame) {
  measured <- frame[is.finite(frame$beta), , drop = FALSE]
  tick <- paste(vapply(seq(0, 1, by = 0.25), function(b) {
    x <- 71 + 407 * b
    sprintf('<line x1="%.1f" x2="%.1f" y1="19" y2="95" stroke="#20340c" stroke-width="1"/>', x, x)
  }, character(1)), collapse = "")
  points <- if (nrow(measured)) paste(vapply(seq_len(nrow(measured)), function(i) {
    row <- measured[i, , drop = FALSE]
    x <- 71 + 407 * row$beta[[1L]]
    cy <- if (row$group[[1L]] == "Tumor") 39 else 76
    cy <- cy + 9 * sin(i * 1.7 + row$beta[[1L]] * 13)
    col <- if (row$group[[1L]] == "Tumor") "#ff4444" else "#baff48"
    sprintf('<circle cx="%.2f" cy="%.2f" r="2.45" fill="%s" fill-opacity=".7"><title>%s | %s | β %.3f</title></circle>',
            x, cy, col, htmltools::htmlEscape(row$sample[[1L]], attribute = FALSE),
            row$group[[1L]], row$beta[[1L]])
  }, character(1)), collapse = "") else ""
  htmltools::HTML(paste0('<svg class="cohort-strip" viewBox="0 0 500 112" preserveAspectRatio="xMidYMid meet" role="img" aria-label="Beta by tissue group">',
    '<rect width="500" height="112" fill="#000"/>', tick,
    '<text x="6" y="42" fill="#ff7878" font-size="13">Tumor</text>',
    '<text x="6" y="78" fill="#c9ff93" font-size="13">Normal</text>',
    '<text x="70" y="108" fill="#8faa5b" font-size="10">0</text>',
    '<text x="474" y="108" fill="#8faa5b" font-size="10">1</text>',
    '<line x1="71" x2="478" y1="94" y2="94" stroke="#476026"/>', points, '</svg>'))
}
