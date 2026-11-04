# Cancer-specific cohort-level CGR analysis: reference-anchored geometry and
# cohort-level summaries of *measured* sample beta values.
# This module does not infer sequence mutations or CpG methylation from DNA.

LUAD_ATLAS_COHORTS <- c(
  "TCGA-LUAD · 450K" = "TCGA_LUAD",
  "GSE32861 · 27K · Selamat" = "GSE32861",
  "GSE62948 · 27K · Mansfield" = "GSE62948",
  "GSE63384 · 27K · Robles" = "GSE63384",
  "GSE32866 · 27K · validation from GSE32867" = "GSE32866"
)

# GSE32861 and GSE32866 are subseries of the same larger study;
# never treat them as independent publications in replication claims.
LUAD_GEO <- unname(LUAD_ATLAS_COHORTS)[-1L]

atlas_key <- function(x) paste(x, collapse = "|")

# Durable, user-owned cache. RDS contains only the *queried gene-window* and
# measured samples/probes; it does not mirror entire TCGA/GEO beta matrices.
# Files persist across project updates, and no arbitrary TTL triggers downloads.
# Override for a portable drive or tests: GFA_DATA_DIR or option(gfa.cache_root).
atlas_cache_root <- function() {
  configured <- getOption("gfa.cache_root", NULL)
  if (is.null(configured) || !nzchar(configured))
    configured <- Sys.getenv("GFA_DATA_DIR", unset =
                              file.path(path.expand("~"), "GeneFractalAtlas", "data"))
  file.path(path.expand(configured), "atlas-v1")
}

atlas_disk_cache <- function(cohort, gene, window_bp, cache_dir = atlas_cache_root()) {
  file.path(cache_dir, sprintf("%s_%s_%s.rds", cohort,
           gsub("[^A-Za-z0-9_-]", "", toupper(gene)), as.integer(window_bp)))
}

atlas_reference_file <- function(cancer, gene, window_bp,
                                 cache_dir = atlas_cache_root()) {
  file.path(cache_dir, sprintf("reference_%s_%s_%d.rds", cancer,
           gsub("[^A-Za-z0-9_-]", "", toupper(gene)), as.integer(window_bp)))
}

atlas_cache_atomic_write <- function(x, file) {
  tryCatch({
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = "atlas-", tmpdir = dirname(file))
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    # On Windows file.rename cannot overwrite an existing destination.
    if (file.exists(file)) unlink(file)
    isTRUE(file.rename(tmp, file))
  }, error = function(e) FALSE)
}

atlas_cached_read <- function(cohort, gene, window_bp, ttl_days = Inf) {
  file <- atlas_disk_cache(cohort, gene, window_bp)
  legacy <- atlas_disk_cache(cohort, gene, window_bp, "cache/atlas-v1")
  for (candidate in unique(c(file, legacy))) {
    if (!file.exists(candidate)) next
    if (is.finite(ttl_days) &&
        as.numeric(difftime(Sys.time(), file.info(candidate)$mtime,
                           units = "days")) > ttl_days) next
    x <- tryCatch(readRDS(candidate), error = function(e) NULL)
    if (!is.list(x) || !identical(x$cohort, cohort) ||
        !is.data.frame(x$measurements) ||
        !all(c("sample", "group", "probe", "beta") %in% names(x$measurements))) next
    if (!identical(candidate, file)) atlas_cache_atomic_write(x, file)
    return(x)
  }
  NULL
}

atlas_cached_write <- function(x, gene, window_bp) {
  atlas_cache_atomic_write(x, atlas_disk_cache(x$cohort, gene, window_bp))
}

atlas_reference_read <- function(cancer, gene, window_bp) {
  file <- atlas_reference_file(cancer, gene, window_bp)
  legacy <- atlas_reference_file(cancer, gene, window_bp, "cache/atlas-v1")
  for (candidate in unique(c(file, legacy))) {
    if (!file.exists(candidate)) next
    cached <- tryCatch(readRDS(candidate), error = function(e) NULL)
    if (is.null(cached) || !is.list(cached$reference) ||
        !is.data.frame(cached$map) ||
        !all(c("probe", "position") %in% names(cached$map)) ||
        !identical(cached$reference$assembly, "GRCh37") ||
        !identical(nchar(cached$reference$sequence), as.integer(window_bp))) next
    if (!identical(candidate, file)) atlas_cache_atomic_write(cached, file)
    return(cached)
  }
  NULL
}

atlas_reference_write <- function(cancer, gene, window_bp, reference, probe_map) {
  atlas_cache_atomic_write(list(reference = reference, map = probe_map),
                           atlas_reference_file(cancer, gene, window_bp))
}

atlas_cache_usage <- function() {
  root <- atlas_cache_root()
  files <- if (dir.exists(root)) list.files(root, pattern = "\\.rds$",
                                             full.names = TRUE, recursive = TRUE) else character()
  size <- if (length(files)) sum(file.info(files)$size, na.rm = TRUE) else 0
  list(path = root, n_files = length(files), bytes = size)
}

# Map probes using the current reference sequence and a small public Xena read.
# Returned positions must have been validated as forward-reference CpG cytosines.
atlas_probe_map <- function(symbol, reference_record, post = NULL, cancer = "LUAD") {
  dataset <- xena_dataset(cancer)
  raw <- xena_post(sprintf("(%s %s)", XENA_SAMPLE_QUERY, xena_literal(dataset)), post = post)
  meta <- all_tcga_sample_info(raw)
  if (!nrow(meta)) stop(sprintf("No Xena TCGA-%s samples available for probe mapping", cancer))
  # A handful of samples helps avoid dropping low-coverage loci from the map.
  use <- unique(c(head(meta$sample[meta$group == "Tumor"], 3L),
                  head(meta$sample[meta$group == "Normal"], 3L)))
  si <- data.frame(sample = use, group = tcga_sample_type(use), paired = FALSE)
  observed <- xena_gene_records(symbol, si, reference_record, post = post, dataset = dataset)
  if (!nrow(observed)) stop("No validated CpG coordinates returned by Xena")
  map <- unique(observed[, c("probe", "position"), drop = FALSE])
  map <- map[grepl("^cg[0-9]{8}$", map$probe) & is.finite(map$position), , drop = FALSE]
  map <- map[!duplicated(map$probe), , drop = FALSE]
  if (!nrow(map)) stop("No validated cg probe IDs in this gene")
  rownames(map) <- NULL
  map
}

# All eligible primary-tumor and solid-tissue-normal TCGA samples for the
# selected gene probes. The Xena gene query runs in bounded sample chunks.
atlas_fetch_tcga <- function(symbol, reference_record, probe_map, post = NULL,
                             chunk_size = 120L, cancer = "LUAD") {
  dataset <- xena_dataset(cancer)
  raw <- xena_post(sprintf("(%s %s)", XENA_SAMPLE_QUERY, xena_literal(dataset)), post = post)
  meta <- all_tcga_sample_info(raw)
  if (!nrow(meta)) stop(sprintf("No eligible TCGA-%s tumor or normal samples", cancer))
  chunks <- split(seq_len(nrow(meta)), ceiling(seq_len(nrow(meta)) / chunk_size))
  frames <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    ids <- chunks[[i]]
    si <- data.frame(sample = meta$sample[ids], group = meta$group[ids], paired = FALSE)
    res <- xena_gene_records(symbol, si, reference_record, post = post, dataset = dataset)
    if (nrow(res)) {
      res <- res[res$probe %in% probe_map$probe, c("sample", "group", "probe", "beta"), drop = FALSE]
      frames[[i]] <- res
    }
  }
  frames <- Filter(function(x) !is.null(x) && nrow(x), frames)
  if (!length(frames)) stop(sprintf("TCGA-%s has no measured CpGs matching the selected reference window", cancer))
  measurements <- do.call(rbind, frames)
  measurements <- measurements[!duplicated(paste(measurements$sample, measurements$probe)), , drop = FALSE]
  if (!all(c("Tumor", "Normal") %in% unique(measurements$group))) {
    stop(sprintf("TCGA-%s: usable probes are not present in both tumor and normal specimens", cancer))
  }
  measurements$cohort <- paste0("TCGA_", cancer)
  list(cohort = paste0("TCGA_", cancer), measurements = measurements,
       available_samples = meta, platform = "450K", source = "UCSC Xena / TCGA")
}

# Strict GEO series-matrix reader. Never assign tissue labels where the
# public sample metadata are ambiguous, and never reinterpret M-values as beta.
geo_matrix_url <- function(accession) {
  study <- strsplit(as.character(accession), "_", fixed = TRUE)[[1L]][[1L]]
  if (!grepl("^GSE[0-9]+$", study)) stop("Invalid GEO study accession")
  prefix <- sub("[0-9]{3}$", "nnn", study)
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
          prefix, study, study)
}

# Keep the legacy fixture interface for offline tests. Live use below is
# streaming so hundreds-of-MB 450K matrices are never read into RAM whole.
geo_download_matrix_lines <- function(accession) {
  url <- geo_matrix_url(accession)
  h <- curl::new_handle(timeout = 120L, connecttimeout = 12L)
  r <- curl::curl_fetch_memory(url, handle = h)
  if (r$status_code != 200L) stop("GEO HTTP ", r$status_code)
  if (length(r$content) > 75 * 1024^2) stop("Too large for legacy in-memory parser; use streaming")
  con <- gzcon(rawConnection(r$content, open = "rb"))
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

geo_meta_values <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
  out <- fields[-1L]
  gsub('^"|"$', "", out)
}

geo_group_from_text <- function(text) {
  s <- tolower(trimws(paste(text, collapse = " ")))
  is_normal <- grepl("adjacent|non[ -]?tumo[u]?r|nontumo[u]?r|normal|non[ -]?malignant|control|benign", s)
  # The token in 'non-tumor' must not itself trigger the tumor class.
  tumor_text <- gsub("non[ -]?tumo[u]?r|nontumo[u]?r", "", s)
  is_tumor <- grepl("adenocarcinoma|tumo[u]?r|carcinoma|cancer|lesion", tumor_text)
  # Ambiguous descriptions can mention BOTH categories (e.g. study-level
  # free text). Never assign them to Normal just because that token comes first.
  if (is_normal && is_tumor) return(NA_character_)
  if (is_normal) return("Normal")
  if (is_tumor) return("Tumor")
  NA_character_
}

geo_sample_info <- function(lines, sample_ids) {
  meta_rows <- lines[grepl("^!Sample_(title|source_name_ch1|characteristics_ch1|description)\t", lines)]
  fields <- lapply(meta_rows, geo_meta_values)
  keys <- sub("\t.*$", "", meta_rows)
  valid <- lengths(fields) == length(sample_ids)
  fields <- fields[valid]
  keys <- keys[valid]
  title <- fields[which(keys == "!Sample_title")]
  if (!length(title)) title <- list(rep("", length(sample_ids)))
  out <- data.frame(sample = sample_ids, title = title[[1L]],
                    group = rep(NA_character_, length(sample_ids)),
                    stringsAsFactors = FALSE)
  for (i in seq_along(sample_ids)) {
    label <- geo_group_from_text(out$title[[i]])
    if (is.na(label)) {
      remainder <- vapply(fields[keys != "!Sample_title"],
                          function(z) z[[i]], character(1))
      label <- geo_group_from_text(remainder)
    }
    out$group[[i]] <- label
  }
  out
}

geo_parse_series_matrix <- function(lines, accession, keep_probes, cohort_id = accession) {
  begin <- which(lines == "!series_matrix_table_begin")
  end <- which(lines == "!series_matrix_table_end")
  if (length(begin) != 1L || length(end) != 1L || end <= begin + 2L) {
    stop(sprintf("%s has no usable processed series matrix", accession))
  }
  header <- geo_meta_values(paste0("header\t", lines[[begin + 1L]]))
  # Above keeps the first ID_REF token. Check it rather than guessing column order.
  if (length(header) < 2L || header[[1L]] != "ID_REF") {
    stop(sprintf("%s matrix is not indexed by ID_REF", accession))
  }
  sample_ids <- header[-1L]
  if (!all(grepl("^GSM[0-9]+$", sample_ids))) {
    stop(sprintf("%s matrix sample IDs not recognized", accession))
  }
  metadata <- geo_sample_info(lines[seq_len(begin - 1L)], sample_ids)
  metadata <- geo_restrict_study_groups(metadata, cohort_id)
  known <- which(!is.na(metadata$group))
  if (!length(known)) stop(sprintf("%s tumor/normal class is ambiguous in GEO metadata", accession))
  keep_probes <- unique(as.character(keep_probes))
  if (!length(keep_probes)) stop("No compatible probe IDs available")

  records <- vector("list", 0L)
  probe_lines <- lines[seq.int(begin + 2L, end - 1L)]
  # Avoid parsing all 27K/450K rows into a potentially very large matrix.
  probe_ids <- gsub('^"|"$', "", sub("\t.*$", "", probe_lines))
  selected <- which(probe_ids %in% keep_probes)
  for (i in selected) {
    vals <- strsplit(probe_lines[[i]], "\t", fixed = TRUE)[[1L]]
    if (length(vals) != length(sample_ids) + 1L) next
    beta <- suppressWarnings(as.numeric(gsub('^"|"$', "", vals[-1L])))
    good <- known[is.finite(beta[known]) & beta[known] >= 0 & beta[known] <= 1]
    if (!length(good)) next
    records[[length(records) + 1L]] <- data.frame(
      sample = sample_ids[good], group = metadata$group[good],
      probe = gsub('^"|"$', "", vals[[1L]]), beta = beta[good], cohort = cohort_id,
      stringsAsFactors = FALSE)
  }
  if (!length(records)) stop(sprintf("%s has no usable beta values for mapped probes (check assay scale)", accession))
  out <- do.call(rbind, records)
  # Strict scale guard: if any nonmissing data in the selected probes are
  # outside beta range, do not silently discard them as if they were NAs.
  for (i in selected) {
    vals <- strsplit(probe_lines[[i]], "\t", fixed = TRUE)[[1L]]
    if (length(vals) != length(sample_ids) + 1L) next
    nums <- suppressWarnings(as.numeric(gsub('^"|"$', "", vals[-1L])))
    nums <- nums[is.finite(nums)]
    if (length(nums) && any(nums < 0 | nums > 1)) {
      stop(sprintf("%s contains values outside [0,1]; M-values cannot be treated as beta", accession))
    }
  }
  out <- out[!duplicated(paste(out$sample, out$probe)), , drop = FALSE]
  present <- unique(out[, c("sample", "group")])
  if (!all(c("Tumor", "Normal") %in% present$group)) {
    stop(sprintf("%s lacks confidently labeled tumor and normal beta measurements", accession))
  }
  list(cohort = cohort_id, measurements = out, available_samples = metadata,
       platform = if (exists("GEO_COHORT_CATALOG")) {
         i <- match(cohort_id, GEO_COHORT_CATALOG$cohort)
         if (is.na(i)) "unknown" else GEO_COHORT_CATALOG$platform[[i]]
       } else "unknown", source = paste("NCBI GEO", accession))
}

# A mixed-organ series must be filtered before it can enter an atlas.
# Ambiguous or disease-unspecific samples are excluded, never guessed.
geo_restrict_study_groups <- function(metadata, cohort_id) {
  id <- as.character(cohort_id)
  title <- tolower(trimws(metadata$title))
  # Explicit public GSM title conventions, not inferred from study-level notes.
  if (grepl("^GSE53051_", id)) {
    organ <- switch(sub("^GSE53051_", "", id),
                    BRCA = "breast", COAD = "colon", THCA = "thyroid", PAAD = "pancreas", "")
    if (!nzchar(organ)) stop("Unknown GSE53051 disease subset")
    metadata$group[] <- NA_character_
    metadata$group[grepl(paste0("^", organ, "_normal_"), title)] <- "Normal"
    tumor_pattern <- if (organ == "thyroid") "^thyroid_ptc_" else paste0("^", organ, "_cancer_")
    metadata$group[grepl(tumor_pattern, title)] <- "Tumor"
  }
  if (grepl("^GSE52955_", id)) {
    organ <- sub("^GSE52955_", "", id)
    # GEO sample titles are NK/RT for kidney, NB/TB for bladder, CP/P for
    # prostate. These are documented GSM titles, not missing organ text.
    metadata$group[] <- NA_character_
    if (organ == "KIRC") {
      metadata$group[grepl("^nk[0-9]+$", title)] <- "Normal"
      metadata$group[grepl("^rt[0-9]+$", title)] <- "Tumor"
    } else if (organ == "BLCA") {
      metadata$group[grepl("^nb[0-9]+$", title)] <- "Normal"
      metadata$group[grepl("^tb[0-9]+$", title)] <- "Tumor"
    } else if (organ == "PRAD") {
      metadata$group[grepl("^cp[0-9]+$", title)] <- "Normal"
      metadata$group[grepl("^p[0-9]+$", title)] <- "Tumor"
    } else stop("Unknown GSE52955 disease subset")
  }
  if (identical(id, "GSE121849")) {
    metadata$group[] <- NA_character_
    metadata$group[grepl("from scc_[0-9]+", title)] <- "Tumor"
    metadata$group[grepl("from lung ", title, fixed = TRUE)] <- "Normal"
  }
  if (identical(id, "GSE61441")) {
    metadata$group[] <- NA_character_
    metadata$group[grepl("ccrcc[0-9]+t$", title)] <- "Tumor"
    metadata$group[grepl("ccrcc[0-9]+n$", title)] <- "Normal"
  }
  if (identical(id, "GSE48684")) {
    metadata$group[] <- NA_character_
    metadata$group[grepl("normal individual", title)] <- "Normal"
    metadata$group[grepl("crc patients", title)] <- "Tumor"
  }
  if (identical(id, "GSE121930")) {
    metadata$group[] <- NA_character_
    metadata$group[grepl("normal esophagus tissues", title)] <- "Normal"
    metadata$group[grepl("escc tissues", title)] <- "Tumor"
  }
  if (identical(id, "GSE87053")) {
    metadata$group[] <- NA_character_
    metadata$group[grepl("adjacent normal", title)] <- "Normal"
    metadata$group[grepl("oscc lesion", title)] <- "Tumor"
  }
  if (identical(id, "GSE56044_LUSC")) {
    # Mixed LUAD/LUSC lung study. Only keep cases explicitly labeled LUSC.
    metadata$group[!is.na(metadata$group) & metadata$group == "Tumor" &
                     !grepl("squamous|lusc", title)] <- NA_character_
  }
  metadata
}

# Parse gzipped public GEO Series Matrix directly as a stream.
# Only probe rows matching the current (validated) reference map are parsed.
geo_parse_series_connection <- function(con, accession, keep_probes, cohort_id = accession,
                                        chunk_size = 500L, max_lines = 650000L) {
  keep_probes <- unique(as.character(keep_probes))
  if (!length(keep_probes)) stop("No compatible CpG probe IDs")
  metadata_lines <- character()
  fields <- NULL
  meta <- NULL
  out <- list()
  n_lines <- 0L
  in_table <- FALSE
  finished <- FALSE
  repeat {
    block <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(block)) break
    for (line in block) {
      n_lines <- n_lines + 1L
      if (n_lines > max_lines) stop("GEO matrix exceeded safe line limit")
      if (identical(line, "!series_matrix_table_end")) { finished <- TRUE; break }
      if (!in_table) {
        if (identical(line, "!series_matrix_table_begin")) { in_table <- TRUE; next }
        if (grepl("^!Sample_(title|source_name_ch1|characteristics_ch1|description)\t", line))
          metadata_lines <- c(metadata_lines, line)
        next
      }
      if (is.null(fields)) {
        fields <- geo_meta_values(paste0("header\t", line))
        if (length(fields) < 3L || fields[[1L]] != "ID_REF" ||
            !all(grepl("^GSM[0-9]+$", fields[-1L])))
          stop("GEO matrix has no recognized ID_REF/sample columns")
        meta <- geo_restrict_study_groups(geo_sample_info(metadata_lines, fields[-1L]), cohort_id)
        next
      }
      probe <- gsub('^"|"$', "", sub("\t.*$", "", line))
      if (!probe %in% keep_probes) next
      values <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(values) != length(fields)) next
      nums <- suppressWarnings(as.numeric(gsub('^"|"$', "", values[-1L])))
      if (any(is.finite(nums) & (nums < 0 | nums > 1)))
        stop(sprintf("%s has non-beta values (outside [0,1])", accession))
      good <- which(!is.na(meta$group) & is.finite(nums))
      if (!length(good)) next
      out[[length(out) + 1L]] <- data.frame(
        sample = meta$sample[good], group = meta$group[good], probe = probe,
        beta = nums[good], cohort = cohort_id, stringsAsFactors = FALSE)
    }
    if (finished) break
  }
  if (!finished) stop(sprintf("%s has no complete beta series-matrix data table", accession))
  if (!length(out)) stop(sprintf("%s: no usable beta records or confidently labeled tissues", accession))
  measurements <- do.call(rbind, out)
  measurements <- measurements[!duplicated(paste(measurements$sample, measurements$probe)), , drop = FALSE]
  if (!all(c("Tumor", "Normal") %in% measurements$group))
    stop(sprintf("%s has no verified tumor-and-normal measurements", cohort_id))
  platform <- if (exists("GEO_COHORT_CATALOG")) {
    i <- match(cohort_id, GEO_COHORT_CATALOG$cohort)
    if (!is.na(i)) GEO_COHORT_CATALOG$platform[[i]] else "unknown"
  } else "unknown"
  list(cohort = cohort_id, measurements = measurements, available_samples = meta,
       platform = platform, source = paste("GEO", accession, "Series Matrix"))
}

atlas_fetch_geo <- function(accession, probe_map, read_lines = NULL) {
  study <- sub("_.*$", "", accession)
  if (!is.null(read_lines))
    return(geo_parse_series_matrix(read_lines(study), study, probe_map$probe, cohort_id = accession))
  # No read-all into memory. No fallback to fabricated beta or to intensity/M values.
  h <- curl::new_handle(timeout = 240L, connecttimeout = 15L,
                        useragent = "gene-fractal-atlas/0.13.0")
  con <- gzcon(curl::curl(geo_matrix_url(study), open = "rb", handle = h))
  on.exit(close(con), add = TRUE)
  geo_parse_series_connection(con, study, probe_map$probe, cohort_id = accession)
}

# One row per sample and measured probe, never per simulated CpG or pseudo-patient.
atlas_summarize <- function(x, probe_map, threshold = 0.6) {
  if (!is.list(x) || !is.data.frame(x$measurements)) stop("Invalid cohort measurements")
  m <- x$measurements
  m <- m[is.finite(m$beta) & m$beta >= 0 & m$beta <= 1 &
         m$group %in% c("Tumor", "Normal") & m$probe %in% probe_map$probe, , drop = FALSE]
  if (!nrow(m)) stop("No compatible measured probes")
  map <- unique(probe_map[, c("probe", "position"), drop = FALSE])
  m <- merge(m, map, by = "probe")
  keys <- split(seq_len(nrow(m)), paste(m$probe, m$group, sep = "|"))
  stats <- do.call(rbind, lapply(keys, function(ix) {
    d <- m[ix, , drop = FALSE]
    data.frame(probe = d$probe[[1L]], position = d$position[[1L]],
               group = d$group[[1L]], mean_beta = mean(d$beta),
               prevalence = mean(d$beta >= threshold),
               n = length(unique(d$sample)), stringsAsFactors = FALSE)
  }))
  rownames(stats) <- NULL
  attr(stats, "cohort") <- x$cohort
  stats
}

atlas_sample_counts <- function(x) {
  m <- x$measurements
  data.frame(cohort = x$cohort,
             tumor = length(unique(m$sample[m$group == "Tumor"])),
             normal = length(unique(m$sample[m$group == "Normal"])),
             n_probes = length(unique(m$probe)),
             platform = x$platform,
             stringsAsFactors = FALSE)
}

atlas_probe_deltas <- function(summary) {
  if (is.null(summary) || !nrow(summary)) return(data.frame())
  t <- summary[summary$group == "Tumor", c("probe", "position", "mean_beta", "prevalence", "n")]
  n <- summary[summary$group == "Normal", c("probe", "position", "mean_beta", "prevalence", "n")]
  if (!nrow(t) || !nrow(n)) return(data.frame())
  out <- merge(t, n, by = c("probe", "position"), suffixes = c("_tumor", "_normal"))
  if (!nrow(out)) return(data.frame())
  out$delta <- out$mean_beta_tumor - out$mean_beta_normal
  out$prev_delta <- out$prevalence_tumor - out$prevalence_normal
  out
}

# Mean-beta range is evaluated against the values displayed in each group
# independently. A normal measurement never acts as a prerequisite for showing
# a tumor methylation mark (or vice versa). No original beta is modified.
atlas_validate_beta_range <- function(range) {
  if (!is.numeric(range) || length(range) != 2L ||
      any(!is.finite(range)) || any(range < 0 | range > 1) ||
      range[[1L]] > range[[2L]])
    stop("Mean beta range must be within [0,1] (minimum <= maximum)")
  invisible(TRUE)
}

atlas_filter_beta_range <- function(summary, range = c(0, 1), min_samples = 2L) {
  atlas_validate_beta_range(range)
  if (length(min_samples) != 1L || !is.finite(min_samples) || min_samples < 1)
    stop("Minimum specimens must be >= 1")
  if (!is.data.frame(summary) ||
      !all(c("probe", "group", "mean_beta", "n") %in% names(summary)))
    stop("Invalid CpG summary")
  summary[summary$group %in% c("Tumor", "Normal") &
            is.finite(summary$mean_beta) & is.finite(summary$n) &
            summary$n >= min_samples &
            summary$mean_beta >= range[[1L]] &
            summary$mean_beta <= range[[2L]], , drop = FALSE]
}

# For differential analysis the same probe must have measurements in BOTH
# groups. Each pair is eligible when EITHER group's mean beta is in range.
# Keep both original rows so delta is never computed from two different probes.
atlas_paired_summary <- function(summary, min_samples = 2L) {
  if (!is.data.frame(summary) ||
      !all(c("probe", "group", "mean_beta", "n") %in% names(summary)))
    stop("Invalid CpG summary")
  eligible <- summary[summary$group %in% c("Tumor", "Normal") &
                        is.finite(summary$mean_beta) & is.finite(summary$n) &
                        summary$n >= min_samples, , drop = FALSE]
  if (!nrow(eligible)) return(eligible)
  t <- eligible$probe[eligible$group == "Tumor"]
  n <- eligible$probe[eligible$group == "Normal"]
  eligible[eligible$probe %in% intersect(t, n), , drop = FALSE]
}

# Descriptive tumor-minus-normal mean beta for the SAME probe. By default
# beta-range eligibility uses either measured group's mean (logical OR).
# The delta cutoff itself remains a separate effect-size filter.
atlas_filter_delta_probes <- function(summary, min_samples = 2L,
                                      cutoff = 0.15, beta_range = c(0, 1)) {
  atlas_validate_beta_range(beta_range)
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff < 0 || cutoff > 1)
    stop("Delta beta cutoff must be in [0,1]")
  if (length(min_samples) != 1L || !is.finite(min_samples) || min_samples < 1)
    stop("Minimum specimens must be >= 1")
  d <- atlas_probe_deltas(summary)
  if (!nrow(d)) return(d)
  in_range <- (is.finite(d$mean_beta_tumor) &
               d$mean_beta_tumor >= beta_range[[1L]] &
               d$mean_beta_tumor <= beta_range[[2L]]) |
              (is.finite(d$mean_beta_normal) &
               d$mean_beta_normal >= beta_range[[1L]] &
               d$mean_beta_normal <= beta_range[[2L]])
  d[is.finite(d$delta) & is.finite(d$n_tumor) & is.finite(d$n_normal) &
      d$n_tumor >= min_samples & d$n_normal >= min_samples &
      in_range & abs(d$delta) >= cutoff, , drop = FALSE]
}

# Keep both groups when either group qualifies; useful for exact-probe
# agreement tables without making the unselected group appear in its own panel.
atlas_filter_beta_pairs <- function(summary, range = c(0, 1), min_samples = 2L) {
  d <- atlas_filter_delta_probes(summary, min_samples, cutoff = 0,
                                 beta_range = range)
  if (!nrow(d)) return(summary[FALSE, , drop = FALSE])
  paired <- atlas_paired_summary(summary, min_samples)
  paired[paired$probe %in% d$probe, , drop = FALSE]
}

# Slightly translucent low-beta marks remain visible at 4-column size while
# high-beta marks have a much stronger red core. Delta uses |effect size|.
atlas_signal_opacity <- function(value, delta = FALSE) {
  strength <- if (isTRUE(delta)) abs(value) / 0.6 else value
  0.52 + 0.46 * pmax(0, pmin(1, strength))
}

# Spatial binning in CGR space. Each probe contributes one mean/prevalence,
# not one dot per sample. Overlapping/similar k-mer histories share cells.
atlas_spatial_cells <- function(probes, ref_coordinates, bins = 28L,
                                value_col = "mean_beta", min_samples = 2L) {
  if (!nrow(probes)) return(data.frame())
  cpg <- ref_coordinates[ref_coordinates$is_cpg & is.finite(ref_coordinates$x) &
                           is.finite(ref_coordinates$y), c("cpg_position", "x", "y"), drop = FALSE]
  x <- merge(probes, cpg, by.x = "position", by.y = "cpg_position")
  if (!nrow(x)) return(data.frame())
  if (!value_col %in% names(x)) stop("Invalid spatial score")
  if ("n" %in% names(x)) x <- x[x$n >= min_samples, , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  x <- x[is.finite(x[[value_col]]), , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  bins <- as.integer(bins)
  x$ix <- pmax(0L, pmin(bins - 1L, floor(x$x * bins)))
  x$iy <- pmax(0L, pmin(bins - 1L, floor(x$y * bins)))
  ix <- split(seq_len(nrow(x)), paste(x$ix, x$iy, sep = ":"))
  ans <- do.call(rbind, lapply(ix, function(k) {
    d <- x[k, , drop = FALSE]
    data.frame(ix = d$ix[[1L]], iy = d$iy[[1L]],
               x = (d$ix[[1L]] + .5) / bins,
               y = (d$iy[[1L]] + .5) / bins,
               value = mean(d[[value_col]]),
               n_probes = length(unique(d$position)), stringsAsFactors = FALSE)
  }))
  rownames(ans) <- NULL
  ans
}

# A consensus needs the *same probe* assayed in >=2 independently labeled
# cohorts, not merely two colored dots in neighboring CGR bins.
atlas_probe_agreement <- function(summary_list, positions, min_cohorts = 2L) {
  rows <- lapply(names(summary_list), function(id) {
    z <- atlas_probe_deltas(summary_list[[id]])
    if (!nrow(z)) return(NULL)
    z$cohort <- id
    z[, c("cohort", "probe", "position", "delta"), drop = FALSE]
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  df <- do.call(rbind, rows)
  spl <- split(seq_len(nrow(df)), df$probe)
  result <- do.call(rbind, lapply(spl, function(ix) {
    d <- df[ix, , drop = FALSE]
    if (length(unique(d$cohort)) < min_cohorts) return(NULL)
    delta <- d$delta[is.finite(d$delta)]
    if (length(delta) < min_cohorts) return(NULL)
    direction <- sign(delta)
    concordance <- max(sum(direction > 0), sum(direction < 0)) / length(delta)
    data.frame(probe = d$probe[[1L]], position = d$position[[1L]],
               delta = stats::median(delta), cohorts = length(unique(d$cohort)),
               concordance = concordance, stringsAsFactors = FALSE)
  }))
  if (is.null(result)) return(data.frame())
  rownames(result) <- NULL
  result
}

# An *illustrative expected movement law*: bulk beta is used as the weight
# between C and mC vertices at assayed CpGs. It is NOT a consensus DNA string,
# a single-molecule methylation path or an empirical average of individual CGRs.
atlas_expected_cgr <- function(record, beta_by_position) {
  if (!nchar(record$sequence)) stop("Empty reference sequence")
  chars <- strsplit(toupper(record$sequence), "", fixed = TRUE)[[1L]]
  n <- length(chars)
  map <- vertex_map("mc_triangle")
  cpg_ref <- make_cgr(record, "mc_triangle")
  v <- map$vertices
  idx <- match(chars, rownames(v))
  if (anyNA(idx)) stop("Expected path requires unambiguous A/C/G/T reference sequence")
  vx <- v[idx, 1L]
  vy <- v[idx, 2L]
  beta <- rep(0, n)
  sites <- which(cpg_ref$is_cpg)
  if (length(sites) && length(beta_by_position)) {
    observed <- beta_by_position[as.character(cpg_ref$cpg_position[sites])]
    found <- is.finite(observed)
    beta[sites[found]] <- observed[found]
  }
  if (any(beta > 0)) {
    vx <- vx + beta * (v["mC", 1L] - v["C", 1L])
    vy <- vy + beta * (v["mC", 2L] - v["C", 2L])
  }
  cpg_ref$x <- as.numeric(stats::filter(vx / 2, filter = 0.5,
                                        method = "recursive", init = map$start[[1L]]))
  cpg_ref$y <- as.numeric(stats::filter(vy / 2, filter = 0.5,
                                        method = "recursive", init = map$start[[2L]]))
  cpg_ref
}

# Draw one group-level landscape. No fake patients or extra CpG measurements.
atlas_svg <- function(reference_coords, probe_stats, geometry,
                      score = "prevalence", delta = FALSE,
                      min_samples = 2L, bins = 28L,
                      show_background = TRUE, agreement = FALSE,
                      min_bin_probes = 1L,
                      click_cohort = "", click_group = "") {
  ref <- reference_coords
  border <- vertex_map(geometry)
  xy <- svg_coords(border$border_x, border$border_y, geometry)
  polygon <- paste(sprintf("%.1f,%.1f", xy$x, xy$y), collapse = " ")
  parts <- c('<svg class="fractal-svg" viewBox="0 0 1000 1000" preserveAspectRatio="xMidYMid meet" role="img" aria-label="CGR cohort methylation landscape">',
             '<rect width="1000" height="1000" fill="#000000"/>',
             sprintf('<polyline points="%s" fill="none" stroke="#39541a" stroke-width="2"/>', polygon))
  # Every A/T/C/G coordinate remains on the same reference path in every
  # cohort panel. Old code rendered only 3400 faint points, visually hiding
  # the sequence behind red methylation cells. This is display-only sampling.
  if (isTRUE(show_background)) {
    stride <- max(1L, ceiling(nrow(ref) / 16000L))
    base <- ref[seq.int(1L, nrow(ref), by = stride), , drop = FALSE]
    colors <- c(A = "#baff48", T = "#e1ff98", C = "#74e660", G = "#a2f779")
    for (nucleotide in names(colors)) {
      d <- base[base$base == nucleotide, , drop = FALSE]
      if (nrow(d)) parts <- c(parts,
        svg_path_points(d, geometry, colors[[nucleotide]], width = 4.4, opacity = .77))
    }
  }
  if (nrow(probe_stats)) {
    value_column <- if (agreement) {
      "delta"
    } else if (delta) {
      if (identical(score, "prevalence")) "prev_delta" else "delta"
    } else {
      score
    }
    cells <- atlas_spatial_cells(probe_stats, ref, bins = bins,
                                 value_col = value_column,
                                 min_samples = min_samples)
    if (nrow(cells)) cells <- cells[cells$n_probes >= min_bin_probes, , drop = FALSE]
    if (nrow(cells)) {
      center <- svg_coords(cells$x, cells$y, geometry)
      marks <- vapply(seq_len(nrow(cells)), function(i) {
        z <- cells$value[[i]]
        if (delta || agreement) {
          color <- if (z >= 0) "#ff3948" else "#44bfff"
          opacity <- atlas_signal_opacity(z, delta = TRUE)
        } else {
          color <- "#ff3545"
          opacity <- atlas_signal_opacity(z)
        }
        radius <- 8.5 + min(7.5, 2 * sqrt(cells$n_probes[[i]]))
        # Compact two-layer marker: dim halo behind an opaque, saturated core.
        # Do not change the measured beta values or spatial-cell calculations.
        core <- sprintf(paste0(
          '<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" fill-opacity="%.3f" pointer-events="none"/>',
          '<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" fill-opacity="%.3f" ',
          'stroke="#300006" stroke-width=".9" pointer-events="none">',
          '<title>%d probes | value %.3f</title></circle>'),
          center$x[[i]], center$y[[i]], radius + 6, color, .08 + .16 * opacity,
          center$x[[i]], center$y[[i]], radius, color, opacity,
          cells$n_probes[[i]], z)
        if (!nzchar(click_cohort) || !nzchar(click_group)) return(core)
        hit <- sprintf('<circle class="atlas-cell-hit" data-bin-ix="%d" data-bin-iy="%d" data-bin-count="%d" data-cohort="%s" data-group="%s" tabindex="0" role="button" aria-label="Inspect spatial cell" cx="%.2f" cy="%.2f" r="%.2f" fill="transparent"/>',
                       cells$ix[[i]], cells$iy[[i]], as.integer(bins),
                       htmltools::htmlEscape(click_cohort, attribute = TRUE),
                       htmltools::htmlEscape(click_group, attribute = TRUE),
                       center$x[[i]], center$y[[i]], max(11, radius + 3))
        paste0(core, hit)
      }, character(1))
      parts <- c(parts, paste(marks, collapse = ""))
    }
  }
  # Draw corner names over the sequence and methylation overlays.
  paste0(paste(c(parts, svg_vertex_labels(geometry)), collapse = ""), "</svg>")
}

# Additional independently processed beta-value cohorts may be loaded from
# user's long-form CSV. A coordinate is optional only when a known validated
# public probe map contains the ID; otherwise assembly/chrom/position are required.
atlas_read_user_csv <- function(path, reference_record, known_map = data.frame(),
                                cancer = NULL) {
  if (is.na(file.info(path)$size) || file.info(path)$size > 100 * 1024^2)
    stop("Cohort CSV exceeds 100 MB")
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- tolower(trimws(names(d)))
  required <- c("cancer", "cohort", "sample", "group", "probe", "beta")
  if (!all(required %in% names(d)))
    stop(sprintf("Missing cohort CSV columns: %s", paste(setdiff(required, names(d)), collapse = ", ")))
  if (!nrow(d)) stop("Empty cohort CSV")
  d$cancer <- toupper(trimws(as.character(d$cancer)))
  if (is.null(cancer) || length(cancer) != 1L || anyNA(d$cancer) ||
      any(d$cancer != cancer)) {
    stop("CSV cancer codes must all match the selected cancer type")
  }
  d$cohort <- gsub("[^A-Za-z0-9_-]", "_", trimws(d$cohort))
  d$sample <- trimws(as.character(d$sample))
  d$probe <- trimws(as.character(d$probe))
  d$group <- tolower(trimws(as.character(d$group)))
  d$group[d$group %in% c("tumour", "tumor", "case")] <- "Tumor"
  d$group[d$group %in% c("normal", "control", "adjacent")] <- "Normal"
  d$beta <- suppressWarnings(as.numeric(d$beta))
  if (anyNA(d[, required]) || any(!nzchar(d$cohort) | !nzchar(d$sample) | !nzchar(d$probe)) ||
      any(!grepl("^cg[0-9]{8}$", d$probe)) ||
      any(!d$group %in% c("Tumor", "Normal")) ||
      any(!is.finite(d$beta) | d$beta < 0 | d$beta > 1)) {
    stop("Invalid cohort/sample/group/probe/beta in uploaded file")
  }
  if (anyDuplicated(paste(d$cohort, d$sample, d$probe, sep = "|")))
    stop("Duplicate cohort/sample/probe rows in uploaded CSV")
  by_sample <- split(d$group, paste(d$cohort, d$sample))
  if (any(vapply(by_sample, function(x) length(unique(x)) != 1L, logical(1))))
    stop("A sample is classified as both Tumor and Normal")

  if (all(c("assembly", "chrom", "position") %in% names(d))) {
    assembly <- trimws(as.character(d$assembly))
    chrom <- sub("^chr", "", trimws(as.character(d$chrom)), ignore.case = TRUE)
    position <- suppressWarnings(as.numeric(d$position))
    if (anyNA(assembly) || any(toupper(assembly) != "GRCH37") ||
        anyNA(chrom) || any(tolower(chrom) != tolower(sub("^chr", "", reference_record$chromosome))) ||
        any(!is.finite(position) | position < 1 | position != floor(position)))
      stop("Uploaded coordinates must be forward CpG C positions on GRCh37 and the reference chromosome")
    ref <- make_cgr(reference_record, "triangle")
    valid_positions <- unique(ref$cpg_position[ref$is_cpg])
    if (any(!position %in% valid_positions))
      stop("Uploaded file contains positions that do not match reference CpG cytosines")
    keyed <- paste(d$probe, position, sep = "|")
    if (any(vapply(split(position, d$probe), function(z) length(unique(z)) > 1L, logical(1))))
      stop("One probe ID has inconsistent genomic positions")
    map <- unique(data.frame(probe = d$probe, position = position))
  } else {
    if (!is.data.frame(known_map) || !all(c("probe", "position") %in% names(known_map)))
      stop("Provide assembly/chrom/position columns when Xena probe map is unavailable")
    if (any(!d$probe %in% known_map$probe))
      stop("Cohort CSV includes unrecognized probe IDs; add validated GRCh37 coordinates")
    map <- known_map[known_map$probe %in% d$probe, c("probe", "position"), drop = FALSE]
  }
  out <- lapply(split(d, d$cohort), function(chunk) {
    id <- paste0("USER_", chunk$cohort[[1L]])
    m <- chunk[, c("sample", "group", "probe", "beta"), drop = FALSE]
    m$cohort <- id
    list(cohort = id, measurements = m,
         available_samples = unique(m[, c("sample", "group")]),
         platform = "uploaded", source = "User CSV (verify independent study provenance)")
  })
  names(out) <- vapply(out, function(x) x$cohort, character(1))
  list(cohorts = out, map = map)
}


# Combine cohort-level probe summaries without pretending pooled mean profiles
# are individual patient sequences. Never mix different cancer types.
# Input: named list of study summaries (each already computed from measured data).
atlas_merge_summaries <- function(summary_list,
                                  mode = c("equal", "weighted", "shared"),
                                  min_samples = 2L) {
  mode <- match.arg(mode)
  if (length(summary_list) < 2L) return(data.frame())
  input_names <- names(summary_list)
  if (is.null(input_names) || anyDuplicated(input_names)) stop("Unique cohort IDs required")
  groups <- lapply(summary_list, function(df) {
    if (!all(c("probe", "position", "group", "mean_beta", "prevalence", "n") %in% names(df))) {
      stop("Incomplete cohort summary")
    }
    d <- df[is.finite(df$mean_beta) & is.finite(df$prevalence) &
              is.finite(df$n) & df$n >= min_samples &
              df$group %in% c("Tumor", "Normal"), , drop = FALSE]
    d
  })
  if (identical(mode, "shared")) {
    # Require both tissue types for each probe in each participating study.
    common <- Reduce(intersect, lapply(groups, function(d) {
      intersect(d$probe[d$group == "Tumor"], d$probe[d$group == "Normal"])
    }))
    groups <- lapply(groups, function(d) d[d$probe %in% common, , drop = FALSE])
  }
  frames <- lapply(seq_along(groups), function(i) {
    d <- groups[[i]]
    d$cohort <- input_names[[i]]
    d
  })
  frames <- Filter(function(z) nrow(z) > 0L, frames)
  if (!length(frames)) return(data.frame())
  all_rows <- do.call(rbind, frames)
  # A probe must not map to two genomic positions across the input cohorts.
  if (any(vapply(split(all_rows$position, all_rows$probe),
                 function(p) length(unique(p)) != 1L, logical(1)))) {
    stop("The same probe has conflicting positions across cohorts")
  }
  key <- paste(all_rows$probe, all_rows$group, sep = "|")
  out <- lapply(split(seq_len(nrow(all_rows)), key), function(ii) {
    d <- all_rows[ii, , drop = FALSE]
    w <- if (mode == "weighted") d$n else rep(1, nrow(d))
    data.frame(probe = d$probe[[1L]], position = d$position[[1L]],
               group = d$group[[1L]],
               mean_beta = stats::weighted.mean(d$mean_beta, w),
               prevalence = stats::weighted.mean(d$prevalence, w),
               n = sum(d$n), n_cohorts = length(unique(d$cohort)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}


# Agreement is based on SAME CpG probe IDs measured in >=2 independent
# cohort labels. The final display groups those validated probes into CGR
# cells; green points remain the common reference sequence.
atlas_spatial_concordance <- function(summary_list, ref_coordinates,
                                      bins = 28L, min_samples = 2L,
                                      delta_cutoff = 0.15, beta_range = c(0, 1)) {
  empty <- list(cells = data.frame(), probes = data.frame(), cohorts = names(summary_list))
  if (length(summary_list) < 2L) return(empty)
  items <- Filter(Negate(is.null), lapply(names(summary_list), function(id) {
    d <- atlas_filter_delta_probes(summary_list[[id]], min_samples, delta_cutoff,
                                  beta_range = beta_range)
    if (!nrow(d)) return(NULL)
    d$cohort <- id
    d[, c("probe", "position", "delta", "cohort"), drop = FALSE]
  }))
  if (length(items) < 2L) return(empty)
  z <- do.call(rbind, items)
  z <- z[!duplicated(paste(z$probe, z$cohort)), , drop = FALSE]
  sample_coverage <- tapply(z$cohort, z$probe, function(x) length(unique(x)))
  shared <- names(sample_coverage)[sample_coverage >= 2L]
  z <- z[z$probe %in% shared, , drop = FALSE]
  if (!nrow(z)) return(empty)
  positions <- ref_coordinates[ref_coordinates$is_cpg &
                               is.finite(ref_coordinates$cpg_position),
                               c("cpg_position", "x", "y"), drop = FALSE]
  positions <- positions[!duplicated(positions$cpg_position), , drop = FALSE]
  z <- merge(z, positions, by.x = "position", by.y = "cpg_position")
  if (!nrow(z)) return(empty)
  bins <- as.integer(bins)
  z$ix <- pmax(0L, pmin(bins - 1L, floor(z$x * bins)))
  z$iy <- pmax(0L, pmin(bins - 1L, floor(z$y * bins)))
  # One datum per cohort per cell: mean delta across exact matching probes.
  groups <- split(seq_len(nrow(z)), paste(z$ix, z$iy, sep = ":"))
  cells <- lapply(groups, function(ix) {
    d <- z[ix, , drop = FALSE]
    by_cohort <- tapply(d$delta, d$cohort, mean)
    if (length(by_cohort) < 2L) return(NULL)
    sign_delta <- sign(by_cohort)
    positive <- sum(sign_delta > 0L)
    negative <- sum(sign_delta < 0L)
    agreeing <- max(positive, negative)
    direction <- if (positive >= negative) "higher" else "lower"
    concordant <- names(by_cohort)[if (direction == "higher") sign_delta > 0L else sign_delta < 0L]
    data.frame(ix = d$ix[[1L]], iy = d$iy[[1L]],
               n_probes = length(unique(d$probe)), n_cohorts = length(by_cohort),
               agreeing = agreeing, fraction = agreeing / length(by_cohort),
               median_delta = stats::median(by_cohort), direction = direction,
               contributors = paste(names(by_cohort), collapse = "|"),
               concordant = paste(concordant, collapse = "|"),
               stringsAsFactors = FALSE)
  })
  cells <- Filter(Negate(is.null), cells)
  if (!length(cells)) return(list(cells = data.frame(), probes = z,
                                 cohorts = names(summary_list)))
  list(cells = do.call(rbind, cells), probes = z, cohorts = names(summary_list))
}

atlas_cohort_colors <- function(ids) {
  palette <- c("#dfff78", "#f8bd53", "#7cd6ff", "#e8a0ff", "#6de2bb",
               "#ff8b98", "#bbb8ff", "#ffd980", "#76e5ed", "#bfffea")
  stats::setNames(rep(palette, length.out = length(ids)), ids)
}

atlas_concordance_svg <- function(reference_coords, geometry, agreement,
                                  bins = 28L, min_bin_probes = 1L) {
  base <- atlas_svg(reference_coords, data.frame(), geometry,
                    show_background = TRUE)
  base <- sub("</svg>$", "", base)
  cells <- agreement$cells
  if (!is.data.frame(cells) || !nrow(cells)) return(paste0(base, "</svg>"))
  cells <- cells[cells$n_probes >= min_bin_probes, , drop = FALSE]
  if (!nrow(cells)) return(paste0(base, "</svg>"))
  # Decimation is display-only; table and CSV retain all actual cells/probes.
  if (nrow(cells) > 1200L)
    cells <- cells[unique(as.integer(round(seq(1, nrow(cells), length.out = 1200L)))), , drop = FALSE]
  colors <- atlas_cohort_colors(agreement$cohorts)
  xy <- svg_coords((cells$ix + .5) / bins, (cells$iy + .5) / bins, geometry)
  parts <- vapply(seq_len(nrow(cells)), function(i) {
    cell <- cells[i, , drop = FALSE]
    agreement_sign <- cell$direction[[1L]]
    color <- if (cell$fraction[[1L]] == 1) {
      if (agreement_sign == "higher") "#ff4859" else "#48c6ff"
    } else "#dba457"
    x <- xy$x[[i]]; y <- xy$y[[i]]
    r <- 8.5 + min(7, sqrt(cell$n_probes[[1L]]) * 2)
    core <- sprintf(paste0(
      '<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" fill-opacity=".22" pointer-events="none"/>',
      '<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" fill-opacity=".99" ',
      'stroke="#280c11" stroke-width=".9" pointer-events="none"/>'),
      x, y, r + 5, color, x, y, r, color)
    ids <- strsplit(cell$contributors[[1L]], "|", fixed = TRUE)[[1L]]
    agree <- strsplit(cell$concordant[[1L]], "|", fixed = TRUE)[[1L]]
    n <- length(ids)
    ticks <- vapply(seq_along(ids), function(k) {
      theta <- (2 * pi * (k - 1L) / n) - pi / 2
      alpha <- if (ids[[k]] %in% agree) 1 else .25
      sprintf('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-opacity="%.2f" stroke-width="2.4" pointer-events="none"/>',
              x + (r + 3) * cos(theta), y + (r + 3) * sin(theta),
              x + (r + 9) * cos(theta), y + (r + 9) * sin(theta),
              colors[[ids[[k]]]], alpha)
    }, character(1))
    hit <- sprintf('<circle class="atlas-cell-hit" data-bin-ix="%d" data-bin-iy="%d" data-bin-count="%d" data-cohort="Agreement" data-group="Difference" tabindex="0" role="button" aria-label="Inspect concordance cell" cx="%.2f" cy="%.2f" r="%.2f" fill="transparent"><title>%d/%d cohorts agree; %d shared CpGs; median delta %.3f</title></circle>',
                   cell$ix[[1L]], cell$iy[[1L]], as.integer(bins), x, y,
                   max(15, r + 10), cell$agreeing[[1L]], cell$n_cohorts[[1L]],
                   cell$n_probes[[1L]], cell$median_delta[[1L]])
    paste0(core, paste(ticks, collapse = ""), hit)
  }, character(1))
  paste0(base, paste(parts, collapse = ""), "</svg>")
}
