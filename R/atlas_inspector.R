# Click-to-inspect support for cohort-level CGR cells.
# Uses only already loaded samples / actual measured beta values. No API calls.

atlas_cell_probe_map <- function(probe_map, ref_coordinates, ix, iy, bins) {
  empty <- data.frame(probe = character(), position = numeric(),
                      x = numeric(), y = numeric(), stringsAsFactors = FALSE)
  if (!is.data.frame(probe_map) || !nrow(probe_map) ||
      !is.data.frame(ref_coordinates) || !nrow(ref_coordinates)) return(empty)
  if (!all(c("probe", "position") %in% names(probe_map)) ||
      !all(c("is_cpg", "cpg_position", "x", "y") %in% names(ref_coordinates))) {
    stop("Incomplete reference coordinates or probe map")
  }
  if (length(bins) != 1L || !is.finite(bins) || bins < 1 ||
      length(ix) != 1L || length(iy) != 1L ||
      !is.finite(ix) || !is.finite(iy) || ix < 0 || iy < 0 ||
      ix >= bins || iy >= bins) stop("Invalid CGR spatial cell")
  ref <- ref_coordinates[ref_coordinates$is_cpg &
                           is.finite(ref_coordinates$cpg_position) &
                           is.finite(ref_coordinates$x) &
                           is.finite(ref_coordinates$y),
                         c("cpg_position", "x", "y"), drop = FALSE]
  if (!nrow(ref)) return(empty)
  names(ref)[[1L]] <- "position"
  mapped <- merge(unique(probe_map[, c("probe", "position"), drop = FALSE]),
                  unique(ref), by = "position")
  if (!nrow(mapped)) return(empty)
  bx <- pmax(0L, pmin(as.integer(bins) - 1L, floor(mapped$x * bins)))
  by <- pmax(0L, pmin(as.integer(bins) - 1L, floor(mapped$y * bins)))
  mapped <- mapped[bx == ix & by == iy, c("probe", "position", "x", "y"), drop = FALSE]
  mapped <- mapped[!duplicated(mapped$probe), , drop = FALSE]
  rownames(mapped) <- NULL
  mapped
}

atlas_cell_probe_table <- function(summaries, probe_ids) {
  empty <- data.frame(cohort = character(), probe = character(), position = numeric(),
                      tumor_beta = numeric(), normal_beta = numeric(),
                      delta = numeric(), tumor_n = integer(), normal_n = integer(),
                      stringsAsFactors = FALSE)
  if (!length(summaries) || !length(probe_ids)) return(empty)
  parts <- lapply(names(summaries), function(id) {
    s <- summaries[[id]]
    s <- s[s$probe %in% probe_ids & s$group %in% c("Tumor", "Normal"), , drop = FALSE]
    if (!nrow(s)) return(NULL)
    tumor <- s[s$group == "Tumor", c("probe", "position", "mean_beta", "n"), drop = FALSE]
    normal <- s[s$group == "Normal", c("probe", "position", "mean_beta", "n"), drop = FALSE]
    names(tumor)[3:4] <- c("tumor_beta", "tumor_n")
    names(normal)[3:4] <- c("normal_beta", "normal_n")
    joined <- merge(tumor, normal, by = c("probe", "position"), all = TRUE)
    if (!nrow(joined)) return(NULL)
    joined$cohort <- id
    joined$delta <- joined$tumor_beta - joined$normal_beta
    joined[, names(empty), drop = FALSE]
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(empty)
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

atlas_cell_sample_rows <- function(cohorts, probe_ids) {
  empty <- data.frame(cohort = character(), sample = character(),
                      group = character(), probe = character(), beta = numeric(),
                      stringsAsFactors = FALSE)
  if (!length(cohorts) || !length(probe_ids)) return(empty)
  parts <- lapply(names(cohorts), function(id) {
    m <- cohorts[[id]]$measurements
    if (is.null(m) || !is.data.frame(m) || !nrow(m)) return(NULL)
    m <- m[m$probe %in% probe_ids & is.finite(m$beta),
           c("sample", "group", "probe", "beta"), drop = FALSE]
    if (!nrow(m)) return(NULL)
    m$cohort <- id
    m[, names(empty), drop = FALSE]
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(empty)
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# Magnify the SAME CGR cell across cohorts; dots are positioned by the
# reference k-mer history, while brightness follows the cohort's mean beta.
# A missing probe is omitted, never painted as zero methylation.
atlas_cell_zoom_svg <- function(cell_map, group_summary, picked) {
  open <- paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 240 240" ',
    'role="img" aria-label="Measured CpG beta values within selected CGR cell">',
    '<rect width="240" height="240" fill="#000000"/>',
    '<rect x="14" y="14" width="212" height="212" fill="none" stroke="#35501d" stroke-width="1"/>')
  if (!nrow(cell_map) || !nrow(group_summary)) return(paste0(open, '</svg>'))
  z <- merge(cell_map[, c("probe", "position", "x", "y"), drop = FALSE],
             group_summary[, c("probe", "position", "mean_beta"), drop = FALSE],
             by = c("probe", "position"))
  z <- z[is.finite(z$mean_beta) & z$mean_beta >= 0 & z$mean_beta <= 1, , drop = FALSE]
  if (!nrow(z)) return(paste0(open, '</svg>'))
  # Stable display-only decimation if a selected cell holds many points.
  if (nrow(z) > 2000L)
    z <- z[unique(as.integer(round(seq(1, nrow(z), length.out = 2000L)))), , drop = FALSE]
  px <- 14 + 212 * pmax(0, pmin(1, z$x * picked$bins - picked$ix))
  py <- 226 - 212 * pmax(0, pmin(1, z$y * picked$bins - picked$iy))
  marks <- vapply(seq_len(nrow(z)), function(i) {
    beta <- z$mean_beta[[i]]
    sprintf('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="#ff3c47" fill-opacity="%.3f"><title>%s | mean beta %.3f</title></circle>',
            px[[i]], py[[i]], 2.5 + 3 * beta, atlas_signal_opacity(beta),
            htmltools::htmlEscape(z$probe[[i]], attribute = FALSE), beta)
  }, character(1))
  paste0(open, paste(marks, collapse = ""), '</svg>')
}
