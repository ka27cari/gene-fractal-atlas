# Descriptive CpG beta differences between public tumor and solid-tissue-normal
# samples. Never infer missing beta values. Never interpret this as a test of
# differential methylation significance.

tcga_delta <- function(samples, methylation) {
  if (!all(c("sample", "group") %in% names(samples)) ||
      !all(c("sample", "position", "beta") %in% names(methylation))) {
    stop("Invalid sample metadata or methylation matrix")
  }
  if (!all(c("Tumor", "Normal") %in% samples$group)) {
    return(list(data = data.frame(), method = "unavailable", n_pairs = 0L))
  }
  joined <- merge(methylation[, c("sample", "position", "beta")],
                  samples[, c("sample", "group")], by = "sample")
  joined <- joined[is.finite(joined$position) & is.finite(joined$beta) &
                     joined$beta >= 0 & joined$beta <= 1, , drop = FALSE]
  if (!nrow(joined)) {
    return(list(data = data.frame(), method = "unavailable", n_pairs = 0L))
  }
  joined$patient <- substr(joined$sample, 1L, 12L)
  tumor_patients <- unique(joined$patient[joined$group == "Tumor"])
  normal_patients <- unique(joined$patient[joined$group == "Normal"])
  patients <- intersect(tumor_patients, normal_patients)
  if (length(patients)) {
    # Matched samples are compared within each patient, then averaged by locus.
    pairs <- lapply(patients, function(patient) {
      x <- joined[joined$patient == patient, , drop = FALSE]
      t <- x[x$group == "Tumor", c("position", "beta"), drop = FALSE]
      n <- x[x$group == "Normal", c("position", "beta"), drop = FALSE]
      # Keep one observation per position per tissue to avoid join expansion.
      t <- t[!duplicated(t$position), , drop = FALSE]
      n <- n[!duplicated(n$position), , drop = FALSE]
      out <- merge(t, n, by = "position", suffixes = c("_tumor", "_normal"))
      if (!nrow(out)) return(NULL)
      out$patient <- patient
      out
    })
    pairs <- Filter(Negate(is.null), pairs)
    if (length(pairs)) {
      both <- do.call(rbind, pairs)
      positions <- sort(unique(both$position))
      out <- do.call(rbind, lapply(positions, function(position) {
        d <- both[both$position == position, , drop = FALSE]
        data.frame(position = position,
                   tumor_beta = mean(d$beta_tumor),
                   normal_beta = mean(d$beta_normal),
                   delta = mean(d$beta_tumor - d$beta_normal),
                   n_pairs = nrow(d))
      }))
      rownames(out) <- NULL
      return(list(data = out, method = "patient-paired", n_pairs = length(pairs)))
    }
  }

  # Unpaired samples: report a descriptive mean difference, not a paired delta.
  t <- joined[joined$group == "Tumor", c("position", "beta"), drop = FALSE]
  n <- joined[joined$group == "Normal", c("position", "beta"), drop = FALSE]
  if (!nrow(t) || !nrow(n)) {
    return(list(data = data.frame(), method = "unavailable", n_pairs = 0L))
  }
  t <- stats::aggregate(beta ~ position, data = t, FUN = mean)
  n <- stats::aggregate(beta ~ position, data = n, FUN = mean)
  names(t)[2L] <- "tumor_beta"
  names(n)[2L] <- "normal_beta"
  out <- merge(t, n, by = "position")
  out$delta <- out$tumor_beta - out$normal_beta
  out$n_pairs <- 0L
  list(data = out, method = "unpaired-group-means", n_pairs = 0L)
}

sequence_statistics <- function(record) {
  bases <- strsplit(toupper(record$sequence), "", fixed = TRUE)[[1L]]
  valid <- bases %in% c("A", "C", "G", "T")
  n <- sum(valid)
  if (!n) return(list(gc = NA_real_, cpg_oe = NA_real_, n_cpg = 0L))
  n_c <- sum(bases == "C")
  n_g <- sum(bases == "G")
  cpg <- if (length(bases) < 2L) 0L else
    sum(head(bases, -1L) == "C" & tail(bases, -1L) == "G")
  expected <- n_c * n_g / n
  list(gc = 100 * (n_c + n_g) / n,
       cpg_oe = if (expected > 0) cpg / expected else NA_real_,
       n_cpg = cpg)
}
