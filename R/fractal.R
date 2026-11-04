# Pure sequence-to-geometry helpers. No external R packages required.

reverse_complement <- function(sequence) {
  sequence <- toupper(gsub("\\s+", "", sequence))
  chars <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  if (!length(chars)) return("")
  complements <- c(A = "T", T = "A", G = "C", C = "G", N = "N")
  result <- unname(complements[rev(chars)])
  result[is.na(result)] <- "N"
  paste0(result, collapse = "")
}

vertex_map <- function(geometry = c("triangle", "square", "mc_triangle")) {
  geometry <- match.arg(geometry)
  if (geometry == "triangle") {
    return(list(
      vertices = rbind(A = c(0.5, sqrt(3) / 2), T = c(0.5, sqrt(3) / 2),
                       C = c(0, 0), G = c(1, 0), mC = c(0, 0)),
      start = c(0.5, sqrt(3) / 6), border_x = c(0, .5, 1, 0),
      border_y = c(0, sqrt(3) / 2, 0, 0), yrange = c(-.04, .96)
    ))
  }
  if (geometry == "square") {
    return(list(
      vertices = rbind(A = c(0, 0), C = c(0, 1), G = c(1, 1), T = c(1, 0), mC = c(0, 1)),
      start = c(0.5, 0.5), border_x = c(0, 0, 1, 1, 0),
      border_y = c(0, 1, 1, 0, 0), yrange = c(-.05, 1.05)
    ))
  }
  list(
    vertices = rbind(A = c(0.5, sqrt(3) / 2), T = c(0.5, sqrt(3) / 2),
                     G = c(0.5, sqrt(3) / 2), C = c(0, 0), mC = c(1, 0)),
    start = c(0.5, sqrt(3) / 6), border_x = c(0, .5, 1, 0),
    border_y = c(0, sqrt(3) / 2, 0, 0), yrange = c(-.06, sqrt(3) / 2 + .06)
  )
}

make_cgr <- function(record, geometry = c("triangle", "square", "mc_triangle"),
                     methylation_sites = data.frame(position = numeric(), beta = numeric()),
                     methylation_threshold = 0.6) {
  geometry <- match.arg(geometry)
  chars <- strsplit(toupper(record$sequence), "", fixed = TRUE)[[1L]]
  n <- length(chars)
  if (n == 0L) stop("Empty DNA sequence")

  layout <- vertex_map(geometry)
  vertices <- layout$vertices
  methyl_positions <- if (nrow(methylation_sites)) methylation_sites$position[methylation_sites$beta >= methylation_threshold] else numeric()

  if (record$strand == 1L) {
    genomic_pos <- record$region_start + seq_len(n) - 1L
  } else {
    genomic_pos <- record$region_end - seq_len(n) + 1L
  }

  is_cpg <- rep(FALSE, n)
  if (n >= 2L) {
    is_cpg[seq_len(n - 1L)] <- chars[seq_len(n - 1L)] == "C" & chars[seq.int(2L, n)] == "G"
  }
  cpg_forward_pos <- genomic_pos - if (record$strand == -1L) 1L else 0L
  is_methyl_c <- is_cpg & cpg_forward_pos %in% methyl_positions

  # Re-label only CpG cytosines with measured beta at/above the threshold.
  # Must exist for all geometries, including the unmodified demo on startup.
  draw_base <- chars
  draw_base[is_methyl_c] <- "mC"

  # For unambiguous DNA the CGR midpoint recurrence can be evaluated in
  # compiled code through stats::filter.  This is much faster than an R loop
  # for 250 kb to 1 Mb windows. Ambiguous nucleotides keep the last state and
  # leave gaps, so they use the precise fallback below.
  if (all(draw_base %in% rownames(vertices))) {
    x <- as.numeric(stats::filter(
      vertices[draw_base, 1L] / 2, filter = 0.5,
      method = "recursive", init = layout$start[[1L]]))
    y <- as.numeric(stats::filter(
      vertices[draw_base, 2L] / 2, filter = 0.5,
      method = "recursive", init = layout$start[[2L]]))
  } else {
    x <- y <- rep(NA_real_, n)
    px <- layout$start[[1L]]
    py <- layout$start[[2L]]
    for (i in seq_len(n)) {
      base <- draw_base[[i]]
      if (base %in% rownames(vertices)) {
        px <- (px + vertices[base, 1L]) / 2
        py <- (py + vertices[base, 2L]) / 2
        x[[i]] <- px
        y[[i]] <- py
      }
    }
  }

  data.frame(
    index = seq_len(n), base = chars, draw_base = draw_base, x = x, y = y,
    genomic_position = genomic_pos, is_cpg = is_cpg, is_methyl_c = is_methyl_c,
    cpg_position = ifelse(is_cpg, cpg_forward_pos, NA_real_),
    stringsAsFactors = FALSE
  )
}

make_demo_records <- function() {
  labels <- c("Synthetic Human", "Synthetic Chimpanzee", "Synthetic Mouse", "Synthetic Zebrafish")
  species <- c("homo_sapiens", "pan_troglodytes", "mus_musculus", "danio_rerio")
  symbols <- c("TP53", "TP53", "Trp53", "tp53")
  assemblies <- c("SYNTHETIC_GRCh38", "SYNTHETIC_PanTro", "SYNTHETIC_GRCm39", "SYNTHETIC_GRCz11")
  weights <- list(c(.30, .20, .21, .29), c(.26, .23, .22, .29),
                  c(.20, .31, .29, .20), c(.28, .24, .24, .24))
  set.seed(31415)
  result <- vector("list", length(labels))
  for (k in seq_along(labels)) {
    sequence <- paste0(sample(c("A", "C", "G", "T"), 12000L, replace = TRUE, prob = weights[[k]]), collapse = "")
    rec <- list(species = species[[k]], species_label = labels[[k]],
                symbol = symbols[[k]], gene_id = paste0("synthetic_", k), assembly = assemblies[[k]],
                chromosome = "demo", region_start = 1L, region_end = 12000L,
                strand = 1L, sequence = sequence, source = "demo", orthology = "synthetic")
    cgr <- make_cgr(rec, "triangle")
    cpg <- cgr[cgr$is_cpg, "cpg_position"]
    keep <- if (length(cpg)) seq.int(1L, length(cpg), by = 3L) else integer()
    if (length(keep)) {
      pos <- cpg[keep]
      rec$embedded_methylation <- data.frame(
        position = pos,
        beta = pmin(.98, pmax(.02, (sin((pos + k * 27) / 89) + 1) / 2))
      )
    } else rec$embedded_methylation <- data.frame(position = numeric(), beta = numeric())
    result[[k]] <- rec
  }
  result
}
