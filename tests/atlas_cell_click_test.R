# Regression: cohort-colored spatial cells must remain clickable and inspectable.
# No external data sources or R Shiny browser are required.
source("R/fractal.R")
source("R/svg.R")
source("R/atlas_cohorts.R")
source("R/atlas_inspector.R")

rec <- list(sequence = paste(rep("AACGATCGGCGTACG", 30L), collapse = ""),
            strand = 1L, region_start = 1001L, region_end = 1450L,
            chromosome = "1", assembly = "GRCh37")
coords <- make_cgr(rec, "triangle")
cpg <- which(coords$is_cpg)
stopifnot(length(cpg) >= 4L)
probe_map <- data.frame(probe = sprintf("cg%08d", seq_len(4L)),
                        position = coords$cpg_position[cpg[1:4]],
                        stringsAsFactors = FALSE)

make_cohort <- function(id, shift = 0) {
  samples <- c("T1", "T2", "N1", "N2")
  values <- c(.90, .80, .20, .15)
  m <- do.call(rbind, lapply(seq_len(nrow(probe_map)), function(i) {
    data.frame(sample = samples, group = c("Tumor", "Tumor", "Normal", "Normal"),
               probe = probe_map$probe[[i]], beta = pmax(0, pmin(1, values + shift)),
               cohort = id, stringsAsFactors = FALSE)
  }))
  list(cohort = id, measurements = m, platform = "450K")
}
cohorts <- list(COHORT_A = make_cohort("COHORT_A"),
                COHORT_B = make_cohort("COHORT_B", .03))
smry <- lapply(cohorts, atlas_summarize, probe_map = probe_map, threshold = .6)
for (geom in c("triangle", "square")) {
  ref <- make_cgr(rec, geom)
  t <- smry$COHORT_A[smry$COHORT_A$group == "Tumor", , drop = FALSE]
  svg <- atlas_svg(ref, t, geom, bins = 28L, click_cohort = "COHORT_A",
                   click_group = "Tumor")
  stopifnot(grepl('class="atlas-cell-hit"', svg, fixed = TRUE),
            grepl('data-cohort="COHORT_A"', svg, fixed = TRUE),
            grepl('data-group="Tumor"', svg, fixed = TRUE),
            grepl('data-bin-ix=', svg, fixed = TRUE))
  fallback <- atlas_svg(ref, t, geom, bins = 28L)
  stopifnot(!grepl('class="atlas-cell-hit"', fallback, fixed = TRUE))

  reference_point <- ref[ref$cpg_position == probe_map$position[[1L]] & ref$is_cpg, , drop = FALSE]
  ix <- min(27L, floor(reference_point$x[[1L]] * 28L))
  iy <- min(27L, floor(reference_point$y[[1L]] * 28L))
  local <- atlas_cell_probe_map(probe_map, ref, ix, iy, 28L)
  stopifnot(nrow(local) >= 1L,
            probe_map$probe[[1L]] %in% local$probe)
  tbl <- atlas_cell_probe_table(smry, local$probe)
  stopifnot(nrow(tbl) >= 2L,
            all(tbl$delta[is.finite(tbl$delta)] > 0),
            all(c("COHORT_A", "COHORT_B") %in% tbl$cohort))
  raw <- atlas_cell_sample_rows(cohorts, local$probe)
  stopifnot(all(c("Tumor", "Normal") %in% raw$group),
            all(raw$probe %in% local$probe),
            length(unique(raw$cohort)) == 2L)
  zoom <- atlas_cell_zoom_svg(local,
                             smry$COHORT_A[smry$COHORT_A$group == "Tumor", , drop = FALSE],
                             list(ix = ix, iy = iy, bins = 28L))
  stopifnot(grepl("<svg", zoom, fixed = TRUE),
            grepl("mean beta", zoom, fixed = TRUE))
}
cat("Cohort cell click targets, exact CpG matching and drill-down passed\n")
