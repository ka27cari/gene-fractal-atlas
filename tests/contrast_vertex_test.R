# Offline regression for visible signals and correctly labeled CGR corners.
source("R/fractal.R")
source("R/svg.R")
source("R/atlas_cohorts.R")

expect_labels <- list(
  triangle = c("A/T", "C", "G"),
  square = c("C", "G", "A", "T"),
  mc_triangle = c("A/T/G", "C", "mC")
)

for (geom in names(expect_labels)) {
  labels <- svg_vertex_labels(geom)
  stopifnot(grepl('class="cgr-vertex-labels"', labels, fixed = TRUE),
            grepl('pointer-events="none"', labels, fixed = TRUE))
  for (label in expect_labels[[geom]]) {
    stopifnot(grepl(paste0(">", label, "</text>"), labels, fixed = TRUE))
  }
  stopifnot(length(gregexpr('<text ', labels, fixed = TRUE)[[1L]]) == length(expect_labels[[geom]]))
}

rec <- list(sequence = paste(rep("ATCGCGTACG", 50), collapse = ""),
            region_start = 1001L, region_end = 1500L,
            strand = 1L, chromosome = "1", assembly = "GRCh37")
for (geom in c("triangle", "square", "mc_triangle")) {
  coords <- make_cgr(rec, geom)
  picks <- which(coords$is_cpg)[1:3]
  # The low-beta fixture reproduces the old near-invisible dot failure.
  z <- data.frame(probe = sprintf("cg%08d", 1:3),
                  position = coords$cpg_position[picks], group = "Tumor",
                  mean_beta = c(.02, .15, .96), prevalence = c(0, 0, 1),
                  n = c(3L, 3L, 3L), stringsAsFactors = FALSE)
  svg <- atlas_svg(coords, z, geom, score = "mean_beta", min_samples = 1L,
                   click_cohort = "test", click_group = "Tumor")
  # Cores have a clear beta-dependent gradient, with a visible minimum.
  cores <- regmatches(svg, gregexpr('<circle[^>]*stroke="#300006"[^>]*>', svg, perl = TRUE))[[1L]]
  stopifnot(length(cores) >= 1L)
  vals <- as.numeric(sub('.*fill-opacity="([0-9.]+)".*', "\\1", cores))
  stopifnot(all(is.finite(vals)), all(vals >= .52), all(vals <= .98))
  stopifnot(grepl('fill-opacity="0.', svg, fixed = TRUE),
            grepl('class="atlas-cell-hit"', svg, fixed = TRUE),
            grepl('class="cgr-vertex-labels"', svg, fixed = TRUE))
  plain <- svg_panel(list(record = rec, coordinates = coords,
                          methylation = list(points = data.frame())),
                     geom, show_cpg = FALSE, show_beta = FALSE)
  stopifnot(grepl('class="cgr-vertex-labels"', plain, fixed = TRUE))
}
cat("CGR corner labels, graded beta contrast and clickable cells passed.\n")
