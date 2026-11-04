# Run from the repository root: Rscript tests/startup_glow_test.R
# Regresses the v0.8.0 draw_base crash before the Shiny server starts.
source("R/fractal.R")
source("R/methylation.R")
source("R/svg.R")

demos <- make_demo_records()
stopifnot(length(demos) == 4L)

# Independent iterative reference: verifies the compiled recurrence as well as
# the previously missing draw_base values for each supported geometry.
reference_coordinates <- function(draw_base, geometry) {
  layout <- vertex_map(geometry)
  x <- y <- rep(NA_real_, length(draw_base))
  px <- layout$start[[1L]]
  py <- layout$start[[2L]]
  for (i in seq_along(draw_base)) {
    base <- draw_base[[i]]
    if (base %in% rownames(layout$vertices)) {
      px <- (px + layout$vertices[base, 1L]) / 2
      py <- (py + layout$vertices[base, 2L]) / 2
      x[[i]] <- px
      y[[i]] <- py
    }
  }
  data.frame(x = x, y = y)
}

for (geometry in c("triangle", "square", "mc_triangle")) {
  for (rec in demos) {
    sites <- methylation_sites_for_record(rec)
    coords <- make_cgr(rec, geometry, methylation_sites = sites,
                       methylation_threshold = 0.6)
    stopifnot(nrow(coords) == nchar(rec$sequence),
              !anyNA(coords$draw_base), all(is.finite(coords$x)),
              all(is.finite(coords$y)))
    ref <- reference_coordinates(coords$draw_base, geometry)
    stopifnot(max(abs(ref$x - coords$x)) < 1e-10,
              max(abs(ref$y - coords$y)) < 1e-10)
  }
}

fragment <- list(sequence = "ACGTCG", region_start = 100L, region_end = 105L,
                 strand = 1L, chromosome = "1", species = "homo_sapiens",
                 assembly = "GRCh38", source = "ensembl")
sites <- data.frame(position = c(101, 104), beta = c(0.95, 0.3))
coords <- make_cgr(fragment, "mc_triangle", methylation_sites = sites,
                   methylation_threshold = 0.6)
stopifnot(identical(coords$draw_base[c(2, 5)], c("mC", "C")),
          identical(which(coords$is_methyl_c), 2L))

# An ambiguous nucleotide must remain a gap, not a fabricated CGR point.
fragment_n <- fragment
fragment_n$sequence <- "ACNTCG"
coords_n <- make_cgr(fragment_n, "square")
stopifnot(is.na(coords_n$x[[3L]]), is.na(coords_n$y[[3L]]),
          identical(coords_n$draw_base[[3L]], "N"))

# The star option affects visuals only: coordinates and methylation values
# must be identical whichever display style is chosen.
rec <- demos[[1L]]
d <- make_cgr(rec, "square", methylation_sites = methylation_sites_for_record(rec))
matched <- match_methylation(rec, d)
stopifnot(matched$matched > 0L)
item <- list(record = rec, coordinates = d, methylation = matched)
shiny_svg <- svg_panel(item, "square", methyl_threshold = 0.6,
                       star_highlights = TRUE)
plain_svg <- svg_panel(item, "square", methyl_threshold = 0.6,
                       star_highlights = FALSE)
stopifnot(grepl('stroke="#ff2138"', shiny_svg, fixed = TRUE),
          grepl('stroke="#fff0c8"', shiny_svg, fixed = TRUE),
          grepl('fill="#ff3030"', shiny_svg, fixed = TRUE),
          !grepl('stroke="#ff2138"', plain_svg, fixed = TRUE),
          grepl('fill="#ff3030"', plain_svg, fixed = TRUE))

item$methylation$points <- data.frame()
no_data_svg <- svg_panel(item, "square", star_highlights = TRUE)
stopifnot(!grepl('stroke="#ff2138"', no_data_svg, fixed = TRUE))

cat("Startup and methylation-star regression checks passed.\n")
