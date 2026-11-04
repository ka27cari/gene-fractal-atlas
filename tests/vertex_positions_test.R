# Corner labels must follow actual vertex_map coordinates, not mirrored plots.
source("R/fractal.R")
source("R/svg.R")
expected <- list(
  triangle = c("A/T", "C", "G"),
  square = c("C", "G", "A", "T"),
  mc_triangle = c("A/T/G", "C", "mC")
)
for (geometry in names(expected)) {
  svg <- svg_vertex_labels(geometry)
  for (lab in expected[[geometry]]) {
    stopifnot(grepl(sprintf(">%s</text>", lab), svg, fixed = TRUE))
  }
  stopifnot(length(gregexpr("<text ", svg, fixed = TRUE)[[1L]]) ==
              length(expected[[geometry]]))
}
# Square semantics: C (upper left), G (upper right), A (lower left), T (lower right).
map <- vertex_map("square")$vertices
stopifnot(map["C",1] < map["G",1], map["C",2] == map["G",2],
          map["A",1] < map["T",1], map["A",2] == map["T",2])
# Shared vertices in the two triangle modes.
map <- vertex_map("triangle")$vertices
stopifnot(all(map["A",] == map["T",]))
map <- vertex_map("mc_triangle")$vertices
stopifnot(all(map["A",] == map["T",]), all(map["T",] == map["G",]))
cat("Vertex mapping semantics passed.\n")
