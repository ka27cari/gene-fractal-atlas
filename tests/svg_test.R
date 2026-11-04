# The plot renderer is a direct SVG, not a dynamically registered Plotly output.
source("R/fractal.R")
source("R/methylation.R")
source("R/svg.R")

rec <- make_demo_records()[[1L]]
for (geometry in c("triangle", "square", "mc_triangle")) {
  sites <- methylation_sites_for_record(rec)
  coords <- make_cgr(rec, geometry, methylation_sites = sites)
  meth <- match_methylation(rec, coords)
  rendered <- svg_panel(list(record = rec, coordinates = coords, methylation = meth),
                        geometry, reveal_fraction = 0.5)
  stopifnot(startsWith(rendered, "<svg "),
            grepl('viewBox="0 0 1000 1000"', rendered, fixed = TRUE),
            grepl('stroke="#bfff47"', rendered, fixed = TRUE),
            grepl('fill="#ff3030"', rendered, fixed = TRUE),
            grepl('</svg>$', rendered))
}

# 1000px SVG viewBox has exact equal scaling for both axes.
triangle <- vertex_map("triangle")
pts <- svg_coords(triangle$border_x[1:3], triangle$border_y[1:3], "triangle")
distxy <- function(i, j) sqrt((pts$x[[i]] - pts$x[[j]])^2 + (pts$y[[i]] - pts$y[[j]])^2)
stopifnot(max(c(distxy(1,2), distxy(2,3), distxy(3,1))) -
            min(c(distxy(1,2), distxy(2,3), distxy(3,1))) < 1e-8)
cat("Native SVG renderer and equilateral geometry checks passed.\n")
