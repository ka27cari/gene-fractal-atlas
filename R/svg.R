# Native SVG rendering: no Plotly/WebGL or dynamically bound widget outputs.
# Every panel is a single SVG, embedded directly by Shiny's renderUI.
# A fixed viewBox + uniform coordinate scaling preserves exact 1:1 geometry.

svg_coords <- function(x, y, geometry) {
  center <- if (identical(geometry, "square")) 840 else 500 + 340 * sqrt(3) / 2
  list(x = 160 + 680 * x, y = center - 680 * y)
}

# Vertex labels are deliberately outside the polygon so that even dense
# methylation cells cannot hide them. The triangle shares an A/T vertex;
# only the square has four independent nucleotide corners.
svg_vertex_labels <- function(geometry) {
  geometry <- match.arg(geometry, c("triangle", "square", "mc_triangle"))
  vertices <- vertex_map(geometry)$vertices
  if (identical(geometry, "square")) {
    labels <- c("C", "G", "A", "T")
    positions <- vertices[c("C", "G", "A", "T"), , drop = FALSE]
    dx <- c(-30, 30, -30, 30)
    dy <- c(-29, -29, 48, 48)
  } else {
    # Select the unique physical corners from the actual coordinate map.
    labels <- if (identical(geometry, "mc_triangle"))
      c("A/T/G", "C", "mC") else c("A/T", "C", "G")
    keys <- if (identical(geometry, "mc_triangle"))
      c("A", "C", "mC") else c("A", "C", "G")
    positions <- vertices[keys, , drop = FALSE]
    dx <- c(0, -25, 25)
    dy <- c(-37, 49, 49)
  }
  xy <- svg_coords(positions[, 1L], positions[, 2L], geometry)
  labels_svg <- vapply(seq_along(labels), function(i) {
    sprintf('<text x="%.2f" y="%.2f" text-anchor="middle">%s</text>',
            xy$x[[i]] + dx[[i]], xy$y[[i]] + dy[[i]], labels[[i]])
  }, character(1))
  paste0('<g class="cgr-vertex-labels" fill="#edffc4" font-size="49" ',
         'font-weight="750" font-family="Arial, sans-serif" ',
         'stroke="#000000" stroke-width="7" paint-order="stroke" ',
         'pointer-events="none">', paste(labels_svg, collapse = ""), '</g>')
}

svg_path_points <- function(d, geometry, color, width, opacity) {
  if (!nrow(d)) return("")
  d <- d[is.finite(d$x) & is.finite(d$y), , drop = FALSE]
  if (!nrow(d)) return("")
  xy <- svg_coords(d$x, d$y, geometry)
  # Tiny horizontal lines with round caps render as point-like circles.
  commands <- paste0(sprintf("M%.1f %.1fh0.1", xy$x, xy$y), collapse = "")
  sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="%s" stroke-opacity="%s" stroke-linecap="round" pointer-events="none"/>',
          commands, color, width, opacity)
}

svg_highlights <- function(d, geometry, color, radius = 5.2, tooltip = NULL) {
  if (!nrow(d)) return("")
  xy <- svg_coords(d$x, d$y, geometry)
  if (is.null(tooltip)) tooltip <- rep("", nrow(d))
  paste(vapply(seq_len(nrow(d)), function(i) {
    if (nzchar(tooltip[[i]])) {
      sprintf('<circle cx="%.2f" cy="%.2f" r="%.1f" fill="%s" fill-opacity=".96"><title>%s</title></circle>',
              xy$x[[i]], xy$y[[i]], radius, color,
              htmltools::htmlEscape(as.character(tooltip[[i]]), attribute = FALSE))
    } else {
      sprintf('<circle cx="%.2f" cy="%.2f" r="%.1f" fill="%s" fill-opacity=".96"/>',
              xy$x[[i]], xy$y[[i]], radius, color)
    }
  }, character(1)), collapse = "")
}

# Efficient deterministic glow. Shared halo paths keep thousands of
# methylation sites compact (rather than thousands of SVG filters or animations).
# Only measured beta sites passed in here receive a halo or star glint.
svg_methylation_stars <- function(d, geometry, tooltip = NULL,
                                  glow = TRUE, methyl_threshold = 0.6) {
  if (!nrow(d)) return("")
  valid <- is.finite(d$x) & is.finite(d$y) & is.finite(d$beta)
  if (!is.null(tooltip)) tooltip <- tooltip[valid]
  d <- d[valid, , drop = FALSE]
  if (!nrow(d)) return("")
  if (!isTRUE(glow)) {
    return(svg_highlights(d, geometry, "#ff3030", radius = 7.6, tooltip = tooltip))
  }

  glow_layers <- c(
    svg_path_points(d, geometry, "#ff2138", width = 34, opacity = .16),
    svg_path_points(d, geometry, "#ff4353", width = 18, opacity = .40),
    svg_highlights(d, geometry, "#ff3030", radius = 7.6, tooltip = tooltip),
    svg_path_points(d, geometry, "#fff0c8", width = 2.6, opacity = .92)
  )
  # A few stronger stars get short, nonanimated four-way glints.
  # Stable selection: no random flicker when a slider is moved.
  strong <- d[d$beta >= max(methyl_threshold, 0.8), , drop = FALSE]
  if (nrow(strong)) {
    if (nrow(strong) > 32L) {
      strong <- strong[unique(round(seq(1L, nrow(strong), length.out = 32L))), , drop = FALSE]
    }
    xy <- svg_coords(strong$x, strong$y, geometry)
    path <- paste0(sprintf("M%.2f %.2fv-9m0 18v-9m-9 0h18",
                           xy$x, xy$y), collapse = "")
    glow_layers <- c(glow_layers, sprintf(
      '<path d="%s" stroke="#ffe3cb" stroke-width="1.8" stroke-opacity=".86" fill="none" stroke-linecap="round" pointer-events="none"/>',
      path))
  }
  paste0(glow_layers, collapse = "")
}

# Transparent hit targets sit above the visual glow; click events are handled
# by one delegated browser listener. A CGR position is not a unique genomic
# coordinate, so only verified methylation marks get selection targets.
svg_click_targets <- function(d, geometry, group = "") {
  if (!nrow(d) || !"cpg_position" %in% names(d)) return("")
  valid <- is.finite(d$x) & is.finite(d$y) & is.finite(d$cpg_position)
  d <- d[valid, , drop = FALSE]
  if (!nrow(d)) return("")
  xy <- svg_coords(d$x, d$y, geometry)
  paste(vapply(seq_len(nrow(d)), function(i) {
    sprintf('<circle class="cpg-hit" data-cpg-position="%.0f" data-cpg-group="%s" cx="%.2f" cy="%.2f" r="11" fill="transparent" tabindex="0" role="button" aria-label="Inspect CpG at position %.0f"><title>Inspect CpG %s</title></circle>',
      d$cpg_position[[i]],
      htmltools::htmlEscape(as.character(group), attribute = TRUE),
      xy$x[[i]], xy$y[[i]], d$cpg_position[[i]], d$cpg_position[[i]])
  }, character(1)), collapse = "")
}

svg_panel <- function(item, geometry, show_cpg = TRUE, show_beta = TRUE,
                      methyl_threshold = 0.6, reveal_fraction = 1,
                      delta = NULL, delta_cutoff = 0.2,
                      star_highlights = TRUE) {
  d <- item$coordinates
  n <- nrow(d)
  last <- max(1L, as.integer(ceiling(n * reveal_fraction)))
  d <- d[d$index <= last, , drop = FALSE]
  stride <- max(1L, as.integer(ceiling(nrow(d) / 6500)))
  background <- d[seq.int(1L, nrow(d), by = stride), , drop = FALSE]
  border <- vertex_map(geometry)
  edge <- svg_coords(border$border_x, border$border_y, geometry)
  polygon <- paste(sprintf("%.2f,%.2f", edge$x, edge$y), collapse = " ")

  layers <- c(
    '<svg xmlns="http://www.w3.org/2000/svg" class="fractal-svg" viewBox="0 0 1000 1000" preserveAspectRatio="xMidYMid meet" role="img" aria-label="DNA chaos game representation">',
    '<rect width="1000" height="1000" fill="#000000"/>',
    sprintf('<polyline points="%s" fill="none" stroke="#39541a" stroke-width="2"/>', polygon),
    svg_path_points(background, geometry, "#bfff47", width = 3.6, opacity = .64)
  )

  if (is.null(delta)) {
    if (isTRUE(show_cpg)) {
      cpg <- d[d$is_cpg, , drop = FALSE]
      # Keep the full CpG table for matching/export; limit only drawn
      # motif points in very large windows to avoid huge SVG payloads.
      if (nrow(cpg) > 6500L) {
        cpg <- cpg[seq.int(1L, nrow(cpg), by = ceiling(nrow(cpg) / 6500)), , drop = FALSE]
      }
      layers <- c(layers, svg_path_points(cpg, geometry, "#e5ff9c", width = 6.2, opacity = .8))
    }
    m <- item$methylation$points
    if (isTRUE(show_beta) && nrow(m)) {
      m <- m[m$index <= last & is.finite(m$beta) & m$beta >= methyl_threshold, , drop = FALSE]
      if (nrow(m)) {
        layers <- c(layers, svg_methylation_stars(
          m, geometry,
          tooltip = sprintf("%s:%s | beta %.3f", item$record$chromosome, m$cpg_position, m$beta),
          glow = star_highlights, methyl_threshold = methyl_threshold
        ))
        if (!is.null(item$record$sample_id))
          layers <- c(layers, svg_click_targets(m, geometry, item$record$sample_id))
      }
    }
  } else {
    if (nrow(delta)) {
      cpg <- d[d$is_cpg, , drop = FALSE]
      joined <- merge(cpg, delta, by.x = "cpg_position", by.y = "position")
      for (direction in c("up", "down")) {
        subset <- if (direction == "up") joined[joined$delta >= delta_cutoff, , drop = FALSE] else
          joined[joined$delta <= -delta_cutoff, , drop = FALSE]
        if (!nrow(subset)) next
        layers <- c(layers, svg_highlights(
          subset, geometry, if (direction == "up") "#ff3030" else "#25c8ff",
          radius = 8.2,
          tooltip = sprintf("%s:%s | tumor %.3f | normal %.3f | delta %+.3f",
                            item$record$chromosome, subset$cpg_position,
                            subset$tumor_beta, subset$normal_beta, subset$delta)
        ), svg_click_targets(subset, geometry, "delta"))
      }
    }
  }
  # Always draw labels last. They never intercept clicks on real CpGs.
  paste0(paste(c(layers, svg_vertex_labels(geometry)), collapse = ""), "</svg>")
}
