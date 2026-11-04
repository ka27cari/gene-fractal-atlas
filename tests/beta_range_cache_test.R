# Offline regression: per-panel mean-beta filtering, paired delta eligibility,
# contrast and persistent RDS. No synthetic patient values used in the app.
source("R/fractal.R")
source("R/svg.R")
source("R/atlas_cohorts.R")

sample_names <- c("T1", "T2", "N1", "N2")
groups <- c("Tumor", "Tumor", "Normal", "Normal")
probes <- sprintf("cg%08d", 1:3)
positions <- c(1001, 1011, 1021)
map <- data.frame(probe = probes, position = positions)
build <- function(normal, tumor, name) {
  values <- Map(function(t, n) c(t, t, n, n), tumor, normal)
  m <- do.call(rbind, lapply(seq_along(probes), function(i) {
    data.frame(sample = sample_names, group = groups,
               probe = probes[[i]], beta = values[[i]], cohort = name)
  }))
  list(cohort = name, measurements = m, platform = "fixture", source = "offline test")
}
a <- build(normal = c(.1, .4, .8), tumor = c(.9, .7, .2), name = "A")
b <- build(normal = c(.4, .7, .8), tumor = c(.8, .8, .1), name = "B")
sa <- atlas_summarize(a, map)
sb <- atlas_summarize(b, map)
stopifnot(nrow(atlas_filter_beta_range(sa, c(0, 1), 2L)) == 6L)
low <- atlas_filter_beta_range(sa, c(0, .2), 2L)
stopifnot(identical(sort(low$probe), sort(probes[c(1, 3)])),
          all(low$group == "Normal" | low$group == "Tumor"),
          identical(low$group[low$probe == probes[[1L]]], "Normal"),
          identical(low$group[low$probe == probes[[3L]]], "Tumor"))
# Tumor high and normal low: both are still compared as a delta.
low_delta <- atlas_filter_delta_probes(sa, 2L, .15, beta_range = c(0, .2))
stopifnot(identical(sort(low_delta$probe), sort(probes[c(1, 3)])),
          all(is.finite(low_delta$mean_beta_tumor)),
          all(is.finite(low_delta$mean_beta_normal)))
# Range where each panel has a different hit: cohort A Tumor cg2, Normal cg2.
mid <- atlas_filter_beta_range(sa, c(.35, .45), 2L)
stopifnot(nrow(mid) == 1L, mid$group[[1L]] == "Normal",
          mid$probe[[1L]] == probes[[2L]])
# An independent range can select a TUMOR probe not eligible in normal.
high <- atlas_filter_beta_range(sa, c(.85, .95), 2L)
stopifnot(nrow(high) == 1L, high$group[[1L]] == "Tumor",
          high$probe[[1L]] == probes[[1L]])
# Delta requires real values in both groups, not merely a qualifying solo group.
few <- sa
few$n[few$probe == probes[[1L]] & few$group == "Normal"] <- 1L
stopifnot(probes[[1L]] %in% atlas_filter_beta_range(few, c(.85, .95), 2L)$probe,
          !probes[[1L]] %in% atlas_filter_delta_probes(
            few, 2L, 0, beta_range = c(.85, .95))$probe)
for (bad in list(c(-.1, 1), c(.5, .1), c(0, 1.1), c(NA_real_, 1))) {
  stopifnot(inherits(try(atlas_filter_beta_range(sa, bad), silent = TRUE), "try-error"),
            inherits(try(atlas_filter_delta_probes(sa, beta_range = bad),
                         silent = TRUE), "try-error"))
}
# Pooled panel selection follows the pooled beta, not pre-selected cohorts.
pooled <- atlas_merge_summaries(list(A = sa, B = sb), mode = "equal")
pooled_vis <- atlas_filter_beta_range(pooled, c(.2, .3), 2L)
stopifnot(nrow(pooled_vis) > 0L,
          all(pooled_vis$mean_beta >= .2 & pooled_vis$mean_beta <= .3))
# Strongly graded red signal, but low beta is not invisible.
stopifnot(atlas_signal_opacity(0) >= .5,
          atlas_signal_opacity(.1) < atlas_signal_opacity(.9),
          atlas_signal_opacity(1) <= 1,
          atlas_signal_opacity(.1, delta = TRUE) <
            atlas_signal_opacity(.5, delta = TRUE))

# Persistent RDS stored outside versioned app directories; no 14-day expiry.
root <- tempfile("gfa-cache-test-")
old <- getOption("gfa.cache_root", NULL)
options(gfa.cache_root = root)
stopifnot(atlas_cached_write(a, "ASXL1", 1000L))
stopifnot(file.exists(atlas_disk_cache("A", "ASXL1", 1000L)))
Sys.setFileTime(atlas_disk_cache("A", "ASXL1", 1000L), Sys.time() - 40 * 86400)
roundtrip <- atlas_cached_read("A", "ASXL1", 1000L)
stopifnot(identical(roundtrip$measurements, a$measurements),
          is.null(atlas_cached_read("A", "TP53", 1000L)))
record <- list(assembly = "GRCh37", sequence = paste(rep("A", 1000L), collapse = ""))
stopifnot(atlas_reference_write("LUAD", "ASXL1", 1000L, record, map))
saved <- atlas_reference_read("LUAD", "ASXL1", 1000L)
stopifnot(identical(saved$map, map), identical(saved$reference, record))
usage <- atlas_cache_usage()
stopifnot(usage$n_files == 2L, usage$bytes > 0)
options(gfa.cache_root = old)
unlink(root, recursive = TRUE)
cat("Per-panel beta range and persistent local cache checks passed\n")
