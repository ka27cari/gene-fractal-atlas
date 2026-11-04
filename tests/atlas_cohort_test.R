# Offline contract checks: no real patient values or web calls required.
source("R/fractal.R")
source("R/svg.R")
source("R/atlas_cohorts.R")

record <- list(
  sequence = paste(rep("AACGATCGGCGTACG", 25L), collapse = ""),
  strand = 1L, region_start = 1001L,
  region_end = 1000L + 15L * 25L,
  chromosome = "1", assembly = "GRCh37")
ref <- make_cgr(record, "triangle")
cpg <- which(ref$is_cpg)
stopifnot(length(cpg) > 3L)
map <- data.frame(probe = sprintf("cg%08d", seq_len(3L)),
                  position = ref$cpg_position[cpg[1:3]], stringsAsFactors = FALSE)

sample <- c("T01", "T02", "N01", "N02")
values <- rbind(c(.9, .8, .1, .2), c(.9, .7, .1, .3), c(.2, .2, .6, .6))
m <- do.call(rbind, lapply(seq_len(3L), function(i) {
  data.frame(sample = sample, group = rep(c("Tumor", "Normal"), each = 2L),
             probe = map$probe[[i]], beta = values[i, ], cohort = "fixture-A")
}))
x <- list(cohort = "fixture-A", measurements = m, platform = "27K")
s <- atlas_summarize(x, map, .6)
stopifnot(nrow(s) == 6L,
          all(s$n == 2L),
          all(s$prevalence >= 0 & s$prevalence <= 1))
d <- atlas_probe_deltas(s)
stopifnot(nrow(d) == 3L, d$delta[[1L]] > 0,
          abs(d$delta[[1L]] - .7) < 1e-6)
for (geom in c("triangle", "square", "mc_triangle")) {
  coords <- make_cgr(record, geom)
  cells <- atlas_spatial_cells(s[s$group == "Tumor", ], coords, bins = 24L)
  stopifnot(nrow(cells) > 0L, all(cells$value >= 0 & cells$value <= 1))
  svg <- atlas_svg(coords, s[s$group == "Tumor", ], geom, bins = 24L)
  stopifnot(grepl("<svg", svg, fixed = TRUE), grepl("</svg>", svg, fixed = TRUE))
  delta_svg <- atlas_svg(coords, d, geom, delta = TRUE,
                         score = "prevalence", min_samples = 1L)
  stopifnot(grepl("<svg", delta_svg, fixed = TRUE))
}

s2 <- s
s2$mean_beta <- pmin(1, pmax(0, s2$mean_beta + .03))
agreement <- atlas_probe_agreement(list(`fixture-A` = s, `fixture-B` = s2), map)
stopifnot(nrow(agreement) == 3L, all(agreement$cohorts == 2L),
          all(agreement$concordance == 1))

# Beta-weighted expected movement is a separate model; reference geometry stays fixed.
base <- make_cgr(record, "mc_triangle")
zero <- atlas_expected_cgr(record, stats::setNames(rep(0, 3L), as.character(map$position)))
stopifnot(max(abs(base$x - zero$x), na.rm = TRUE) < 1e-8)
one <- atlas_expected_cgr(record, stats::setNames(rep(1, 3L), as.character(map$position)))
stopifnot(any(abs(one$x - zero$x) > 1e-6))

# Minimal but realistic processed GEO Series Matrix fixture.
fake <- c(
  '!Sample_title\t"Lung adenocarcinoma A"\t"Adjacent non-tumor lung B"\t"Lung adenocarcinoma C"\t"Adjacent non-tumor lung D"',
  '!Sample_source_name_ch1\t"lung"\t"lung"\t"lung"\t"lung"',
  '!series_matrix_table_begin',
  '"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"',
  sprintf('"%s"\t"0.90"\t"0.10"\t"0.80"\t"0.20"', map$probe[[1L]]),
  sprintf('"%s"\t"0.60"\t"0.70"\t"0.40"\t"0.50"', map$probe[[2L]]),
  '!series_matrix_table_end'
)
f <- geo_parse_series_matrix(fake, "GSE32861", map$probe)
stopifnot(nrow(f$measurements) == 8L,
          identical(sort(unique(f$measurements$group)), c("Normal", "Tumor")))
bad <- fake
bad[[5L]] <- sub('"0.90"', '"-1.90"', bad[[5L]], fixed = TRUE)
stopifnot(inherits(try(geo_parse_series_matrix(bad, "GSE32861", map$probe), silent = TRUE), "try-error"))
cat("LUAD cohort atlas offline contracts passed\n")

# Import only real-position-compatible user-supplied cohorts.
csv_file <- tempfile(fileext = ".csv")
local <- data.frame(
  cancer = "LUAD", cohort = "OTHER_LUAD", sample = c("U01", "U02"),
  group = c("Tumor", "Normal"), probe = rep(map$probe[[1L]], 2L),
  beta = c(.81, .23), assembly = rep("GRCh37", 2L),
  chrom = rep("1", 2L), position = rep(map$position[[1L]], 2L))
utils::write.csv(local, csv_file, row.names = FALSE)
u <- atlas_read_user_csv(csv_file, record, map, cancer = "LUAD")
stopifnot(length(u$cohorts) == 1L,
          identical(names(u$cohorts), "USER_OTHER_LUAD"))
local$position <- local$position + 8L
utils::write.csv(local, csv_file, row.names = FALSE)
stopifnot(inherits(try(atlas_read_user_csv(csv_file, record, map, cancer = "LUAD"), silent = TRUE), "try-error"))
unlink(csv_file)
