# Run from repository root: Rscript tests/compare_test.R
source("R/fractal.R")
source("R/compare.R")

coords <- vertex_map("mc_triangle")$vertices
triangle <- coords[c("A", "C", "mC"), , drop = FALSE]
dist_xy <- function(a, b) sqrt(sum((a - b)^2))
sides <- c(dist_xy(triangle[1, ], triangle[2, ]),
           dist_xy(triangle[2, ], triangle[3, ]),
           dist_xy(triangle[1, ], triangle[3, ]))
stopifnot(max(sides) - min(sides) < 1e-10)

sq <- vertex_map("square")$vertices[c("A", "C", "G", "T"), , drop = FALSE]
stopifnot(all(abs(c(dist_xy(sq[1,],sq[2,]), dist_xy(sq[2,],sq[3,]),
                    dist_xy(sq[3,],sq[4,]), dist_xy(sq[4,],sq[1,])) - 1) < 1e-10))

samples <- data.frame(sample=c("TCGA-AA-0001-01", "TCGA-AA-0001-11",
                               "TCGA-AA-0002-01", "TCGA-AA-0002-11"),
                      group=c("Tumor", "Normal", "Tumor", "Normal"))
methyl <- data.frame(sample = rep(samples$sample, each=2),
                     position = rep(c(100, 104), times=4),
                     beta = c(.8, .1, .2, .4, .9, .3, .3, .2))
result <- tcga_delta(samples, methyl)
stopifnot(result$method == "patient-paired", result$n_pairs == 2L,
          nrow(result$data) == 2L,
          isTRUE(all.equal(result$data$delta, c(.6, -.1), tolerance=1e-12)))

# No matched patient: differences use separate group means, labeled as such.
samples$sample[3:4] <- c("TCGA-AA-0003-01", "TCGA-AA-0004-11")
methyl$sample[5:8] <- rep(samples$sample[3:4], each=2)
res2 <- tcga_delta(samples, methyl)
stopifnot(res2$method == "patient-paired", res2$n_pairs == 1L)
# Single completely unpaired control set.
samples$sample[c(2, 4)] <- c("TCGA-AA-0005-11", "TCGA-AA-0006-11")
methyl$sample[3:4] <- samples$sample[2]
methyl$sample[7:8] <- samples$sample[4]
res3 <- tcga_delta(samples, methyl)
stopifnot(res3$method == "unpaired-group-means", all(is.finite(res3$data$delta)))

st <- sequence_statistics(list(sequence = "CCGGAT"))
stopifnot(isTRUE(all.equal(st$gc, 100 * 4/6)),
          isTRUE(all.equal(st$cpg_oe, 1.5)), st$n_cpg == 1L)
cat("Geometry and tumor/normal comparison tests passed.\n")
