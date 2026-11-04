# Offline cohort composition: cancer filters and pooling on measured data only.
source("R/fractal.R")
source("R/svg.R")
source("R/xena.R")
source("R/atlas_cohorts.R")
source("R/data_catalog.R")

stopifnot(all(c("LUAD", "BRCA", "LUSC", "PAAD", "CHOL") %in% unname(TCGA_COHORTS)))
stopifnot("GSE32861" %in% unname(catalog_cohorts("LUAD")))
stopifnot(!"GSE32861" %in% unname(catalog_cohorts("BRCA")))
stopifnot(!"TCGA_LUAD" %in% unname(catalog_cohorts("LUSC")))

# Both studies share one probe; study A has 2 measured samples/group,
# study B has 8. Equal means should differ from specimen-weighted means.
a <- data.frame(probe = rep("cg00000001", 2), position = 103,
                group = c("Tumor", "Normal"),
                mean_beta = c(.9, .2), prevalence = c(1, 0), n = c(2L, 2L))
b <- data.frame(probe = rep("cg00000001", 2), position = 103,
                group = c("Tumor", "Normal"),
                mean_beta = c(.3, .4), prevalence = c(.25, .5), n = c(8L, 8L))
inputs <- list(`TCGA_LUAD` = a, `GSE32861` = b)
equal <- atlas_merge_summaries(inputs, mode = "equal")
weighted <- atlas_merge_summaries(inputs, mode = "weighted")
shared <- atlas_merge_summaries(inputs, mode = "shared")
stopifnot(nrow(equal) == 2L, nrow(weighted) == 2L, nrow(shared) == 2L)
t <- equal[equal$group == "Tumor", ]
tw <- weighted[weighted$group == "Tumor", ]
stopifnot(abs(t$mean_beta - .6) < 1e-12, abs(tw$mean_beta - .42) < 1e-12)
stopifnot(identical(t$n_cohorts, 2L), tw$n == 10L)
# Shared mode: probe not measured in both cohorts must be discarded.
a2 <- rbind(a, data.frame(probe = "cg00000002", position = 108,
                          group = "Tumor", mean_beta = .8, prevalence = 1, n = 2L))
shared <- atlas_merge_summaries(list(A = a2, B = b), mode = "shared")
stopifnot(nrow(shared) == 2L, all(shared$probe == "cg00000001"))

record <- list(sequence = "AACGATCGGCGTACG", strand = 1L, region_start = 101L,
               region_end = 115L, chromosome = "1", assembly = "GRCh37")
coords <- make_cgr(record, "triangle")
site <- coords$cpg_position[which(coords$is_cpg)[[1L]]]
tmp <- tempfile(fileext = ".csv")
x <- data.frame(cancer = "BRCA", cohort = "EXTERNAL_BRCA", sample = c("T1", "N1"),
                group = c("Tumor", "Normal"), probe = "cg00000001", beta = c(.9, .2),
                assembly = "GRCh37", chrom = "1", position = site)
utils::write.csv(x, tmp, row.names = FALSE)
stopifnot(inherits(try(atlas_read_user_csv(tmp, record, cancer = "LUAD"), silent = TRUE), "try-error"))
loaded <- atlas_read_user_csv(tmp, record, cancer = "BRCA")
stopifnot(length(loaded$cohorts) == 1L)
unlink(tmp)
cat("Cancer filter and cohort pooling checks passed\n")
