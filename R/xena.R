# Read a small, real sample subset from selectable TCGA methylation matrices.
# Public UCSC Xena legacy TCGA hub, HM450 probe map coordinates = hg19.
# No guessed beta values, sample identifiers or genomic coordinates.

XENA_HUB <- "https://tcga.xenahubs.net"
# Selectable TCGA 450K cancer cohorts with tumor (01) and potential matched
# solid-tissue-normal (11) measurements; sample availability is validated live.
TCGA_COHORTS <- c(
  "Lung adenocarcinoma (LUAD)" = "LUAD",
  "Lung squamous carcinoma (LUSC)" = "LUSC",
  "Breast carcinoma (BRCA)" = "BRCA",
  "Colon adenocarcinoma (COAD)" = "COAD",
  "Kidney clear cell carcinoma (KIRC)" = "KIRC",
  "Liver carcinoma (LIHC)" = "LIHC",
  "Head and neck carcinoma (HNSC)" = "HNSC",
  "Prostate adenocarcinoma (PRAD)" = "PRAD",
  "Thyroid carcinoma (THCA)" = "THCA",
  "Uterine endometrial carcinoma (UCEC)" = "UCEC",
  "Bladder urothelial carcinoma (BLCA)" = "BLCA",
  "Cervical squamous/endocervical (CESC)" = "CESC",
  "Stomach adenocarcinoma (STAD)" = "STAD",
  "Cholangiocarcinoma (CHOL)" = "CHOL",
  "Esophageal carcinoma (ESCA)" = "ESCA",
  "Kidney papillary carcinoma (KIRP)" = "KIRP",
  "Pancreatic adenocarcinoma (PAAD)" = "PAAD",
  "Rectal adenocarcinoma (READ)" = "READ"
)
XENA_DATASET <- "TCGA.LUAD.sampleMap/HumanMethylation450"

xena_dataset <- function(cancer) {
  if (length(cancer) != 1L || is.na(cancer) ||
      !cancer %in% unname(TCGA_COHORTS)) stop("Unsupported TCGA cohort")
  sprintf("TCGA.%s.sampleMap/HumanMethylation450", cancer)
}

xena_cancer_label <- function(cancer) {
  names(TCGA_COHORTS)[match(cancer, unname(TCGA_COHORTS))]
}

XENA_SAMPLE_QUERY <- paste(
  '(fn [dataset]',
  ' (map :value',
  '  (query {:select [:value]',
  '          :from [:dataset]',
  '          :join [:field [:= :dataset.id :dataset_id]',
  '                 :code [:= :field.id :field_id]]',
  '          :where [:and',
  '                  [:= :dataset.name dataset]',
  '                  [:= :field.name "sampleID"]]})))',
  sep = "\n"
)

XENA_GENE_QUERY <- paste(
  '(fn [dataset samples genes]',
  ' (let [probemap (:probemap (car (query {:select [:probemap]',
  '              :from [:dataset]',
  '              :where [:= :name dataset]})))',
  '       position (xena-query {:select ["name" "position"]',
  '                :from [probemap]',
  '                :where [:in :any "genes" genes]})',
  '       probes (position "name")]',
  '  [position (fetch [{:table dataset :samples samples :columns probes}])]))',
  sep = "\n"
)

xena_literal <- function(value) {
  # JSON double-quoted strings and arrays are EDN-compatible for these values.
  jsonlite::toJSON(value, auto_unbox = TRUE, pretty = FALSE)
}

xena_post <- function(query, post = NULL) {
  if (!is.null(post)) return(post(query))
  h <- curl::new_handle()
  curl::handle_setopt(h, post = TRUE, postfields = query,
                      timeout = 35L, connecttimeout = 8L,
                      useragent = "gene-fractal-atlas/0.7.0")
  curl::handle_setheaders(h, `Content-Type` = "text/plain", Accept = "application/json")
  r <- curl::curl_fetch_memory(paste0(XENA_HUB, "/data/"), handle = h)
  if (r$status_code >= 400L) stop(sprintf("UCSC Xena HTTP %d", r$status_code))
  jsonlite::fromJSON(rawToChar(r$content), simplifyVector = FALSE)
}

# Sample type is encoded in TCGA barcode characters 14-15:
# 01 = Primary Solid Tumor, 11 = Solid Tissue Normal.
tcga_sample_type <- function(barcode) {
  out <- rep(NA_character_, length(barcode))
  type <- substr(barcode, 14L, 15L)
  good <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}", barcode)
  out[good & type == "01"] <- "Tumor"
  out[good & type == "11"] <- "Normal"
  out
}

choose_tcga_samples <- function(all_samples, per_group = 2L) {
  samples <- unique(as.character(unlist(all_samples, use.names = FALSE)))
  types <- tcga_sample_type(samples)
  samples <- samples[!is.na(types)]
  types <- types[!is.na(types)]
  tumors <- sort(samples[types == "Tumor"])
  normals <- sort(samples[types == "Normal"])
  if (!length(tumors) || !length(normals)) {
    stop("Selected TCGA dataset lacks tumor or solid-tissue-normal sample IDs")
  }
  # Prefer patient pairs; otherwise choose unpaired tissue samples and label them.
  patient <- function(z) substr(z, 1L, 12L)
  matched <- sort(intersect(patient(tumors), patient(normals)))
  if (length(matched)) {
    matched <- head(matched, per_group)
    sel_t <- vapply(matched, function(id) sort(tumors[patient(tumors) == id])[[1L]], character(1))
    sel_n <- vapply(matched, function(id) sort(normals[patient(normals) == id])[[1L]], character(1))
    selection <- as.vector(rbind(sel_t, sel_n))
    return(data.frame(sample = selection,
                      group = tcga_sample_type(selection),
                      paired = TRUE, stringsAsFactors = FALSE))
  }
  selection <- c(head(tumors, per_group), head(normals, per_group))
  data.frame(sample = selection, group = tcga_sample_type(selection),
             paired = FALSE, stringsAsFactors = FALSE)
}

# Convert Xena probe position to candidate forward-reference C positions.
# Xena probe maps generally use zero-based, half-open coordinates, but
# confirm candidates against CpGs found in the actual GRCh37 sequence.
parse_probe_position <- function(position) {
  if (is.null(position)) return(list(chrom = NA_character_, start = NA_real_))
  v <- unlist(position, use.names = TRUE)
  named <- names(v)
  if (!is.null(named) && any(c("chrom", "chr", "chromosome") %in% named)) {
    chrom_key <- intersect(c("chrom", "chr", "chromosome"), named)[1L]
    start_key <- intersect(c("chromStart", "chromstart", "start", "position"), named)
    if (length(start_key)) {
      return(list(chrom = as.character(v[[chrom_key]]),
                  start = suppressWarnings(as.numeric(v[[start_key[[1L]]]]))))
    }
  }
  if (length(v) >= 2L && grepl("^(chr)?([0-9]+|X|Y|M|MT)$", as.character(v[[1L]]), ignore.case = TRUE)) {
    return(list(chrom = as.character(v[[1L]]),
                start = suppressWarnings(as.numeric(v[[2L]]))))
  }
  if (length(v) >= 1L && grepl("^(chr)?[[:alnum:]_]+:[0-9]+(-[0-9]+)?$", as.character(v[[1L]]))) {
    txt <- as.character(v[[1L]])
    pieces <- strsplit(txt, ":", fixed = TRUE)[[1L]]
    return(list(chrom = pieces[[1L]],
                start = as.numeric(strsplit(pieces[[2L]], "-", fixed = TRUE)[[1L]][[1L]])))
  }
  list(chrom = NA_character_, start = NA_real_)
}

xena_gene_records <- function(symbol, sample_info, reference_record, post = NULL,
                              dataset = XENA_DATASET) {
  payload <- sprintf("(%s %s %s %s)", XENA_GENE_QUERY,
                     xena_literal(dataset), xena_literal(sample_info$sample),
                     paste0("[", xena_literal(symbol), "]"))
  result <- xena_post(payload, post = post)
  if (!is.list(result) || length(result) < 2L || is.null(result[[1L]])) {
    stop("Xena response did not contain probe map and values")
  }
  position <- result[[1L]]
  if (is.null(position$name) || is.null(position$position)) {
    stop("Xena probe map missing names or positions")
  }
  ids <- as.character(unlist(position$name, use.names = FALSE))
  maps <- position$position
  values <- result[[2L]]
  if (!length(ids)) stop(sprintf("No HM450 probes annotated to %s", symbol))
  if (length(ids) != length(maps) || length(values) != length(ids)) {
    stop("Xena probe map and methylation matrix lengths disagree")
  }

  # Both the input and reference must be assembly-consistent.
  if (!identical(reference_record$assembly, "GRCh37")) {
    stop("TCGA legacy probe coordinates require a GRCh37 reference sequence")
  }
  reference_cgr <- make_cgr(reference_record, "triangle")
  cpg <- sort(unique(reference_cgr$cpg_position[reference_cgr$is_cpg]))
  chrom <- sub("^chr", "", reference_record$chromosome, ignore.case = TRUE)

  frames <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    loc <- parse_probe_position(maps[[i]])
    if (!is.finite(loc$start) || is.na(loc$chrom) ||
        !identical(tolower(sub("^chr", "", loc$chrom, ignore.case = TRUE)), tolower(chrom))) next
    # Probe mapping may use 0- or 1-based starts. Keep only an observed CpG C.
    candidate <- c(loc$start + 1, loc$start)
    valid <- candidate[candidate %in% cpg]
    if (!length(valid)) next
    betas <- suppressWarnings(as.numeric(unlist(values[[i]], use.names = FALSE)))
    if (length(betas) != nrow(sample_info)) {
      stop(sprintf("Xena values for probe %s do not match sample count", ids[[i]]))
    }
    good <- is.finite(betas) & betas >= 0 & betas <= 1
    if (!any(good)) next
    frames[[i]] <- data.frame(sample = sample_info$sample[good],
                             group = sample_info$group[good],
                             paired = sample_info$paired[good],
                             species = "homo_sapiens", assembly = "GRCh37",
                             chrom = chrom, position = valid[[1L]],
                             beta = betas[good], probe = ids[[i]],
                             stringsAsFactors = FALSE)
  }
  mapped <- Filter(Negate(is.null), frames)
  if (!length(mapped)) {
    stop("No Xena HM450 probes aligned to validated CpGs in the GRCh37 sequence window")
  }
  out <- do.call(rbind, mapped)
  # A rare duplicate probe position is summarized by the first observed probe.
  out <- out[!duplicated(paste(out$sample, out$position, sep = "|")), , drop = FALSE]
  rownames(out) <- NULL
  out
}

fetch_tcga <- function(symbol, cancer = "LUAD", window_bp = 20000L, post = NULL,
                       reference_fetch = fetch_grch37_human) {
  gene <- trimws(symbol)
  if (!valid_gene_symbol(gene)) stop("Invalid human gene symbol")
  dataset <- xena_dataset(cancer)
  raw_samples <- xena_post(sprintf("(%s %s)", XENA_SAMPLE_QUERY,
                                   xena_literal(dataset)), post = post)
  sample_info <- choose_tcga_samples(raw_samples)
  if (!nrow(sample_info)) stop("No qualifying TCGA samples")
  rec <- reference_fetch(gene, window_bp = window_bp)
  values <- xena_gene_records(gene, sample_info, rec, post = post, dataset = dataset)
  supported <- sample_info$sample %in% unique(values$sample)
  sample_info <- sample_info[supported, , drop = FALSE]
  if (!all(c("Tumor", "Normal") %in% sample_info$group)) {
    stop("Available TCGA probe values do not include both tissue classes")
  }
  list(record = rec, samples = sample_info, methylation = values,
       cohort = cancer, cancer_label = xena_cancer_label(cancer),
       platform = "HumanMethylation450", dataset = dataset, hub = XENA_HUB,
       note = "Tumor: TCGA type 01; normal tissue: TCGA type 11; hg19/GRCh37")
}

# Preserve the old entry point for tests and existing scripts.
fetch_tcga_luad <- function(symbol, window_bp = 20000L, post = NULL,
                            reference_fetch = fetch_grch37_human) {
  fetch_tcga(symbol, cancer = "LUAD", window_bp = window_bp,
             post = post, reference_fetch = reference_fetch)
}
