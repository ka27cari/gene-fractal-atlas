# Independent sequence fallback. Gene coordinates: MyGene.info; DNA: UCSC.
# No methylation values are synthesized here. Never silently mix assemblies.
# UCSC query: 0-based half-open; records exposed to app: 1-based inclusive.
# Supported assembly mapping follows MyGene.info's documented genomic_pos
# reference versions. Other species are deliberately not guessed.

UCSC_API_BASE <- "https://api.genome.ucsc.edu"
MYGENE_API_BASE <- "https://mygene.info/v3"

UCSC_ASSEMBLIES <- list(
  homo_sapiens = list(taxid = 9606L, genome = "hg38", assembly = "GRCh38"),
  mus_musculus = list(taxid = 10090L, genome = "mm10", assembly = "GRCm38"),
  rattus_norvegicus = list(taxid = 10116L, genome = "rn6", assembly = "Rnor_6.0"),
  danio_rerio = list(taxid = 7955L, genome = "danRer10", assembly = "GRCz10")
)

# Retrieve several independent HTTPS URLs in a single concurrent batch.
# A remote failure is returned per request and never throws through Shiny.
public_json_batch <- function(urls, timeout = 9L) {
  n <- length(urls)
  out <- vector("list", n)
  errors <- rep("", n)
  if (!n) return(list(results = out, errors = errors))
  pool <- curl::new_pool(total_con = min(4L, n), host_con = min(4L, n))
  for (i in seq_len(n)) {
    local({
      index <- i
      handle <- curl::new_handle(url = urls[[index]])
      curl::handle_setopt(handle, timeout = as.integer(timeout),
                          connecttimeout = 4L, followlocation = TRUE,
                          useragent = "gene-fractal-atlas/0.7.1")
      curl::handle_setheaders(handle, Accept = "application/json")
      curl::multi_add(
        handle = handle, pool = pool,
        done = function(result) {
          if (result$status_code >= 400L) {
            errors[index] <<- sprintf("HTTP %d", result$status_code)
            return(invisible(NULL))
          }
          parsed <- tryCatch(
            jsonlite::fromJSON(rawToChar(result$content), simplifyVector = FALSE),
            error = function(e) {
              errors[index] <<- paste("Invalid JSON:", conditionMessage(e))
              NULL
            }
          )
          out[index] <<- list(parsed)
        },
        fail = function(message) { errors[index] <<- as.character(message) }
      )
    })
  }
  tryCatch(curl::multi_run(pool = pool),
           error = function(e) errors[] <<- conditionMessage(e))
  for (i in seq_len(n)) {
    if (is.null(out[[i]]) && !nzchar(errors[[i]])) errors[[i]] <- "No response"
  }
  list(results = out, errors = errors)
}

mygene_query_url <- function(symbol, species) {
  setting <- UCSC_ASSEMBLIES[[species]]
  if (is.null(setting)) stop("No verified assembly mapping for species")
  sym <- preferred_symbol(symbol, species)
  paste0(MYGENE_API_BASE, "/query?q=",
         utils::URLencode(paste0("symbol:", sym), reserved = TRUE),
         "&species=", setting$taxid,
         "&fields=symbol,taxid,entrezgene,ensembl,genomic_pos,genomic_pos_hg19",
         "&size=10")
}

# Choose only exact symbol/taxon matches on a primary chromosome.
# Many hits can be transcript isoforms or mappings on alternate contigs.
mygene_pick_gene <- function(payload, symbol, species, assembly = "current") {
  setting <- UCSC_ASSEMBLIES[[species]]
  if (is.null(setting) || is.null(payload$hits)) return(NULL)
  possible <- unique(tolower(symbol_candidates(symbol, species)))
  key <- if (identical(assembly, "GRCh37")) "genomic_pos_hg19" else "genomic_pos"
  for (hit in payload$hits) {
    if (is.null(hit$symbol) || is.null(hit$taxid) ||
        !(tolower(as.character(hit$symbol)) %in% possible) ||
        !identical(as.integer(hit$taxid), setting$taxid)) next
    positions <- hit[[key]]
    if (is.null(positions)) next
    if (!is.null(positions$chr)) positions <- list(positions)
    for (pos in positions) {
      chrom <- if (!is.null(pos$chr)) as.character(pos$chr) else ""
      begin <- suppressWarnings(as.numeric(pos$start))
      finish <- suppressWarnings(as.numeric(pos$end))
      strand <- suppressWarnings(as.integer(pos$strand))
      if (!grepl("^([0-9]{1,2}|X|Y|MT|M)$", chrom, ignore.case = TRUE)) next
      if (length(begin) != 1L || length(finish) != 1L ||
          !is.finite(begin) || !is.finite(finish) || begin < 1 ||
          finish < begin || length(strand) != 1L || !is.finite(strand) || !strand %in% c(-1L, 1L)) next
      ensembl_id <- if (!is.null(hit$ensembl)) {
        x <- hit$ensembl
        if (!is.null(x$gene)) as.character(x$gene) else
          if (is.list(x) && length(x) && !is.null(x[[1L]]$gene)) as.character(x[[1L]]$gene) else ""
      } else ""
      id <- if (nzchar(ensembl_id)) ensembl_id else
        if (!is.null(hit$entrezgene)) paste0("NCBI:", hit$entrezgene) else "mygene_annotation"
      return(list(chr = chrom, start = as.integer(begin), end = as.integer(finish),
                  strand = strand, gene_id = id, symbol = as.character(hit$symbol)))
    }
  }
  NULL
}

# Construct a 1-based inclusive window near the 5' gene boundary.
# The window is then converted to UCSC's 0-based, half-open indexing.
ucsc_window <- function(gene, window_bp) {
  anchor <- if (gene$strand == 1L) gene$start else gene$end
  begin <- as.integer(max(1, anchor - window_bp %/% 2L))
  list(start = begin, end = begin + window_bp - 1L)
}

ucsc_sequence_url <- function(species, gene, window_bp, assembly = "current") {
  settings <- UCSC_ASSEMBLIES[[species]]
  genome <- if (identical(assembly, "GRCh37")) "hg19" else settings$genome
  interval <- ucsc_window(gene, window_bp)
  sprintf("%s/getData/sequence?genome=%s&chrom=chr%s&start=%d&end=%d",
          UCSC_API_BASE, genome, gene$chr, interval$start - 1L, interval$end)
}

ucsc_record <- function(payload, gene, species, window_bp, assembly = "current") {
  setting <- UCSC_ASSEMBLIES[[species]]
  expected <- if (identical(assembly, "GRCh37")) "GRCh37" else setting$assembly
  genome <- if (identical(assembly, "GRCh37")) "hg19" else setting$genome
  dna <- if (!is.null(payload$dna)) toupper(gsub("\\s+", "", as.character(payload$dna))) else ""
  if (!nzchar(dna) || nchar(dna) != window_bp || grepl("[^ACGTN]", dna)) {
    stop("UCSC returned a missing or invalid reference DNA window")
  }
  if (!is.null(payload$genome) && !identical(as.character(payload$genome), genome)) {
    stop("UCSC sequence assembly differs from annotation assembly")
  }
  interval <- ucsc_window(gene, window_bp)
  list(species = species, species_label = species_label(species),
       symbol = gene$symbol, gene_id = gene$gene_id,
       assembly = expected, chromosome = gene$chr,
       region_start = interval$start, region_end = interval$end,
       strand = gene$strand,
       sequence = if (gene$strand == -1L) reverse_complement(dna) else dna,
       source = "ucsc", orthology = "symbol_match_unverified")
}

# api parameter is injectable so both provider failures and positive paths
# can be tested without touching a network.
fetch_ucsc_fallback <- function(symbol, selected, window_bp,
                                assembly = "current", api = public_json_batch) {
  selected <- unique(selected)
  allowed <- intersect(selected, names(UCSC_ASSEMBLIES))
  if (identical(assembly, "GRCh37")) allowed <- intersect(allowed, "homo_sapiens")
  records <- list()
  errors <- character()
  if (!length(allowed)) return(list(records = records, errors = errors))
  queries <- vapply(allowed, function(s) mygene_query_url(symbol, s), character(1))
  lookups <- api(queries)
  genes <- list()
  for (i in seq_along(allowed)) {
    s <- allowed[[i]]
    if (is.null(lookups$results[[i]])) {
      errors <- c(errors, sprintf("%s MyGene: %s", species_label(s), lookups$errors[[i]]))
      next
    }
    g <- mygene_pick_gene(lookups$results[[i]], symbol, s, assembly)
    if (is.null(g)) {
      errors <- c(errors, sprintf("%s MyGene: no exact symbol and coordinate match", species_label(s)))
      next
    }
    genes[[s]] <- g
  }
  ready <- names(genes)
  if (!length(ready)) return(list(records = records, errors = errors))
  urls <- vapply(ready, function(s) ucsc_sequence_url(s, genes[[s]], window_bp, assembly), character(1))
  sequences <- api(urls)
  for (i in seq_along(ready)) {
    s <- ready[[i]]
    if (is.null(sequences$results[[i]])) {
      errors <- c(errors, sprintf("%s UCSC sequence: %s", species_label(s), sequences$errors[[i]]))
      next
    }
    rec <- tryCatch(ucsc_record(sequences$results[[i]], genes[[s]], s, window_bp, assembly),
                    error = function(e) {
                      errors <<- c(errors, sprintf("%s UCSC: %s", species_label(s), conditionMessage(e)))
                      NULL
                    })
    if (!is.null(rec)) records[[s]] <- rec
  }
  list(records = records, errors = errors)
}

fetch_ucsc_grch37_human <- function(symbol, window_bp = 20000L, api = public_json_batch) {
  result <- fetch_ucsc_fallback(symbol, "homo_sapiens", as.integer(window_bp),
                                assembly = "GRCh37", api = api)
  if (!length(result$records)) stop(paste(result$errors, collapse = "; "))
  result$records[[1L]]
}
