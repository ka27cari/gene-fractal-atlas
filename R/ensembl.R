# Ensembl REST client and human-anchored ortholog discovery.
# Reference docs: https://rest.ensembl.org/documentation/info/homology_id
#                 https://rest.ensembl.org/documentation/info/sequence_region

ENSEMBL_BASE <- "https://rest.ensembl.org"

SPECIES <- c(
  "Human" = "homo_sapiens",
  "Chimpanzee" = "pan_troglodytes",
  "Mouse" = "mus_musculus",
  "Rat" = "rattus_norvegicus",
  "Dog" = "canis_lupus_familiaris",
  "Chicken" = "gallus_gallus",
  "Zebrafish" = "danio_rerio"
)

# Small pragmatic fallback map for common symbol casing / historical aliases.
# It is used only if Compara or same-symbol lookup fails.
LOCAL_SYMBOL_FALLBACKS <- list(
  TP53 = c(
    homo_sapiens = "TP53",
    pan_troglodytes = "TP53",
    mus_musculus = "Trp53",
    rattus_norvegicus = "Tp53",
    canis_lupus_familiaris = "TP53",
    gallus_gallus = "TP53",
    danio_rerio = "tp53"
  ),
  BRCA1 = c(
    homo_sapiens = "BRCA1",
    pan_troglodytes = "BRCA1",
    mus_musculus = "Brca1",
    rattus_norvegicus = "Brca1",
    canis_lupus_familiaris = "BRCA1",
    gallus_gallus = "BRCA1",
    danio_rerio = "brca1"
  ),
  TET2 = c(
    homo_sapiens = "TET2",
    pan_troglodytes = "TET2",
    mus_musculus = "Tet2",
    rattus_norvegicus = "Tet2",
    canis_lupus_familiaris = "TET2",
    gallus_gallus = "TET2",
    danio_rerio = "tet2"
  )
)

ensembl_handle <- function(endpoint) {
  handle <- curl::new_handle(url = paste0(ENSEMBL_BASE, endpoint))
  curl::handle_setopt(handle, timeout = 7L, connecttimeout = 4L,
                      useragent = "gene-fractal-atlas/0.8.0")
  curl::handle_setheaders(handle, Accept = "application/json")
  handle
}

ensembl_json_get <- function(endpoint) {
  out <- tryCatch({
    result <- curl::curl_fetch_memory(url = paste0(ENSEMBL_BASE, endpoint), handle = ensembl_handle(endpoint))
    code <- result$status_code
    if (code >= 400L) {
      msg <- if (code == 404L) "Not found (HTTP 404)" else
        if (code == 429L) "Ensembl rate limit (HTTP 429)" else
          sprintf("Ensembl HTTP %s", code)
      return(list(result = NULL, error = msg))
    }
    parsed <- jsonlite::fromJSON(rawToChar(result$content), simplifyVector = FALSE)
    list(result = parsed, error = "")
  }, error = function(e) list(result = NULL, error = conditionMessage(e)))
  out
}

# Concurrent reads reduce latency when several species are selected.
# Every result is retained in its original request order; failures stay local.
ensembl_json_batch <- function(endpoints) {
  n <- length(endpoints)
  responses <- vector("list", n)
  errors <- rep("", n)
  if (n == 0L) return(list(results = responses, errors = errors))

  pool <- curl::new_pool(total_con = min(5L, n), host_con = min(5L, n))
  for (i in seq_len(n)) {
    local({
      index <- i
      curl::multi_add(
        handle = ensembl_handle(endpoints[[index]]), pool = pool,
        done = function(result) {
          code <- result$status_code
          if (code >= 400L) {
            errors[index] <<- if (code == 404L) "Not found (HTTP 404)" else
              if (code == 429L) "Ensembl rate limit (HTTP 429)" else
                sprintf("Ensembl HTTP %s", code)
            return(invisible(NULL))
          }
          parsed <- tryCatch(
            jsonlite::fromJSON(rawToChar(result$content), simplifyVector = FALSE),
            error = function(e) {
              errors[index] <<- sprintf("Invalid Ensembl JSON: %s", conditionMessage(e))
              NULL
            }
          )
          responses[index] <<- list(parsed)
        },
        fail = function(message) {
          errors[index] <<- as.character(message)
        }
      )
    })
  }
  curl::multi_run(pool = pool)
  for (i in seq_len(n)) {
    if (is.null(responses[[i]]) && !nzchar(errors[[i]])) errors[[i]] <- "No response from Ensembl"
  }
  list(results = responses, errors = errors)
}

valid_gene_symbol <- function(symbol) {
  length(symbol) == 1L && !is.na(symbol) &&
    grepl("^[A-Za-z0-9_.-]{1,40}$", symbol)
}

species_label <- function(species) names(SPECIES)[match(species, unname(SPECIES))]

# Reject ambiguous user shorthand before requesting a remote API.
# BRCA names a TCGA cancer cohort, but not an unambiguous gene symbol.
gene_input_issue <- function(symbol) {
  if (length(symbol) != 1L || is.na(symbol) || !nzchar(trimws(symbol))) {
    return("Enter a gene symbol such as TP53, BRCA1 or BRCA2.")
  }
  if (identical(toupper(trimws(symbol)), "BRCA")) {
    return("BRCA is ambiguous. Choose BRCA1 or BRCA2; BRCA is a cancer cohort code, not one gene.")
  }
  if (!valid_gene_symbol(trimws(symbol))) {
    return("Invalid gene symbol. Use letters, digits, underscore, period or hyphen.")
  }
  NULL
}

symbol_candidates <- function(symbol, species) {
  base <- trimws(as.character(symbol))
  upper <- toupper(base)
  candidates <- c(base, upper, tolower(base),
                  paste0(toupper(substr(base, 1L, 1L)), tolower(substr(base, 2L, nchar(base)))))
  fallback <- LOCAL_SYMBOL_FALLBACKS[[upper]]
  if (!is.null(fallback) && species %in% names(fallback)) {
    candidates <- c(candidates, fallback[[species]])
  }
  unique(candidates[nzchar(candidates)])
}

# Select one Ensembl Compara ortholog per species, preferring one-to-one and
# then high-confidence/high-sequence-identity hits. Never guess by gene name.
resolve_orthologs <- function(payload, target_species) {
  ids <- stats::setNames(rep(NA_character_, length(target_species)), target_species)
  types <- stats::setNames(rep(NA_character_, length(target_species)), target_species)
  warnings <- character()
  if (is.null(payload$data) || !length(payload$data)) {
    return(list(ids = ids, types = types, warnings = warnings))
  }
  homologies <- payload$data[[1L]]$homologies
  if (is.null(homologies)) homologies <- list()
  for (species in target_species) {
    hits <- Filter(function(h) {
      !is.null(h$target$species) &&
        identical(tolower(as.character(h$target$species)), species) &&
        !is.null(h$target$id) && grepl("^ortholog", h$type)
    }, homologies)
    if (!length(hits)) next
    one_to_one <- vapply(hits, function(h) identical(h$type, "ortholog_one2one"), logical(1))
    confidence <- vapply(hits, function(h) {
      z <- suppressWarnings(as.numeric(h$confidence))
      if (length(z) != 1L || !is.finite(z)) 0 else z
    }, numeric(1))
    identity <- vapply(hits, function(h) {
      z <- suppressWarnings(as.numeric(h$target$perc_id))
      if (length(z) != 1L || !is.finite(z)) 0 else z
    }, numeric(1))
    gene_ids <- vapply(hits, function(h) as.character(h$target$id), character(1))
    chosen <- order(-as.integer(one_to_one), -confidence, -identity, gene_ids)[1L]
    ids[[species]] <- gene_ids[[chosen]]
    types[[species]] <- hits[[chosen]]$type
    if (length(unique(gene_ids)) > 1L) {
      warnings <- c(warnings, sprintf(
        "%s: %d ortholog candidates; selected %s (%s)",
        species_label(species), length(unique(gene_ids)), ids[[species]], types[[species]]))
    }
  }
  list(ids = ids, types = types, warnings = warnings)
}

# Build one region request from a gene metadata response. 5-prime gene boundary
# is not necessarily a transcript-specific transcription start site.
gene_interval <- function(gene, window_bp) {
  if (is.null(gene) || !identical(gene$object_type, "Gene")) {
    stop("Ensembl did not return a gene")
  }
  strand <- as.integer(gene$strand)
  if (length(strand) != 1L || !strand %in% c(-1L, 1L)) stop("Invalid gene strand")
  anchor <- if (strand == 1L) as.integer(gene$start) else as.integer(gene$end)
  if (length(anchor) != 1L || !is.finite(anchor)) stop("Invalid gene coordinate")
  chrom <- as.character(gene$seq_region_name)
  if (length(chrom) != 1L || !nzchar(chrom)) stop("Missing chromosome")
  start <- max(1L, anchor - window_bp %/% 2L)
  end <- start + window_bp - 1L
  region <- sprintf("%s:%d..%d:1", chrom, start, end)
  list(start = start, end = end, chromosome = chrom, strand = strand, region = region)
}

record_from_sequence <- function(species, gene, interval, seq_result, window_bp,
                                 orthology_type = "reference") {
  if (is.null(seq_result$seq)) stop("No sequence returned")
  dna <- toupper(gsub("\\s+", "", seq_result$seq))
  if (nchar(dna) != window_bp) stop("Ensembl sequence length differs from requested window")
  # The gene metadata is authoritative for its reference assembly. Parsing
  # the sequence identifier alone previously attached GRCh38 to non-human
  # species. Reject conflicting assemblies rather than silently mislabel DNA.
  id_fields <- strsplit(if (is.null(seq_result$id)) "" else seq_result$id,
                        ":", fixed = TRUE)[[1L]]
  assembly_from_sequence <- if (length(id_fields) >= 3L) id_fields[[2L]] else ""
  assembly_from_gene <- if (!is.null(gene$assembly_name))
    as.character(gene$assembly_name) else ""
  if (length(assembly_from_gene) != 1L) assembly_from_gene <- ""
  if (nzchar(assembly_from_gene) && nzchar(assembly_from_sequence) &&
      !identical(assembly_from_gene, assembly_from_sequence)) {
    stop(sprintf("Genome assembly conflict: gene=%s; sequence=%s",
                 assembly_from_gene, assembly_from_sequence))
  }
  assembly <- if (nzchar(assembly_from_gene)) assembly_from_gene else assembly_from_sequence
  if (!nzchar(assembly)) stop("Unknown genome assembly; sequence rejected")
  display_name <- if (!is.null(gene$display_name) && nzchar(gene$display_name))
    as.character(gene$display_name) else as.character(gene$id)
  list(species = species, species_label = species_label(species),
       symbol = display_name, gene_id = as.character(gene$id),
       assembly = assembly, chromosome = interval$chromosome,
       region_start = interval$start, region_end = interval$end,
       strand = interval$strand,
       sequence = if (interval$strand == -1L) reverse_complement(dna) else dna,
       source = "ensembl", orthology = orthology_type)
}

lookup_symbol_once <- function(species, symbol) {
  endpoint <- sprintf("/lookup/symbol/%s/%s?expand=0",
                      species, utils::URLencode(symbol, reserved = TRUE))
  ensembl_json_get(endpoint)
}

lookup_symbol_fallback <- function(species, symbol, api_batch = ensembl_json_batch) {
  candidates <- symbol_candidates(symbol, species)
  endpoints <- vapply(candidates, function(candidate) {
    sprintf("/lookup/symbol/%s/%s?expand=0", species,
            utils::URLencode(candidate, reserved = TRUE))
  }, character(1))
  lookup <- api_batch(endpoints)
  available <- which(!vapply(lookup$results, is.null, logical(1)))
  if (length(available)) {
    chosen <- available[[1L]]
    return(list(result = lookup$results[[chosen]], error = "", candidate = candidates[[chosen]]))
  }
  list(result = NULL,
       error = sprintf("No direct symbol match for %s (tried %s)",
                       species_label(species), paste(candidates, collapse = ", ")),
       candidate = NA_character_)
}

# Human gene symbol is the anchor. Other species are resolved by orthology IDs.
# If that fails, species-specific symbol lookup fallbacks are attempted.
fetch_comparative_regions <- function(symbol, selected, window_bp = 20000L,
                                      cache = NULL, api_batch = ensembl_json_batch) {
  if (!valid_gene_symbol(symbol)) stop("Invalid human gene symbol")
  if (!length(selected) || anyNA(selected) || any(!selected %in% unname(SPECIES))) {
    stop("Unsupported or missing species")
  }
  if (length(window_bp) != 1L || !is.finite(window_bp) || window_bp < 1000 || window_bp > 1000000) {
    stop("Invalid window size")
  }
  symbol <- trimws(symbol)
  window_bp <- as.integer(window_bp)
  selected <- unique(as.character(selected))
  records <- stats::setNames(vector("list", length(selected)), selected)
  problems <- character()
  new_species <- character()
  for (species in selected) {
    key <- paste(species, toupper(symbol), window_bp, sep = "|")
    if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
      records[species] <- list(get(key, envir = cache, inherits = FALSE))
    } else new_species <- c(new_species, species)
  }
  if (!length(new_species)) return(list(records = unname(records), errors = problems))

  human <- api_batch(sprintf("/lookup/symbol/homo_sapiens/%s?expand=0",
                             utils::URLencode(symbol, reserved = TRUE)))
  if (is.null(human$results[[1L]])) {
    return(list(records = unname(Filter(Negate(is.null), records)),
                errors = c(problems, sprintf("Human gene %s: %s", symbol, human$errors[[1L]]))))
  }
  reference <- human$results[[1L]]

  gene_meta <- stats::setNames(vector("list", length(new_species)), new_species)
  orthology <- stats::setNames(rep("reference", length(new_species)), new_species)
  if ("homo_sapiens" %in% new_species) gene_meta[["homo_sapiens"]] <- reference

  others <- setdiff(new_species, "homo_sapiens")
  if (length(others)) {
    endpoint <- sprintf("/homology/id/homo_sapiens/%s?type=orthologues&sequence=none&aligned=0",
                        utils::URLencode(as.character(reference$id), reserved = TRUE))
    homology <- api_batch(endpoint)
    mapped_ids <- stats::setNames(rep(NA_character_, length(others)), others)
    mapped_types <- stats::setNames(rep(NA_character_, length(others)), others)

    if (is.null(homology$results[[1L]])) {
      problems <- c(problems, sprintf("Ortholog lookup: %s", homology$errors[[1L]]))
    } else {
      mapped <- resolve_orthologs(homology$results[[1L]], others)
      problems <- c(problems, mapped$warnings)
      mapped_ids[names(mapped$ids)] <- mapped$ids
      mapped_types[names(mapped$types)] <- mapped$types
    }

    available <- others[!is.na(mapped_ids[others])]
    if (length(available)) {
      lookup <- api_batch(sprintf("/lookup/id/%s?expand=0", unname(mapped_ids[available])))
      for (i in seq_along(available)) {
        species <- available[[i]]
        if (is.null(lookup$results[[i]])) {
          problems <- c(problems, sprintf("%s: ortholog-id lookup failed: %s",
                                          species_label(species), lookup$errors[[i]]))
        } else {
          gene_meta[[species]] <- lookup$results[[i]]
          orthology[[species]] <- mapped_types[[species]]
        }
      }
    }

    unresolved <- others[vapply(gene_meta[others], is.null, logical(1))]
    if (length(unresolved)) {
      for (species in unresolved) {
        fallback <- lookup_symbol_fallback(species, symbol, api_batch = api_batch)
        if (!is.null(fallback$result)) {
          gene_meta[[species]] <- fallback$result
          orthology[[species]] <- paste0("symbol_fallback:", fallback$candidate)
          problems <- c(problems, sprintf("%s: used direct symbol fallback (%s)",
                                          species_label(species), fallback$candidate))
        } else {
          problems <- c(problems, fallback$error)
        }
      }
    }
  }

  ready <- names(gene_meta)[!vapply(gene_meta, is.null, logical(1))]
  intervals <- vector("list", length(ready))
  keep <- logical(length(ready))
  for (i in seq_along(ready)) {
    species <- ready[[i]]
    intervals[[i]] <- tryCatch(gene_interval(gene_meta[[species]], window_bp),
                               error = function(e) {
                                 problems <<- c(problems, sprintf("%s: %s", species_label(species), conditionMessage(e)))
                                 NULL
                               })
    keep[[i]] <- !is.null(intervals[[i]])
  }
  ready <- ready[keep]
  intervals <- intervals[keep]

  if (length(ready)) {
    endpoints <- vapply(seq_along(ready), function(i) {
      sprintf("/sequence/region/%s/%s", ready[[i]],
              utils::URLencode(intervals[[i]]$region, reserved = FALSE))
    }, character(1))
    sequences <- api_batch(endpoints)
    for (i in seq_along(ready)) {
      species <- ready[[i]]
      if (is.null(sequences$results[[i]])) {
        problems <- c(problems, sprintf("%s: %s", species_label(species), sequences$errors[[i]]))
        next
      }
      rec <- tryCatch(record_from_sequence(species, gene_meta[[species]], intervals[[i]],
                                           sequences$results[[i]], window_bp, orthology[[species]]),
                      error = function(e) {
                        problems <<- c(problems, sprintf("%s: %s", species_label(species), conditionMessage(e)))
                        NULL
                      })
      if (!is.null(rec)) {
        records[[species]] <- rec
        if (!is.null(cache)) {
          key <- paste(species, toupper(symbol), window_bp, sep = "|")
          assign(key, rec, envir = cache)
        }
      }
    }
  }

  ordered <- Filter(Negate(is.null), records[selected])
  list(records = unname(ordered), errors = unique(problems[nzchar(problems)]))
}

# Backward-compatible single-species entry point.
fetch_gene_region <- function(species, species_label, symbol, window_bp = 20000L) {
  result <- fetch_comparative_regions(symbol, selected = species, window_bp = window_bp)
  if (!length(result$records)) stop(paste(result$errors, collapse = "; "))
  result$records[[1L]]
}

# Human genomic window in the GRCh37 assembly, for comparison with
# legacy TCGA HumanMethylation450 probe coordinates (hg19).
fetch_grch37_human <- function(symbol, window_bp = 20000L) {
  if (!valid_gene_symbol(symbol)) stop("Invalid human gene symbol")
  base <- "https://grch37.rest.ensembl.org"
  path <- sprintf("/lookup/symbol/homo_sapiens/%s?expand=0",
                  utils::URLencode(symbol, reserved = TRUE))
  fetch_json <- function(endpoint) {
    url <- paste0(base, endpoint)
    handle <- curl::new_handle()
    curl::handle_setopt(handle, timeout = 7L, connecttimeout = 4L,
                        useragent = "gene-fractal-atlas/0.8.0")
    curl::handle_setheaders(handle, Accept = "application/json")
    response <- curl::curl_fetch_memory(url, handle = handle)
    if (response$status_code >= 400L) {
      stop(sprintf("Ensembl GRCh37 returned HTTP %d", response$status_code))
    }
    jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
  }
  primary <- tryCatch({
    gene <- fetch_json(path)
    interval <- gene_interval(gene, as.integer(window_bp))
    seq <- fetch_json(sprintf("/sequence/region/homo_sapiens/%s",
                              utils::URLencode(interval$region, reserved = FALSE)))
    rec <- record_from_sequence("homo_sapiens", gene, interval, seq,
                                as.integer(window_bp), orthology_type = "reference")
    rec$assembly <- "GRCh37"
    rec
  }, error = function(e) NULL)
  if (!is.null(primary)) return(primary)
  # TCGA 450K is mapped to GRCh37; only a matching hg19 window is acceptable.
  fetch_ucsc_grch37_human(symbol, as.integer(window_bp))
}
