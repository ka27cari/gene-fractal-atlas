# Regression test for v0.2.0 curl::multi_add(url=...) crash.
# Run from repository root: Rscript tests/api_contract_test.R
# Uses libcurl's local file:// protocol so no Ensembl/internet is required.
source("R/ensembl.R")

stopifnot(requireNamespace("curl", quietly = TRUE))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))

# curl::multi_add does not accept a URL argument. URLs belong on handles.
stopifnot("handle" %in% names(formals(curl::multi_add)))
stopifnot(!"url" %in% names(formals(curl::multi_add)))

sandbox <- tempfile("ensembl-curl-contract-")
dir.create(sandbox)
writeLines('{"result":"first","value":1}', file.path(sandbox, "a.json"))
writeLines('{"result":"second","value":2}', file.path(sandbox, "b.json"))
old_base <- ENSEMBL_BASE
ENSEMBL_BASE <- paste0("file://", normalizePath(sandbox, winslash = "/"))

responses <- ensembl_json_batch(c("/b.json", "/a.json", "/missing.json"))
stopifnot(length(responses$results) == 3L)
stopifnot(identical(responses$results[[1L]]$result, "second"))
stopifnot(identical(responses$results[[2L]]$result, "first"))
stopifnot(is.null(responses$results[[3L]]))
stopifnot(nzchar(responses$errors[[3L]]))
stopifnot(!nzchar(responses$errors[[1L]]), !nzchar(responses$errors[[2L]]))

ENSEMBL_BASE <- old_base
unlink(sandbox, recursive = TRUE)
cat("Curl multi_add contract test passed (parallel local requests, order and failures).\n")
