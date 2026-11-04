# Gene Fractal Atlas | R Shiny
# From this directory: shiny::runApp(".")

options(shiny.maxRequestSize = 15 * 1024^2)
library(shiny)

source("R/fractal.R", local = TRUE)
source("R/ensembl.R", local = TRUE)
source("R/fast.R", local = TRUE)
source("R/ucsc.R", local = TRUE)
source("R/methylation.R", local = TRUE)
source("R/xena.R", local = TRUE)
source("R/compare.R", local = TRUE)
source("R/svg.R", local = TRUE)
source("R/cpg_cohorts.R", local = TRUE)
source("R/atlas_cohorts.R", local = TRUE)
source("R/atlas_inspector.R", local = TRUE)
source("R/data_catalog.R", local = TRUE)

sequence_cache <- new.env(parent = emptyenv())

WINDOW_OPTIONS <- c("5 kb" = 5000, "10 kb" = 10000, "20 kb" = 20000,
                    "50 kb" = 50000, "100 kb" = 100000,
                    "250 kb" = 250000, "500 kb" = 500000, "1 Mb" = 1000000)
COMMON_GENES <- CANDIDATE_GENES
GEOMETRIES <- c("Triangle · AT / C / G" = "triangle",
                "Square · A / C / G / T" = "square",
                "Triangle · ATG / C / methyl-C" = "mc_triangle")

status_banner <- function(status) {
  if (is.null(status)) return(NULL)
  tags$div(class = paste("load-banner", status$kind), status$text)
}

plot_grid <- function(items, geometry, show_cpg, show_beta, threshold, reveal_pct,
                      star_highlights = TRUE) {
  tags$div(class = "fractal-grid",
    lapply(items, function(item) {
      rec <- item$record
      label <- if (is.null(rec$sample_id)) rec$species_label else rec$group
      second <- if (is.null(rec$sample_id)) rec$symbol else rec$sample_id
      tags$div(class = "fractal-panel",
        tags$div(class = "fractal-label", tags$strong(label), tags$span(second)),
        tags$div(class = "fractal-square",
          HTML(svg_panel(item, geometry, show_cpg, show_beta,
                         methyl_threshold = threshold,
                         reveal_fraction = reveal_pct / 100,
                         star_highlights = star_highlights)))
      )
    })
  )
}

metric <- function(label, value, detail) {
  tags$div(class = "metric", tags$div(class = "metric-label", label),
           tags$div(class = "metric-value", value),
           tags$div(class = "metric-detail", detail))
}

nfmt <- function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)

source_footer <- function() {
  tags$footer(class = "site-footer",
    tags$div("Sources: ",
      tags$a(href = "https://rest.ensembl.org", target = "_blank", "Ensembl"), " · ",
      tags$a(href = "https://mygene.info", target = "_blank", "MyGene.info"), " · ",
      tags$a(href = "https://api.genome.ucsc.edu", target = "_blank", "UCSC Genome"), " · ",
      tags$a(href = "https://xenabrowser.net/datapages/", target = "_blank", "UCSC Xena / TCGA 450K"), " · ",
      tags$a(href = "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes",
             target = "_blank", "TCGA sample codes"), " · ",
      tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/", target = "_blank",
             "NCBI GEO (selected cancer-specific series)")),
    tags$span("TCGA 'normal' = solid tissue normal, not healthy controls"))
}

ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    shiny::includeCSS("www/styles.css"),
    tags$script(HTML('
      (function(){
        function selectCpG(el){
          if (!el || !window.Shiny || !Shiny.setInputValue) return;
          const pos = Number(el.getAttribute("data-cpg-position"));
          if (!Number.isSafeInteger(pos) || pos < 1) return;
          Shiny.setInputValue("picked_cpg", {
            position: pos,
            source: el.getAttribute("data-cpg-group") || "",
            nonce: Date.now()
          }, {priority: "event"});
        }
        function selectAtlasCell(el){
          if (!el || !window.Shiny || !Shiny.setInputValue) return;
          const ix = Number(el.getAttribute("data-bin-ix"));
          const iy = Number(el.getAttribute("data-bin-iy"));
          const bins = Number(el.getAttribute("data-bin-count"));
          const cohort = el.getAttribute("data-cohort") || "";
          const group = el.getAttribute("data-group") || "";
          if (!Number.isSafeInteger(ix) || !Number.isSafeInteger(iy) ||
              !Number.isSafeInteger(bins) || bins < 1 ||
              ix < 0 || iy < 0 || ix >= bins || iy >= bins) return;
          Shiny.setInputValue("picked_atlas_cell", {
            ix: ix, iy: iy, bins: bins, cohort: cohort, group: group,
            nonce: Date.now()
          }, {priority: "event"});
          document.querySelectorAll(".atlas-cell-hit.is-selected").forEach(function(node){
            node.classList.remove("is-selected");
          });
          el.classList.add("is-selected");
        }
        document.addEventListener("click", function(e){
          const el = e.target.closest ? e.target.closest(".atlas-cell-hit") : null;
          if (el) selectAtlasCell(el);
          else selectCpG(e.target.closest ? e.target.closest(".cpg-hit") : null);
        });
        document.addEventListener("keydown", function(e){
          if (e.key !== "Enter" && e.key !== " ") return;
          const cell = e.target.closest ? e.target.closest(".atlas-cell-hit") : null;
          if (cell) { e.preventDefault(); selectAtlasCell(cell); return; }
          const hit = e.target.closest ? e.target.closest(".cpg-hit") : null;
          if (hit) { e.preventDefault(); selectCpG(hit); }
        });
        document.addEventListener("shiny:connected", function(){
          Shiny.addCustomMessageHandler("focus-atlas-cell", function(){
            window.setTimeout(function(){
              const box = document.getElementById("atlas_cell_inspector");
              if (box) box.scrollIntoView({behavior:"smooth",block:"start"});
            }, 80);
          });
          Shiny.addCustomMessageHandler("focus-cpg", function(){
            window.setTimeout(function(){
              const box = document.getElementById("cpg_browser_card");
              if (box) box.scrollIntoView({behavior:"smooth",block:"start"});
            }, 250);
          });
        });
      })();
    '))
  ),
  tags$div(class = "app-shell", style = "background:#000;color:#f0f5df;min-height:100vh;",
    tags$header(class = "masthead",
      tags$div(class = "mark", "✦"),
      tags$div(class = "brand",
        tags$div(class = "eyebrow", "DNA SEQUENCE VISUALIZATION"),
        tags$h1("Gene Fractal ", tags$span("Atlas"))
      )
    ),
    tags$div(class = "workspace-tabs",
      radioButtons("workspace", NULL,
                   choices = c("Gene comparison" = "gene", "Cancer / normal" = "cancer"),
                   selected = "gene", inline = TRUE)
    ),
    tags$div(class = "layout-grid",
      tags$main(class = "content",
        conditionalPanel(condition = "input.workspace == 'gene'",
          tags$div(class = "overview",
            tags$div(tags$div(class = "section-kicker", "GENE COMPARISON"),
                     tags$h2(textOutput("gene_headline", inline = TRUE)),
                     tags$p(class = "subheading", textOutput("gene_subtitle", inline = TRUE))),
            tags$div(class = "source-tag", textOutput("gene_source_tag", inline = TRUE))
          ),
          tags$div(class = "plot-card hero-card",
            tags$div(class = "card-head",
              tags$div(tags$div(class = "section-kicker", "CHAOS GAME REPRESENTATION"),
                       tags$strong("Species")),
              downloadLink("gene_csv", "↓ Export point data", class = "action-link")
            ),
            uiOutput("gene_status"),
            uiOutput("gene_atlas"),
            tags$div(class = "plot-footnote",
              "Green: sequence / CpG · Red: β above threshold. ",
              "Synthetic examples are not reference or patient data.")
          ),
          uiOutput("gene_metrics"),
          tags$div(class = "detail-card",
            tags$div(class = "card-head", tags$strong("Reference coverage")),
            tableOutput("gene_coverage"),
            tags$p(class = "helper", textOutput("gene_match_note", inline = TRUE)),
            tags$details(class = "technical-details",
              tags$summary("Lookup details"),
              tags$pre(class = "error-text", textOutput("gene_errors", inline = TRUE)))
          )
        ),
        conditionalPanel(condition = "input.workspace == 'cancer'",
          tags$div(class = "overview",
            tags$div(tags$div(class = "section-kicker", "COHORT METHYLATION"),
                     tags$h2(textOutput("cancer_headline", inline = TRUE)),
                     tags$p(class = "subheading", textOutput("cancer_subtitle", inline = TRUE))),
            tags$div(class = "source-tag", textOutput("cancer_source_tag", inline = TRUE))
          ),
          tags$div(class = "plot-card hero-card",
            tags$div(class = "card-head",
              tags$div(tags$div(class = "section-kicker", "TUMOR / NORMAL"),
                       tags$strong("Cohort landscapes")),
              tags$div(class = "card-actions",
                downloadLink("atlas_export", "↓ Export probe summaries", class = "action-link"),
                tags$span(" · "),
                downloadLink("atlas_delta_export", "↓ Export Δβ", class = "action-link"))
            ),
            uiOutput("atlas_status"),
            tags$div(class = "atlas-sequence-legend",
              tags$span(tags$i(style="background:#baff48;"), "A"),
              tags$span(tags$i(style="background:#e1ff98;"), "T"),
              tags$span(tags$i(style="background:#74e660;"), "C"),
              tags$span(tags$i(style="background:#a2f779;"), "G"),
              tags$span(tags$i(style="background:#ff3948;"), "Measured β / Δβ")),
            uiOutput("atlas_landscapes"),
            uiOutput("atlas_agreement"),
            tags$div(class = "plot-footnote",
              "Green A/T/C/G dots: the shared reference DNA sequence. Red: higher mean beta; ",
              "blue/red difference: Normal/Tumor comparison. Click a cell to inspect measured CpGs. ",
              "GRCh37 reference coordinates; no imputation.")
          ),
          uiOutput("atlas_cell_inspector"),
          tags$div(class = "detail-card",
            tags$div(class = "card-head", tags$strong("Data coverage")),
            tableOutput("atlas_coverage"),
            conditionalPanel(condition = "input.atlas_show_delta",
              tags$div(class = "card-head", tags$strong("Tumor − Normal β differences")),
              tableOutput("atlas_delta_table")),
            tags$details(class = "technical-details", tags$summary("Data sources and errors"),
              tableOutput("dataset_inventory"),
              tags$pre(class = "error-text", textOutput("atlas_errors", inline = TRUE)))
          )
        )
      ),
      tags$aside(class = "controls",
        conditionalPanel(condition = "input.workspace == 'gene'",
          tags$div(class = "section-kicker", "01 / GENE"),
          tags$h3("Sequence comparison"),
          selectInput("gene_preset", "Candidate genes",
                      choices = c("Custom symbol" = "custom", stats::setNames(COMMON_GENES, COMMON_GENES)),
                      selected = "TP53"),
          textInput("gene_symbol", "Human gene symbol", value = "TP53",
                    placeholder = "TP53, BRCA1, BRCA2 ..."),
          uiOutput("gene_input_hint"),
          radioButtons("gene_mode", NULL,
                       choices = c("Live sequence" = "live", "Offline example" = "demo"),
                       selected = "live", inline = TRUE),
          selectizeInput("gene_species", "Species", choices = SPECIES,
                         selected = c("homo_sapiens", "pan_troglodytes",
                                      "mus_musculus", "rattus_norvegicus"),
                         multiple = TRUE, options = list(maxItems = length(SPECIES))),
          selectInput("gene_window", "Window around 5′ gene boundary",
                      choices = WINDOW_OPTIONS, selected = 50000),
          tags$p(class = "helper", "100–250 kb shows more structure. 500 kb–1 Mb needs more data and memory."),
          actionButton("gene_load", "Load species", class = "btn-primary"),
          tags$div(class = "separator"),
          tags$details(class = "control-section", tags$summary("Display settings"),
          selectInput("gene_geometry", "Mode", choices = GEOMETRIES, selected = "square"),
          sliderInput("gene_methyl_threshold", "Methylated C: β threshold",
                      min = 0, max = 1, value = 0.6, step = 0.05),
          tags$p(class = "control-note",
                 "β ranges from 0 to 1. CpGs with measured β ≥ threshold are red; in ATG/C/mC mode they use the mC vertex. Reference DNA alone does not provide β."),
          checkboxInput("gene_show_cpg", "Show CpG motifs", TRUE),
          checkboxInput("gene_show_beta", "Show methylation points", TRUE),
          checkboxInput("gene_stars", "Star highlights", TRUE),
          sliderInput("gene_reveal", "Sequence shown (%)", min = 10, max = 100,
                      value = 100, step = 10)
          ),
          tags$div(class = "separator"),
          tags$div(class = "section-kicker", "03 / METHYLATION CSV"),
          fileInput("gene_file", "Upload CSV", accept = ".csv",
                    buttonLabel = "Choose CSV", placeholder = "No file selected"),
          selectInput("gene_sample", "Sample", choices = "sample_1"),
          tags$p(class = "helper",
                 "species, assembly, chrom, position, beta; optional: sample. Forward-reference CpG C coordinates."),
          downloadLink("template", "↓ Download CSV template", class = "quiet-link")
        ),
        conditionalPanel(condition = "input.workspace == 'cancer'",
          tags$div(class = "section-kicker", "CANCER DATA"),
          selectInput("cancer_type", "Cancer type", choices = TCGA_COHORTS, selected = "LUAD"),
          selectizeInput("cancer_gene_choice", "Gene", choices = catalog_gene_choices(),
                         selected = "TP53", options = list(create = TRUE, maxItems = 1)),
          uiOutput("cancer_input_hint"),
          selectizeInput("atlas_cohorts", "Cohorts", choices = catalog_cohorts("LUAD"),
                         selected = c("TCGA_LUAD", "GSE32861", "GSE62948", "GSE63384"),
                         multiple = TRUE),
          uiOutput("atlas_external_hint"),
          selectInput("atlas_window", "Reference window", choices = WINDOW_OPTIONS,
                      selected = 100000),
          actionButton("atlas_load", "Load cohort measurements", class = "btn-primary"),
          tags$p(class = "helper", "After the first successful download, the selected gene-window's measured CpGs and specimens are stored on this computer and reused offline."),
          selectInput("atlas_data_mode", "Data access",
                      choices = c("Saved first · download missing" = "cache",
                                  "Offline only · skip missing" = "offline",
                                  "Refresh selected datasets" = "refresh"),
                      selected = "cache"),
          textOutput("atlas_cache_status"),
          tags$div(class = "separator"),
          selectInput("atlas_geometry", "Geometry",
                      choices = c("Triangle · AT / C / G" = "triangle",
                                  "Square · A / C / G / T" = "square"),
                      selected = "triangle"),
          tags$p(class = "control-note",
                 "Red opacity follows the displayed mean β (lighter = lower, brighter = higher). There is no high-methylation cutoff. Green A/T/C/G dots are unfiltered reference DNA."),
          sliderInput("atlas_beta_range", "Mean β range (per Tumor / Normal panel)",
                      min = 0, max = 1, value = c(0, 1), step = 0.025),
          tags$p(class = "control-note", "Tumor and Normal panels each show CpGs whose OWN measured mean β is inside this range. For Δβ and cross-cohort agreement, include a probe when either group is in range, but require both groups to be measured. Green A/T/C/G reference points stay unchanged."),
          checkboxInput("atlas_show_delta", "Show Tumor − Normal Δβ panels", TRUE),
          sliderInput("atlas_delta_cutoff", "Minimum absolute Δβ per probe", min = 0,
                      max = 0.5, value = 0.15, step = 0.025),
          tags$p(class = "control-note",
                 "Δβ = mean β(Tumor) − mean β(Normal). Red: tumor higher; blue: tumor lower. This filters Δβ panels AND cross-cohort agreement, never mean β. Not a p-value."),
          selectInput("atlas_merge", "Combine cohorts",
                      choices = c("Separate cohorts" = "separate",
                                  "Also show equal-cohort pooled" = "equal",
                                  "Also show specimen-weighted pooled" = "weighted",
                                  "Also show shared-probe pooled" = "shared"),
                      selected = "separate"),
          tags$p(class = "control-note",
                 "Equal: studies weigh equally; weighted: sample counts; shared: only probes measured in both groups across all studies. Agreement uses the same probe IDs before spatial binning."),
          tags$details(class = "control-section", tags$summary("Advanced / import"),
            sliderInput("atlas_bins", "Spatial detail", min = 12, max = 48, value = 28, step = 4),
            sliderInput("atlas_min_samples", "Minimum specimens per probe/group",
                        min = 1, max = 15, value = 2, step = 1),
            sliderInput("atlas_min_bin_probes", "Minimum probes per spatial cell",
                        min = 1, max = 10, value = 1, step = 1),
            fileInput("atlas_files", "Cohort β CSV (optional)", accept = ".csv", multiple = TRUE),
            tags$p(class = "helper", "CSV requires cancer, cohort, sample, group, probe, beta. Unmapped probes also require assembly, chrom, position. Same cancer as selected is mandatory. Some public GEO datasets expose only raw intensities and require offline conversion."),
            downloadLink("atlas_csv_template", "↓ CSV template", class = "quiet-link")
          )
        )
      )
    ),
    source_footer()
  )
)

server <- function(input, output, session) {
  gene_records <- reactiveVal(make_demo_records())
  gene_status <- reactiveVal(NULL)
  gene_errors <- reactiveVal(character())
  gene_loaded_symbol <- reactiveVal("SYNTHETIC")
  cancer_result <- reactiveVal(NULL)
  cancer_status <- reactiveVal(NULL)
  cancer_errors <- reactiveVal(character())
  cancer_cache <- new.env(parent = emptyenv())
  probe_cache <- new.env(parent = emptyenv())
  selected_cpg <- reactiveVal(NULL)
  cohort_rows <- reactiveVal(list())
  cohort_failures <- reactiveVal(character())
  scan_meta <- reactiveVal(NULL)
  atlas_bundle <- reactiveVal(NULL)
  atlas_failures <- reactiveVal(character())
  atlas_notice <- reactiveVal("Choose a cancer, gene and measured datasets.")
  atlas_coord_cache <- new.env(parent = emptyenv())
  picked_atlas_cell <- reactiveVal(NULL)

  observeEvent(input$gene_preset, {
    if (!identical(input$gene_preset, "custom")) {
      updateTextInput(session, "gene_symbol", value = input$gene_preset)
    }
  })
  observeEvent(input$cancer_type, {
    choices <- catalog_cohorts(input$cancer_type)
    selected <- if (identical(input$cancer_type, "LUAD"))
      intersect(c("TCGA_LUAD", "GSE32861", "GSE62948", "GSE63384"), unname(choices))
      else head(unname(choices), 3L)
    updateSelectizeInput(session, "atlas_cohorts", choices = choices,
                         selected = selected, server = FALSE)
    # Never leave results from a different cancer under the new cancer's heading.
    atlas_bundle(NULL)
    picked_atlas_cell(NULL)
    cancer_result(NULL)
    atlas_notice(sprintf("Choose datasets for %s, then load measurements.", input$cancer_type))
    atlas_failures(character())
    cancer_status(NULL)
  }, ignoreInit = TRUE)

  output$atlas_external_hint <- renderUI({
    extra <- catalog_external_cohorts(input$cancer_type)
    if (!nrow(extra)) return(NULL)
    tags$details(class = "technical-details",
      tags$summary(sprintf("Other public studies (%d · CSV import)", nrow(extra))),
      tags$p(class = "helper", "Measured data are public, but these large GEO matrices are not imported automatically."),
      tags$div(class = "helper", lapply(seq_len(nrow(extra)), function(i) {
        tags$div(tags$a(href = extra$url[[i]], target = "_blank", extra$cohort[[i]]),
                 paste0(" · ", extra$note[[i]]))
      })))
  })

  output$dataset_inventory <- renderTable({
    b <- atlas_bundle()
    loaded <- if (!is.null(b) && identical(b$cancer, input$cancer_type)) names(b$cohorts) else character()
    counts <- if (length(loaded)) stats::setNames(
      vapply(b$cohorts, function(x) length(unique(x$measurements$probe)), integer(1)), loaded
    ) else integer()
    catalog_rows(input$cancer_type, loaded = loaded, loaded_probes = counts)
  }, striped = FALSE, bordered = FALSE, spacing = "s", rownames = FALSE)
  output$gene_input_hint <- renderUI({
    issue <- gene_input_issue(input$gene_symbol)
    if (is.null(issue)) return(NULL)
    tags$p(class = "input-warning", issue)
  })
  output$cancer_input_hint <- renderUI({
    issue <- gene_input_issue(input$cancer_gene_choice)
    if (is.null(issue)) return(NULL)
    tags$p(class = "input-warning", issue)
  })

  gene_uploaded <- reactive({
    if (is.null(input$gene_file)) return(list(data = NULL, error = NULL))
    tryCatch(list(data = read_methylation(input$gene_file$datapath), error = NULL),
             error = function(e) list(data = NULL, error = conditionMessage(e)))
  })
  observeEvent(gene_uploaded(), {
    up <- gene_uploaded()
    if (!is.null(up$error)) showNotification(up$error, type = "error", duration = 12)
    else if (!is.null(up$data)) updateSelectInput(session, "gene_sample",
      choices = unique(up$data$sample), selected = unique(up$data$sample)[[1L]])
  }, ignoreInit = TRUE)

  observeEvent(input$gene_load, {
    gene <- trimws(input$gene_symbol)
    if (identical(input$gene_mode, "demo")) {
      gene_records(make_demo_records())
      gene_status(NULL)
      gene_errors(character())
      gene_loaded_symbol("SYNTHETIC")
      return(invisible(NULL))
    }
    issue <- gene_input_issue(gene)
    if (!is.null(issue)) {
      gene_status(list(kind = "fail", text = issue))
      gene_errors(issue)
      return(invisible(NULL))
    }
    if (!length(input$gene_species)) {
      gene_status(list(kind = "fail", text = "Choose at least one species."))
      return(invisible(NULL))
    }
    result <- tryCatch(
      withProgress(message = sprintf("Loading %s", gene), value = 0,
        fetch_comparative_regions_fast(gene, input$gene_species,
                                       as.integer(input$gene_window), cache = sequence_cache)),
      error = function(e) list(records = list(), errors = conditionMessage(e)))
    gene_errors(result$errors)
    if (length(result$records)) {
      gene_records(result$records)
      gene_loaded_symbol(toupper(gene))
      count <- length(result$records)
      total <- length(unique(input$gene_species))
      gene_status(list(kind = if (count == total) "ok" else "partial",
                       text = sprintf("%s · %d/%d species loaded%s", gene, count, total,
                                      if (count < total) " · unavailable species omitted" else "")))
    } else {
      previous <- if (identical(gene_records()[[1L]]$source, "demo"))
        "The visible panels are synthetic examples." else
        paste("The visible panels still show", gene_loaded_symbol(), "from an earlier load.")
      gene_status(list(kind = "fail", text = paste("Sequence lookup failed.", previous)))
      showNotification("Sequence lookup failed. See lookup details.", type = "error", duration = 10)
    }
  })

  gene_data <- reactive({
    methylation <- gene_uploaded()$data
    lapply(gene_records(), function(rec) {
      sites <- methylation_sites_for_record(rec, methylation, input$gene_sample)
      threshold_for_geometry <- if (identical(input$gene_geometry, "mc_triangle"))
        input$gene_methyl_threshold else 0.6
      coords <- make_cgr(rec, input$gene_geometry,
                         methylation_sites = sites,
                         methylation_threshold = threshold_for_geometry)
      list(record = rec, coordinates = coords,
           methylation = match_methylation(rec, coords, methylation, input$gene_sample))
    })
  })
  output$gene_status <- renderUI(status_banner(gene_status()))
  output$gene_headline <- renderText({
    if (identical(gene_records()[[1L]]$source, "demo")) "Synthetic sequences" else
      sprintf("%s · %d species", gene_loaded_symbol(), length(gene_records()))
  })
  output$gene_subtitle <- renderText({
    if (identical(gene_records()[[1L]]$source, "demo")) "Offline example with synthetic DNA and β values" else
      sprintf("Reference DNA · %s bases per species", nfmt(nchar(gene_records()[[1L]]$sequence)))
  })
  output$gene_source_tag <- renderText({
    if (identical(gene_records()[[1L]]$source, "demo")) "SYNTHETIC" else
      paste(unique(vapply(gene_records(), `[[`, character(1), "source")), collapse = " / ")
  })
  output$gene_atlas <- renderUI({
    plot_grid(gene_data(), input$gene_geometry, input$gene_show_cpg,
              input$gene_show_beta, input$gene_methyl_threshold, input$gene_reveal,
              star_highlights = input$gene_stars)
  })
  output$gene_metrics <- renderUI({
    items <- gene_data()
    metric_values <- c(length(items), sum(vapply(items, function(x) nrow(x$coordinates), numeric(1))),
                       sum(vapply(items, function(x) sum(x$coordinates$is_cpg), numeric(1))),
                       sum(vapply(items, function(x) x$methylation$matched, numeric(1))))
    tags$div(class = "metrics",
      metric("PANELS", nfmt(metric_values[[1L]]), "species / demo"),
      metric("BASES", nfmt(metric_values[[2L]]), "processed"),
      metric("CpGs", nfmt(metric_values[[3L]]), "sequence motifs"),
      metric("β MATCHED", nfmt(metric_values[[4L]]), "measured sites"))
  })
  output$gene_coverage <- renderTable({
    do.call(rbind, lapply(gene_data(), function(x) {
      r <- x$record
      data.frame(Species = r$species_label, Gene = r$symbol, Assembly = r$assembly,
                 Source = r$source, Lookup = if (is.null(r$orthology)) "n/a" else r$orthology,
                 `Gene ID` = r$gene_id,
                 Bases = nfmt(nchar(r$sequence)), CpGs = sum(x$coordinates$is_cpg),
                 `GC %` = round(sequence_statistics(r)$gc, 1),
                 `β matched` = x$methylation$matched, check.names = FALSE)
    }))
  }, striped = FALSE, bordered = FALSE, spacing = "s", rownames = FALSE)
  output$gene_match_note <- renderText({
    if (identical(gene_records()[[1L]]$source, "demo"))
      return("Synthetic demonstration only. No patient or reference methylation was used.")
    up <- gene_uploaded()
    if (!is.null(up$error)) return(paste("CSV issue:", up$error))
    if (is.null(up$data)) return("No methylation CSV uploaded. CpG motif locations are not methylation measurements.")
    sprintf("Uploaded sample: %s. Matched CpGs require the same species, reference assembly, chromosome and coordinate.", input$gene_sample)
  })
  output$gene_errors <- renderText(paste(unique(gene_errors()), collapse = "\n"))

  observeEvent(input$cancer_load, {
    gene <- trimws(input$cancer_gene_choice)
    issue <- gene_input_issue(gene)
    if (!is.null(issue)) {
      cancer_status(list(kind = "fail", text = issue))
      cancer_errors(issue)
      return(invisible(NULL))
    }
    cache_key <- paste(input$cancer_type, toupper(gene), as.integer(input$cancer_window), sep = "|")
    result <- if (exists(cache_key, envir = cancer_cache, inherits = FALSE)) {
      get(cache_key, envir = cancer_cache, inherits = FALSE)
    } else tryCatch(withProgress(message = sprintf("Loading TCGA-%s", input$cancer_type), value = 0,
                   fetch_tcga(gene, cancer = input$cancer_type,
                              window_bp = as.integer(input$cancer_window))),
                  error = function(e) e)
    if (inherits(result, "error")) {
      cancer_status(list(kind = "fail", text = "TCGA request failed. No new sample data loaded."))
      cancer_errors(conditionMessage(result))
      return(invisible(NULL))
    }
    assign(cache_key, result, envir = cancer_cache)
    cancer_result(result)
    selected_cpg(NULL)
    cohort_rows(list())
    cohort_failures(character())
    scan_meta(NULL)
    cancer_status(list(kind = "ok", text = sprintf("TCGA-%s · %d tissues loaded",
                                                     result$cohort, nrow(result$samples))))
    cancer_errors(character())
  })

  cancer_data <- reactive({
    t <- cancer_result()
    if (is.null(t)) return(list())
    lapply(seq_len(nrow(t$samples)), function(i) {
      record <- t$record
      record$sample_id <- t$samples$sample[[i]]
      record$group <- t$samples$group[[i]]
      sites <- methylation_sites_for_record(record, t$methylation, record$sample_id)
      threshold_for_geometry <- if (identical(input$cancer_geometry, "mc_triangle"))
        input$cancer_methyl_threshold else 0.6
      coords <- make_cgr(record, input$cancer_geometry,
                         methylation_sites = sites,
                         methylation_threshold = threshold_for_geometry)
      list(record = record, coordinates = coords,
           methylation = match_methylation(record, coords, t$methylation, record$sample_id))
    })
  })
  cancer_delta <- reactive({
    t <- cancer_result()
    req(!is.null(t))
    tcga_delta(t$samples, t$methylation)
  })
  # Memoize reference coordinates. Adjusting the |Delta beta| display cutoff
  # must never recompute a 500 kb / 1 Mb chaos-game trajectory.
  cancer_reference <- reactive({
    t <- cancer_result()
    if (is.null(t)) return(NULL)
    list(record = t$record,
         coordinates = make_cgr(t$record, input$cancer_geometry),
         methylation = list(points = data.frame()))
  })
  output$cancer_headline <- renderText({
    b <- atlas_bundle()
    if (is.null(b)) sprintf("%s · cohort comparison", input$cancer_type) else
      sprintf("%s · %s cohorts", b$gene, b$cancer)
  })
  output$cancer_subtitle <- renderText({
    "Tumor / normal · measured array beta values · GRCh37"
  })
  output$cancer_source_tag <- renderText({
    if (is.null(atlas_bundle())) "NOT LOADED" else "MEASURED DATA"
  })
  output$cancer_status <- renderUI(status_banner(cancer_status()))
  output$cancer_atlas <- renderUI({
    if (is.null(cancer_result()))
      return(tags$p(class = "empty-message", "Choose a cohort and click Load tumor / normal."))
    plot_grid(cancer_data(), input$cancer_geometry, input$cancer_show_cpg,
              input$cancer_show_beta, input$cancer_methyl_threshold, input$cancer_reveal,
              star_highlights = input$cancer_stars)
  })
  output$cancer_delta_container <- renderUI({
    t <- cancer_result()
    if (is.null(t)) return(NULL)
    difference <- cancer_delta()
    reference <- cancer_reference()
    tags$div(class = "plot-card delta-card",
      tags$div(class = "card-head",
        tags$strong("Tumor − normal · Δβ"),
        downloadLink("delta_csv", "↓ Export Δβ", class = "action-link")),
      if (!nrow(difference$data)) {
        tags$p(class = "empty-message", "No shared measured CpG sites in this window.")
      } else {
        tags$div(class = "delta-layout",
          tags$div(class = "delta-figure",
            tags$div(class = "fractal-square",
              HTML(svg_panel(reference, input$cancer_geometry, FALSE, FALSE,
                             reveal_fraction = input$cancer_reveal / 100,
                             delta = difference$data, delta_cutoff = input$cancer_delta_cutoff)))),
          tags$div(class = "delta-summary",
            tags$p(class = "helper", sprintf("%s shared CpGs · |Δβ| ≥ %.2f · %s",
                                           nfmt(nrow(difference$data)), input$cancer_delta_cutoff,
                                           if (identical(difference$method, "patient-paired")) "Paired" else "Unpaired")),
            tags$p(class = "helper", "Red: β higher in tumor. Blue: β lower in tumor. No statistical significance implied."),
            tableOutput("delta_table"))
        )
      }
    )
  })
  output$delta_table <- renderTable({
    d <- cancer_delta()$data
    d <- d[abs(d$delta) >= input$cancer_delta_cutoff, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    d <- head(d[order(-abs(d$delta)), , drop = FALSE], 8L)
    data.frame(Position = nfmt(d$position), `Tumor β` = round(d$tumor_beta, 3),
               `Normal β` = round(d$normal_beta, 3), `Δβ` = round(d$delta, 3),
               check.names = FALSE)
  }, rownames = FALSE, striped = FALSE, bordered = FALSE, spacing = "s")
  output$cancer_metrics <- renderUI({
    t <- cancer_result()
    if (is.null(t)) return(NULL)
    d <- cancer_data()
    tags$div(class = "metrics",
      metric("TISSUES", nfmt(nrow(t$samples)), "actual TCGA samples"),
      metric("BASES", nfmt(nchar(t$record$sequence)), "reference window"),
      metric("MEASUREMENTS", nfmt(nrow(t$methylation)), "β rows"),
      metric("β MATCHED", nfmt(sum(vapply(d, function(x) x$methylation$matched, numeric(1)))), "CpG sites"))
  })
  output$cancer_coverage <- renderTable({
    t <- cancer_result()
    if (is.null(t)) return(NULL)
    do.call(rbind, lapply(cancer_data(), function(item) {
      rec <- item$record
      data.frame(Group = rec$group, Sample = rec$sample_id,
                 Gene = rec$symbol, Assembly = rec$assembly,
                 `β matched` = item$methylation$matched, check.names = FALSE)
    }))
  }, striped = FALSE, bordered = FALSE, spacing = "s", rownames = FALSE)
  output$cancer_match_note <- renderText({
    t <- cancer_result()
    if (is.null(t)) return("TCGA data has not been loaded. No demonstration β values are substituted.")
    sprintf("TCGA-%s · type 01 primary tumor; type 11 solid tissue normal. CpG positions checked against GRCh37. %s.",
            t$cohort, if (all(t$samples$paired)) "Paired tissues" else "Unpaired tissues")
  })
  output$cancer_errors <- renderText(paste(unique(cancer_errors()), collapse = "\n"))

  observeEvent(input$cpg_all_cohorts, {
    updateSelectizeInput(session, "cpg_cohorts", selected = unname(TCGA_COHORTS))
  })

  # Clicking a visible CpG is an exact genomic selection, never an image pixel
  # lookup. The current reference/probe map has already been checked against DNA.
  observeEvent(input$picked_cpg, {
    t <- cancer_result()
    if (is.null(t) || !identical(input$workspace, "cancer")) return(invisible(NULL))
    pos <- suppressWarnings(as.numeric(input$picked_cpg$position))
    probe <- cpg_probe_for_position(t$methylation, pos)
    if (is.null(probe)) {
      showNotification("Selected position has no validated HM450 probe identifier.",
                       type = "warning", duration = 8)
      return(invisible(NULL))
    }
    selected_cpg(list(position = pos, probe = probe, chrom = t$record$chromosome,
                      gene = t$record$symbol, assembly = t$record$assembly))
    cohort_rows(list())
    cohort_failures(character())
    scan_meta(NULL)
    if (!length(input$cpg_cohorts)) {
      showNotification("Choose at least one cohort on the right.", type = "warning")
      return(invisible(NULL))
    }
    scan_cohorts(probe, input$cpg_cohorts)
    session$onFlushed(function() {
      session$sendCustomMessage("focus-cpg", list())
    }, once = TRUE)
  })

  scan_cohorts <- function(probe, cohorts) {
    cohorts <- unique(intersect(as.character(cohorts), unname(TCGA_COHORTS)))
    if (!length(cohorts)) return(invisible(NULL))
    found <- list()
    errors <- character()
    withProgress(message = paste("CpG", probe, "· TCGA samples"), value = 0, {
      for (i in seq_along(cohorts)) {
        cohort <- cohorts[[i]]
        incProgress(1 / length(cohorts), detail = sprintf("%s · %d/%d", cohort, i, length(cohorts)))
        key <- paste(cohort, probe, sep = "|")
        result <- if (exists(key, envir = probe_cache, inherits = FALSE)) {
          get(key, envir = probe_cache, inherits = FALSE)
        } else {
          tryCatch(fetch_tcga_probe_cohort(probe, cohort), error = function(e) e)
        }
        if (inherits(result, "error")) {
          errors <- c(errors, sprintf("TCGA-%s: %s", cohort, conditionMessage(result)))
        } else {
          assign(key, result, envir = probe_cache)
          found[[cohort]] <- result
        }
      }
    })
    cohort_rows(found)
    cohort_failures(errors)
    scan_meta(list(probe = probe, cohorts = cohorts))
  }

  observeEvent(input$cpg_rescan, {
    selected <- selected_cpg()
    if (!is.null(selected)) scan_cohorts(selected$probe, input$cpg_cohorts)
  })

  output$cpg_side_status <- renderUI({
    picked <- selected_cpg()
    if (is.null(picked)) return(NULL)
    tags$div(class = "cpg-side-note", tags$strong(paste("Selected:", picked$probe)),
      tags$div(paste(picked$chrom, nfmt(picked$position), sep = ":")),
      actionButton("cpg_rescan", "Refresh selected cohorts", class = "btn-secondary"))
  })

  output$cpg_browser <- renderUI({
    t <- cancer_result()
    if (is.null(t)) return(NULL)
    picked <- selected_cpg()
    if (is.null(picked)) return(tags$div(class = "plot-card cpg-explorer",
      tags$div(class = "card-head", tags$strong("CpG cohort explorer")),
      tags$p(class = "empty-message", "Select a highlighted CpG in a tumor/normal panel or the Δβ map to inspect all assayed samples across TCGA cohorts.")))
    rows <- cohort_rows()
    failures <- cohort_failures()
    tags$div(class = "plot-card cpg-explorer", id = "cpg_browser_card",
      tags$div(class = "card-head", tags$div(tags$strong("CpG cohort explorer"),
        tags$div(class = "selected-probe", sprintf("%s · %s:%s · %s", picked$probe, picked$chrom,
                                                       nfmt(picked$position), picked$assembly))),
        if (length(rows)) downloadLink("cpg_download", "↓ Download all samples", class = "action-link")),
      tags$p(class = "helper", style = "padding:0 12px;",
        sprintf("%d of %d requested cohorts loaded · all eligible 01/11 samples per loaded cohort; no subsampling.",
                length(rows), if (is.null(scan_meta())) length(input$cpg_cohorts) else length(scan_meta()$cohorts))),
      if (length(rows)) tags$div(class = "cohort-results-grid",
        lapply(names(rows), function(cohort) {
          frame <- rows[[cohort]]
          stat <- cpg_cohort_summary(frame)
          tval <- if (is.finite(stat$tumor_median)) sprintf("%.3f", stat$tumor_median) else "n/a"
          nval <- if (is.finite(stat$normal_median)) sprintf("%.3f", stat$normal_median) else "n/a"
          delta <- if (is.finite(stat$mean_delta)) sprintf("%+.3f", stat$mean_delta) else "n/a"
          paired_delta <- if (is.finite(stat$paired_median_delta))
            sprintf("%+.3f", stat$paired_median_delta) else "n/a"
          tags$div(class = "cohort-result",
            tags$div(class = "cohort-title", tags$strong(paste("TCGA-", cohort, sep = "")),
                     tags$span(sprintf("%s paired patients", stat$paired_patients))),
            cpg_cohort_strip(frame),
            tags$div(class = "cohort-stat",
              sprintf("Tumor %s samples (%s patients) · median β %s", stat$tumor_samples, stat$tumor_patients, tval)),
            tags$div(class = "cohort-stat",
              sprintf("Normal %s samples (%s patients) · median β %s", stat$normal_samples, stat$normal_patients, nval)),
            tags$div(class = "cohort-stat", sprintf("Mean tumor − normal β: %s · paired median Δβ: %s · missing %s", delta, paired_delta, stat$missing_beta)))
        })
      ),
      if (length(failures)) tags$details(class = "technical-details", open = "open",
        tags$summary(sprintf("Unavailable cohorts (%d)", length(failures))),
        tags$pre(class = "error-text", paste(failures, collapse = "\n"))),
      tags$p(class = "cpg-caveat", "Points represent all nonmissing assayed samples in each loaded cohort. Tumor = primary tumor (01); Normal = solid tissue normal (11), not population healthy control. The mean difference is descriptive and not a significance test."
      )
    )
  })
  output$cpg_download <- downloadHandler(
    filename = function() {
      picked <- selected_cpg()
      req(!is.null(picked))
      sprintf("TCGA_%s_all_assayed_samples.csv", picked$probe)
    },
    content = function(file) {
      rows <- cohort_rows()
      req(length(rows) > 0L)
      utils::write.csv(do.call(rbind, unname(rows)), file, row.names = FALSE, na = "")
    }
  )


  # ---- Cancer-specific multicohort spatial atlas ----
  observeEvent(input$atlas_load, {
    picked_atlas_cell(NULL)
    symbol <- toupper(trimws(input$cancer_gene_choice))
    issue <- gene_input_issue(symbol)
    if (!is.null(issue)) {
      atlas_failures(issue)
      atlas_notice("Nothing loaded: enter an exact human gene symbol.")
      return(invisible(NULL))
    }
    cancer <- input$cancer_type
    selected <- unique(intersect(as.character(input$atlas_cohorts), unname(catalog_cohorts(cancer))))
    if (!length(selected) && (is.null(input$atlas_files) || !nrow(input$atlas_files))) {
      atlas_failures("Select a dataset for the selected cancer or upload a cohort CSV")
      return(invisible(NULL))
    }
    window_bp <- as.integer(input$atlas_window)
    offline_only <- identical(input$atlas_data_mode, "offline")
    force_refresh <- identical(input$atlas_data_mode, "refresh")
    base_key <- paste(cancer, symbol, window_bp, sep = "|")
    withProgress(message = paste(cancer, "cohorts ·", symbol), value = 0, {
      reference <- NULL
      probe_map <- NULL
      old <- isolate(atlas_bundle())
      if (!force_refresh && !is.null(old) && identical(old$key, base_key)) {
        reference <- old$reference
        probe_map <- old$map
      } else {
        cached <- if (force_refresh) NULL else
          atlas_reference_read(cancer, symbol, window_bp)
        if (!is.null(cached)) {
          reference <- cached$reference
          probe_map <- cached$map
        }
      }
      errors <- character()
      if (offline_only && (is.null(reference) || is.null(probe_map))) {
        atlas_failures("Reference + probe mapping not saved for this cancer / gene / window. Load it online once.")
        atlas_notice("Offline mode: no local reference. Previous results kept.")
        return(invisible(NULL))
      }
      if (is.null(reference) || is.null(probe_map)) {
        fetched <- tryCatch({
          rec <- fetch_grch37_human(symbol, window_bp = window_bp)
          map <- tryCatch(atlas_probe_map(symbol, rec, cancer = cancer), error = function(e) e)
          if (inherits(map, "error") && !is.null(input$atlas_files) && nrow(input$atlas_files)) {
            imported <- atlas_read_user_csv(input$atlas_files$datapath[[1L]], rec, cancer = cancer)
            map <- imported$map
          }
          if (inherits(map, "error")) stop(conditionMessage(map))
          list(record = rec, map = map)
        }, error = function(e) e)
        if (inherits(fetched, "error")) {
          atlas_failures(sprintf("Reference / probe-map request: %s", conditionMessage(fetched)))
          atlas_notice("Reference lookup failed; previous cohort atlas kept.")
          return(invisible(NULL))
        }
        reference <- fetched$record
        probe_map <- fetched$map
        atlas_reference_write(cancer, symbol, window_bp, reference, probe_map)
      }
      incProgress(0.10, detail = sprintf("Reference + %d validated probes", nrow(probe_map)))
      found <- list()
      reused <- 0L
      fetched_count <- 0L
      for (cohort in selected) {
        cached <- if (force_refresh) NULL else
          atlas_cached_read(cohort, symbol, window_bp)
        result <- if (!is.null(cached)) cached else if (offline_only) {
          simpleError(sprintf("%s is not saved locally; skipped without a network request", cohort))
        } else tryCatch(
          if (identical(cohort, paste0("TCGA_", cancer))) {
            atlas_fetch_tcga(symbol, reference, probe_map, cancer = cancer)
          } else {
            atlas_fetch_geo(cohort, probe_map)
          }, error = function(e) e)
        if (inherits(result, "error")) {
          errors <- c(errors, sprintf("%s: %s", cohort, conditionMessage(result)))
        } else {
          found[[cohort]] <- result
          if (is.null(cached)) {
            atlas_cached_write(result, symbol, window_bp)
            fetched_count <- fetched_count + 1L
          } else reused <- reused + 1L
        }
        incProgress(.90 / max(1L, length(selected)), detail = cohort)
      }
      if (!is.null(input$atlas_files) && nrow(input$atlas_files)) {
        for (f in input$atlas_files$datapath) {
          incoming <- tryCatch(atlas_read_user_csv(f, reference, probe_map, cancer = cancer), error = function(e) e)
          if (inherits(incoming, "error")) {
            errors <- c(errors, sprintf("Uploaded cohort: %s", conditionMessage(incoming)))
          } else {
            common_ids <- intersect(probe_map$probe, incoming$map$probe)
            check <- merge(probe_map[probe_map$probe %in% common_ids, , drop = FALSE],
                           incoming$map[incoming$map$probe %in% common_ids, , drop = FALSE],
                           by = "probe", suffixes = c("_reference", "_upload"))
            if (nrow(check) && any(check$position_reference != check$position_upload)) {
              errors <- c(errors, "Uploaded cohort: probe positions conflict with the validated GRCh37 map")
            } else {
              for (id in names(incoming$cohorts)) found[[id]] <- incoming$cohorts[[id]]
              probe_map <- rbind(probe_map,
                                 incoming$map[!incoming$map$probe %in% probe_map$probe,
                                              , drop = FALSE])
            }
          }
        }
      }
      atlas_failures(errors)
      if (!length(found)) {
        atlas_bundle(NULL)
        atlas_notice(sprintf("No measured %s cohorts loaded. No synthetic results are shown.", cancer))
      } else {
        atlas_bundle(list(key = base_key, gene = symbol, cancer = cancer, reference = reference,
                          map = probe_map, cohorts = found, requested = selected))
        atlas_notice(sprintf("%d %s cohorts loaded · %d reused locally · %d downloaded · %d uploaded CSV. Measured β only.",
                             length(found), cancer, reused, fetched_count,
                             if (is.null(input$atlas_files)) 0L else nrow(input$atlas_files)))
        rm(list = ls(atlas_coord_cache), envir = atlas_coord_cache)
      }
    })
  })

  output$atlas_cache_status <- renderText({
    # This is the cache on the user's own machine, not a repository download.
    invisible(atlas_notice())
    info <- atlas_cache_usage()
    sprintf("Saved data: %d RDS file(s), %s · %s",
            info$n_files, sprintf("%.2f MB", info$bytes / 1024^2),
            info$path)
  })

  output$atlas_status <- renderUI({
    b <- atlas_bundle()
    tags$div(class = "load-banner", style = "background:#071800;color:#c7f586;",
      atlas_notice(),
      if (!is.null(b)) tags$span(sprintf("  ·  %s · %s validated probes", b$gene,
                                          nfmt(nrow(b$map))))
    )
  })

  atlas_summaries <- reactive({
    b <- atlas_bundle()
    if (is.null(b)) return(list())
    stats::setNames(lapply(b$cohorts, function(x) atlas_summarize(x, b$map,
                                                      threshold = 0.5)), names(b$cohorts))
  })

  atlas_visible_summaries <- reactive({
    raw <- atlas_summaries()
    if (!length(raw)) return(raw)
    beta_limits <- input$atlas_beta_range
    required_n <- input$atlas_min_samples
    stats::setNames(lapply(raw, function(s) atlas_filter_beta_range(
      s, range = beta_limits, min_samples = required_n)), names(raw))
  })

  # Pool UNFILTERED measured cohort summaries first, then apply the beta range
  # to each pooled Tumor / Normal value. This avoids averaging only preselected
  # studies and keeps the pooled estimate interpretable.
  atlas_pooled_summaries <- reactive({
    raw <- atlas_summaries()
    mode <- input$atlas_merge
    if (length(raw) < 2L || is.null(mode) || identical(mode, "separate")) return(NULL)
    if ("GSE32861" %in% names(raw) && "GSE32866" %in% names(raw))
      raw[["GSE32866"]] <- NULL
    if (length(raw) < 2L) return(NULL)
    min_n <- input$atlas_min_samples
    merged_mean <- atlas_merge_summaries(raw, mode = mode, min_samples = min_n)
    pair_rows <- lapply(raw, atlas_paired_summary, min_samples = min_n)
    merged_pair <- atlas_merge_summaries(pair_rows, mode = mode,
                                         min_samples = min_n)
    list(mean = if (nrow(merged_mean))
           atlas_filter_beta_range(merged_mean, input$atlas_beta_range, min_n)
         else merged_mean,
         paired = merged_pair, n_studies = length(raw))
  })

  atlas_core_probes <- reactive({
    z <- atlas_visible_summaries()
    if (!length(z)) return(character())
    Reduce(intersect, lapply(z, function(d) unique(d$probe)))
  })

  atlas_coords <- reactive({
    b <- atlas_bundle()
    if (is.null(b)) return(NULL)
    key <- paste(b$key, input$atlas_geometry, sep = "|")
    if (!exists(key, envir = atlas_coord_cache, inherits = FALSE)) {
      assign(key, make_cgr(b$reference, input$atlas_geometry), envir = atlas_coord_cache)
    }
    get(key, envir = atlas_coord_cache, inherits = FALSE)
  })

  atlas_panel <- function(cohort, group, summary, reference, coordinates, score,
                          geometry, threshold, shared_probes, bins, min_n,
                          background, style, counts) {
    id <- cohort
    d <- summary
    if (length(shared_probes)) d <- d[d$probe %in% shared_probes, , drop = FALSE]
    if (identical(group, "Difference")) {
      z <- atlas_probe_deltas(d)
      if (nrow(z)) z <- z[z$n_tumor >= min_n & z$n_normal >= min_n, , drop = FALSE]
    } else {
      z <- d[d$group == group & d$n >= min_n, , drop = FALSE]
    }
    coords <- coordinates
    if (identical(style, "expected") && identical(geometry, "mc_triangle") &&
        !identical(group, "Difference") && nrow(z)) {
      profile <- stats::setNames(z$mean_beta, as.character(z$position))
      coords <- tryCatch(atlas_expected_cgr(reference, profile), error = function(e) coordinates)
    }
    visual <- atlas_svg(coords, z, geometry,
                        score = score, delta = identical(group, "Difference"),
                        min_samples = min_n, bins = bins,
                        show_background = background,
                        min_bin_probes = input$atlas_min_bin_probes)
    n_assays <- if (group == "Tumor") counts$tumor[[1L]] else
      if (group == "Normal") counts$normal[[1L]] else min(counts$tumor[[1L]], counts$normal[[1L]])
    tags$div(class = "fractal-panel atlas-panel",
      tags$div(class = "fractal-label",
        tags$strong(paste(id, "·", group)),
        tags$span(sprintf("%s probes · %s assays", nrow(z), nfmt(n_assays)))),
      tags$div(class = "fractal-square", htmltools::HTML(visual))
    )
  }

  output$atlas_landscapes <- renderUI({
    b <- atlas_bundle()
    if (is.null(b)) return(tags$p(class = "empty-message", "Choose a cancer, cohorts and gene, then load measured data."))
    summaries <- atlas_visible_summaries()
    raw_summaries <- atlas_summaries()
    geom <- input$atlas_geometry
    score <- "mean_beta"
    coords <- atlas_coords()
    mode <- input$atlas_merge
    if (is.null(mode)) mode <- "separate"
    if (!mode %in% c("separate", "equal", "weighted", "shared")) mode <- "separate"
    min_n <- input$atlas_min_samples
    panel_for <- function(id, group, smry, counts, pooled = FALSE) {
      z <- smry[smry$group == group & smry$n >= min_n, , drop = FALSE]
      svg <- atlas_svg(coords, z, geom, score = score,
                       min_samples = min_n, bins = input$atlas_bins,
                       min_bin_probes = input$atlas_min_bin_probes,
                       show_background = TRUE, click_cohort = id,
                       click_group = group)
      n_assays <- if (group == "Tumor") counts$tumor[[1L]] else counts$normal[[1L]]
      tags$div(class = "fractal-panel atlas-panel",
        tags$div(class = "fractal-label", tags$strong(paste(id, "·", group)),
          tags$span(sprintf("%s probes · %s assays", nrow(z), nfmt(n_assays)))),
        tags$div(class = "fractal-square", htmltools::HTML(svg)))
    }
    study_panels <- unlist(lapply(names(summaries), function(id) {
      counts <- atlas_sample_counts(b$cohorts[[id]])
      list(panel_for(id, "Tumor", summaries[[id]], counts),
           panel_for(id, "Normal", summaries[[id]], counts))
    }), recursive = FALSE)
    # Separate row of one delta map per cohort; keeps the normal/tumor pairs
    # contiguous in the first grid and avoids pretending each point is a patient.
    diff_panels <- NULL
    if (isTRUE(input$atlas_show_delta)) {
      diff_panels <- tags$div(class = "cohort-atlas-stack",
        tags$div(class = "section-kicker",
                 sprintf("Tumor − Normal · |Δβ| ≥ %.3f", input$atlas_delta_cutoff)),
        tags$div(class = "cohort-atlas-row",
          lapply(names(summaries), function(id) {
            d <- atlas_filter_delta_probes(raw_summaries[[id]], min_samples = min_n,
                                           cutoff = input$atlas_delta_cutoff,
                                           beta_range = input$atlas_beta_range)
            svg <- atlas_svg(coords, d, geom, delta = TRUE,
                             min_samples = min_n, bins = input$atlas_bins,
                             min_bin_probes = input$atlas_min_bin_probes,
                             show_background = TRUE, click_cohort = id,
                             click_group = "Difference")
            tags$div(class = "fractal-panel atlas-panel",
              tags$div(class = "fractal-label", tags$strong(paste(id, "· Δβ")),
                       tags$span(sprintf("%s passing CpGs · descriptive", nrow(d)))),
              tags$div(class = "fractal-square", htmltools::HTML(svg)))
          })))
    }
    pooled <- NULL
    pool_result <- atlas_pooled_summaries()
    if (!is.null(pool_result)) {
      pooled_summary <- pool_result$mean
      paired_pooled_summary <- pool_result$paired
      if (nrow(pooled_summary) || nrow(paired_pooled_summary)) {
        pooled <- tags$div(class = "cohort-atlas-stack",
          tags$div(class = "section-kicker", sprintf("Combined · %s · %d studies",
                           mode, pool_result$n_studies)),
          if ("GSE32861" %in% names(summaries) && "GSE32866" %in% names(summaries))
            tags$p(class = "helper", "GSE32866 excluded from pooled summary to avoid double-counting the GSE32861 parent study."),
          tags$div(class = "cohort-atlas-row",
            lapply(c("Tumor", "Normal"), function(group) {
              z <- pooled_summary[pooled_summary$group == group, , drop = FALSE]
              svg <- atlas_svg(coords, z, geom, score = score,
                               min_samples = 1L, bins = input$atlas_bins,
                               min_bin_probes = input$atlas_min_bin_probes,
                               show_background = TRUE, click_cohort = "Pooled",
                               click_group = group)
              tags$div(class = "fractal-panel atlas-panel",
                tags$div(class = "fractal-label", tags$strong(paste("Pooled ·", group)),
                  tags$span(sprintf("%d probes · %d contributing studies", nrow(z), pool_result$n_studies))),
                tags$div(class = "fractal-square", htmltools::HTML(svg)))
            })
          ),
          if (isTRUE(input$atlas_show_delta)) {
            d <- atlas_filter_delta_probes(paired_pooled_summary,
                                           min_samples = min_n,
                                           cutoff = input$atlas_delta_cutoff,
                                           beta_range = input$atlas_beta_range)
            tags$div(class = "cohort-atlas-row",
              tags$div(class = "fractal-panel atlas-panel",
                tags$div(class = "fractal-label", tags$strong("Pooled · Δβ"),
                         tags$span(sprintf("%s passing CpGs", nrow(d)))),
                tags$div(class = "fractal-square", htmltools::HTML(
                  atlas_svg(coords, d, geom, delta = TRUE,
                            min_samples = min_n, bins = input$atlas_bins,
                            min_bin_probes = input$atlas_min_bin_probes,
                            show_background = TRUE, click_cohort = "Pooled",
                            click_group = "Difference")))))
          }
        )
      }
    }
    tags$div(class = "cohort-atlas-stack",
      tags$div(class = "section-kicker", sprintf("%s · %d loaded cohorts · Mean β %.3f–%.3f",
                 b$cancer, length(summaries), input$atlas_beta_range[[1L]],
                 input$atlas_beta_range[[2L]])),
      tags$div(class = "cohort-atlas-row", study_panels),
      diff_panels,
      if (mode != "separate" && length(summaries) < 2L)
        tags$p(class = "input-warning", "Pooling requires at least two loaded cohorts."),
      if (mode != "separate" && length(summaries) >= 2L && is.null(pooled))
        tags$p(class = "input-warning", "No probes satisfy this pooling mode."),
      pooled
    )
  })


  # Colored CGR cells in the cohort atlas are selectable again.
  # Drill-down uses the already loaded measurements for the SAME cancer.
  observeEvent(input$picked_atlas_cell, {
    b <- atlas_bundle()
    event <- input$picked_atlas_cell
    if (is.null(b) || !identical(input$workspace, "cancer")) return(invisible(NULL))
    ix <- suppressWarnings(as.integer(event$ix))
    iy <- suppressWarnings(as.integer(event$iy))
    bins <- suppressWarnings(as.integer(event$bins))
    cohort <- as.character(event$cohort)
    group <- as.character(event$group)
    if (length(ix) != 1L || length(iy) != 1L || length(bins) != 1L ||
        anyNA(c(ix, iy, bins)) || bins != input$atlas_bins ||
        ix < 0L || iy < 0L || ix >= bins || iy >= bins ||
        length(cohort) != 1L || length(group) != 1L ||
        !cohort %in% c(names(b$cohorts), "Pooled", "Agreement") ||
        !group %in% c("Tumor", "Normal", "Difference")) return(invisible(NULL))
    if (identical(cohort, "Pooled") &&
        (identical(input$atlas_merge, "separate") || length(b$cohorts) < 2L))
      return(invisible(NULL))
    if (identical(cohort, "Agreement") && length(b$cohorts) < 2L)
      return(invisible(NULL))
    picked_atlas_cell(list(ix = ix, iy = iy, bins = bins, cohort = cohort,
                           group = group, key = b$key))
    session$onFlushed(function() session$sendCustomMessage("focus-atlas-cell", list()),
                      once = TRUE)
  })

  # Cutoff/bin/geometry changes redefine the selected area; avoid stale results.
  observeEvent(list(input$atlas_geometry, input$atlas_bins,
                    input$atlas_min_samples, input$atlas_min_bin_probes,
                    input$atlas_delta_cutoff, input$atlas_merge,
                    input$atlas_beta_range), {
    picked_atlas_cell(NULL)
  }, ignoreInit = TRUE)

  atlas_cell_details <- reactive({
    b <- atlas_bundle()
    picked <- picked_atlas_cell()
    if (is.null(b) || is.null(picked) || !identical(picked$key, b$key)) return(NULL)
    coords <- atlas_coords()
    map <- atlas_cell_probe_map(b$map, coords, picked$ix, picked$iy, picked$bins)
    if (!nrow(map)) return(list(picked = picked, map = map, summaries = list(),
                                table = data.frame(), samples = data.frame()))
    summaries <- atlas_visible_summaries()
    # Use the same source as the clicked image, not a newly filtered proxy.
    source_summary <- if (identical(picked$cohort, "Agreement")) {
      NULL
    } else if (identical(picked$cohort, "Pooled")) {
      pooled <- atlas_pooled_summaries()
      if (is.null(pooled)) NULL else
        if (identical(picked$group, "Difference")) pooled$paired else pooled$mean
    } else if (identical(picked$group, "Difference")) {
      atlas_summaries()[[picked$cohort]]
    } else summaries[[picked$cohort]]
    if (identical(picked$cohort, "Agreement")) {
      agrees <- atlas_agreement_data()$probes
      if (!nrow(agrees)) return(NULL)
      source_rows <- agrees[agrees$ix == picked$ix & agrees$iy == picked$iy, , drop = FALSE]
    } else {
      if (is.null(source_summary) || !nrow(source_summary)) return(NULL)
      source_rows <- if (identical(picked$group, "Difference")) {
      atlas_filter_delta_probes(source_summary, input$atlas_min_samples,
                                input$atlas_delta_cutoff,
                                beta_range = input$atlas_beta_range)
    } else {
      min_n <- if (identical(picked$cohort, "Pooled")) 1L else input$atlas_min_samples
      source_summary[source_summary$group == picked$group &
                       source_summary$n >= min_n, , drop = FALSE]
      }
    }
    map <- map[map$probe %in% source_rows$probe, , drop = FALSE]
    if (!nrow(map)) return(NULL)
    # The zoom panels honor the per-group beta range; the numerical detail
    # table retains BOTH original measured group means for the selected CpG.
    # This prevents a real delta from becoming NA when only one mean qualifies.
    list(picked = picked, map = map, summaries = summaries,
         table = atlas_cell_probe_table(atlas_summaries(), map$probe),
         samples = atlas_cell_sample_rows(b$cohorts, map$probe))
  })

  output$atlas_cell_inspector <- renderUI({
    b <- atlas_bundle()
    if (is.null(b)) return(NULL)
    d <- atlas_cell_details()
    if (is.null(d)) return(tags$div(class = "detail-card cell-help",
      tags$div(class = "card-head", tags$strong("Pattern details")),
      tags$p(class = "helper", "Click any colored cell in the cohort landscape to inspect its measured CpGs, beta levels and tumor–normal differences across the loaded cohorts.")))
    picked <- d$picked
    tags$div(class = "detail-card atlas-cell-inspector", id = "atlas_cell_inspector",
      tags$div(class = "card-head",
        tags$div(tags$strong(sprintf("Pattern details · %s / %s", picked$cohort, picked$group)),
                 tags$div(class = "selected-probe",
                          sprintf("Cell (%d, %d) · %d bins · %d measured CpGs",
                                  picked$ix + 1L, picked$iy + 1L, picked$bins, nrow(d$map)))),
        downloadLink("atlas_cell_download", "↓ Export specimens", class = "action-link")),
      tags$p(class = "helper", style = "padding:0 12px;",
             "The selected cell is compared across cohorts. Each Tumor / Normal panel filters its own mean β. Difference and agreement panels require both tissue groups; either may qualify for the β range. Nearby CGR points need not represent the same genomic locus."),
      tags$div(class = "atlas-cell-cohorts", lapply(names(d$summaries), function(id) {
        s <- d$summaries[[id]]
        tags$div(class = "atlas-cell-cohort",
          tags$div(class = "atlas-cell-cohort-label", id),
          tags$div(class = "atlas-cell-mini-pair",
            lapply(c("Tumor", "Normal"), function(group) {
              z <- s[s$group == group & s$probe %in% d$map$probe, , drop = FALSE]
              tags$div(class = "atlas-cell-mini",
                tags$div(class = "atlas-cell-mini-title",
                         sprintf("%s · %d CpGs", group, nrow(z))),
                htmltools::HTML(atlas_cell_zoom_svg(d$map, z, picked)))
            }))
        )
      })),
      tags$div(class = "card-head", tags$strong("All specimens at one CpG")),
      tags$div(class = "atlas-probe-controls",
        selectizeInput("atlas_focus_probe", "Measured CpG in selected cell",
                       choices = sort(unique(d$map$probe)),
                       selected = sort(unique(d$map$probe))[[1L]],
                       options = list(maxItems = 1)),
        downloadLink("atlas_probe_download", "↓ Export selected CpG", class = "action-link")),
      uiOutput("atlas_probe_distributions"),
      tags$div(class = "card-head", tags$strong("Measured CpGs (largest |Δβ| first)")),
      tableOutput("atlas_cell_table")
    )
  })


  output$atlas_probe_distributions <- renderUI({
    d <- atlas_cell_details()
    if (is.null(d)) return(NULL)
    probe <- input$atlas_focus_probe
    if (length(probe) != 1L || is.na(probe) || !probe %in% d$map$probe) return(NULL)
    rows <- d$samples[d$samples$probe == probe, , drop = FALSE]
    if (!nrow(rows)) return(tags$p(class = "helper", "No measured specimens for this CpG."))
    tags$div(class = "atlas-probe-strips", lapply(unique(rows$cohort), function(id) {
      x <- rows[rows$cohort == id, , drop = FALSE]
      tags$div(class = "atlas-probe-strip",
        tags$div(class = "atlas-cell-cohort-label",
                 sprintf("%s · %d Tumor / %d Normal specimens", id,
                         sum(x$group == "Tumor"), sum(x$group == "Normal"))),
        cpg_cohort_strip(x))
    }))
  })

  output$atlas_probe_download <- downloadHandler(
    filename = function() {
      d <- atlas_cell_details()
      req(!is.null(d), length(input$atlas_focus_probe) == 1L)
      sprintf("%s_%s_specimens.csv", input$cancer_type, input$atlas_focus_probe)
    },
    content = function(file) {
      d <- atlas_cell_details()
      req(!is.null(d), length(input$atlas_focus_probe) == 1L,
          input$atlas_focus_probe %in% d$map$probe)
      rows <- d$samples[d$samples$probe == input$atlas_focus_probe, , drop = FALSE]
      utils::write.csv(rows, file, row.names = FALSE, na = "")
    }
  )

  output$atlas_cell_table <- renderTable({
    d <- atlas_cell_details()
    if (is.null(d) || !nrow(d$table)) return(NULL)
    x <- d$table
    x <- x[order(-abs(x$delta), na.last = TRUE), , drop = FALSE]
    x <- head(x, 24L)
    data.frame(Cohort = x$cohort, CpG = x$probe, Position = x$position,
               `Tumor β` = round(x$tumor_beta, 3), `Normal β` = round(x$normal_beta, 3),
               `Δβ` = round(x$delta, 3), `Tumor n` = x$tumor_n,
               `Normal n` = x$normal_n, check.names = FALSE)
  }, rownames = FALSE, striped = FALSE, bordered = FALSE, spacing = "s")

  output$atlas_cell_download <- downloadHandler(
    filename = function() {
      d <- atlas_cell_details()
      req(!is.null(d))
      sprintf("%s_%s_cell_%d_%d_samples.csv", input$cancer_type,
              gsub("[^A-Za-z0-9_-]", "", input$cancer_gene_choice),
              d$picked$ix + 1L, d$picked$iy + 1L)
    },
    content = function(file) {
      d <- atlas_cell_details()
      req(!is.null(d), nrow(d$samples) > 0L)
      utils::write.csv(d$samples, file, row.names = FALSE, na = "")
    }
  )

  # Same-probe comparison first, spatial cell aggregation second.
  atlas_agreement_data <- reactive({
    b <- atlas_bundle()
    if (is.null(b) || length(b$cohorts) < 2L)
      return(list(cells = data.frame(), probes = data.frame(), cohorts = character()))
    summaries <- atlas_summaries()
    # Subseries of the same Selamat study are not two independent studies.
    if ("GSE32861" %in% names(summaries) && "GSE32866" %in% names(summaries))
      summaries[["GSE32866"]] <- NULL
    atlas_spatial_concordance(summaries, atlas_coords(),
                              bins = input$atlas_bins,
                              min_samples = input$atlas_min_samples,
                              delta_cutoff = input$atlas_delta_cutoff,
                              beta_range = input$atlas_beta_range)
  })

  output$atlas_agreement <- renderUI({
    b <- atlas_bundle()
    if (is.null(b)) return(NULL)
    if (length(b$cohorts) < 2L)
      return(tags$p(class = "helper", style = "padding:9px;",
                    "Cross-cohort agreement needs at least two successfully loaded independent studies."))
    z <- atlas_agreement_data()
    cells <- z$cells
    palette <- atlas_cohort_colors(z$cohorts)
    legend <- tags$div(class = "atlas-cohort-legend",
      lapply(names(palette), function(id) tags$span(
        tags$i(style = paste0("background:", palette[[id]], ";")), id)))
    if (!nrow(cells)) return(tags$div(class = "detail-card",
      tags$div(class = "card-head", tags$strong("Cross-cohort pattern agreement")),
      tags$p(class = "helper", style = "padding:14px;",
             "No matching CpG IDs passed the minimum sample count and Δβ cutoff in two independent studies."),
      legend))
    cells <- cells[cells$n_probes >= input$atlas_min_bin_probes, , drop = FALSE]
    svg <- atlas_concordance_svg(atlas_coords(), input$atlas_geometry, z,
                                 bins = input$atlas_bins,
                                 min_bin_probes = input$atlas_min_bin_probes)
    tags$div(class = "plot-card atlas-agreement-card",
      tags$div(class = "card-head", tags$strong("Cross-cohort pattern agreement"),
               tags$span(sprintf("%d spatial cells", nrow(cells))),
               downloadLink("atlas_agreement_export", "↓ Export agreement", class = "action-link")),
      tags$div(class = "agreement-layout",
        tags$div(class = "fractal-panel agreement-panel",
          tags$div(class = "fractal-square", htmltools::HTML(svg))),
        tags$div(class = "agreement-description",
          tags$p("Only CpGs measured in ≥2 independent cohorts enter this map. Spatial bins then show whether cohort-level mean Δβ directions agree."),
          tags$p("Red: unanimous tumor hypermethylation; blue: unanimous hypomethylation; amber: discordant studies. Colored radial ticks identify contributing cohorts; dim ticks disagree with the majority."),
          legend,
          tableOutput("atlas_agreement_table")))
    )
  })

  output$atlas_agreement_table <- renderTable({
    z <- atlas_agreement_data()$cells
    if (!nrow(z)) return(NULL)
    z <- z[order(-z$fraction, -z$n_cohorts, -z$n_probes), , drop = FALSE]
    head(data.frame(Cell = paste(z$ix + 1L, z$iy + 1L, sep = ","),
                    `Same-probe CpGs` = z$n_probes, `Cohorts` = z$n_cohorts,
                    `Agreeing` = z$agreeing,
                    `Agreement %` = round(z$fraction * 100),
                    `Median Δβ` = round(z$median_delta, 3),
                    check.names = FALSE), 12L)
  }, striped = FALSE, bordered = FALSE, spacing = "s", rownames = FALSE)

  output$atlas_agreement_export <- downloadHandler(
    filename = function() {
      b <- atlas_bundle(); req(!is.null(b))
      sprintf("%s_%s_spatial_agreement.csv", b$cancer, b$gene)
    },
    content = function(file) {
      z <- atlas_agreement_data()$cells
      utils::write.csv(z, file, row.names = FALSE, na = "")
    })

  output$atlas_top_probes <- renderTable({
    b <- atlas_bundle()
    req(!is.null(b))
    summary <- lapply(atlas_summaries(), atlas_filter_beta_pairs,
                      range = input$atlas_beta_range,
                      min_samples = input$atlas_min_samples)
    if ("GSE32861" %in% names(summary) && "GSE32866" %in% names(summary))
      summary[["GSE32866"]] <- NULL
    p <- atlas_probe_agreement(summary, b$map, 2L)
    if (!nrow(p)) return(NULL)
    p <- p[order(-p$concordance, -abs(p$delta)), , drop = FALSE]
    head(p, 8L)
  }, rownames = FALSE, striped = FALSE)

  output$atlas_delta_table <- renderTable({
    b <- atlas_bundle()
    if (is.null(b) || !isTRUE(input$atlas_show_delta)) return(NULL)
    d <- atlas_summaries()
    do.call(rbind, lapply(names(d), function(cohort) {
      complete <- atlas_filter_delta_probes(d[[cohort]], input$atlas_min_samples, 0,
                                            beta_range = input$atlas_beta_range)
      passing <- atlas_filter_delta_probes(d[[cohort]], input$atlas_min_samples,
                                           input$atlas_delta_cutoff,
                                           beta_range = input$atlas_beta_range)
      data.frame(Cohort = cohort,
                 `Compared CpGs` = nrow(complete),
                 `Passing |Δβ|` = nrow(passing),
                 `Median Δβ (passing)` = if (nrow(passing)) round(stats::median(passing$delta), 3) else NA_real_,
                 check.names = FALSE)
    }))
  }, rownames = FALSE, striped = FALSE, bordered = FALSE)

  output$atlas_delta_export <- downloadHandler(
    filename = function() {
      b <- atlas_bundle()
      req(!is.null(b))
      sprintf("%s_%s_delta_beta.csv", b$cancer, b$gene)
    },
    content = function(file) {
      d <- atlas_summaries()
      req(length(d))
      frames <- Filter(Negate(is.null), lapply(names(d), function(cohort) {
        z <- atlas_filter_delta_probes(d[[cohort]], input$atlas_min_samples,
                                       input$atlas_delta_cutoff,
                                       beta_range = input$atlas_beta_range)
        if (!nrow(z)) return(NULL)
        z$cohort <- rep(cohort, nrow(z))
        z
      }))
      if (!length(frames)) {
        utils::write.csv(data.frame(cohort = character(), probe = character(),
                                    position = numeric(), delta = numeric()),
                         file, row.names = FALSE)
      } else {
        utils::write.csv(do.call(rbind, frames), file, row.names = FALSE, na = "")
      }
    }
  )

  output$atlas_coverage <- renderTable({
    b <- atlas_bundle()
    if (is.null(b)) return(NULL)
    visible <- atlas_visible_summaries()
    do.call(rbind, lapply(names(b$cohorts), function(id) {
      x <- atlas_sample_counts(b$cohorts[[id]])
      x$selected_cpgs <- length(unique(visible[[id]]$probe))
      x
    }))
  }, striped = FALSE, bordered = FALSE, rownames = FALSE)

  output$atlas_errors <- renderText({
    if (!length(atlas_failures())) "No recorded source errors" else
      paste(atlas_failures(), collapse = "\n")
  })

  output$atlas_csv_template <- downloadHandler(
    filename = function() "cohort_long_template.csv",
    content = function(file) file.copy("data/cohort_long_template.csv", file, overwrite = TRUE)
  )

  output$atlas_export <- downloadHandler(
    filename = function() {
      b <- atlas_bundle()
      req(!is.null(b))
      sprintf("%s_cohort_atlas_%s.csv", b$cancer, b$gene)
    },
    content = function(file) {
      s <- atlas_visible_summaries()
      req(length(s))
      x <- do.call(rbind, lapply(names(s), function(id) {
        d <- s[[id]]
        d$cohort <- id
        d
      }))
      utils::write.csv(x, file, row.names = FALSE)
    }
  )

  output$gene_csv <- downloadHandler(
    filename = function() sprintf("gene_fractal_%s_%s.csv", gene_loaded_symbol(), input$gene_geometry),
    content = function(file) {
      frames <- lapply(gene_data(), function(item) {
        r <- item$record
        d <- item$coordinates
        d$beta <- NA_real_
        if (item$methylation$matched) {
          ix <- match(item$methylation$points$index, d$index)
          d$beta[ix] <- item$methylation$points$beta
        }
        data.frame(source = r$source, species = r$species, gene = r$symbol,
                   assembly = r$assembly, chromosome = r$chromosome, d,
                   stringsAsFactors = FALSE)
      })
      utils::write.csv(do.call(rbind, frames), file, row.names = FALSE, na = "")
    }
  )
  output$cancer_csv <- downloadHandler(
    filename = function() {
      t <- cancer_result()
      req(!is.null(t))
      sprintf("TCGA_%s_%s_samples.csv", t$cohort, t$record$symbol)
    },
    content = function(file) {
      req(!is.null(cancer_result()))
      frames <- lapply(cancer_data(), function(item) {
        rec <- item$record
        d <- item$coordinates
        d$beta <- NA_real_
        if (item$methylation$matched) {
          ix <- match(item$methylation$points$index, d$index)
          d$beta[ix] <- item$methylation$points$beta
        }
        data.frame(sample = rec$sample_id, group = rec$group,
                   gene = rec$symbol, assembly = rec$assembly,
                   chromosome = rec$chromosome, d, stringsAsFactors = FALSE)
      })
      utils::write.csv(do.call(rbind, frames), file, row.names = FALSE, na = "")
    }
  )
  output$delta_csv <- downloadHandler(
    filename = function() {
      t <- cancer_result()
      req(!is.null(t))
      sprintf("TCGA_%s_%s_delta.csv", t$cohort, t$record$symbol)
    },
    content = function(file) utils::write.csv(cancer_delta()$data, file, row.names = FALSE)
  )
  output$template <- downloadHandler(
    filename = function() "methylation_template.csv",
    content = function(file) file.copy("data/methylation_template.csv", file, overwrite = TRUE)
  )
}

shinyApp(ui = ui, server = server)
