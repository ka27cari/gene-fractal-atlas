# Technical details

# Gene Fractal Atlas

R Shiny application for DNA sequence Chaos Game Representations (CGR) and **measured** tumor/normal DNA methylation spatial patterns. Each selected cancer type lists its own public studies; no LUAD-only workflow.

## Run

```r
install.packages(c("shiny", "curl", "jsonlite", "htmltools"))
shiny::runApp(".")
```

Run commands from the extracted `gene-fractal-atlas` directory. Live modes require an internet connection. No patient measurements are bundled or synthesized as real observations. Real, queried gene-window CpG/sample measurements are kept as compressed `.rds` files in the **user-local permanent cache**, default `~/GeneFractalAtlas/data/atlas-v1/` (Windows: under your home folder). The release ZIP includes **no patient data**. Once a gene-window/cohort has successfully loaded, later loads reuse that saved copy even without an internet connection. Data does not expire after 14 days. Choose **Saved first** for normal use, **Offline only** to skip missing cohorts without waiting for network timeouts, or **Refresh selected datasets** to replace the copies. The UI reports the storage path and disk use. To put the cache on another drive, set the environment variable `GFA_DATA_DIR` before starting R (the app appends `atlas-v1`). Existing `cache/atlas-v1/` RDS files from older releases are migrated on use when present; copy that directory into the updated project first if needed.

## Analysis workspaces

**Gene comparison:** reference DNA sequences for selected orthologs. Uses triangular or square CGR and the optional experimental ATG/C/mC map. The separate beta cutoff in this workspace is needed **only** to label a C as methylated when the user supplies beta measurements or explicitly selects a synthetic demo. `ASXL1` is available among gene candidates.

**Cancer / normal:** select one of 18 TCGA cancer codes, human gene, compatible public GEO datasets, and a reference window. Every successfully loaded cohort contributes an adjacent **Tumor / Normal** pair, arranged in a four-column, row-major grid. Mean measured beta values are drawn in red; the reference DNA A/T/C/G path is rendered in visible green shades in **every panel**. No high-methylation beta cutoff is applied to mean-beta panels.

- **Shape:** 3-vertex triangle (AT/C/G) or square (A/C/G/T). A/T share the top vertex in the default triangle. Both use the same GRCh37 reference for comparable spatial coordinates.
- **Select by each panel's mean β:** the two-ended slider defaults to **0–1** (all eligible values). A Tumor panel includes a CpG only if *its tumor-group mean β* lies inside the interval, and a Normal panel independently checks *its normal-group mean β*. A probe may appear in just one panel. A Δβ or cross-cohort agreement calculation requires measured means in **both** groups but accepts a probe when **either** group is in the chosen interval; the Δβ cutoff is then applied. Pooled mean panels filter the pooled mean after pooling the source cohorts. Pooled Δβ uses paired probe contributions and the pooled group mean to determine eligibility. Spatial click selection and exports follow these exact definitions. Green reference DNA remains unfiltered.
- **Δβ:** optionally draw a map for every loaded cohort below its tumor/normal facets. Red = tumor-higher; blue = tumor-lower. The adjustable per-probe `|mean tumor beta - mean normal beta|` cutoff defaults to 0.15 and is applied before binning. It is **not** a statistical significance test.
- **Cross-cohort agreement:** available with at least two independently sourced, **successfully loaded** cohorts from the same cancer. Only identical measured CpG probe IDs satisfying per-cohort sample and Δβ criteria are considered. They are then aggregated into fixed CGR bins. Red = unanimous positive tumor-normal difference; blue = unanimous negative difference; amber = mixed directions. A ring of **cohort-colored ticks** identifies participants; dim ticks disagree with the bin's majority. This map is descriptive, not evidence of biological replication by itself. GSE32861 and GSE32866 share a parent study and are not counted as independent from each other.
- **Click a colored bin:** display same-region Tumor / Normal patterns from all loaded cohorts, exact CpG IDs, probe-level summary and sample measurements, with CSV download. CGR-near points are not necessarily genomically adjacent.
- **Combine cohorts:** separate; equal-cohort means; specimen-count-weighted means; or only probes measured in both groups of every selected cohort. Original cohort facets remain visible. This is descriptive pooling, not a batch-corrected meta-analysis. Mismatched cancer types are rejected.
- **Advanced:** spatial-bin resolution, minimum per-probe specimen count, minimum CpGs per cell, validated user-provided cohort beta CSV.

## Cancer-specific public GEO catalog

The selector also offers public GEO series. An entry in the catalog means **a GEO study and its tumor/normal design were located**, **not that the external series matrix has been downloaded and parsed successfully here**. The adapter tries its public `.txt.gz` GEO Series Matrix as a stream, retaining only mapped `cg` probe IDs and samples whose tumor/normal labels can be established. Any missing processed matrix, ambiguous labels, incompatible probe IDs or non-beta assay values cause an explicit error. No fabricated or proxy values are substituted. Large 450K matrices may take substantial time or require the offline CSV import route.

| Cancer type | Additional GEO study IDs offered (beyond TCGA) |
|:--|:--|
| LUAD | GSE32861, GSE62948, GSE63384, GSE32866 |
| LUSC | GSE121849; LUSC-labeled subset of GSE56044 |
| BRCA | GSE69914; breast-only subset of GSE53051 |
| COAD | GSE131013, GSE48684; colon-only subset of GSE53051 |
| KIRC | GSE61441; kidney-only subset of GSE52955 |
| LIHC | GSE54503 |
| HNSC | GSE87053 (oral squamous cancer) |
| PRAD | GSE76938; prostate-only subset of GSE52955 |
| THCA | thyroid PTC-only subset of GSE53051 |
| BLCA | bladder-only subset of GSE52955 |
| CHOL | GSE49656 |
| PAAD | GSE181740; pancreas-only subset of GSE53051 |
| CESC | GSE211668 |
| STAD | GSE164988 |
| ESCA | GSE121930 (esophageal squamous cancer) |
| UCEC, KIRP, READ | TCGA; no additional tumor/normal GEO adapter was responsibly established |

**Do not treat a selected but failed dataset as loaded.** The app's **Data coverage** and **Data sources and errors** sections report the actual specimen and probe counts and provider failures. Detailed public-study links and documented GEO platform labels are kept in [`research/geo_cohort_catalog.csv`](../research/geo_cohort_catalog.csv). Some series (e.g. GSE53051, GSE52955) contain mixed cancers; these are explicitly restricted to matching disease samples using confirmed GEO sample naming rules.

Adjacent non-tumor tissue and healthy donor controls are scientifically different groups. This release uses a broad `Normal` category for descriptive displays, so downstream analysis must check and, if necessary, stratify control source, patient pairing, age, tissue origin and processing batch before drawing biological conclusions.

## Optional user cohort CSV

Download the template inside **Cancer / normal**. Required fields:

```csv
cancer,cohort,sample,group,probe,beta
```

Cancer must exactly match the chosen code; `group` is `Tumor` or `Normal`; beta in [0,1]; CpGs use `cg` IDs. When not already in the verified reference mapping, use additional `assembly,chrom,position` fields (GRCh37, forward genomic CpG cytosine position). Conflicting reference coordinates, invalid beta values, duplicate observations and tumor/normal cancer-code mismatches are rejected. Imported studies are labeled as user-supplied data, not as independently verified GEO records.

## Scientific scope and limitations

- A 450K/27K methylation array provides measured beta values at assayed CpGs, **not** patient DNA sequences. The green dots are from the same **reference DNA** for all tumor/normal panels. Tumor-specific sequence variation cannot be inferred from them.
- CGR patterns largely reflect k-mer context and reference composition. A visually unusual motif or agreement bin does not establish a clinically significant marker without platform-coverage, batch, tissue-composition and probe-density controls and independent replication.
- Probe IDs must be identical across studies to enter the agreement computation, but separate arrays can differ in technical coverage. Cross-platform missingness is never filled by artificial values.
- Pairing of patient samples is not assumed in mean-beta and Δβ calculations. This release performs no formal hypothesis tests and no WGS/WGBS import.
- A substantial fraction of publicly accessible GEO studies distribute raw intensities or supplementary tables instead of the Series Matrix format expected by the online adapter. Those records may still require offline beta-matrix processing and CSV import. The public study list should **not** be mistaken for a list of live-verified loaders.

## Check locally with R

```r
source("tests/startup_glow_test.R")
source("tests/atlas_cohort_test.R")
source("tests/beta_range_cache_test.R")
source("tests/cancer_merge_test.R")
source("tests/multicancer_concordance_test.R")
shiny::runApp(".")
```

GitHub Actions uses R 4.4 for source parsing and the regression test suite. The release ZIP was built in an environment without R; **actual R execution, external provider downloads, and successful online loading of all catalog entries are not verified**. Source parsing and offline fixture tests are intended to run in CI or on the user's R installation.

## Public sequence research (not yet an app feature)

The `research/` folder contains audits of publicly available reference DNA, open tumor/normal WGS and WGS-derived mutation tables. Public variant calls are distinct from fully downloadable raw WGS. No controlled-access clinical data should be committed to this repository.

## Primary data services

- [UCSC Xena TCGA 450K](https://xenabrowser.net/datapages/)
- [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)
- [UCSC Genome API](https://api.genome.ucsc.edu/), [Ensembl REST](https://rest.ensembl.org/), [MyGene.info](https://mygene.info/)

### CGR corner labels and color contrast (v0.13.1)

- The standard **triangle** uses three vertices: `A/T` together at the top, `C` at bottom-left and `G` at bottom-right. It is *not* a four-distinct-base triangle.
- The **square** uses `C` upper-left, `G` upper-right, `A` lower-left, `T` lower-right.
- The experimental **methyl-C triangle** uses `A/T/G` together at the top, `C` bottom-left and `mC` bottom-right.
- Saturated red/blue cohort spatial markers are always visibly opaque on the black background. A faint surrounding halo helps distinguish them from green reference sequence points. Mean beta and delta-beta values remain available in SVG tooltips; marker opacity is deliberately **not** a reliable linear quantitative color scale.
