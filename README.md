# Gene Fractal Atlas

**Gene Fractal Atlas** is an R Shiny app for exploring DNA sequence patterns with the **Chaos Game Representation (CGR)** and visualizing measured DNA methylation in cancer and non-tumor tissues.

![preview]("docs/preview.png")

The aim is to compare **spatial patterns**, not only individual CpG sites. Cohorts can be viewed separately or summarized together, with reference sequence points and tumor–normal methylation signals shown in the same coordinate system.

## What you can explore

- **DNA geometry:** triangular or square CGR of reference A/T/C/G sequences, with labeled vertices.
- **Cancer cohorts:** tumor and normal panels for each successfully loaded dataset, organized by cancer type.
- **Methylation:** mean β, tumor-minus-normal **Δβ**, and an adjustable **per-panel mean β range**.
- **Reproducibility:** cross-cohort agreement for shared measured CpG probes and selectable cohort pooling methods.
- **Pattern inspection:** click a colored spatial bin to inspect contributing CpGs and cohort measurements.
- **Local storage:** save queried gene-window measurements in a persistent local cache for faster reuse and offline exploration.

## Get started

You need R and an internet connection for the first live data download. In R, from the repository directory:

```r
install.packages(c("shiny", "curl", "jsonlite", "htmltools"))
shiny::runApp(".")
```

Choose **Gene comparison** for reference sequence exploration or **Cancer / normal** for methylation cohort comparisons. Successfully fetched measurements are saved by default to `~/GeneFractalAtlas/data/atlas-v1/` (outside this repository). Use **Offline only** for data already saved locally.

## Data and interpretation

The current cancer analysis uses **public TCGA/GEO array methylation measurements**, where available. The selected cancer type determines which cohorts are offered. **A listed GEO study is not a guarantee that its β matrix can be downloaded and parsed automatically**; the app reports loading failures and also supports an import template.

The green sequence points are from **reference DNA**. The array data do **not** provide patient-specific tumor sequences, and the red methylation marks are not sequence mutations. CGR patterns are exploratory and do not establish statistical or clinical significance on their own. WGS/WGBS integration is *not* implemented in this release.

See [technical details, data coverage, and limitations](docs/TECHNICAL.md) for cohort IDs, input schema, pooling definitions, and scientific cautions. The public sequence/WGS resource investigation is in [`research/`](research/).

## Files

| Location | Contents |
|---|---|
| `app.R` | Shiny application |
| `R/` | Sequence, CGR, methylation, cohort, and SVG calculations |
| `www/` | Interface styling |
| `data/` | Example CSV templates (no patient records) |
| `tests/` | Source and regression tests |
| `research/` | Public resource inventory and source notes |
| `docs/preview.png` | README preview image, replaceable with a current screenshot |

## Verification

This GitHub-ready package contains **analysis code updated in v0.14.2** to apply β ranges per tissue panel and restore β-dependent red intensity. The GitHub preview remains a historical screenshot. Source and archive integrity checks were performed, but a complete **R/Shiny runtime and live data retrieval have not been verified** here. A [GitHub Actions workflow](.github/workflows/r-smoke.yml) is included for automated R tests.

Data providers: [UCSC Xena](https://xenabrowser.net/datapages/) · [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) · [Ensembl REST](https://rest.ensembl.org/) · [UCSC Genome API](https://api.genome.ucsc.edu/)
