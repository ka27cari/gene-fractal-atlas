# v0.14.2

- Filter Tumor and Normal CpG means independently against the same β range instead of requiring Normal as a gate.
- Δβ and cross-cohort agreement include a probe when either group mean is in range, requiring measured values from both groups.
- Recompute pooled group means before range selection; only matched Tumor/Normal cohort contributions enter pooled Δβ.
- Restore noticeable red opacity gradients (0.52–0.98) while preserving green A/T/C/G reference dots and clickable spatial cells.
- Replace normal-baseline regression checks with group-level β filter/contrast checks.

## 0.14.1

- GitHub landing page with a screenshot from an earlier app build, explicitly labeled as such.
- Short README with extended technical details moved to `docs/TECHNICAL.md`.
- Git ignore rules expanded to reduce accidental upload of raw data and user-local caches.
- No changes to the v0.14.0 R analysis code or dataset loading behavior.

# Changelog

## v0.14.0 (2026-11-04)

- Store only successfully retrieved, real gene-window measurements and reference/probe maps as compressed RDS in the persistent user data directory (`~/GeneFractalAtlas/data/atlas-v1` by default), shared across app updates. Set `GFA_DATA_DIR` for another location. Existing v0.13.x `cache/atlas-v1` entries can be migrated on use. No 14-day cache expiry; Data access has saved-first, strict offline-only and explicit refresh modes. Cached data can be used offline once a requested gene-window and its cohorts have been downloaded.
- Show the actual local cache folder, compressed file count and disk use. Do **not** package proprietary, controlled or synthesized sample measurements in release ZIPs.
- Add inclusive per-cohort **normal-group mean β range** slider. Require at least the selected minimum sample count in both groups for displayed CpGs. Apply it to cohort Tumor/Normal facets, Δβ, pooled summaries, shared-cohort agreement, cell inspection and summary/Δβ exports. Leave green A/T/C/G reference points unchanged.
- Add offline regression tests for normal β selection, pooled agreement input filtering and persistent-cache roundtrip.

## v0.13.1 (2026-11-04)

- Boosted cohort-level tumor/normal, difference and agreement CGR markers: bright saturated red/blue cores, translucent halos and larger minimum radius. No methylation values, bins or cohort merges changed; tooltips still show the exact values.
- Added outside-corner labels to all SVG plots. Square: A/C/G/T; standard triangle: shared A/T apex with C and G bottom vertices; methyl-C triangle: shared A/T/G apex with C and mC bottom vertices.
- Kept green A/T/C/G reference traces and spatial-bin click targets; text labels do not intercept pointer events.
- Added offline regression assertions for label positions, clear markers and click overlays. Actual Shiny execution remains unverified without an R runtime.


## v0.13.0 (2026-11-04)

- Expanded the catalog from LUAD-only external GEO series to 25 publicly registered tumor/normal cancer-specific GEO entries across 15/18 TCGA cancer codes; TCGA remains a per-cancer selectable baseline. Online adapter availability is **unverified**, not labeled as successfully loaded without real beta data.
- Added a streaming gzip GEO Series Matrix reader for gene-window probes; preserved the small-memory fixture parser. Invalid assay scales, ambiguous tissue classes and incomplete public matrices fail explicitly.
- Added study-specific sample-label rules (GSE121849, GSE52955, GSE53051, GSE61441, GSE48684, GSE87053, GSE121930), including disease-specific filtering for mixed-cancer series.
- Removed the redundant high-methylation beta cutoff from Cancer / normal. Retained the meaningful `|Δβ|` filter and restored clear per-cohort tumor/normal/difference views.
- Strengthened always-visible green reference A/T/C/G CGR dots; added per-base legend. No patient sequence was inferred.
- Added same-probe, cross-cohort spatial-direction agreement map using red/blue/amber centers and cohort-colored radial ticks. Clicking a bin still opens measured CpG and specimen summaries.
- Maintained four-column facets, ASXL1, triangle/square options, and existing cohort merging and CSV exports.
- Added multi-cancer catalog, sample-label, streamed gzip, green-background and concordance offline regression fixtures to GitHub Actions.
- No local R runtime and no external GEO/Xena API access: static/package checks only; actual Shiny execution and every online data adapter are **not verified**.

## v0.12.1

- Restored click-to-inspect interactions for selected cohort-level spatial bins, measured CpGs, tumor/normal patterns and specimen-level CSV export.

## v0.12.0

- Added ASXL1 and measured Delta-beta maps with adjustable probe cutoff; catalogued public WGS research resources.

## v0.11.0

- Defined selected-cancer cohort layouts and separate/equal/weighted/shared-probe pooling. Added public-first WGS audit documentation.
