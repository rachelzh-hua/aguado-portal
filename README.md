# Aguado Lab Aortic Valve Portal

Interactive single-cell + spatial transcriptomics portal for **GEO GSE306288**
("Spatial transcriptomic profiling of the human aortic valve reveals cellular sex
differences near sites of calcification", Baddour et al., Aguado Lab, UC San Diego).

The link opens on **Live** — a cinematic walkthrough of the analysis with a
streaming console (see below). The **Portal** segment is the full interactive
portal described here, one click away and unchanged.

Built in the **Explorer** chrome — an indigo gradient header with
`Live / Portal / Insights / Report` mode segments and a live dataset-scope pill; a
**Workspace** sample rail on the left; a coordinated **multi-pane Vitessce grid**
in the center (Spatial + UMAP on top; Cell Sets + Gene List + Heatmap below); and
a **Copilot** panel on the right (dataset context, suggested actions, and a
plain-language command box that routes gene names / keywords to the analysis
panels). Gene search, lasso A/B compare, per-sample inspector, integrated Atlas,
Pathways and **Spatial** (squidpy statistics) remain as Copilot tabs, and a
**▦ Grid** stage plots every Visium section side by side. Fully static — the
Copilot routes commands locally (no LLM backend).

## Contents
```
app.html                 self-contained portal (Plotly + marked from CDN)
vitessce.html            coordinated multi-view explorer (Vitessce, loaded from esm.sh CDN)
data/
  vitessce/              per-dataset AnnData-Zarr stores + Vitessce configs + index.json
  manifest.json          {samples:[...], meta:{...}, pathways:[{key,name,category,genes}], signatures:[...]}
  precomputed.json       per-sample counts, marker genes, signature + pathway stats
  <sid>.json             {modality, spatial:[[x,y]], umap:[[x,y]], obs:{cell_type, *_score, <pathway keys>}}
  <sid>.expr.json/.bin   quantized (uint8, gene-major) log-norm expression, 700 genes
  atlas.json             integrated scRNA UMAP + obs (subtype, lineage, sample, sex, disease)
  atlas.expr.json/.bin   compact key-gene expression for the atlas
  atlas_meta.json        subtype composition, markers, ligand-receptor axes, Moran's I
  atlas_pathways.json    subtype×pathway activity matrix (z) + sex/disease Δ-activity contrasts
  sq_stats.json          squidpy spatial statistics per Visium section (Spatial tab)
  graph/<sid>.json       6-NN spatial-graph edges in portal spot indexing (Grid edges)
  tissue/<sid>.jpg       cropped, display-enhanced H&E background + index.json geometry
```

## Data
17 samples: **8 Visium spatial** (tissue map) + **9 scRNA-seq** (UMAP). All human
aortic valve leaflet. Sex/disease per sample come from the GEO series matrix.
Patients **1** and **23** are profiled by both modalities (cross-modality anchor).

- Colors: **Cell type**, a signature (**Osteogenic/Calcification**, **VIC
  Activation**, **Inflammatory**, **Endothelial**, **Proliferation**), or one of
  **14 pathway-activity** scores (see below) — so you can see *where* a program is
  active on the tissue.
- Whatever the map is coloured by, the **Gene** tab shows that channel's summary
  stats, distribution and **per-cell-type violins** (genes, signatures, pathways,
  H&E stains, deconvolution fractions alike). The Composition and Top Markers
  panels do *not* follow the colour channel: both are properties of the fixed
  cell-type labels (counts, and per-label 1-vs-rest Wilcoxon markers), so they
  only change with the sample or the A/B selection.
- Cell types are computational (Leiden + relative-enrichment marker scoring). On
  Visium (spot-level, ECM-rich), labels are the *dominant program* per spot.
  Types: VIC, activated VIC (myofibroblast), endothelial, macrophage, T/NK,
  B/plasma, mast. (The ACTA2/TAGLN program in diseased valve is activated VICs,
  not vascular smooth muscle; FABP4/TREM2 cells are foam macrophages.)
- **Paper headline genes are guaranteed searchable in every sample** (bypass the
  700-gene HVG cap): `COMP`, `FOS`/`FOSB`/`JUN`/`JUNB`/`EGR1` (AP-1),
  `CD44`, `SPP1`, plus foam-mac (`TREM2`,`GPNMB`,`APOE`) and CAVD/osteogenic genes.
  Where a gene is undetected in a sample it is carried as a flat-zero column
  ("not detected") so cross-sample search stays consistent.

## Provenance & caveats — read before interpreting

This portal is an **independent Scanpy re-analysis, not a reproduction** of the
paper. Key differences from Baddour et al.:

| | This portal | Paper |
|---|---|---|
| Toolkit | Scanpy, log1p | Seurat v4.4.0 + SCTransform v2 |
| Scope | per-sample, **no integration** | 12-patient integration (incl. Xu et al. 2020 data **not in this GEO deposit**) |
| Clustering | Leiden 0.6–1.3 | resolution 0.01 / 0.1 |
| Cell types | 7 marker-scored (per-sample) | 4 (VIC, macrophage, VEC, T cell) + VIC1–4 / Mac1–3 subclusters |
| Integration / subtypes / CellChat / Moran's I | **approximated in the Atlas tab** (heuristics) | full methods with statistical tests |

### Integrated Atlas (the "Atlas" tab)

`pipeline/integrate.py` builds a cross-sample integrated view the per-sample
portal can't — the closest we get to the paper's structural analyses:

- **Harmony** integration of all 9 scRNA samples (batch = sample), re-clustered
  (Leiden) and re-embedded (UMAP) → `atlas.json` (47,996 cells).
- **VIC / macrophage subtypes** by hierarchical marker scoring: VIC quiescent /
  myofibroblast / osteogenic / inflammatory, and macrophage resident / inflammatory
  / foam.
- **Subtype composition** by sex and disease.
- **Ligand→receptor axes** — a CellChat-*style* **co-expression heuristic** (mean
  ligand in the top sender subtype × mean receptor in the top receiver), not a
  permutation-tested CellChat run.
- **Moran's I** spatial autocorrelation (squidpy) averaged across the 8 Visium
  samples → top spatially-variable genes.

These are **independent heuristics**, useful for orientation, not restatements of
the paper's integrated, statistically-tested results.

### Pathway analysis (the "Pathways" tab)

Both pipelines score **14 curated CAVD-relevant pathway gene sets** per cell/spot
with scanpy `score_genes` (a control-gene-corrected mean, i.e. AddModuleScore-style
activity — **not** permutation-tested GSEA/ORA). Panels: TGF-β/VIC activation,
EMT/fibrosis, BMP/osteogenic, WNT, Notch, NF-κB/TNF, interferon, complement,
MHC-II antigen presentation, ECM/MMP remodeling, angiogenesis, hypoxia, lipid/foam-
cell, senescence/SASP. Members are canonical Hallmark/KEGG/Reactome genes pruned
toward genes detectable in valve tissue (`pipeline/process.py::PATHWAYS`).

- **Per-sample:** each pathway is a **Color dropdown** option — color the tissue /
  UMAP by pathway activity, and it feeds the Inspector and Compare A/B panels.
- **Atlas (`integrate.py`):** a **subtype × pathway** activity heatmap (z-scored
  across atlas cells) plus **male-vs-female** and **diseased-vs-healthy** Δ-activity
  contrasts → `atlas_pathways.json`, rendered in the **Pathways** tab.

Δ-activity contrasts are **descriptive** — no statistical test — and the disease
contrast is **confounded** (all 3 healthy samples are male; see below).

**Cohort limitation:** this deposit has only **3 healthy samples, all male** — there
is **no healthy female**. Sex×disease comparisons are confounded; treat any
cross-group difference as descriptive, not causal.

The portal is best used for **spatial exploration and gene lookup** (where does
`COMP` sit relative to calcification? where are macrophages?), not as a restatement
of the paper's integrated, statistically-tested findings.

## Live mode (the **● Live** segment) — cinematic preview + console

The default landing view once the abstract modal is dismissed. It is an overlay across
the workspace rows, so the Portal underneath is never rebuilt or disturbed by the switch;
`setMode('portal')` just hides it.

**The preview.** Every cell/spot in the payload is one particle on a 2-D canvas. A *layout*
assigns each particle a target position, colour and alpha; the render loop eases the current
state toward the target (`k = .075`) with a slow sinusoidal drift. Switching layouts therefore
*moves the same cells* rather than redrawing a new chart — the atlas UMAP folds into the real
tissue coordinates of a leaflet, and back. Particles with no role in a layout park on a faded
ring at α .045 instead of vanishing.

Payload: an even subsample of `atlas.json` (9,000 of 47,996 cells, so every donor is
represented) plus the full spot sets of `P1A` (male, diseased) and `P6A` (female, diseased).
Layouts: atlas UMAP by lineage / subtype / sex / disease, one Visium section by cell type,
one Visium section by any score or gene, and a **paired** layout that puts both leaflets on
one canvas under a *shared* colour scale — the sex contrast, made honest by construction.

**The console.** A terminal strip streams the pipeline as it runs (`▸` step, `✓` result,
`◆` finding, `dim` chatter). Every number in it is computed at runtime from the same
`data/*.json` the Portal renders — spot and cell counts from the manifest, subtype counts
from `atlas_meta.json`, sex means from `sex_stats.json`, Moran's I and the ligand–receptor
axis from `atlas_meta.json` — so the console and the maps can never disagree. Once the story
ends the console keeps narrating: `switchSample`, `selectGene`, `switchTab`, `onColorChange`,
`clearSelection` and `clearGene` are wrapped so any action in the Portal prints a line.

**The story** (8 scenes, scrub / pause / restart): ingest → annotate 7 lineages → the VIC
disease axis (quiescent 18,743 → activated/myofibroblast 16,347; osteogenic VICs are rare in
the scRNA atlas and read out spatially via COMP) → into the P1A leaflet at real (x,y) → paint
the osteogenic signature → paint `COMP` (highest Moran's I, 0.324) → male vs female side by
side → molecule/programme/niche/axis. (Macrophages resolve into resident / inflammatory /
foam-lipid (FABP5/MSR1, 683 cells) after the 2026-07-31 annotation fix.)

The sex scene reports what the data actually shows, which is *not* the intuitive story:
across the 8 Visium sections the osteogenic signature is near parity (male 0.51 vs female
0.49) while **VIC activation runs higher in the female leaflets** (0.40 vs 0.73, Δ −0.33).
These are section-level means over n=4 per sex — a descriptive contrast, not a powered test,
and the Inspector says so on screen.

Autostart is gated on the abstract modal being closed, so the story never plays out behind
it. Deferred console lines carry an epoch token and drop if the viewer has already moved on.

## Multi-section grid (the **▦ Grid** button)

The portal's `sq.pl.spatial_scatter`: any number of Visium sections × one or two colour
channels, in one grid. It carries the arguments that plot has, as live controls:

- **Sections**: chips select which libraries are plotted; the grid re-lays itself out.
- **Channel 1 / Channel 2**: any two of cell type, signature, pathway, deconvolution
  fraction, H&E stain, or the active gene. Each channel gets **one colour scale across
  every section**, so a colour means the same value in all panels.
- **⇄ Sections first / Channels first**: `library_first`, i.e. whether a row is a channel
  (sections across) or a section (channels across).
- **⬡ Edges**: draws the 6-NN spatial graph (`connectivity_key`) beneath the spots. This
  is the same graph the neighborhood enrichment, centrality and Moran's I are computed on,
  so the analytics and the picture are the same object.
- **⛓ Link crop**: zoom or drag any panel and the identical `crop_coord` is applied to
  every panel. Squidpy's semantics exactly, including the consequence that one section's
  crop coordinates can fall outside a smaller section. **◫ Crop to A** takes the crop from
  lasso bucket A's bounding box; **⤢ Reset** clears it.
- **Size**: global spot size, plus per-panel `−` / `+` (squidpy's per-library `size` list).
- **◍ H&E** background, **▬ µm** 500 µm scale bar, and **⤓** per-panel PNG export at 3×
  for figures.

The scale bar and every µm figure in the portal are calibrated from the Visium array
itself (two spots two `array_col` apart are 100 µm centre-to-centre), **not** from
`spot_diameter_fullres` — Space Ranger reports that as the 65 µm render diameter, not the
55 µm capture spot, which would shrink every distance by 18%.

## Spatial statistics (the **Spatial** tab)

`pipeline/squidpy_stats.py` precomputes the rest of the squidpy toolbox per Visium section
into `data/sq_stats.json`; the tab plots it for the section you are on:

- **Co-occurrence** (`sq.gr.co_occurrence`): p(type | a spot of the chosen type within r) ÷
  p(type), as a function of radius in µm. Answers *at what distance* a lineage is enriched
  around another, which the pairwise neighborhood z-score cannot.
- **Ripley's L** (`sq.gr.ripley`, 100 permutations): clustered vs indistinguishable-from-
  random, with the permuted 95% envelope drawn. Visium labels sit on a fixed lattice, so
  this is clustering of the *labels*, not of physical cells.
- **Graph hubs** (`sq.gr.centrality_scores`): degree / closeness / clustering per lineage
  over the spatial graph.
- **Adjacency** (`sq.gr.interaction_matrix`): the observed fraction of graph edges between
  each pair of lineages — the untested counterpart of the ⬡ Neighbors permutation z-score.
- **SVG: sepal vs Moran** (`sq.gr.sepal` on the hex-grid graph vs `sq.gr.spatial_autocorr`
  Moran's I with 100 permutations, over the 700 exported genes). The two rank genes
  differently on purpose: Moran rewards smooth gradients, sepal rewards sparse patchy
  structure. Click any point to colour the map by that gene.
- **H&E ↔ expression** (`sq.im.calculate_image_features`, summary + Haralick texture, plus
  skimage `rgb2hed` colour deconvolution): per-spot histology features correlated
  (Spearman) with the expression programmes at the same spots. The three stain densities
  (**hematoxylin** = nuclei, **eosin** = ECM/collagen, and hematoxylin SD = texture) are
  written back into `data/<sid>.json` and appear in **Color by → Histology**, so the
  histology itself can be put on the map next to the transcriptomic programmes.

Single-cell-type sections (P20V, P23V are all-VIC) have no cell-type graph structure, so
co-occurrence / centrality / adjacency are skipped there and the panel says so.

## Vitessce coordinated multi-view (the **Vitessce** tab)

A [Vitessce](http://vitessce.io) linked-brushing explorer — where **spatial ↔ UMAP
↔ heatmap ↔ cell-sets ↔ gene-list** are all coordinated (select a cell set or gene
and every view updates together) — is embedded **directly in `app.html`** as the
default **Explorer** center view. It renders as a **multi-pane grid** by default
(Visium: Spatial + Scatterplot(UMAP) on top; Cell Sets + Gene List + Heatmap on the
bottom row — the Cell Sets / Gene List panes are injected client-side in
`vtApplyLayout` by cloning the heatmap's coordination scopes, so no config
regeneration is needed). The `UMAP` / `Heatmap` toggles in the Copilot panel drop
back to a single-hero pane. The
engine (React + Vitessce, from the esm.sh CDN) is **lazy-loaded on first open** so it
adds nothing to initial page load; a dataset dropdown in the side panel switches
between the integrated atlas and any of the 17 samples. The CommonJS interop shims
the Vitessce bundle needs are set only at that first open, so they can't interfere
with Plotly's / marked's UMD detection.

The standalone `vitessce.html` (same viewer, full-window) is retained as an optional
fallback but is no longer linked from the portal.

- **Data** (`pipeline/vitessce_export.py`): reconstructs a compact AnnData for each
  dataset directly from the already-exported `data/*.json` + `*.expr.bin` (no scanpy
  re-run, so the Vitessce view is byte-consistent with the main portal), restricted
  to the curated gene set, and writes an **AnnData-Zarr (v2)** store +
  `<sid>.config.json` under `data/vitessce/`.
- **Engine (self-hosted, no CDN)**: React + Vitessce are **vendored onto this origin**
  under `data/vitessce/vendor/` (importmap points there), so the viewer needs **no
  third-party CDN at runtime** — it works on locked-down institutional networks that
  block esm.sh. The bundle's blosc zarr codec (`dist/blosc-*.mjs`, wasm inlined) is
  vendored too, so Zarr chunks decode offline; heavy unused view chunks (GeoTIFF
  image codecs, HiGlass, Neuroglancer, 3D-text, SpatialWrapper) are left as
  never-reached esm.sh fallbacks to keep the vendored set to ~9.5 MB. Re-vendor with
  `python3 pipeline/vendor_esm.py portal/data/vitessce/vendor`. The page still
  resolves each config's relative store URLs against the config file's own location,
  so it stays fully static and deploy-portable.
- **Static-host note (important):** Zarr v2 metadata are dotfiles (`.zarray`,
  `.zattrs`, `.zgroup`, `.zmetadata`) and **Vercel 404s dotfiles**. The export
  therefore also writes a non-dot copy of each (`.zarray`→`zarray` …) and
  `vercel.json` rewrites `…/.zarray` → `…/zarray`. `python3 -m http.server` serves
  the dotfiles directly, so local dev needs no rewrite. (Verified end-to-end against
  a Vercel-mimicking server that 404s dotfiles.)

## Run locally
```bash
cd portal && python3 -m http.server 8899   # then open http://localhost:8899/app.html
```
Fully static — no backend, no API keys, no analytics. Deploy the `portal/` folder
to any static host (e.g. Vercel) as-is; `app.html` is the entry point, and the
`vitessce.html` explorer + esm.sh CDN both work from the same static deploy.

## Regenerate data
```bash
python3 ../pipeline/process.py           # per-sample data/ from ../raw/ (scanpy)
python3 ../pipeline/integrate.py         # integrated atlas.* (Harmony + squidpy)
python3 ../pipeline/squidpy_stats.py     # data/sq_stats.json + data/graph/* + HE_* obs (~13 min)
python3 ../pipeline/vitessce_export.py   # data/vitessce/* (AnnData-Zarr + Vitessce configs)
```
`squidpy_stats.py` needs squidpy ≥1.6 with the image extras (scikit-image, tifffile, dask);
`/home/ubuntu/analysis-env/bin/python` has them. It reads `data/h5ad/<sid>.h5ad` and asserts
each h5ad's row order matches `data/<sid>.json` positionally before exporting graph edges,
since the edge indices address portal spots.
Bump `DATA_VERSION` in `app.html` after regenerating to bust browser cache
(`vitessce_export.py` reads the exported data/ files, so run it last).
