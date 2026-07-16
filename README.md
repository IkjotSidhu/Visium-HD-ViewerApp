# Visium HD Viewer

An R Shiny app for interactively exploring **10x Genomics Visium HD** spatial transcriptomics data stored in Seurat objects — no coding required.

Load a Seurat `.RDS` file and explore it through spatial plots, UMAP embeddings, gene expression comparisons, cell-type composition, and metadata summaries. Every plot is publication-ready and exports to PDF or PNG.

<!-- Add a screenshot here once you have one:
![Visium HD Viewer](docs/screenshots/overview.png)
-->

---

## Features

| Tab | What it does |
|---|---|
| **Spatial** | Plot clusters or gene expression over tissue images. Choose sample/image, point size, alpha, and colour scale. |
| **UMAP / Reduction** | `DimPlot` and `FeaturePlot` on any dimensionality reduction, with optional cluster labels. |
| **Feature Expression** | Violin, box, and dot plots for one or many genes, grouped by any metadata column — with optional statistics. |
| **Composition** | Stacked/grouped bar charts showing cell-type or cluster makeup per sample. |
| **Metadata** | Bar charts, histograms, density plots, and a searchable data table for any metadata column. |

**Across the whole app:**

- **Correct cluster ordering** — clusters sort numerically (0, 1, 2 … 10, 11) rather than alphabetically (0, 1, 10, 11, 2 …).
- **25 colour palettes** — a built-in 60-colour palette, 14 GraphPad Prism palettes via [ggprism](https://csdaw.github.io/ggprism/), and 5 ColorBrewer sets. Palettes extend automatically if you have more clusters than colours.
- **Per-tab themes** — Classic, Prism, Minimal, or Black & White, chosen independently for each plot type.
- **Prism-style statistics** — Wilcoxon or t-test with significance brackets on violin and box plots.
- **Biologist-friendly labels** — the app reports *bins* (Visium HD), *spots* (standard Visium), or *cells* (scRNA-seq) as appropriate, and says *genes* rather than *features*. Bin size (8 µm / 16 µm) is auto-detected.
- **Export anything** — every plot saves as PDF (vector) or PNG at a width and height you specify.

---

## Requirements

- **R ≥ 4.2** (developed on R 4.5.2)
- Enough RAM to hold your Seurat object. Visium HD objects are large — a 16 µm binned object may need ~8–16 GB, and an 8 µm object can exceed 24 GB.

---

## Installation

```bash
git clone https://github.com/IkjotSidhu/Visium-HD-ViewerApp.git
cd Visium-HD-ViewerApp
Rscript install_packages.R
```

The installer pulls everything from CRAN and skips anything already present.

<details>
<summary><b>Seurat installation trouble?</b></summary>

Seurat depends on system libraries that may need installing first.

**macOS** (with [Homebrew](https://brew.sh)):
```bash
brew install hdf5 gdal geos proj
```

**Ubuntu / Debian:**
```bash
sudo apt-get install libhdf5-dev libgdal-dev libgeos-dev libproj-dev
```

Then re-run `Rscript install_packages.R`.
</details>

---

## Usage

Launch from a terminal:

```bash
Rscript launch_app.R
```

Or from an R console / RStudio:

```r
shiny::runApp("app.R", launch.browser = TRUE)
```

Your browser will open the app. Then:

1. Click **Browse for RDS File** to pick your Seurat object with a native file dialog — or paste a full path into the text box.
2. Click **Load Object**. Large objects take a minute or two.
3. Explore the tabs. Adjust the palette and font size in the sidebar; each tab has its own theme selector.
4. Set a width/height and click **Save** to export any plot as PDF or PNG.

> **Note:** the app reads your file directly from disk — nothing is uploaded or copied anywhere. It runs entirely on your own machine.

---

## Supported data

Built for **Visium HD** Seurat objects, but works with anything Seurat-shaped:

| Data type | Spatial tab | UMAP tab | Expression / Composition / Metadata |
|---|:---:|:---:|:---:|
| Visium HD (8 µm / 16 µm bins) | ✅ | ✅ | ✅ |
| Standard Visium (spots) | ✅ | ✅ | ✅ |
| scRNA-seq (no images) | — | ✅ | ✅ |

If your object has multiple assays (e.g. `Spatial.008um` and `sketch`), switch between them with the **Active Assay** selector in the sidebar.

---

## Statistics

Optional significance testing on **violin** and **box** plots. Enable **Add statistics** in the Feature Expression tab.

- **Test** — Wilcoxon rank-sum (non-parametric, the default and the right choice for expression data) or t-test (parametric).
- **Comparisons** — all pairwise, or every group against one reference group.
- **Display** — significance stars (`*`, `**`, `***`), adjusted p-value, or raw p-value.
- **Multiple-testing correction** — Benjamini-Hochberg (FDR; the genomics standard and the default), Bonferroni (conservative, family-wise error rate), or none.

Brackets are drawn Prism-style and stack automatically. In multi-gene views each facet gets its own independent comparisons.

---

## Troubleshooting

**"vector memory limit reached" / R crashes on load**
Your object is bigger than available RAM. Try the 16 µm binned object instead of 8 µm, or run on a machine with more memory. On macOS you can raise R's limit by adding `R_MAX_VSIZE=32Gb` to `~/.Renviron`.

**Spatial tab says "No spatial images in this object"**
The object has no image data — expected for scRNA-seq. Use the UMAP tab instead.

**"Statistics error" notification**
Usually too few cells in one group for the chosen test. The plot still renders without brackets.

**Plot is slow with many bins**
Expected with hundreds of thousands of bins. Reduce point size, or switch to the `sketch` assay if your object has one.

---

## Built with

[Shiny](https://shiny.posit.co/) ·
[Seurat](https://satijalab.org/seurat/) ·
[ggplot2](https://ggplot2.tidyverse.org/) ·
[ggprism](https://csdaw.github.io/ggprism/) ·
[rstatix](https://rpkgs.datanovia.com/rstatix/) ·
[bslib](https://rstudio.github.io/bslib/)

## License

Released under the [MIT License](LICENSE) — free to use, modify, and distribute.
