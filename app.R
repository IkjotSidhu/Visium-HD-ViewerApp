library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)
library(ggprism)
library(patchwork)
library(DT)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(dplyr)
library(tidyr)
library(rstatix)

# ============================================================
# PALETTE DEFINITIONS  (top-level so server can access them)
# ============================================================

# 60-color project palette from the main analysis workflow
PROJECT_CUSTOM <- c(
  "#4876FF","#CD853F","#8B4513","#6B8E23","#708090","#8B008B","#2F4F4F",
  "#8B0000","#483D8B","#556B2F","#8B4726","#4682B4","#6A5ACD","#A0522D",
  "#5F9EA0","#9370DB","#BC8F8F","#B8860B","#3CB371","#7B68EE","#CD5C5C",
  "#4169E1","#D2691E","#6495ED","#DC143C","#008B8B","#B22222","#228B22",
  "#DAA520","#808000","#DB7093","#48D1CC","#C71585","#191970","#800000",
  "#BA55D3","#9932CC","#FF8C00","#E9967A","#8FBC8F","#00CED1","#9400D3",
  "#696969","#1E90FF","#B0C4DE","#FF6347","#40E0D0","#EE82EE","#FF4500",
  "#DA70D6","#87CEEB","#6A5ACD","#00FA9A","#87CEFA","#FFA07A","#32CD32",
  "#00FF7F","#7FFF00","#ADFF2F","#CAF178"
)

.gp <- function(nm) tryCatch(ggprism_data$colour_palettes[[nm]], error = function(e) NULL)

PALETTES <- list(
  "Project Custom (60)"  = PROJECT_CUSTOM,
  "Prism Dark"           = .gp("prism_dark"),
  "Prism Light"          = .gp("prism_light"),
  "Prism Dark 2"         = .gp("prism_dark2"),
  "Prism Light 2"        = .gp("prism_light2"),
  "Candy Bright"         = .gp("candy_bright"),
  "Candy Soft"           = .gp("candy_soft"),
  "Pastels"              = .gp("pastels"),
  "Warm Pastels"         = .gp("warm_pastels"),
  "Ocean"                = .gp("ocean"),
  "Flames"               = .gp("flames"),
  "Floral"               = .gp("floral"),
  "Muted Rainbow"        = .gp("muted_rainbow"),
  "Colorblind Safe"      = .gp("colorblind_safe"),
  "Autumn Leaves"        = .gp("autumn_leaves"),
  "Stained Glass"        = .gp("stained_glass"),
  "Warm & Sunny"         = .gp("warm_and_sunny"),
  "Starry"               = .gp("starry"),
  "Spring"               = .gp("spring"),
  "Summer"               = .gp("summer"),
  "Brewer Set1"          = brewer.pal(9,  "Set1"),
  "Brewer Set2"          = brewer.pal(8,  "Set2"),
  "Brewer Set3"          = brewer.pal(12, "Set3"),
  "Brewer Dark2"         = brewer.pal(8,  "Dark2"),
  "Brewer Paired"        = brewer.pal(12, "Paired")
)
PALETTES <- PALETTES[!sapply(PALETTES, is.null)]   # drop any that failed


# ============================================================
# HELPER FUNCTIONS
# ============================================================

# Sort factor levels numerically when all levels look like numbers (0,1,2…10,11)
# Otherwise falls back to alphabetical sort.
order_factor <- function(x) {
  x  <- as.character(x)
  lvls <- unique(x)
  nums <- suppressWarnings(as.numeric(lvls))
  sorted_lvls <- if (!any(is.na(nums))) lvls[order(nums)] else sort(lvls)
  factor(x, levels = sorted_lvls)
}

# Return n named colours from a palette, extending via interpolation if needed.
get_cat_colors <- function(pal_name, levels_vec) {
  n    <- length(levels_vec)
  cols <- PALETTES[[pal_name]]
  if (is.null(cols)) cols <- hue_pal()(n)
  out  <- if (n <= length(cols)) cols[seq_len(n)] else colorRampPalette(cols)(n)
  setNames(out, levels_vec)
}

# Return a ggplot2 theme object (used with patchwork & operator).
get_theme_obj <- function(theme_name, base_size = 12) {
  switch(theme_name,
    prism   = theme_prism(base_size = base_size),
    classic = theme_classic(base_size = base_size),
    minimal = theme_minimal(base_size = base_size),
    bw      = theme_bw(base_size = base_size),
    theme()   # identity / no change
  )
}

# Add a ggplot2 theme to a single plot.
apply_theme <- function(p, theme_name, base_size = 12) {
  p + get_theme_obj(theme_name, base_size)
}

# Continuous fill scale for spatial / feature plots
cont_fill_scale <- function(scale_name) {
  switch(scale_name,
    viridis = scale_fill_viridis_c(option = "viridis"),
    plasma  = scale_fill_viridis_c(option = "plasma"),
    YlOrRd  = scale_fill_distiller(palette = "YlOrRd", direction = 1),
    Blues   = scale_fill_distiller(palette = "Blues",  direction = 1),
    scale_fill_viridis_c()
  )
}

# Continuous color scale for feature plots on reductions
cont_color_scale <- function(scale_name) {
  switch(scale_name,
    viridis = scale_color_viridis_c(option = "viridis"),
    plasma  = scale_color_viridis_c(option = "plasma"),
    YlOrRd  = scale_color_distiller(palette = "YlOrRd", direction = 1),
    Blues   = scale_color_distiller(palette = "Blues",  direction = 1),
    scale_color_viridis_c()
  )
}

error_plot <- function(msg) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, color = "firebrick") +
    theme_void()
}

save_plot <- function(p, fmt, w, h, file) {
  if (fmt == "PDF") {
    pdf(file, width = w, height = h); print(p); dev.off()
  } else {
    png(file, width = w, height = h, units = "in", res = 150)
    print(p); dev.off()
  }
}

# Compute pairwise or vs-reference stats for a long-format expression data frame.
# Returns a stat_res data frame annotated with y positions for ggprism::add_pvalue(),
# or NULL on failure (a notification is shown to the user).
compute_stats <- function(df_long, lvls, input, session) {
  tryCatch({
    # Build comparison list (NULL = all pairwise)
    comps <- if (isTRUE(input$fe_stat_compare == "ref")) {
      ref <- input$fe_stat_ref
      lapply(setdiff(lvls, ref), function(g) c(ref, g))
    } else NULL

    gdf <- df_long %>% group_by(gene)

    raw <- switch(input$fe_stat_test,
      wilcox = gdf %>%
        wilcox_test(expr ~ group, comparisons = comps) %>%
        adjust_pvalue(method = input$fe_stat_padj),

      ttest  = gdf %>%
        t_test(expr ~ group, comparisons = comps) %>%
        adjust_pvalue(method = input$fe_stat_padj)
    )

    stat_res <- raw %>%
      add_significance() %>%
      add_xy_position(x = "group", scales = "free", step.increase = 0.1)

    # Pre-format p-value columns for readable on-plot labels.
    label_col <- input$fe_stat_label

    if (label_col %in% c("p", "p.adj")) {
      stat_res$label_fmt <- vapply(
        stat_res[[label_col]],
        function(v) if (is.na(v)) "ns" else scales::pvalue(v, accuracy = 0.001, add_p = TRUE),
        character(1)
      )
      attr(stat_res, "label_col") <- "label_fmt"
    } else {
      attr(stat_res, "label_col") <- "p.adj.signif"
    }

    stat_res

  }, error = function(e) {
    showNotification(paste("Statistics error:", conditionMessage(e)),
                     type = "warning", duration = 6, session = session)
    NULL
  })
}

# Reusable per-tab theme selector widget.
# Each tab gets its own input ID so themes can differ across plots.
theme_picker_ui <- function(id, selected = "classic") {
  selectInput(id, "Theme",
              choices  = c("Classic"       = "classic",
                           "Prism"         = "prism",
                           "Minimal"       = "minimal",
                           "Black & White" = "bw"),
              selected = selected)
}


# ============================================================
# UI
# ============================================================
ui <- page_sidebar(
  title = "Visium HD Viewer",
  theme = bs_theme(
    bootswatch = "flatly",
    primary    = "#2C3E50",
    secondary  = "#18BC9C"
  ),

  sidebar = sidebar(
    width = 310,

    # ── Object loading ────────────────────────────────────────
    card(
      card_header("Load Seurat Object"),
      # Opens the native OS file-picker; path goes into the text box below.
      # No file is copied — only the path is captured.
      actionButton("browse_btn", "Browse for RDS File",
                   icon  = icon("folder-open"),
                   class = "btn-outline-primary w-100"),
      br(), br(),
      textInput("rds_path", label = NULL,
                placeholder = "…or paste full path here"),
      actionButton("load_btn", "Load Object",
                   class = "btn-primary w-100"),
      br(),
      verbatimTextOutput("obj_info", placeholder = TRUE)
    ),

    # ── Active assay ──────────────────────────────────────────
    card(
      card_header("Active Assay"),
      selectInput("active_assay", label = NULL, choices = NULL)
    ),

    # ── Global plot style ─────────────────────────────────────
    card(
      card_header("Plot Style"),
      selectInput("cat_palette", "Colour Palette",
                  choices  = names(PALETTES),
                  selected = "Project Custom (60)"),
      numericInput("base_size", "Base Font Size", 12, 8, 24, 1)
    )
  ),

  navset_card_tab(
    id = "main_tabs",

    # ── SPATIAL ──────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("map"), "Spatial"),

      layout_sidebar(
        sidebar = sidebar(
          open = TRUE,
          selectInput("sp_type", "Plot Type",
                      choices = c("Clusters / Metadata" = "dim",
                                  "Gene Expression"     = "feature")),

          conditionalPanel("input.sp_type == 'dim'",
            selectInput("sp_color_by", "Color By", choices = NULL),
            checkboxInput("sp_label", "Show Labels", TRUE),
            numericInput("sp_label_size", "Label Size", 3, 1, 10, 0.5)
          ),

          conditionalPanel("input.sp_type == 'feature'",
            selectizeInput("sp_gene", "Gene", choices = NULL,
                           options = list(placeholder = "Type gene name...")),
            selectInput("sp_color_scale", "Color Scale",
                        choices = c("Viridis" = "viridis",
                                    "Plasma"  = "plasma",
                                    "YlOrRd"  = "YlOrRd",
                                    "Blues"   = "Blues"))
          ),

          selectInput("sp_image", "Image / Sample", choices = NULL),
          numericInput("sp_pt",    "Point Size", 1,   0.1, 10, 0.1),
          numericInput("sp_alpha", "Alpha",      1.0, 0.1,  1, 0.05),
          theme_picker_ui("sp_theme"),
          hr(),
          h6("Save Plot"),
          fluidRow(
            column(6, selectInput("sp_fmt", NULL, choices = c("PNG", "PDF"))),
            column(3, numericInput("sp_w", "W", 8, 2, 30, 1)),
            column(3, numericInput("sp_h", "H", 7, 2, 30, 1))
          ),
          downloadButton("sp_dl", "Save", class = "btn-success w-100")
        ),
        plotOutput("sp_plot", height = "620px")
      )
    ),

    # ── REDUCTION / UMAP ─────────────────────────────────────
    nav_panel(
      title = tagList(icon("circle-dot"), "UMAP / Reduction"),

      layout_sidebar(
        sidebar = sidebar(
          open = TRUE,
          selectInput("dr_red", "Reduction", choices = NULL),
          selectInput("dr_type", "Plot Type",
                      choices = c("Clusters / Metadata" = "dim",
                                  "Gene Expression"     = "feature")),

          conditionalPanel("input.dr_type == 'dim'",
            selectInput("dr_color_by", "Color By", choices = NULL),
            checkboxInput("dr_label",  "Show Labels",  TRUE),
            checkboxInput("dr_repel",  "Repel Labels", TRUE)
          ),

          conditionalPanel("input.dr_type == 'feature'",
            selectizeInput("dr_gene", "Gene", choices = NULL,
                           options = list(placeholder = "Type gene name...")),
            selectInput("dr_color_scale", "Color Scale",
                        choices = c("Viridis" = "viridis",
                                    "Plasma"  = "plasma",
                                    "YlOrRd"  = "YlOrRd",
                                    "Blues"   = "Blues"))
          ),

          numericInput("dr_pt", "Point Size", 0.5, 0.1, 5, 0.1),
          theme_picker_ui("dr_theme"),
          hr(),
          h6("Save Plot"),
          fluidRow(
            column(6, selectInput("dr_fmt", NULL, choices = c("PNG", "PDF"))),
            column(3, numericInput("dr_w", "W", 8, 2, 30, 1)),
            column(3, numericInput("dr_h", "H", 7, 2, 30, 1))
          ),
          downloadButton("dr_dl", "Save", class = "btn-success w-100")
        ),
        plotOutput("dr_plot", height = "620px")
      )
    ),

    # ── FEATURE EXPRESSION ───────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-bar"), "Feature Expression"),

      layout_sidebar(
        sidebar = sidebar(
          open = TRUE,
          selectizeInput("fe_genes", "Genes (one or more)", choices = NULL,
                         multiple = TRUE,
                         options = list(placeholder = "Type gene names...")),
          selectInput("fe_type", "Plot Type",
                      choices = c("Violin"   = "violin",
                                  "Dot Plot" = "dot",
                                  "Box Plot" = "box")),
          selectInput("fe_group", "Group By", choices = NULL),

          # ── Statistics controls (violin / box only) ───────────
          hr(),
          checkboxInput("fe_show_stats", "Add statistics", FALSE),
          conditionalPanel(
            "input.fe_show_stats && input.fe_type != 'dot'",
            selectInput("fe_stat_test", "Test",
                        choices = c("Wilcoxon (non-parametric)" = "wilcox",
                                    "t-test (parametric)"       = "ttest"),
                        selected = "wilcox"),
            selectInput("fe_stat_label", "Show as",
                        choices = c("Significance stars" = "p.adj.signif",
                                    "Adjusted p-value"   = "p.adj",
                                    "Raw p-value"        = "p"),
                        selected = "p.adj.signif"),
            selectInput("fe_stat_compare", "Comparisons",
                        choices = c("All pairwise"       = "pairwise",
                                    "vs reference group" = "ref"),
                        selected = "pairwise"),
            conditionalPanel("input.fe_stat_compare == 'ref'",
              selectInput("fe_stat_ref", "Reference group", choices = NULL)
            ),
            selectInput("fe_stat_padj", "P-value adjustment",
                        choices = c("Benjamini-Hochberg" = "BH",
                                    "Bonferroni"         = "bonferroni",
                                    "None"               = "none"),
                        selected = "BH")
          ),

          theme_picker_ui("fe_theme"),
          hr(),
          h6("Save Plot"),
          fluidRow(
            column(6, selectInput("fe_fmt", NULL, choices = c("PNG", "PDF"))),
            column(3, numericInput("fe_w", "W", 10, 2, 30, 1)),
            column(3, numericInput("fe_h", "H",  7, 2, 30, 1))
          ),
          downloadButton("fe_dl", "Save", class = "btn-success w-100")
        ),
        plotOutput("fe_plot", height = "620px")
      )
    ),

    # ── COMPOSITION ──────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-pie"), "Composition"),

      layout_sidebar(
        sidebar = sidebar(
          open = TRUE,
          selectInput("co_x",    "X Axis (group by)",    choices = NULL),
          selectInput("co_fill", "Fill (cell type / cluster)", choices = NULL),
          selectInput("co_type", "Plot Type",
                      choices = c("Proportion (stacked)" = "prop",
                                  "Count (stacked)"      = "count_stack",
                                  "Count (grouped)"      = "count_group")),
          selectInput("co_sort", "Sort X Axis",
                      choices = c("As-is"                 = "none",
                                  "Alphabetical"          = "alpha",
                                  "Total cells (desc)"    = "total_desc",
                                  "Total cells (asc)"     = "total_asc")),
          checkboxInput("co_flip",  "Flip Coordinates", FALSE),
          checkboxInput("co_angle", "Rotate X Labels",  TRUE),
          theme_picker_ui("co_theme"),
          hr(),
          h6("Save Plot"),
          fluidRow(
            column(6, selectInput("co_fmt", NULL, choices = c("PNG", "PDF"))),
            column(3, numericInput("co_w", "W", 9, 2, 30, 1)),
            column(3, numericInput("co_h", "H", 6, 2, 30, 1))
          ),
          downloadButton("co_dl", "Save", class = "btn-success w-100")
        ),
        plotOutput("co_plot", height = "580px")
      )
    ),

    # ── METADATA ─────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("table"), "Metadata"),

      layout_sidebar(
        sidebar = sidebar(
          open = TRUE,
          selectInput("me_col", "Metadata Column", choices = NULL),
          selectInput("me_type", "View",
                      choices = c("Bar Chart"  = "bar",
                                  "Histogram"  = "hist",
                                  "Density"    = "density",
                                  "Data Table" = "table")),
          selectInput("me_group", "Group By (optional)",
                      choices = c("None" = "none")),
          theme_picker_ui("me_theme"),
          hr(),
          h6("Save Plot"),
          fluidRow(
            column(6, selectInput("me_fmt", NULL, choices = c("PNG", "PDF"))),
            column(3, numericInput("me_w", "W", 8, 2, 30, 1)),
            column(3, numericInput("me_h", "H", 6, 2, 30, 1))
          ),
          downloadButton("me_dl", "Save", class = "btn-success w-100")
        ),
        uiOutput("me_content")
      )
    )
  )
)


# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  rv <- reactiveValues(obj = NULL)

  # ── Browse button — opens the native OS file picker ────────
  # file.choose() is a blocking call that returns the selected path.
  # Because the Shiny process runs locally, this opens the macOS / Windows
  # file dialog on the user's own machine — no file is uploaded or copied.
  observeEvent(input$browse_btn, {
    path <- tryCatch(
      file.choose(),          # native OS dialog
      error = function(e) ""  # "" if user cancels or dialog unavailable
    )
    if (nchar(path) > 0) {
      updateTextInput(session, "rds_path", value = path)
    }
  })

  # ── Load Object ────────────────────────────────────────────
  observeEvent(input$load_btn, {
    path <- trimws(input$rds_path)
    if (nchar(path) == 0) { showNotification("Please enter a file path.", type = "warning"); return() }
    if (!file.exists(path)) { showNotification("File not found.", type = "error"); return() }

    id <- showNotification("Loading… large objects may take a minute.", duration = NULL, type = "message")
    on.exit(removeNotification(id))

    tryCatch({
      o <- readRDS(path)
      if (!inherits(o, "Seurat")) { showNotification("Not a Seurat object.", type = "error"); return() }
      rv$obj <- o
      populate_controls(o)
      showNotification("Loaded successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)), type = "error")
    })
  })

  # ── Populate UI Controls ───────────────────────────────────
  populate_controls <- function(o) {
    assays     <- Assays(o)
    meta_cols  <- colnames(o@meta.data)
    images     <- names(o@images)
    reductions <- names(o@reductions)
    features   <- rownames(o)

    def_assay <- DefaultAssay(o)
    def_clust <- if ("seurat_clusters" %in% meta_cols) "seurat_clusters" else meta_cols[1]
    def_red   <- if ("umap.sketch" %in% reductions) "umap.sketch"
                 else if ("umap"   %in% reductions) "umap"
                 else if (length(reductions) > 0)   reductions[1]
                 else NULL
    def_ident <- if ("orig.ident" %in% meta_cols) "orig.ident" else meta_cols[1]

    updateSelectInput(session, "active_assay", choices = assays, selected = def_assay)

    # Spatial
    updateSelectInput(session, "sp_color_by", choices = meta_cols, selected = def_clust)
    updateSelectInput(session, "sp_image",
                      choices = if (length(images) > 0) images else c("(no images)" = "NONE"))
    updateSelectizeInput(session, "sp_gene", choices = features, server = TRUE)

    # Reduction
    updateSelectInput(session, "dr_red",
                      choices = if (length(reductions) > 0) reductions else c("(none)" = "NONE"),
                      selected = def_red)
    updateSelectInput(session, "dr_color_by", choices = meta_cols, selected = def_clust)
    updateSelectizeInput(session, "dr_gene",  choices = features, server = TRUE)

    # Feature expression
    updateSelectizeInput(session, "fe_genes", choices = features, server = TRUE)
    updateSelectInput(session, "fe_group",    choices = meta_cols, selected = def_clust)
    # Initialise reference-group choices from the default cluster column
    def_clust_lvls <- levels(order_factor(o@meta.data[[def_clust]]))
    updateSelectInput(session, "fe_stat_ref", choices = def_clust_lvls,
                      selected = def_clust_lvls[1])

    # Composition
    updateSelectInput(session, "co_x",    choices = meta_cols, selected = def_ident)
    updateSelectInput(session, "co_fill", choices = meta_cols, selected = def_clust)

    # Metadata
    updateSelectInput(session, "me_col",   choices = meta_cols, selected = meta_cols[1])
    updateSelectInput(session, "me_group", choices = c("None" = "none", meta_cols))
  }

  observeEvent(input$active_assay, {
    o <- rv$obj
    req(o, input$active_assay, input$active_assay %in% Assays(o))
    DefaultAssay(rv$obj) <- input$active_assay
  })

  # Keep reference-group choices in sync with the Group By selector
  observeEvent(input$fe_group, {
    o <- rv$obj
    req(o, input$fe_group, input$fe_group %in% colnames(o@meta.data))
    lvls <- levels(order_factor(o@meta.data[[input$fe_group]]))
    updateSelectInput(session, "fe_stat_ref", choices = lvls, selected = lvls[1])
  })

  # ── Object Info ────────────────────────────────────────────
  output$obj_info <- renderText({
    o <- rv$obj
    if (is.null(o)) return("No object loaded.")

    assays     <- Assays(o)
    has_images <- length(o@images) > 0

    # ── Detect Visium HD bin size from assay names ──────────
    # Assay names like "Spatial.008um" or "Spatial.016um"
    bin_sizes <- character(0)
    for (a in assays) {
      m <- regmatches(a, regexpr("(?i)(?<=\\.)0*(\\d+)um", a, perl = TRUE))
      if (length(m) > 0) {
        num <- as.integer(sub("(?i)um$", "", m[[1]], perl = TRUE))
        bin_sizes <- c(bin_sizes, paste0(num, " µm"))
      }
    }
    bin_sizes    <- unique(bin_sizes)
    is_visium_hd <- length(bin_sizes) > 0

    # ── Choose correct terminology ───────────────────────────
    # Visium HD  → "Bins"   (binned square grid)
    # Visium std → "Spots"  (circular capture spots)
    # scRNA-seq  → "Cells"
    spot_label <- if (is_visium_hd) "Bins" else if (has_images) "Spots" else "Cells"

    # ── Fixed-width label helper (left-aligned, 12 chars) ───
    lbl <- function(x) formatC(x, width = -12, flag = "-")

    # ── Build output lines ───────────────────────────────────
    lines <- c(
      if (is_visium_hd) paste0(lbl("Binning:"), paste(bin_sizes, collapse = ", ")),
      paste0(lbl(paste0(spot_label, ":")),  format(ncol(o), big.mark = ",")),
      paste0(lbl("Genes:"),       format(nrow(o), big.mark = ",")),
      paste0(lbl("Assays:"),      paste(assays, collapse = ", ")),
      paste0(lbl("Images:"),      if (has_images)           paste(names(o@images),     collapse = ", ") else "none"),
      paste0(lbl("Reductions:"),  if (length(o@reductions)) paste(names(o@reductions), collapse = ", ") else "none")
    )

    paste(lines, collapse = "\n")
  })


  # ============================================================
  # SPATIAL TAB
  # ============================================================
  sp_plot_r <- reactive({
    o <- rv$obj
    req(o)

    if (length(o@images) == 0 || isTRUE(input$sp_image == "NONE"))
      return(error_plot("No spatial images in this object."))

    img <- input$sp_image
    pt  <- input$sp_pt
    al  <- input$sp_alpha

    if (input$sp_type == "dim") {
      req(input$sp_color_by)
      color_by <- input$sp_color_by

      tryCatch({
        # Apply numeric-aware factor ordering
        o@meta.data[[color_by]] <- order_factor(o@meta.data[[color_by]])
        lvls <- levels(o@meta.data[[color_by]])
        cols <- get_cat_colors(input$cat_palette, lvls)

        p <- suppressWarnings(
          SpatialDimPlot(o,
                         group.by       = color_by,
                         images         = img,
                         label          = input$sp_label,
                         label.size     = input$sp_label_size,
                         pt.size.factor = pt,
                         alpha          = al) +
            scale_color_manual(values = cols) +
            scale_fill_manual(values  = cols) +
            theme(legend.position = "right") +
            ggtitle(color_by)
        )
        apply_theme(p, input$sp_theme, input$base_size)

      }, error = function(e) error_plot(conditionMessage(e)))

    } else {
      gene <- input$sp_gene
      req(gene, nchar(gene) > 0)

      tryCatch({
        p <- SpatialFeaturePlot(o, features = gene, images = img,
                                pt.size.factor = pt, alpha = al) +
          cont_fill_scale(input$sp_color_scale)
        apply_theme(p, input$sp_theme, input$base_size)
      }, error = function(e) error_plot(conditionMessage(e)))
    }
  })

  output$sp_plot <- renderPlot({ sp_plot_r() })
  output$sp_dl   <- downloadHandler(
    filename = function() paste0("spatial_plot.", tolower(input$sp_fmt)),
    content  = function(file) save_plot(sp_plot_r(), input$sp_fmt, input$sp_w, input$sp_h, file)
  )


  # ============================================================
  # REDUCTION / UMAP TAB
  # ============================================================
  dr_plot_r <- reactive({
    o <- rv$obj
    req(o)

    if (length(o@reductions) == 0 || isTRUE(input$dr_red == "NONE"))
      return(error_plot("No reductions available in this object."))

    red <- input$dr_red
    pt  <- input$dr_pt

    if (input$dr_type == "dim") {
      req(input$dr_color_by)
      color_by <- input$dr_color_by

      tryCatch({
        o@meta.data[[color_by]] <- order_factor(o@meta.data[[color_by]])
        lvls <- levels(o@meta.data[[color_by]])
        cols <- get_cat_colors(input$cat_palette, lvls)

        p <- suppressWarnings(
          DimPlot(o,
                  reduction = red,
                  group.by  = color_by,
                  label     = input$dr_label,
                  repel     = input$dr_repel,
                  pt.size   = pt,
                  raster    = FALSE) +
            scale_color_manual(values = cols) +
            ggtitle(paste(toupper(red), "–", color_by))
        )
        apply_theme(p, input$dr_theme, input$base_size)

      }, error = function(e) error_plot(conditionMessage(e)))

    } else {
      gene <- input$dr_gene
      req(gene, nchar(gene) > 0)

      tryCatch({
        p <- FeaturePlot(o, features = gene, reduction = red,
                         pt.size = pt, raster = FALSE) +
          cont_color_scale(input$dr_color_scale)
        apply_theme(p, input$dr_theme, input$base_size)
      }, error = function(e) error_plot(conditionMessage(e)))
    }
  })

  output$dr_plot <- renderPlot({ dr_plot_r() })
  output$dr_dl   <- downloadHandler(
    filename = function() paste0("reduction_plot.", tolower(input$dr_fmt)),
    content  = function(file) save_plot(dr_plot_r(), input$dr_fmt, input$dr_w, input$dr_h, file)
  )


  # ============================================================
  # FEATURE EXPRESSION TAB
  # ============================================================
  fe_plot_r <- reactive({
    o <- rv$obj
    req(o, input$fe_genes, length(input$fe_genes) > 0, input$fe_group)

    genes <- input$fe_genes
    grp   <- input$fe_group

    tryCatch({
      # Pre-sort group factor
      o@meta.data[[grp]] <- order_factor(o@meta.data[[grp]])
      lvls <- levels(o@meta.data[[grp]])
      cols <- get_cat_colors(input$cat_palette, lvls)

      # ── Helper: extract a long data frame for violin / box ──
      get_expr_long <- function() {
        ad <- GetAssayData(o, layer = "data")
        ok <- genes[genes %in% rownames(ad)]
        if (length(ok) == 0) stop("None of the requested genes were found in the active assay.")
        df <- data.frame(
          group = o@meta.data[[grp]],
          t(as.matrix(ad[ok, , drop = FALSE])),
          check.names = FALSE
        )
        df_long <- pivot_longer(df, cols = -group, names_to = "gene", values_to = "expr")
        df_long$group <- factor(df_long$group, levels = lvls)
        df_long
      }

      if (input$fe_type == "violin") {
        # Build as pure ggplot2 violin so stats brackets attach cleanly to facets
        df_long <- get_expr_long()

        p <- ggplot(df_long, aes(x = group, y = expr, fill = group)) +
          geom_violin(trim = FALSE, scale = "width") +
          geom_boxplot(width = 0.08, fill = "white",
                       outlier.size = 0.5, outlier.alpha = 0.4) +
          scale_fill_manual(values = cols) +
          facet_wrap(~gene, scales = "free_y") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none") +
          labs(x = grp, y = "Expression")

        if (isTRUE(input$fe_show_stats)) {
          stat_res <- compute_stats(df_long, lvls, input, session)
          if (!is.null(stat_res))
            p <- p + add_pvalue(stat_res,
                                label        = attr(stat_res, "label_col"),
                                tip.length   = 0.01,
                                bracket.size = 0.4,
                                label.size   = 3.2,
                                inherit.aes  = FALSE)
        }

        apply_theme(p, input$fe_theme, input$base_size)

      } else if (input$fe_type == "dot") {
        p <- DotPlot(o, features = genes, group.by = grp) +
          scale_color_viridis_c() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          coord_flip()
        apply_theme(p, input$fe_theme, input$base_size)

      } else if (input$fe_type == "box") {
        df_long <- get_expr_long()

        p <- ggplot(df_long, aes(x = group, y = expr, fill = group)) +
          geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3) +
          scale_fill_manual(values = cols) +
          facet_wrap(~gene, scales = "free_y") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none") +
          labs(x = grp, y = "Expression")

        if (isTRUE(input$fe_show_stats)) {
          stat_res <- compute_stats(df_long, lvls, input, session)
          if (!is.null(stat_res))
            p <- p + add_pvalue(stat_res,
                                label        = attr(stat_res, "label_col"),
                                tip.length   = 0.01,
                                bracket.size = 0.4,
                                label.size   = 3.2,
                                inherit.aes  = FALSE)
        }

        apply_theme(p, input$fe_theme, input$base_size)
      }

    }, error = function(e) error_plot(conditionMessage(e)))
  })

  output$fe_plot <- renderPlot({ fe_plot_r() })
  output$fe_dl   <- downloadHandler(
    filename = function() paste0("expression_plot.", tolower(input$fe_fmt)),
    content  = function(file) save_plot(fe_plot_r(), input$fe_fmt, input$fe_w, input$fe_h, file)
  )


  # ============================================================
  # COMPOSITION TAB
  # ============================================================
  co_plot_r <- reactive({
    o <- rv$obj
    req(o, input$co_x, input$co_fill)

    x_var    <- input$co_x
    fill_var <- input$co_fill
    meta     <- o@meta.data

    tryCatch({
      # Build counts
      df <- meta %>%
        mutate(
          .x    = as.character(.data[[x_var]]),
          .fill = as.character(.data[[fill_var]])
        ) %>%
        group_by(.x, .fill) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(.x) %>%
        mutate(prop = n / sum(n)) %>%
        ungroup()

      # Sort X axis
      x_order <- switch(input$co_sort,
        alpha      = sort(unique(df$.x)),
        total_desc = df %>% group_by(.x) %>% summarise(tot = sum(n), .groups="drop") %>%
                       arrange(desc(tot)) %>% pull(.x),
        total_asc  = df %>% group_by(.x) %>% summarise(tot = sum(n), .groups="drop") %>%
                       arrange(tot) %>% pull(.x),
        unique(df$.x)   # none / as-is
      )
      df$.x <- factor(df$.x, levels = x_order)

      # Numeric-aware ordering for fill
      fill_lvls <- order_factor(df$.fill) |> levels()
      df$.fill  <- factor(df$.fill, levels = fill_lvls)
      cols      <- get_cat_colors(input$cat_palette, fill_lvls)

      # Build plot
      if (input$co_type == "prop") {
        p <- ggplot(df, aes(x = .x, y = prop, fill = .fill)) +
          geom_col(position = "stack", width = 0.8) +
          scale_y_continuous(labels = percent_format(accuracy = 1)) +
          labs(x = x_var, y = "Proportion", fill = fill_var)

      } else if (input$co_type == "count_stack") {
        p <- ggplot(df, aes(x = .x, y = n, fill = .fill)) +
          geom_col(position = "stack", width = 0.8) +
          labs(x = x_var, y = "Cell Count", fill = fill_var)

      } else {
        p <- ggplot(df, aes(x = .x, y = n, fill = .fill)) +
          geom_col(position = position_dodge(width = 0.85), width = 0.8) +
          labs(x = x_var, y = "Cell Count", fill = fill_var)
      }

      p <- p + scale_fill_manual(values = cols)

      if (input$co_angle)
        p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

      if (input$co_flip)
        p <- p + coord_flip()

      apply_theme(p, input$co_theme, input$base_size)

    }, error = function(e) error_plot(conditionMessage(e)))
  })

  output$co_plot <- renderPlot({ co_plot_r() })
  output$co_dl   <- downloadHandler(
    filename = function() paste0("composition_plot.", tolower(input$co_fmt)),
    content  = function(file) save_plot(co_plot_r(), input$co_fmt, input$co_w, input$co_h, file)
  )


  # ============================================================
  # METADATA TAB
  # ============================================================
  me_plot_r <- reactive({
    o <- rv$obj
    req(o, input$me_col, input$me_type != "table")

    meta   <- o@meta.data
    col    <- input$me_col
    grp    <- if (input$me_group != "none") input$me_group else NULL
    vals   <- meta[[col]]
    is_num <- is.numeric(vals)

    tryCatch({
      if (input$me_type == "bar") {
        if (!is.null(grp)) {
          df <- meta %>%
            group_by(across(all_of(c(grp, col)))) %>%
            summarise(n = n(), .groups = "drop")
          df[[grp]] <- order_factor(df[[grp]])
          df[[col]] <- order_factor(df[[col]])
          fill_lvls <- levels(df[[col]])
          cols      <- get_cat_colors(input$cat_palette, fill_lvls)

          p <- ggplot(df, aes(x = .data[[grp]], y = n,
                              fill = factor(.data[[col]], levels = fill_lvls))) +
            geom_col(position = "fill") +
            scale_y_continuous(labels = percent_format()) +
            scale_fill_manual(values = cols) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = grp, y = "Proportion", fill = col)
        } else {
          vals_f    <- order_factor(vals)
          fill_lvls <- levels(vals_f)
          cols      <- get_cat_colors(input$cat_palette, fill_lvls)
          df        <- as.data.frame(table(vals_f))
          colnames(df) <- c("value", "count")
          df$value  <- factor(df$value, levels = fill_lvls)

          p <- ggplot(df, aes(x = value, y = count, fill = value)) +
            geom_col(show.legend = FALSE) +
            scale_fill_manual(values = cols) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = col, y = "Count")
        }

      } else if (input$me_type == "hist") {
        if (!is_num) {
          vals_f    <- order_factor(vals)
          fill_lvls <- levels(vals_f)
          cols      <- get_cat_colors(input$cat_palette, fill_lvls)
          df        <- as.data.frame(table(vals_f))
          colnames(df) <- c("value", "count")
          df$value  <- factor(df$value, levels = fill_lvls)
          p <- ggplot(df, aes(x = value, y = count, fill = value)) +
            geom_col(show.legend = FALSE) +
            scale_fill_manual(values = cols) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = col, y = "Count")
        } else {
          p <- ggplot(meta, aes(x = .data[[col]]))
          if (!is.null(grp)) {
            grp_lvls <- levels(order_factor(meta[[grp]]))
            cols     <- get_cat_colors(input$cat_palette, grp_lvls)
            meta[[grp]] <- factor(meta[[grp]], levels = grp_lvls)
            p <- p + geom_histogram(aes(fill = .data[[grp]]),
                                    position = "identity", alpha = 0.6, bins = 40) +
              scale_fill_manual(values = cols) + labs(fill = grp)
          } else {
            p <- p + geom_histogram(fill = "#2C3E50", bins = 40)
          }
          p <- p + labs(x = col, y = "Count")
        }

      } else if (input$me_type == "density") {
        if (!is_num)
          return(error_plot("Density plot requires a numeric column."))
        p <- ggplot(meta, aes(x = .data[[col]]))
        if (!is.null(grp)) {
          grp_lvls <- levels(order_factor(meta[[grp]]))
          cols     <- get_cat_colors(input$cat_palette, grp_lvls)
          meta[[grp]] <- factor(meta[[grp]], levels = grp_lvls)
          p <- p + geom_density(aes(fill = .data[[grp]]), alpha = 0.5) +
            scale_fill_manual(values = cols) + labs(fill = grp)
        } else {
          p <- p + geom_density(fill = "#18BC9C", alpha = 0.7)
        }
        p <- p + labs(x = col)
      }

      apply_theme(p, input$me_theme, input$base_size)

    }, error = function(e) error_plot(conditionMessage(e)))
  })

  output$me_content <- renderUI({
    if (input$me_type == "table") DTOutput("me_table")
    else                          plotOutput("me_plot", height = "560px")
  })

  output$me_plot  <- renderPlot({ me_plot_r() })

  output$me_table <- renderDT({
    o <- rv$obj; req(o)
    datatable(o@meta.data,
              options = list(pageLength = 25, scrollX = TRUE),
              rownames = TRUE)
  })

  output$me_dl <- downloadHandler(
    filename = function() paste0("metadata_plot.", tolower(input$me_fmt)),
    content  = function(file) {
      p <- me_plot_r(); req(p)
      save_plot(p, input$me_fmt, input$me_w, input$me_h, file)
    }
  )
}

shinyApp(ui, server)
