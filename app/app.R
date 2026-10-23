# Maize & rice leaf developmental gradient explorer
# Data: Wang et al. (2014) Nat Biotechnol 32:1158 (GEO GSE54274).
# Start with:  shiny::runApp("app")   (from the project root)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(DT)
  library(ggplot2)
})

DEFAULT_GENES <- c("GRMZM2G083841", "GRMZM2G306345", "GRMZM2G121878")   # pep1, pdk1, cah1

inline <- function(...) div(class = "d-flex flex-wrap gap-3 align-items-end", ...)
small_note <- function(...) div(class = "text-muted small", ...)
dt <- function(df, digits = 3, pageLength = 10, selection = "none", ...) {
  for (cl in grep("description", names(df), value = TRUE))      # keep rows one line high
    df[[cl]] <- ifelse(nchar(df[[cl]]) > 70, paste0(substr(df[[cl]], 1, 68), "..."), df[[cl]])
  num <- names(df)[vapply(df, is.double, TRUE)]
  big <- nrow(df) > pageLength
  w <- datatable(df, rownames = FALSE, selection = selection, filter = if (big) "top" else "none",
                 class = "compact stripe hover nowrap",
                 extensions = "Buttons",
                 options = list(pageLength = pageLength, scrollX = TRUE, dom = if (big) "Bfrtip" else "Bt",
                                buttons = c("copy", "csv")), ...)
  if (length(num)) w <- formatSignif(w, num, digits)
  w
}

ui <- page_navbar(
  title = "MaizeExpressG",
  id = "nav",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2E7D32",
                   "navbar-bg" = "#1B5E20", base_font = font_collection("system-ui", "sans-serif")),
  fillable = FALSE,
  footer = citation_footer(),
  sidebar = sidebar(
    width = 330,
    radioButtons("species", "Species", SPECIES, inline = TRUE),
    radioButtons("idtype", "Maize gene IDs", c("B73 v5 (Zm00001eb)" = "v5", "B73 v3 (GRMZM)" = "v3"),
                 inline = TRUE),
    selectInput("dataset", "Expression data", DATASETS),
    selectizeInput("genes", "Genes", choices = NULL, multiple = TRUE,
                   options = list(placeholder = "ID, symbol or description...", maxOptions = 60)),
    textAreaInput("paste_ids", NULL, rows = 2,
                  placeholder = "...or paste IDs/symbols (GRMZM, Zm00001d, Zm00001eb, LOC_Os, pep1)"),
    div(class = "d-flex gap-2",
        actionButton("add_pasted", "Add pasted", class = "btn-sm btn-outline-primary"),
        actionButton("clear_genes", "Clear", class = "btn-sm btn-outline-secondary")),
    uiOutput("gene_msg"),
    hr(class = "my-1"),
    radioButtons("scale", "Scale", c("log2(x + 1)" = "log2", "Linear (FPKM / RPKM)" = "linear"),
                 inline = TRUE),
    selectInput("xaxis", "Gradient axis", XAXES),
    conditionalPanel("input.dataset == 'geo'",
      sliderInput("min_reads", "Exclude libraries with fewer mapped reads than",
                  min = 0, max = 2e6, value = DEFAULT_MIN_READS, step = 50000, sep = ",")),
    uiOutput("lib_note")
  ),

  # ---- Profiles -------------------------------------------------------------
  nav_panel("Profiles", icon = icon("chart-line"),
    card(full_screen = TRUE,
      card_header(inline(
        checkboxInput("show_reps", "Show replicate libraries", TRUE),
        radioButtons("err", "Error bars", c(SE = "se", SD = "sd", None = "none"), inline = TRUE),
        checkboxInput("zscore", "z-score each gene (compare shapes)", FALSE),
        checkboxInput("show_ligule", "Mark leaf-2 ligule (maize)", TRUE))),
      plotlyOutput("profile_plot", height = "480px"),
      uiOutput("missing_note")),
    card(card_header("Gene information"), DTOutput("gene_info")),
    card(card_header(div(class = "d-flex justify-content-between",
                         "Section means (current scale)",
                         downloadButton("dl_profiles", "Download library-level values", class = "btn-sm"))),
         DTOutput("means_table"))
  ),

  # ---- Model fitting --------------------------------------------------------
  nav_panel("Fit models", icon = icon("wave-square"),
    card(
      card_header(inline(
        selectInput("fit_gene", "Gene", choices = NULL, width = "260px"),
        checkboxGroupInput("fit_models", "Models", MODELS, selected = c("linear", "quadratic", "segmented"),
                           inline = TRUE),
        selectInput("fit_focus", "Focus model", MODELS, width = "220px"))),
      card_body(inline(
        sliderInput("fit_range", "Fit range (gradient axis units)", min = 0, max = 15, value = c(0, 15),
                    step = 0.5, width = "320px"),
        radioButtons("fit_on", "Fit to", c("Replicate libraries" = "reps", "Section means" = "means"),
                     inline = TRUE),
        numericInput("spline_df", "Spline df", 3, min = 2, max = 8, width = "100px"),
        numericInput("ci_level", "CI level", 0.95, min = 0.5, max = 0.999, step = 0.01, width = "100px"))),
      plotlyOutput("fit_plot", height = "430px"),
      uiOutput("fit_msg")),
    layout_columns(col_widths = c(7, 5),
      card(card_header("Model comparison (lower AIC/BIC is better)"), DTOutput("fit_compare"),
           small_note("Lack-of-fit p compares each curve with the section-means model using the ",
                      "replicate libraries (pure error). Small values mean the curve misses real ",
                      "section-to-section structure.")),
      card(card_header("Focus model"), uiOutput("fit_derived"), DTOutput("fit_coef"),
           div(class = "d-flex gap-2 mt-2",
               downloadButton("dl_fit_coef", "Coefficients", class = "btn-sm"),
               downloadButton("dl_fit_data", "Data + fitted", class = "btn-sm")))),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Residuals of focus model"), plotlyOutput("fit_resid", height = "300px")),
      card(card_header("Nested tests"), DTOutput("fit_nested"))),
    card(card_header("Do the selected genes share one slope? (ANCOVA, linear in x, same range)"),
         uiOutput("ancova_msg"), DTOutput("ancova_table"))
  ),

  # ---- Genome-wide screen ---------------------------------------------------
  nav_panel("Gradient screen", icon = icon("magnifying-glass-chart"),
    card(
      card_header("Fit a model to every expressed gene and rank by significance"),
      card_body(inline(
        selectInput("scr_model", "Model", MODELS[1:4], width = "180px"),
        sliderInput("scr_range", "Range (gradient axis units)", min = 0, max = 15, value = c(0, 15),
                    step = 0.5, width = "300px"),
        numericInput("scr_minexpr", "Min. max section mean (FPKM)", 1, min = 0, step = 0.5, width = "170px"),
        numericInput("scr_fdr", "FDR cutoff", 0.05, min = 0, max = 1, step = 0.01, width = "110px"),
        actionButton("scr_run", "Run screen", class = "btn-primary", icon = icon("play")))),
      uiOutput("scr_summary")),
    layout_columns(col_widths = c(6, 6),
      card(full_screen = TRUE, card_header("Click a point to plot the gene"),
           plotlyOutput("scr_volcano", height = "420px")),
      card(full_screen = TRUE, card_header("Selected genes"), plotlyOutput("scr_gene_plot", height = "420px"))),
    card(card_header(div(class = "d-flex justify-content-between align-items-center",
                         "Results (select rows to plot them)",
                         div(class = "d-flex gap-2",
                             actionButton("scr_add", "Add selected to gene list", class = "btn-sm btn-outline-primary"),
                             downloadButton("dl_screen", "Download all results", class = "btn-sm")))),
         DTOutput("scr_table"))
  ),

  # ---- Base vs tip patterns ------------------------------------------------
  nav_panel("Base vs tip", icon = icon("arrows-up-down"),
    card(
      card_header("Genes high at the base that decline to the tip, and low at the base that rise to it"),
      card_body(inline(
        numericInput("pf_fold", "Min. fold change, base vs tip", 4, min = 1, step = 0.5, width = "190px"),
        numericInput("pf_fdr", "FDR, linear trend", 0.01, min = 0, max = 1, step = 0.005, width = "130px"),
        numericInput("pf_minexpr", "Min. max section mean (FPKM)", 5, min = 0, step = 1, width = "190px"),
        numericInput("pf_tipfold", "Tip within x-fold of other extreme", 2, min = 1, step = 0.5, width = "210px"),
        sliderInput("pf_end", "Base / tip window (% of sections)", min = 10, max = 40, value = 20, step = 5,
                    width = "230px"),
        actionButton("pf_run", "Update", class = "btn-primary", icon = icon("rotate")))),
      uiOutput("pf_summary")),
    layout_columns(col_widths = c(6, 6),
      card(full_screen = TRUE, card_header("All genes in each pattern (z-scored section means; line = median, band = 10-90%)"),
           plotOutput("pf_overview", height = "380px")),
      card(full_screen = TRUE, card_header("Selected genes (top 6 of the list if none selected)"),
           plotlyOutput("pf_gene_plot", height = "380px"))),
    card(card_header(div(class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
                         radioButtons("pf_which", NULL, PATTERNS, inline = TRUE),
                         div(class = "d-flex gap-2",
                             actionButton("pf_add", "Add selected to gene list", class = "btn-sm btn-outline-primary"),
                             downloadButton("dl_pf", "Download this list", class = "btn-sm")))),
         DTOutput("pf_table"),
         small_note("Declining: highest section in the basal third, tip within the set fold of the lowest section. ",
                    "Rising: lowest section in the basal third, tip within the set fold of the highest section ",
                    "(rise-then-plateau genes count; genes peaking mid-leaf do not). Both need the fold change ",
                    "between base and tip windows and a linear trend (Gradient screen model) at the FDR. ",
                    "shape: 'gradual' when |Spearman rho| >= 0.8; half_change_cm: where half the change is reached."))
  ),

  # ---- Co-expression --------------------------------------------------------
  nav_panel("Co-expression", icon = icon("diagram-project"),
    card(
      card_header("Find genes whose gradient profile matches a query"),
      card_body(inline(
        radioButtons("cx_type", "Query", c("A gene" = "gene", "A metabolite / enzyme trait" = "metab",
                                          "A template shape" = "template"), inline = TRUE),
        conditionalPanel("input.cx_type == 'gene'", selectInput("cx_gene", "Query gene", NULL, width = "260px")),
        conditionalPanel("input.cx_type == 'metab'", selectizeInput("cx_metab", "Trait", NULL, width = "300px")),
        conditionalPanel("input.cx_type == 'template'",
          inline(selectInput("cx_template", "Shape", c("increasing", "decreasing", "peak", "dip",
                                                       "step_up", "step_down"), width = "150px"),
                 sliderInput("cx_at", "Peak / step at section", 1, 15, 8, step = 1, width = "220px"))),
        conditionalPanel("input.dataset == 'udm'",
          radioButtons("cx_target", "Search in", SPECIES, inline = TRUE)))),
      card_body(inline(
        radioButtons("cx_method", "Correlation", c(Pearson = "pearson", Spearman = "spearman"), inline = TRUE),
        radioButtons("cx_dir", "Direction", c(Positive = "pos", Negative = "neg", Both = "both"), inline = TRUE),
        numericInput("cx_minexpr", "Min. max section mean", 1, min = 0, step = 0.5, width = "150px"),
        numericInput("cx_n", "Show top", 100, min = 10, max = 5000, step = 10, width = "110px"),
        actionButton("cx_run", "Search", class = "btn-primary", icon = icon("play")))),
      small_note("Correlations use section means (15 maize / 11 rice / 15 UDM points); p-values are ",
                 "indicative only. Metabolite traits come from separate plants sampled in the same sections. ",
                 "UDM profiles are fitted smooth curves, so their correlations run very high - rank by r, ignore p.")),
    layout_columns(col_widths = c(6, 6),
      card(full_screen = TRUE, card_header("Query (black) and best matches (z-scored section means)"),
           plotlyOutput("cx_plot", height = "420px")),
      card(card_header(div(class = "d-flex justify-content-between align-items-center", "Matches",
                           div(class = "d-flex gap-2",
                               actionButton("cx_add", "Add selected to gene list", class = "btn-sm btn-outline-primary"),
                               downloadButton("dl_cx", "Download", class = "btn-sm")))),
           DTOutput("cx_table")))
  ),

  # ---- Heatmap --------------------------------------------------------------
  nav_panel("Heatmap & clusters", icon = icon("table-cells"),
    card(
      card_body(inline(
        radioButtons("hm_source", "Genes", c("Current gene list" = "list", "Annotation search" = "keyword",
                                            "C4 shuttle genes (Supp. Table 7)" = "c4",
                                            "Significant genes from last screen" = "screen"), inline = TRUE),
        conditionalPanel("input.hm_source == 'keyword'",
          textInput("hm_keyword", "Regex on symbol/description", "photosystem|rubisco|chlorophyll a-b binding",
                    width = "330px")),
        numericInput("hm_max", "Max genes", 300, min = 2, max = 3000, width = "110px"),
        numericInput("hm_minexpr", "Min. max section mean", 1, min = 0, step = 0.5, width = "150px"),
        checkboxInput("hm_z", "z-score rows", TRUE),
        numericInput("hm_k", "Clusters (k)", 4, min = 1, max = 20, width = "100px"),
        actionButton("hm_run", "Draw", class = "btn-primary", icon = icon("play")))),
      uiOutput("hm_msg")),
    layout_columns(col_widths = c(7, 5),
      card(full_screen = TRUE, card_header("Section means, rows clustered (1 - Pearson r, average linkage)"),
           plotlyOutput("hm_plot", height = "640px")),
      card(full_screen = TRUE, card_header("Cluster mean profiles"), plotOutput("hm_clusters", height = "640px"))),
    card(card_header(div(class = "d-flex justify-content-between align-items-center", "Cluster membership",
                         downloadButton("dl_hm", "Download", class = "btn-sm"))),
         DTOutput("hm_table"))
  ),

  # ---- Maize vs rice --------------------------------------------------------
  nav_panel("Maize <-> Rice", icon = icon("right-left"),
    card(card_header(inline(
           "Selected genes and their orthologs (Supp. Table 3)",
           checkboxInput("or_z", "z-score", TRUE))),
         layout_columns(col_widths = c(6, 6),
           div(h6("Unified gradient (paper UDM)"), plotlyOutput("or_udm_plot", height = "420px")),
           div(h6("Native gradients (GEO section means)"), plotlyOutput("or_native_plot", height = "420px"))),
         small_note("Left: the paper's unified developmental model (UDM) puts both species on one ",
                    "15-point axis. Right: GEO replicate data on each species' own axis (maize solid, rice dashed).")),
    layout_columns(col_widths = c(8, 4),
      card(card_header(div(class = "d-flex justify-content-between align-items-center",
                           "All ortholog pairs: correlation of UDM profiles (log2)",
                           downloadButton("dl_pairs", "Download", class = "btn-sm"))),
           DTOutput("pairs_table")),
      card(card_header("Selected pair"), plotlyOutput("pair_plot", height = "300px"),
           plotlyOutput("pairs_hist", height = "220px")))
  ),

  # ---- Metabolites ----------------------------------------------------------
  nav_panel("Metabolites", icon = icon("flask"),
    card(card_header(inline(
           selectizeInput("met_traits", "Traits", NULL, multiple = TRUE, width = "520px"),
           checkboxInput("met_z", "z-score (traits have different units)", TRUE),
           actionButton("met_to_cx", "Find genes correlated with the first trait", class = "btn-sm btn-outline-primary"))),
         plotlyOutput("met_plot", height = "460px"),
         small_note("Supplementary Tables 1-2: mean +/- SD per section (maize n = 6, rice n = 4 replicates). ",
                    "Maize sections -4..-1 and 1..11 in the original table are shown as sections 1-15.")),
    card(card_header("Values"), DTOutput("met_table"))
  ),

  # ---- QC -------------------------------------------------------------------
  nav_panel("Libraries & QC", icon = icon("vial"),
    layout_columns(col_widths = c(6, 6),
      card(full_screen = TRUE, card_header("PCA of libraries (log2 FPKM, 2,000 most variable genes)"),
           plotlyOutput("qc_pca", height = "420px")),
      card(full_screen = TRUE, card_header("Mapped reads per library (Supp. Table 17)"),
           plotlyOutput("qc_depth", height = "420px"))),
    card(card_header("Libraries"), DTOutput("qc_table"),
         small_note("r_vs_other_reps: Pearson correlation (log2 FPKM, genes with mean FPKM > 5) between a ",
                    "library and the mean of the other replicates of its section."))
  ),

  nav_panel("About", icon = icon("circle-info"), card(card_body(about_page()))),
  nav_spacer(),
  nav_item(tags$a(href = "https://doi.org/10.1038/nbt.3019", target = "_blank", "Wang et al. 2014"))
)

server <- function(input, output, session) {

  # ---- shared state -------------------------------------------------------
  sp <- reactive(input$species)
  ds <- reactive(input$dataset)
  min_reads <- reactive(if (ds() == "geo") input$min_reads else 0)
  xaxis <- reactive({
    xa <- input$xaxis
    if (xa == "pos_ligule" && sp() != "maize") "pos_cm" else xa
  })
  xlab <- reactive(x_label(xaxis(), ds()))
  xunit <- reactive(if (ds() == "udm") "UDM position" else switch(xaxis(),
    pos_cm = "cm", pos_ligule = "cm", section = "section", sec_label = "section label"))
  ylab <- reactive(if (input$scale == "log2") "log2(expression + 1)" else
    if (ds() == "geo") "FPKM" else "RPKM")
  ligule_x <- reactive({
    if (!isTRUE(input$show_ligule) || sp() != "maize" || ds() == "udm") return(NULL)
    switch(xaxis(), pos_cm = 4, pos_ligule = 0, section = 4.5, sec_label = 5, NULL)
  })
  genes <- reactive(input$genes %||% character())
  last_species <- reactiveVal("maize")

  idtype <- reactive(input$idtype %||% "v5")
  sn <- function(ids) short_name(ids, idtype())
  set_genes <- function(ids, species = sp()) {
    updateSelectizeInput(session, "genes", choices = gene_choices(species, isolate(idtype())),
                         selected = head(unique(ids), 40), server = TRUE)
  }
  set_genes(DEFAULT_GENES, "maize")
  observeEvent(input$idtype, set_genes(isolate(genes())), ignoreInit = TRUE)

  # Switching species carries the gene list over to orthologs.
  observeEvent(input$species, {
    old <- isolate(input$genes); prev <- last_species(); last_species(sp())
    if (prev == sp()) return()
    o <- LG$orthologs
    new <- if (sp() == "rice") unique(o$rice[o$maize %in% old]) else unique(o$maize[o$rice %in% old])
    new <- head(new, 12)
    set_genes(new)
    output$gene_msg <- renderUI(if (length(old)) small_note(sprintf(
      "Switched to %d %s ortholog(s) of the previous genes.", length(new), sp())))
  }, ignoreInit = TRUE)

  observeEvent(input$add_pasted, {
    r <- resolve_ids(input$paste_ids)
    other <- r$ids[gene_species(r$ids) != sp()]
    ids <- setdiff(r$ids, other)
    set_genes(c(genes(), ids))
    updateTextAreaInput(session, "paste_ids", value = "")
    multi <- r$hits[r$hits$id %in% ids, ]
    multi <- names(which(table(multi$alias) > 1))
    msg <- c(if (length(ids)) sprintf("Added %d gene(s).", length(ids)),
             if (length(multi)) sprintf("%s match several v3 genes (merged/complex in v5); all were added.",
                                        paste(head(multi, 5), collapse = ", ")),
             if (length(other)) sprintf("%d belong to the other species (switch species to use them): %s",
                                        length(other), paste(head(other, 5), collapse = ", ")),
             if (length(r$missing)) sprintf("Not found: %s", paste(head(r$missing, 10), collapse = ", ")))
    output$gene_msg <- renderUI(small_note(HTML(paste(msg, collapse = "<br>"))))
  })
  observeEvent(input$clear_genes, set_genes(character()))
  add_genes <- function(ids) {
    set_genes(c(genes(), ids))
    showNotification(sprintf("Added %d gene(s) to the gene list.", length(ids)), duration = 3)
  }

  output$lib_note <- renderUI({
    if (ds() != "geo") return(small_note(if (ds() == "paper")
      "Section means published with the paper (no replicates)." else
      "UDM: both species interpolated onto 15 shared developmental positions."))
    s <- LG$samples[LG$samples$species == sp(), ]
    kept <- s$mapped_reads >= min_reads()
    per <- tapply(kept, s$section, sum)
    small_note(sprintf("%d of %d libraries kept.%s", sum(kept), nrow(s),
      if (any(per < 2)) paste0(" Sections with < 2 libraries: ",
                               paste(names(per)[per < 2], collapse = ", "), ".") else ""))
  })

  # Axis range for the range sliders
  observe({
    s <- dataset_samples(sp(), ds())
    xa <- if (ds() == "udm") "section" else xaxis()
    r <- range(s[[xa]])
    step <- if (xa %in% c("pos_cm", "pos_ligule")) 0.5 else 1
    for (id in c("fit_range", "scr_range"))
      updateSliderInput(session, id, min = r[1], max = r[2], value = r, step = step)
    k <- length(unique(s$section))
    updateSliderInput(session, "cx_at", max = k, value = min(isolate(input$cx_at), k))
  })

  # ---- Profiles -----------------------------------------------------------
  prof <- reactive({
    req(length(genes()) > 0)
    d <- profile_long(sp(), ds(), genes(), input$scale, min_reads(), xaxis())
    validate(need(!is.null(d), "None of the selected genes are in this dataset."))
    m <- summarise_sections(d)
    if (isTRUE(input$zscore)) { z <- zscore_long(d, m); d <- z$d; m <- z$means }
    list(d = d, means = m)
  })

  output$profile_plot <- renderPlotly({
    validate(need(length(genes()) > 0, "Pick one or more genes in the sidebar."))
    p <- prof()
    profile_plot(p$d, p$means, xlab(), if (isTRUE(input$zscore)) "z-score" else ylab(),
                 show_reps = input$show_reps, err = input$err, ligule_x = ligule_x(), namer = sn)
  })

  output$missing_note <- renderUI({
    miss <- setdiff(genes(), rownames(dataset_matrix(sp(), ds())))
    if (length(miss)) small_note(sprintf("Not in this dataset: %s", paste(sn(miss), collapse = ", ")))
  })

  output$gene_info <- renderDT({
    req(length(genes()) > 0)
    g <- gene_table(genes())
    id <- g$`ID (v3)`
    g$`ID (v3)` <- ifelse(g$Species == "maize",
      sprintf('<a href="https://www.maizegdb.org/gene_center/gene/%s" target="_blank">%s</a>',
              ifelse(g$`B73 v4` != "" & !grepl(",", g$`B73 v4`), g$`B73 v4`, id), id),
      sprintf('<a href="https://rice.uga.edu/cgi-bin/ORF_infopage.cgi?orf=%s" target="_blank">%s</a>', id, id))
    g$`B73 v5` <- vapply(strsplit(g$`B73 v5`, ", "), function(v) paste(sprintf(
      '<a href="https://www.maizegdb.org/gene_center/gene/%s" target="_blank">%s</a>', v, v), collapse = ", "), "")
    if (all(g$Species == "maize")) g <- g[, c(7, 1, 8, 3, 4, 5, 6, 9)] else g$`v3 -> v5` <- NULL
    datatable(g, escape = FALSE, rownames = FALSE, selection = "none",
              options = list(dom = "t", pageLength = 50, scrollX = TRUE))
  })

  output$means_table <- renderDT({
    m <- prof()$means
    w <- reshape(m[, c("gene", "x", "mean")], idvar = "gene", timevar = "x", direction = "wide")
    names(w) <- sub("^mean\\.", "x = ", names(w))
    w <- cbind(Symbol = GENES$symbol[match(w$gene, GENES$id)], add_v5(w, "gene"))
    dt(w, pageLength = 25)
  })

  output$dl_profiles <- download_csv(function() {
    d <- profile_long(sp(), ds(), genes(), input$scale, min_reads(), xaxis())
    names(d)[names(d) == "x"] <- xaxis(); names(d)[names(d) == "raw"] <- "expression"
    names(d)[names(d) == "value"] <- paste0("value_", input$scale)
    add_v5(d, "gene")
  }, "leaf_gradient_profiles")

  # ---- Fit models ---------------------------------------------------------
  observe({
    g <- genes()
    updateSelectInput(session, "fit_gene", choices = setNames(g, sn(g)),
                      selected = if (isolate(input$fit_gene) %in% g) isolate(input$fit_gene) else g[1])
    updateSelectInput(session, "cx_gene", choices = setNames(g, sn(g)),
                      selected = if (isolate(input$cx_gene) %in% g) isolate(input$cx_gene) else g[1])
  })

  fit_data <- reactive({
    req(input$fit_gene, input$fit_range)
    d <- profile_long(sp(), ds(), input$fit_gene, input$scale, min_reads(), xaxis())
    validate(need(!is.null(d), "Gene not in this dataset."))
    d <- d[d$x >= input$fit_range[1] - 1e-9 & d$x <= input$fit_range[2] + 1e-9 & is.finite(d$value), ]
    if (input$fit_on == "means" || ds() != "geo") {
      d <- do.call(rbind, lapply(split(d, d$section), function(z)
        data.frame(gene = z$gene[1], sample = paste("section", z$section[1]), section = z$section[1],
                   x = z$x[1], value = mean(z$value))))
    }
    d$y <- d$value
    validate(need(length(unique(d$x)) >= 4, "Need at least 4 positions in the fit range."))
    d
  })

  fits <- reactive({
    d <- fit_data()
    mods <- union(input$fit_models, input$fit_focus)
    f <- lapply(setNames(mods, mods), function(m) fit_one(d, m, input$spline_df))
    f[!vapply(f, is.null, TRUE)]
  })

  output$fit_msg <- renderUI({
    failed <- setdiff(union(input$fit_models, input$fit_focus), names(fits()))
    d <- fit_data()
    small_note(sprintf("Fitting %s on %d %s. Formula of focus model: %s", input$fit_gene, nrow(d),
                       if (any(duplicated(d$section))) "replicate libraries" else "section means",
                       model_formula_text(input$fit_focus, input$spline_df)),
               if (length(failed)) tags$span(class = "text-danger",
                 sprintf(" Did not converge / not estimable: %s.",
                         paste(names(MODELS)[MODELS %in% failed], collapse = ", "))))
  })

  output$fit_plot <- renderPlotly({
    f <- fits(); d <- fit_data()
    validate(need(length(f) > 0, "No model could be fitted."))
    curves <- do.call(rbind, lapply(f, predict_curve, d = d, level = input$ci_level))
    focus <- names(MODELS)[MODELS == input$fit_focus]
    fit_plot(d, curves, focus, xlab(), ylab(), ligule_x())
  })

  fit_stats <- reactive({
    f <- fits(); d <- fit_data()
    st <- do.call(rbind, lapply(f, model_stats, d = d))
    st$dAIC <- st$AIC - min(st$AIC)
    st[order(st$AIC), c("Model", "n", "Parameters", "R2", "Adj_R2", "Sigma", "AIC", "dAIC", "BIC",
                        "F_p", "LackOfFit_p")]
  })
  output$fit_compare <- renderDT(dt(fit_stats(), pageLength = 10))

  focus_fit <- reactive({
    f <- fits()[[input$fit_focus]]
    validate(need(!is.null(f), "Focus model could not be fitted."))
    f
  })
  output$fit_derived <- renderUI({
    lines <- derived_quantities(focus_fit(), fit_data(), input$scale, xunit())
    tags$ul(lapply(lines, tags$li))
  })
  output$fit_coef <- renderDT(dt(coef_table(focus_fit(), input$ci_level)))
  output$fit_resid <- renderPlotly({
    f <- focus_fit(); d <- fit_data()
    rd <- data.frame(fitted = fitted(f$fit), resid = resid(f$fit), x = d$x, sample = d$sample)
    plot_ly(rd, x = ~fitted, y = ~resid, type = "scatter", mode = "markers",
            text = ~sprintf("%s<br>x = %s", sample, x), hoverinfo = "text+y",
            marker = list(color = "#2E7D32", opacity = 0.7)) |>
      base_layout("Fitted", "Residual") |>
      layout(shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 0, y1 = 0,
                                line = list(dash = "dot", color = "#888"))))
  })
  output$fit_nested <- renderDT(dt(nested_tests(fit_data()), pageLength = 5))

  output$dl_fit_coef <- download_csv(function() {
    add_v5(cbind(gene = input$fit_gene, model = input$fit_focus, coef_table(focus_fit(), input$ci_level)), "gene")
  }, "fit_coefficients")
  output$dl_fit_data <- download_csv(function() {
    d <- fit_data(); out <- d[, c("gene", "sample", "section", "x", "y")]
    for (m in names(fits())) out[[paste0("fitted_", m)]] <- fitted(fits()[[m]]$fit)
    add_v5(out, "gene")
  }, "fit_data")

  ancova <- reactive({
    g <- genes()
    if (length(g) < 2) return(NULL)
    d <- profile_long(sp(), ds(), g, input$scale, min_reads(), xaxis())
    if (is.null(d)) return(NULL)
    d <- d[d$x >= input$fit_range[1] - 1e-9 & d$x <= input$fit_range[2] + 1e-9, ]
    if (input$fit_on == "means" || ds() != "geo")
      d <- aggregate(value ~ gene + section + x, data = d, FUN = mean)
    d$y <- d$value
    slope_comparison(d)
  })
  output$ancova_msg <- renderUI({
    a <- ancova()
    if (is.null(a)) return(small_note("Select two or more genes to compare slopes."))
    tags$p(sprintf("Test of equal slopes (gene x position interaction): F = %.3g, p = %.3g.", a$F, a$p))
  })
  output$ancova_table <- renderDT({ a <- ancova(); req(a); dt(add_v5(a$slopes, "Gene")) })

  # ---- Gradient screen ----------------------------------------------------
  screen <- eventReactive(input$scr_run, {
    s <- dataset_samples(sp(), ds(), min_reads())
    xa <- if (ds() == "udm") "section" else xaxis()
    keep <- s[[xa]] >= input$scr_range[1] - 1e-9 & s[[xa]] <= input$scr_range[2] + 1e-9
    s <- s[keep, ]
    validate(need(length(unique(s$section)) >= 4, "Need at least 4 sections in the range."))
    m <- dataset_matrix(sp(), ds())[, s$sample, drop = FALSE]
    secs <- sort(unique(s$section))
    A <- sapply(secs, function(k) (s$section == k) / sum(s$section == k))
    max_mean <- apply(m %*% A, 1, max)
    m <- m[is.finite(max_mean) & max_mean >= input$scr_minexpr, , drop = FALSE]
    validate(need(nrow(m) > 0, "No genes pass the expression filter."))
    withProgress(message = "Fitting all genes...", {
      res <- screen_genes(transform_values(m, input$scale), s, s[[xa]], input$scr_model, input$spline_df)
    })
    res$max_section_mean <- max_mean[res$id]
    res$symbol <- GENES$symbol[match(res$id, GENES$id)]
    res$description <- GENES$description[match(res$id, GENES$id)]
    lead <- c("id", "symbol", if (input$scr_model == "linear") c("slope", "slope_se"), "R2", "p", "q")
    res <- res[order(res$p), c(lead, setdiff(names(res), c(lead, "description")), "description")]
    res <- add_v5(res)
    attr(res, "settings") <- list(model = input$scr_model, species = sp(), dataset = ds(),
                                  scale = input$scale, range = input$scr_range, xaxis = xa)
    res
  })

  output$scr_summary <- renderUI({
    r <- screen(); st <- attr(r, "settings"); sig <- r$q <= input$scr_fdr
    parts <- sprintf("%s genes tested; %s with FDR <= %s", format(nrow(r), big.mark = ","),
                     format(sum(sig, na.rm = TRUE), big.mark = ","), input$scr_fdr)
    if (st$model == "linear") parts <- paste0(parts, sprintf(" (%s increasing, %s decreasing along the range)",
      format(sum(sig & r$slope > 0, na.rm = TRUE), big.mark = ","),
      format(sum(sig & r$slope < 0, na.rm = TRUE), big.mark = ",")))
    tags$p(class = "mb-0", strong(parts), small_note(sprintf(
      "%s, %s, %s scale, range %s-%s. slope = change per %s; peak_x = position of the highest section mean; section_p = one-way ANOVA across sections; lackfit_p = curve vs section means.",
      names(MODELS)[MODELS == st$model], names(DATASETS)[DATASETS == st$dataset], st$scale,
      st$range[1], st$range[2], xunit())))
  })

  scr_clicked <- reactiveVal(NULL)
  output$scr_volcano <- renderPlotly({
    r <- screen(); st <- attr(r, "settings")
    r$nlq <- -log10(pmax(r$q, 1e-300)); r$sig <- ifelse(r$q <= input$scr_fdr, "significant", "n.s.")
    xv <- if (st$model == "linear") ~slope else ~peak_x
    xt <- if (st$model == "linear") paste("Slope per", xunit()) else paste("Peak position (", xunit(), ")")
    plot_ly(r, x = xv, y = ~nlq, color = ~sig, colors = c("n.s." = "#BBBBBB", significant = "#2E7D32"),
            type = "scattergl", mode = "markers", marker = list(size = 4, opacity = 0.6),
            customdata = ~id, source = "volcano",
            text = ~paste0(sn(id), "<br>R2 = ", signif(R2, 3),
                           "<br>q = ", signif(q, 3)), hoverinfo = "text") |>
      base_layout(xt, "-log10(FDR q)") |> event_register("plotly_click")
  })
  observeEvent(input[["plotly_click-volcano"]], {   # input, so nothing is read before the plot exists
    scr_clicked(event_data("plotly_click", source = "volcano")$customdata)
  })

  scr_selected_ids <- reactive({
    r <- screen(); rows <- input$scr_table_rows_selected
    ids <- if (length(rows)) r$id[rows] else character()
    unique(c(scr_clicked(), ids))
  })

  output$scr_gene_plot <- renderPlotly({
    ids <- head(scr_selected_ids(), 12)
    validate(need(length(ids) > 0, "Click a point or select table rows."))
    st <- attr(screen(), "settings")
    d <- profile_long(st$species, st$dataset, ids, st$scale, min_reads(), st$xaxis)
    m <- summarise_sections(d)
    profile_plot(d, m, xlab(), ylab(), show_reps = TRUE, err = "se", ligule_x = ligule_x(), namer = sn)
  })

  output$scr_table <- renderDT(dt(screen(), selection = "multiple"))
  observeEvent(input$scr_add, {
    ids <- scr_selected_ids()
    if (length(ids)) add_genes(ids)
  })
  output$dl_screen <- download_csv(function() screen(), "gradient_screen")

  # ---- Base vs tip patterns -----------------------------------------------
  patterns <- eventReactive(list(input$pf_run, sp(), ds(), min_reads()), {
    validate(need(input$pf_fold >= 1, "Fold change must be at least 1."),
             need(input$pf_tipfold >= 1, "Tip tolerance must be at least 1."))
    withProgress(message = "Classifying genes...",
      find_patterns(sp(), ds(), min_reads(), input$pf_minexpr, input$pf_fold, input$pf_fdr,
                    input$pf_tipfold, input$pf_end / 100))
  })
  pf_list <- reactive({ r <- patterns(); r[r$pattern %in% input$pf_which, ] })

  output$pf_summary <- renderUI({
    r <- patterns(); st <- attr(r, "settings")
    n <- function(p, sh = NULL) {
      ok <- if (is.null(sh)) TRUE else r$shape %in% sh
      format(sum(r$pattern %in% p & ok), big.mark = ",")
    }
    tags$p(class = "mb-0",
      strong(sprintf("%s high at base & declining (%s gradual, %s early switch); %s low at base & rising (%s gradual, %s early switch).",
        n(PATTERNS[1]), n(PATTERNS[1], "gradual"), n(PATTERNS[1], "early switch, then plateau"),
        n(PATTERNS[2]), n(PATTERNS[2], "gradual"), n(PATTERNS[2], "early switch, then plateau"))),
      small_note(sprintf("%s: %s genes expressed (some section mean >= %g FPKM); %s; %d sections, base = first %d, tip = last %d.",
        tools::toTitleCase(st$species), format(nrow(r), big.mark = ","), st$min_fpkm,
        names(DATASETS)[DATASETS == st$dataset], st$n_sections, st$n_base_tip, st$n_base_tip)))
  })

  output$pf_overview <- renderPlot({
    r <- patterns()
    validate(need(any(!is.na(r$pattern)), "No genes match these thresholds."))
    pattern_figure(r, xlab = if (ds() == "udm") "UDM position" else "Distance from leaf base (cm)",
                   ligule_x = if (sp() == "maize" && ds() != "udm") 4)
  })

  pf_ids <- reactive({
    l <- pf_list(); rows <- input$pf_table_rows_selected
    if (length(rows)) pattern_ids(l)[rows] else head(pattern_ids(l), 6)
  })
  output$pf_gene_plot <- renderPlotly({
    ids <- head(pf_ids(), 12)
    validate(need(length(ids) > 0, "No genes in this list."))
    d <- profile_long(sp(), ds(), ids, "log2", min_reads(), xaxis())
    profile_plot(d, summarise_sections(d), xlab(), "log2(expression + 1)", show_reps = TRUE, err = "se",
                 ligule_x = ligule_x(), namer = sn)
  })

  output$pf_table <- renderDT({
    l <- pf_list()
    cols <- intersect(c("v5_id", "v3_id", "gene_id", "v3_to_v5", "symbol", "shape", "base_FPKM", "tip_FPKM",
                        "log2FC_tip_vs_base", "slope_log2_per_cm", "q_linear", "spearman_rho",
                        "half_change_cm", "description"), names(l))
    dt(l[, cols], selection = "multiple")
  })
  observeEvent(input$pf_add, {
    rows <- input$pf_table_rows_selected
    if (!length(rows)) return(showNotification("Select rows first.", type = "warning"))
    add_genes(pattern_ids(pf_list())[rows])
  })
  output$dl_pf <- downloadHandler(
    filename = function() sprintf("%s_%s_%s.csv", sp(), input$pf_which, format(Sys.Date())),
    content = function(file) write.csv(pf_list(), file, row.names = FALSE, na = ""))

  # ---- Co-expression ------------------------------------------------------
  observe({
    met <- LG$metabolites[[sp()]]$info
    updateSelectizeInput(session, "cx_metab", choices = met$trait, selected = met$trait[1])
    pick <- intersect(c("Chlorophyll a [metabolite platform]", "Starch [metabolite platform]",
                        "PEP-Carboxylase [enzyme activity platform]",
                        "Rubisco (maximal) [enzyme activity platform]"), met$trait)
    updateSelectizeInput(session, "met_traits", choices = met$trait,
                         selected = if (length(pick)) pick else met$trait[1:2])
    updateRadioButtons(session, "cx_target", selected = sp())
  })

  coexp <- eventReactive(input$cx_run, {
    target <- if (ds() == "udm") input$cx_target else sp()
    M <- section_mean_matrix(target, ds(), input$scale, min_reads())
    raw <- section_mean_matrix(target, ds(), "linear", min_reads())
    M <- M[apply(raw, 1, max) >= input$cx_minexpr, , drop = FALSE]
    secs <- colnames(M)
    q <- switch(input$cx_type,
      gene = {
        req(input$cx_gene)
        qm <- section_mean_matrix(sp(), ds(), input$scale, min_reads(), ids = input$cx_gene)
        validate(need(nrow(qm) == 1, "Query gene is not in this dataset."))
        qm[1, secs]
      },
      metab = {
        validate(need(ds() != "udm", "Metabolite traits are on the native section axis; choose GEO or paper data."))
        LG$metabolites[[sp()]]$mean[input$cx_metab, secs]
      },
      template = setNames(template_profile(input$cx_template, length(secs), input$cx_at), secs))
    validate(need(sum(is.finite(q)) >= 4, "Query has fewer than 4 finite values."))
    res <- correlate_profiles(M, q, input$cx_method)
    if (input$cx_type == "gene") res <- res[res$id != input$cx_gene, ]
    res <- res[is.finite(res$r), ]
    res <- switch(input$cx_dir, pos = res[order(-res$r), ], neg = res[order(res$r), ],
                  both = res[order(-abs(res$r)), ])
    res <- head(res, input$cx_n)
    res$symbol <- GENES$symbol[match(res$id, GENES$id)]
    res$description <- GENES$description[match(res$id, GENES$id)]
    qlab <- switch(input$cx_type, gene = sn(input$cx_gene), metab = input$cx_metab,
                   template = paste("template:", input$cx_template))
    list(res = add_v5(res[, c("id", "symbol", "r", "p", "q", "description")]), q = q, M = M, qlab = qlab,
         target = target)
  })

  output$cx_table <- renderDT(dt(coexp()$res, digits = 4, selection = "multiple"))
  output$cx_plot <- renderPlotly({
    cx <- coexp(); rows <- input$cx_table_rows_selected
    ids <- if (length(rows)) cx$res$id[rows] else head(cx$res$id, 10)
    z <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    secs <- as.numeric(names(cx$q))
    p <- plot_ly()
    cols <- gene_cols(length(ids))
    for (i in seq_along(ids))
      p <- add_lines(p, x = secs, y = z(cx$M[ids[i], ]), name = sn(ids[i]),
                     line = list(color = cols[i], width = 1.5))
    p <- add_trace(p, x = secs, y = z(cx$q), name = paste("query:", cx$qlab), type = "scatter",
                   mode = "lines+markers", line = list(color = "black", width = 4),
                   marker = list(color = "black"))
    base_layout(p, if (ds() == "udm") "UDM position" else "Section (base = 1)", "z-score")
  })
  observeEvent(input$cx_add, {
    cx <- coexp(); rows <- input$cx_table_rows_selected
    if (!length(rows)) return(showNotification("Select rows first.", type = "warning"))
    if (cx$target != sp()) return(showNotification("Those genes are from the other species.", type = "warning"))
    add_genes(cx$res$id[rows])
  })
  output$dl_cx <- download_csv(function() coexp()$res, "coexpression")

  # ---- Heatmap ------------------------------------------------------------
  hm <- eventReactive(input$hm_run, {
    ids <- switch(input$hm_source,
      list = genes(),
      c4 = LG$c4$id,
      keyword = {
        g <- GENES[GENES$species == sp(), ]
        hit <- grepl(input$hm_keyword, paste(g$symbol, g$description), ignore.case = TRUE)
        g$id[hit]
      },
      screen = {
        r <- tryCatch(screen(), error = function(e) NULL)
        validate(need(!is.null(r), "Run the gradient screen first."))
        r$id[r$q <= input$scr_fdr]
      })
    validate(need(length(ids) > 0, "No genes to show."))
    raw <- section_mean_matrix(sp(), ds(), "linear", min_reads(), ids = ids)
    M <- section_mean_matrix(sp(), ds(), input$scale, min_reads(), ids = ids)
    keep <- apply(raw, 1, max) >= input$hm_minexpr & apply(M, 1, sd) > 0
    M <- M[keep, , drop = FALSE]
    n_found <- nrow(M)
    if (nrow(M) > input$hm_max) M <- M[order(-apply(M, 1, sd))[seq_len(input$hm_max)], , drop = FALSE]
    validate(need(nrow(M) >= 2, "Fewer than 2 expressed genes in this set."))
    Z <- if (input$hm_z) t(scale(t(M))) else M
    hc <- hclust(as.dist(1 - cor(t(M))), method = "average")
    k <- min(input$hm_k, nrow(M))
    cl <- cutree(hc, k)
    list(Z = Z[hc$order, , drop = FALSE], cl = cl, order = hc$order, n_in = length(ids),
         n_found = n_found)
  })

  output$hm_msg <- renderUI({
    h <- hm()
    small_note(sprintf("%d genes requested, %d expressed above the cutoff, %d shown (most variable first if capped).",
                       h$n_in, h$n_found, nrow(h$Z)))
  })
  output$hm_plot <- renderPlotly({
    h <- hm(); Z <- h$Z
    rownames(Z) <- paste0("[", h$cl[rownames(Z)], "] ", sn(rownames(Z)))
    heatmap_plot(Z, if (ds() == "udm") "UDM position" else "Section (base = 1)",
                 row_labels = nrow(Z) <= 150, zlab = if (input$hm_z) "z-score" else ylab())
  })
  output$hm_clusters <- renderPlot({
    h <- hm(); Z <- h$Z
    long <- data.frame(gene = rep(rownames(Z), ncol(Z)), section = rep(as.numeric(colnames(Z)), each = nrow(Z)),
                       value = as.vector(Z))
    long$cluster <- factor(h$cl[long$gene])
    cm <- aggregate(value ~ cluster + section, long, mean)
    sizes <- table(h$cl)
    lab <- setNames(sprintf("Cluster %s (n = %d)", names(sizes), as.integer(sizes)), names(sizes))
    ggplot(long, aes(section, value, group = gene)) +
      geom_line(colour = "grey75", linewidth = 0.3) +
      geom_line(data = cm, aes(group = cluster), colour = "#2E7D32", linewidth = 1.3) +
      facet_wrap(~cluster, labeller = labeller(cluster = lab), ncol = 2) +
      labs(x = "Section (base = 1)", y = if (input$hm_z) "z-score" else ylab()) +
      theme_minimal(base_size = 13)
  })
  hm_table <- reactive({
    h <- hm(); ids <- rownames(h$Z)
    add_v5(data.frame(id = ids, cluster = h$cl[ids], symbol = GENES$symbol[match(ids, GENES$id)],
                      description = GENES$description[match(ids, GENES$id)]))
  })
  output$hm_table <- renderDT(dt(hm_table(), pageLength = 15))
  output$dl_hm <- download_csv(hm_table, "heatmap_clusters")

  # ---- Maize vs rice ------------------------------------------------------
  orth_set <- reactive({
    g <- head(genes(), 6)
    req(length(g) > 0)
    o <- LG$orthologs
    if (sp() == "maize") list(maize = g, rice = unique(o$rice[o$maize %in% g])) else
      list(maize = unique(o$maize[o$rice %in% g]), rice = g)
  })

  species_lines <- function(p, M, sp_name, dash, z, cols) {
    for (i in seq_len(nrow(M))) {
      y <- log2(M[i, ] + 1); if (z) y <- (y - mean(y)) / sd(y)
      p <- add_trace(p, x = as.numeric(colnames(M)), y = y, type = "scatter", mode = "lines+markers",
                     name = paste0(sn(rownames(M)[i]), " [", sp_name, "]"),
                     line = list(color = cols[i], dash = dash, width = 2), marker = list(color = cols[i], size = 5))
    }
    p
  }

  output$or_udm_plot <- renderPlotly({
    os <- orth_set()
    mz <- LG$udm$maize[intersect(os$maize, rownames(LG$udm$maize)), , drop = FALSE]
    rc <- LG$udm$rice[intersect(os$rice, rownames(LG$udm$rice)), , drop = FALSE]
    validate(need(nrow(mz) + nrow(rc) > 0, "No UDM values for these genes."))
    p <- plot_ly()
    p <- species_lines(p, mz, "maize", "solid", input$or_z, gene_cols(nrow(mz)))
    p <- species_lines(p, rc, "rice", "dash", input$or_z, rev(gene_cols(12))[seq_len(nrow(rc))])
    base_layout(p, "UDM unified gradient position", if (input$or_z) "z-score of log2" else "log2(RPKM + 1)")
  })

  output$or_native_plot <- renderPlotly({
    os <- orth_set()
    one <- function(spn, ids) {
      ids <- intersect(ids, rownames(LG$expr[[spn]]))
      if (!length(ids)) return(NULL)
      M <- section_mean_matrix(spn, "geo", "linear", input$min_reads, ids = ids)
      x <- LG$sections$pos_cm[LG$sections$species == spn][order(LG$sections$section[LG$sections$species == spn])]
      colnames(M) <- x[as.numeric(colnames(M))]
      M
    }
    mz <- one("maize", os$maize); rc <- one("rice", os$rice)
    p <- plot_ly()
    if (!is.null(mz)) p <- species_lines(p, mz, "maize", "solid", input$or_z, gene_cols(nrow(mz)))
    if (!is.null(rc)) p <- species_lines(p, rc, "rice", "dash", input$or_z, rev(gene_cols(12))[seq_len(nrow(rc))])
    base_layout(p, "Distance from leaf base (cm)", if (input$or_z) "z-score of log2" else "log2(FPKM + 1)")
  })

  pairs_v5 <- reactive(add_v5(ortholog_pairs(), "maize", "maize_"))
  output$pairs_table <- renderDT(dt(pairs_v5(), selection = "single", pageLength = 10))
  output$dl_pairs <- download_csv(pairs_v5, "ortholog_pairs_UDM")
  output$pair_plot <- renderPlotly({
    row <- input$pairs_table_rows_selected
    validate(need(length(row) == 1, "Select a pair in the table."))
    pr <- ortholog_pairs()[row, ]
    p <- plot_ly()
    p <- species_lines(p, LG$udm$maize[pr$maize, , drop = FALSE], "maize", "solid", input$or_z, PAL[1])
    p <- species_lines(p, LG$udm$rice[pr$rice, , drop = FALSE], "rice", "dash", input$or_z, PAL[2])
    base_layout(p, "UDM position", if (input$or_z) "z-score" else "log2(RPKM + 1)")
  })
  output$pairs_hist <- renderPlotly({
    plot_ly(ortholog_pairs(), x = ~r_UDM, type = "histogram", nbinsx = 60,
            marker = list(color = "#2E7D32")) |>
      base_layout("Pearson r of UDM profiles", "Ortholog pairs")
  })

  # ---- Metabolites --------------------------------------------------------
  output$met_plot <- renderPlotly({
    tr <- input$met_traits
    validate(need(length(tr) > 0, "Pick one or more traits."))
    met <- LG$metabolites[[sp()]]
    sec <- LG$sections[LG$sections$species == sp(), ]
    sec <- sec[order(sec$section), ]
    xa <- if (ds() == "udm") "section" else xaxis()
    p <- plot_ly(); cols <- gene_cols(length(tr))
    for (i in seq_along(tr)) {
      y <- met$mean[tr[i], ]; e <- met$sd[tr[i], ]
      if (input$met_z) { e <- e / sd(y, na.rm = TRUE); y <- (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE) }
      p <- add_trace(p, x = sec[[xa]], y = y, type = "scatter", mode = "lines+markers", name = tr[i],
                     error_y = list(array = e, color = cols[i]),
                     line = list(color = cols[i]), marker = list(color = cols[i]))
    }
    p <- base_layout(p, if (ds() == "udm") "Section (base = 1)" else xlab(),
                     if (input$met_z) "z-score" else "Trait value (units in table)")
    if (!is.null(ligule_x()) && ds() != "udm") p <- layout(p, shapes = ligule_shape(ligule_x()))
    p
  })
  output$met_table <- renderDT({
    met <- LG$metabolites[[sp()]]
    v <- as.data.frame(met$mean); names(v) <- paste("section", names(v))
    dt(cbind(met$info[, c("name", "unit", "platform")], v), pageLength = 10)
  })
  observeEvent(input$met_to_cx, {
    req(input$met_traits)
    if (ds() == "udm") updateSelectInput(session, "dataset", selected = "geo")
    updateRadioButtons(session, "cx_type", selected = "metab")
    updateSelectizeInput(session, "cx_metab", selected = input$met_traits[1])
    nav_select("nav", "Co-expression")
    showNotification("Query set - press Search.", duration = 4)
  })

  # ---- QC -----------------------------------------------------------------
  qc <- reactive(library_qc(sp()))
  output$qc_pca <- renderPlotly({
    pc <- library_pca(sp(), 0)
    sc <- pc$scores
    sc$kept <- ifelse(sc$mapped_reads >= input$min_reads, "kept", "excluded (low depth)")
    plot_ly(sc, x = ~PC1, y = ~PC2, color = ~section, symbol = ~kept, symbols = c("x", "circle"),
            type = "scatter", mode = "markers", marker = list(size = 10),
            text = ~sprintf("%s<br>%s mapped reads", sample, format(mapped_reads, big.mark = ",")),
            hoverinfo = "text") |>
      colorbar(title = "Section") |>
      base_layout(sprintf("PC1 (%.1f%%)", 100 * pc$var[1]), sprintf("PC2 (%.1f%%)", 100 * pc$var[2]))
  })
  output$qc_depth <- renderPlotly({
    s <- qc()
    s$kept <- s$mapped_reads >= input$min_reads
    plot_ly(s, x = ~sample, y = ~mapped_reads, type = "bar",
            marker = list(color = ifelse(s$kept, "#2E7D32", "#D55E00")),
            text = ~sprintf("%s: %s", sample, format(mapped_reads, big.mark = ",")), hoverinfo = "text") |>
      base_layout("", "Mapped reads") |>
      layout(xaxis = list(categoryorder = "array", categoryarray = s$sample, tickfont = list(size = 8)),
             yaxis = list(type = "log"),
             shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper",
                                y0 = max(input$min_reads, 1), y1 = max(input$min_reads, 1),
                                line = list(dash = "dot", color = "#888"))))
  })
  output$qc_table <- renderDT({
    s <- qc()
    dt(s[, c("sample", "rep", "section", "sec_label", "pos_cm", "total_reads", "mapped_reads",
             "r_vs_other_reps")], pageLength = 15)
  })
}

shinyApp(ui, server)
