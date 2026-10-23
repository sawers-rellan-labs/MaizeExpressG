# Source publication of every expression, ortholog and metabolite table in the app.
PAPER_AUTHORS <- paste(
  "Wang L, Czedik-Eysenberg A, Mertz RA, Si Y, Tohge T, Nunes-Nesi A, Arrivault S, Dedow LK,",
  "Bryant DW, Zhou W, Xu J, Weissmann S, Studer A, Li P, Zhang C, LaRue T, Shao Y, Ding Z, Sun Q,",
  "Patel RV, Turgeon R, Zhu X, Provart NJ, Mockler TC, Fernie AR, Stitt M, Liu P, Brutnell TP")
PAPER_TITLE <- "Comparative analyses of C4 and C3 photosynthesis in developing leaves of maize and rice"
PAPER_BIBTEX <- "@article{Wang2014,
  author  = {Wang, Lin and Czedik-Eysenberg, Angelika and Mertz, Rachel A. and Si, Yaqing and
             Tohge, Takayuki and Nunes-Nesi, Adriano and Arrivault, Stephanie and Dedow, Lauren K. and
             Bryant, Douglas W. and Zhou, Wen and Xu, Jiajia and Weissmann, Sarit and Studer, Anthony and
             Li, Pinghua and Zhang, Cankui and LaRue, Therese and Shao, Ying and Ding, Zehong and Sun, Qi and
             Patel, Rohan V. and Turgeon, Robert and Zhu, Xinguang and Provart, Nicholas J. and
             Mockler, Todd C. and Fernie, Alisdair R. and Stitt, Mark and Liu, Peng and Brutnell, Thomas P.},
  title   = {Comparative analyses of {C4} and {C3} photosynthesis in developing leaves of maize and rice},
  journal = {Nature Biotechnology},
  volume  = {32},
  number  = {11},
  pages   = {1158--1165},
  year    = {2014},
  doi     = {10.1038/nbt.3019}
}"

# Citation pieces are pasted into one HTML string so no stray spaces appear before punctuation.
ext_link <- function(text, href) as.character(a(text, href = href, target = "_blank"))
doi_link <- function() ext_link("doi:10.1038/nbt.3019", "https://doi.org/10.1038/nbt.3019")
geo_link <- function() ext_link("GEO GSE54274", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54274")

# Shown at the bottom of every page.
citation_footer <- function() div(
  class = "text-muted small border-top mt-3 pt-2 pb-3",
  HTML(paste0("Data: Wang L, Czedik-Eysenberg A, Mertz RA, et al. (2014) ", PAPER_TITLE,
              ". <em>Nature Biotechnology</em> 32:1158-1165. ", doi_link(), " | ", geo_link(),
              ". Please cite this paper when you use these data (full citation under About).")))

about_page <- function() {
  s <- LG$samples
  tagList(
    h4("How to cite"),
    p("All expression, ortholog and metabolite data in MaizeExpressG were downloaded from this paper ",
      "(supplementary tables) and its GEO series. Please cite it when you use them:"),
    tags$blockquote(class = "border-start border-3 ps-3", HTML(paste0(
      PAPER_AUTHORS, " (2014) ", PAPER_TITLE, ". <em>Nature Biotechnology</em> 32(11):1158-1165. ",
      doi_link(), ". PMID ", ext_link("25306245", "https://pubmed.ncbi.nlm.nih.gov/25306245/"), ". ",
      "Data: ", geo_link(), "; reads: ",
      ext_link("SRA SRP018823", "https://www.ncbi.nlm.nih.gov/sra?term=SRP018823"), "."))),
    tags$details(tags$summary("BibTeX"), tags$pre(PAPER_BIBTEX)),
    p(HTML(paste0("Other resources used: gene annotation, v3-to-v4/v5 cross-references and gene models from ",
      ext_link("MaizeGDB", "https://www.maizegdb.org"), "; assembly chain files from ",
      ext_link("Ensembl Plants", "https://plants.ensembl.org"), "; rice annotation from the ",
      ext_link("Rice Genome Annotation Project (MSU7)", "https://rice.uga.edu"), "; ProtNLM protein names from ",
      ext_link("UniProt", "https://www.uniprot.org"), "."))),
    h4("Data"),
    tags$ul(
      tags$li(strong("GEO replicates (default). "),
        a("GSE54274", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54274", target = "_blank"),
        sprintf(" processed FPKM tables (TopHat 2.0.10 + cuffdiff 2.1.1; maize B73 AGPv3, rice MSU7). Maize: third leaf of 9-day-old B73 seedlings, fifteen 1-cm segments, %d libraries (replicates 1, 2, 3, 4, 6). Rice: third leaf, eleven 2-cm segments, %d libraries (4 replicates).",
                sum(s$species == "maize"), sum(s$species == "rice"))),
      tags$li(strong("Paper Supp. Table 13. "), "RPKM section means published with the paper (M1-M15, R1-R11). ",
              "These come from the authors' own pipeline and differ in scale from the GEO FPKM values."),
      tags$li(strong("Paper UDM. "), "The 'Unified' sheet of Supp. Table 13: maize and rice expression ",
              "interpolated onto 15 shared developmental positions by the unified developmental model."),
      tags$li(strong("Orthologs: "), "Supp. Table 3. ", strong("Metabolites and enzyme activities: "),
              "Supp. Tables 1-2. ", strong("Library read counts: "), "Supp. Table 17."),
      tags$li(strong("Annotation: "), "MaizeGDB AGPv3 gene data (symbols, names, products), MaizeGDB v3-to-v4 ",
              "cross-reference (ambiguous mappings dropped), v3-to-v5 (all matches, see below), UniProt ProtNLM names for v5 models, ",
              "and MSU7 rice locus annotation. Maize genes without a description take their rice ortholog's ",
              "(marked '[rice ortholog]').")),
    h4("B73 v5 gene IDs"),
    p("Expression was quantified on B73 v3 gene models. v5 IDs come from MaizeGDB's v3-to-v5 cross-reference ",
      "(", a("B73_gene_xref", href = "https://download.maizegdb.org/Pan-genes/B73_gene_xref/", target = "_blank"),
      "); the 'v3 -> v5' column says whether that match is 1:1, split (one v3 to several v5, shown as 'Zm... +n'), ",
      "merged (several v3 to one v5) or complex. Genes the cross-reference misses were matched by genome position: ",
      "their exons were lifted AGPv3 -> AGPv4 -> v5 with the Ensembl Plants chain files and intersected with v5 exons ",
      "('liftover (strong)': >= 50% of the v3 exons and >= 20% of the v5 exons overlap; 'liftover (partial)': >= 20% ",
      "of the v3 exons). On the cross-reference's 1:1 genes the liftover agrees for 99.6% of strong and 97.2% of ",
      "partial matches. A v5 gene claimed by several v3 genes is shown with the v3 ID in brackets. Each v3 gene keeps ",
      "its own values. Run scripts/03_export_v5_tables.R to write the maize tables with v5 IDs to exports/."),
    h4("Positions along the leaf"),
    tags$ul(
      tags$li("GEO and Supp. Table 17 label the maize sections 1-4 and 6-16 (no 'sec5'). The app numbers them ",
              "1-15 from the base, which matches columns M1-M15 of Supp. Table 13 (checked by correlation). ",
              "Distance from base = (section - 0.5) cm for maize and (2 x section - 1) cm for rice (segment midpoints). ",
              "Use 'GEO section label' as the axis if you prefer the original labels."),
      tags$li("Supp. Table 1 numbers the maize sections -4..-1 and 1..11, i.e. relative to the ligule of leaf 2; ",
              "the app therefore places that ligule at 4 cm from the base (dotted line). Treat it as approximate.")),
    h4("Statistics"),
    tags$ul(
      tags$li("log2(x + 1) is recommended for model fitting: slopes become log2 fold-change per unit of the axis."),
      tags$li("Libraries below the mapped-read cutoff (sidebar; default 250,000) are excluded everywhere. ",
              "Several maize libraries have very few reads (section 6, i.e. GEO 'sec7', keeps only one replicate at the default)."),
      tags$li("Lack-of-fit tests compare a fitted curve with the section-means (one-way ANOVA) model using ",
              "between-replicate variation; they need replicate libraries."),
      tags$li("The segmented model's breakpoint is found by grid search; its confidence band is conditional on the ",
              "breakpoint, and the breakpoint interval is an F-based profile interval."),
      tags$li("Base vs tip: a gene is 'declining' when its highest section is in the basal third, the tip is within ",
              "the set fold of its lowest section, the base-vs-tip change reaches the fold cutoff and the linear ",
              "trend passes the FDR; 'rising' is the mirror image. Genes that rise and then plateau count, genes that ",
              "peak mid-leaf do not."),
      tags$li("Genome-wide screens use ordinary least squares on transformed FPKM for every gene with ",
              "Benjamini-Hochberg FDR. This is fast and transparent but is not a count-based model (edgeR/DESeq2/limma-voom); ",
              "use it to explore and rank genes.")),
    h4("Source code"),
    p(HTML(paste0(ext_link("github.com/sawers-rellan-labs/MaizeExpressG",
                           "https://github.com/sawers-rellan-labs/MaizeExpressG"), " (MIT license for the code; ",
                  "data remain under their original sources, cited above)."))),
    h4("Rebuilding the data"),
    tags$pre("Rscript scripts/01_download_data.R\nRscript scripts/02_build_app_data.R"),
    p(class = "text-muted small", "Data file built ", LG$info$built, ".")
  )
}
