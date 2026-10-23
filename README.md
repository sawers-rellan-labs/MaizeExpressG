# MaizeExpressG

An R Shiny app for exploring gene expression along the maize and rice leaf
developmental gradients from

> Wang L, Czedik-Eysenberg A, Mertz RA, *et al.* (2014) Comparative analyses of C4 and C3
> photosynthesis in developing leaves of maize and rice. *Nature Biotechnology* 32:1158–1165.
> [doi:10.1038/nbt.3019](https://doi.org/10.1038/nbt.3019) · GEO [GSE54274](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54274)

## Citation

All expression, ortholog and metabolite data in this app were downloaded from the supplementary
tables and GEO series of this paper. Please cite it when you use them:

> Wang L, Czedik-Eysenberg A, Mertz RA, Si Y, Tohge T, Nunes-Nesi A, Arrivault S, Dedow LK, Bryant DW,
> Zhou W, Xu J, Weissmann S, Studer A, Li P, Zhang C, LaRue T, Shao Y, Ding Z, Sun Q, Patel RV, Turgeon R,
> Zhu X, Provart NJ, Mockler TC, Fernie AR, Stitt M, Liu P, Brutnell TP (2014) Comparative analyses of C4
> and C3 photosynthesis in developing leaves of maize and rice. *Nature Biotechnology* 32(11):1158–1165.
> [doi:10.1038/nbt.3019](https://doi.org/10.1038/nbt.3019). PMID 25306245.
> Data: [GEO GSE54274](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54274); reads: SRA SRP018823.

<details><summary>BibTeX</summary>

```bibtex
@article{Wang2014,
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
}
```

</details>

Other resources used: gene annotation, v3→v4/v5 cross-references and gene models from
[MaizeGDB](https://www.maizegdb.org); assembly chain files from [Ensembl Plants](https://plants.ensembl.org);
rice annotation from the [Rice Genome Annotation Project (MSU7)](https://rice.uga.edu); ProtNLM protein
names from [UniProt](https://www.uniprot.org).

The app shows a short citation at the bottom of every page and the full one under **About**.

## Run it

```bash
Rscript run_app.R
```

Or, from R/RStudio in this folder: `shiny::runApp("app")`.
You need the R packages shiny, bslib, plotly, DT and ggplot2 (`run_app.R` installs any that are missing).
The processed data (`app/data/leafgrad.rds`, 21 MB) is already built. To rebuild it from the original sources:

```bash
Rscript scripts/01_download_data.R
Rscript scripts/04_liftover_v3_to_v5.R  # position-based v3 -> v5 matches (optional, ~1 min)
Rscript scripts/02_build_app_data.R
Rscript scripts/03_export_v5_tables.R   # maize tables with B73 v5 IDs -> exports/
```

## B73 v5 gene IDs

The expression data were quantified on the B73 **v3** gene models (GRMZM…). The app converts them to
**v5** (Zm00001eb…) and shows v5 IDs by default (sidebar switch: *Maize gene IDs*). Every results
table and download gets `v5_id` and `v3_to_v5` columns, and you can search or paste v3, v4 or v5 IDs.

Two sources, in order:

1. **MaizeGDB's v3→v5 cross-reference** (`data-raw/B73v3_to_B73v5.tsv`), which covers the v3 filtered gene set.
2. **Genome position, for genes the cross-reference misses** (`scripts/04_liftover_v3_to_v5.R`). The exons of
   each v3 gene model are lifted AGPv3 → AGPv4 → v5 with the Ensembl Plants assembly chain files and
   intersected (same strand) with the exons of the v5 gene models (protein-coding + non-coding). The v5 gene
   with the most overlapping exon bases is kept if it covers ≥ 50% of the v3 exons and ≥ 20% of its own
   exons (**strong**), or ≥ 20% of the v3 exons (**partial**). Checked against the cross-reference's 1:1
   genes, the liftover picks the same v5 gene for 99.6% of strong matches (n = 22,827) and 97.2% of partial
   ones (n = 1,612).

| v3 → v5 mapping | GEO genes | Supp. Table 13 genes | Meaning |
|---|---|---|---|
| 1:1 | 22,621 | 17,001 | one v3 gene ↔ one v5 gene |
| split | 1,830 | 1,348 | one v3 gene → several v5 genes (all listed; shown as `Zm… +n`) |
| merged | 1,960 | 1,287 | several v3 genes → one v5 gene |
| complex | 1,269 | 891 | several ↔ several (mostly tandem gene families) |
| liftover (strong) | 2,769 | 1,711 | matched by position, strong exon overlap |
| liftover (partial) | 1,919 | 1,125 | matched by position, partial exon overlap |
| no v5 match | 16,560 | 4,890 | mostly v3 low-confidence models where v5 has no gene model (15,369 of the GEO ones lift over fine but hit no v5 gene) |

When a v5 gene is claimed by more than one v3 gene (merged genes, and 1,982 of the 5,165 liftover
matches, often fragments next to a gene the cross-reference already maps), the app shows it as
`Zm… [GRMZM…]` so the v3 genes stay distinguishable.

Expression is kept per v3 gene: values are never summed or split across v5 models. `exports/` holds:

- `maize_v3_to_v5_mapping.csv`: every maize gene in the data, with v5, v4, symbol and mapping type, plus
  the liftover result for every gene (`liftover_v5`, `liftover_class`, exon overlap fractions) to cross-check
- `maize_GEO_FPKM_by_library_v5.csv`: 75 libraries; column names give library, section and position
- `maize_SuppTable13_RPKM_v5.csv`, `maize_UDM_v5.csv`: the paper's tables
- `CITATION.txt`: where each file comes from and how to cite it

For a table keyed by v5 ID, keep the rows with `v3_to_v5 == "1:1"` (add `liftover (strong)` if you accept
position-based matches).

## Data in the app

| Dataset (sidebar) | What it is | Genes |
|---|---|---|
| GEO replicates (default) | Per-library FPKM from GEO (TopHat/cuffdiff; maize AGPv3, rice MSU7). Maize: fifteen 1-cm segments × 5 replicates; rice: eleven 2-cm segments × 4 replicates | 48,928 maize / 33,872 rice |
| Paper Supp. Table 13 | RPKM section means as published (M1–M15, R1–R11) | 28,253 / 19,520 |
| Paper UDM | The "Unified" sheet: both species on 15 shared developmental positions | 27,365 / 21,228 |

Also included: maize–rice orthologs (Supp. Table 3, 84k pairs), metabolites and enzyme
activities (Supp. Tables 1–2), C4 shuttle genes (Supp. Table 7), reads per library
(Supp. Table 17), and gene annotation (MaizeGDB AGPv3 symbols/names, v3→v4/v5 IDs, ProtNLM
names, MSU7 rice annotation).

## What you can do

- **Profiles**: plot any genes along the gradient (replicates, mean ± SE/SD, log2 or linear,
  z-scored to compare shapes). Search by ID, symbol or description, or paste v3/v4/v5 maize IDs,
  symbols or LOC_Os loci. Switching species turns the gene list into its orthologs.
- **Fit models**: linear, quadratic, cubic, natural spline, segmented (breakpoint with a
  profile CI) and 4-parameter logistic (half-way transition position). Restrict the fit to any
  part of the leaf. Compare models (R², AIC/BIC, F tests, lack-of-fit against section means),
  get plain-language slopes (fold change per cm), residuals, nested tests, and an ANCOVA asking
  whether several genes share one slope.
- **Gradient screen**: fit a model to every expressed gene (≈0.2 s), with BH FDR, slope, peak
  position, one-way ANOVA and lack-of-fit p-values. Volcano plot (click to plot a gene), CSV export.
- **Base vs tip**: list the genes that are high at the base and decline to the tip, and those that are
  low at the base and rise to it, with adjustable fold change, FDR, expression cutoff, tip tolerance and
  base/tip window. Overview of every gene in each pattern, profiles of selected genes, add to the gene
  list, CSV download. Works for maize and rice and for all three datasets.
- **Co-expression**: rank genes by correlation with a query gene, a metabolite/enzyme trait, or a
  template shape (rising, falling, peak, step). On the UDM axis you can query maize and search rice.
- **Heatmap & clusters**: gene list, annotation regex (e.g. `photosystem|rubisco`), C4 genes, or
  screen hits → clustered heatmap, k cluster mean profiles, membership table.
- **Maize ↔ Rice**: genes and their orthologs on the unified axis and on each species' own axis;
  table of all ortholog pairs with UDM profile correlations (find conserved or diverged pairs).
- **Metabolites**, **Libraries & QC** (PCA, read depth, replicate agreement), **About**.

## Opposite base-to-tip patterns

The app's **Base vs tip** tab and `scripts/05_gradient_patterns.R` (run `Rscript scripts/05_gradient_patterns.R maize`,
or `rice`) share one function (`find_patterns()` in `app/R/25_patterns.R`). They find genes that are high at
the leaf base and decline towards the tip, and genes that are low at the base and rise towards the tip.
Defaults (all adjustable in the tab), on the GEO replicate data (libraries ≥ 250,000 mapped reads,
log2(FPKM + 1)):

- expressed: some section mean ≥ 5 FPKM
- ≥ 4-fold difference between the base (first 20% of sections) and the tip (last 20%)
- declining: highest section in the basal third, tip within 2-fold of the lowest section;
  rising: lowest section in the basal third, tip within 2-fold of the highest section
  (genes that rise and then plateau count; genes that peak mid-leaf and fall again do not)
- linear trend along the leaf with FDR < 0.01 (same model as the app's *Gradient screen*)

Results go to `results/`: one CSV per pattern (v5/v3 IDs, annotation, base and tip FPKM, fold change,
slope, FDR, `shape` = gradual vs early switch then plateau, `half_change_cm` = where half of the change is
reached, and the section means) plus `<species>_gradient_patterns.png`. Thresholds are at the top of the script.
With the defaults: maize 3,470 declining / 1,288 rising; rice 2,452 declining / 1,199 rising.

## Things to know

- **Section positions.** GEO labels the maize sections 1–4 and 6–16 (there is no "sec5"). The app
  numbers them 1–15 from the base; this matches columns M1–M15 of Supp. Table 13. Distance from the
  base is the segment midpoint: (section − 0.5) cm for maize, (2 × section − 1) cm for rice. Supp.
  Table 1 numbers maize sections −4…−1, 1…11, so the leaf-2 ligule is drawn at 4 cm (approximate).
- **Low-depth libraries.** Several maize libraries have very few mapped reads (e.g. `maize_rep2_sec7`:
  7,274 reads) and correlate poorly with their replicates. Libraries below 250,000 mapped reads are
  excluded by default (slider in the sidebar). At that cutoff maize section 6 (GEO `sec7`) keeps a
  single library.
- **GEO vs paper values.** The GEO FPKM (AGPv3, cuffdiff) and the paper's RPKM (Supp. Table 13) come
  from different pipelines and differ in scale for some genes. Use GEO when you want replicates.
- **Statistics.** Models are ordinary least squares (or nls) on log2(FPKM + 1) or FPKM. The screen is
  meant for exploring and ranking genes; for formal differential expression use a count-based
  method (edgeR, DESeq2, limma-voom) on the SRA reads (SRP018823).

## Layout

```
app/app.R              UI and server
app/R/00_data.R        data access, gene lookup
app/R/10_models.R      single-gene model fitting
app/R/20_genomewide.R  vectorised screen, co-expression, orthologs, QC
app/R/25_patterns.R    base-vs-tip pattern finder (app tab + scripts/05)
app/R/30_plots.R       plotly figures
app/R/40_about.R       About tab
app/data/leafgrad.rds  processed data
scripts/               download, build, liftover and v5-export scripts
data-derived/          liftover results (v3_to_v5_liftover.tsv)
exports/               maize expression tables with v5 IDs
data-raw/              original downloads (GEO, supplementary tables, annotation)
```
