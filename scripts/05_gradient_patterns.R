# Genes with opposite base-to-tip patterns along the leaf gradient:
#   base_high_declining: high at the leaf base, falling towards the tip
#   base_low_rising:     low at the leaf base, climbing towards the tip
# Same rules as the app's "Base vs tip" tab (find_patterns() in app/R/25_patterns.R).
#
# Data: GEO GSE54274 replicate FPKM (Wang et al. 2014, Nat Biotechnol 32:1158),
# log2(FPKM + 1), libraries below MIN_READS mapped reads excluded (as in the app).
# Run from the project root:  Rscript scripts/05_gradient_patterns.R [maize|rice]

suppressPackageStartupMessages(library(ggplot2))
args      <- commandArgs(TRUE)
SP        <- if (length(args)) args[1] else "maize"   # (SPECIES is taken by the app code)
MIN_READS <- 250000   # library depth cutoff
MIN_FPKM  <- 5        # gene must reach this section mean somewhere
MIN_FOLD  <- 4        # base vs tip
FDR       <- 0.01     # linear trend along the leaf
TIP_FOLD  <- 2        # tip within this fold of the opposite extreme
END_FRAC  <- 0.2      # base / tip = first / last 20% of sections

owd <- setwd("app"); for (f in list.files("R", full.names = TRUE)) source(f); setwd(owd)

res <- find_patterns(SP, "geo", MIN_READS, MIN_FPKM, MIN_FOLD, FDR, TIP_FOLD, END_FRAC)

dir.create("results", showWarnings = FALSE)
for (p in PATTERNS)
  write.csv(res[res$pattern %in% p, ], sprintf("results/%s_%s.csv", SP, p), row.names = FALSE, na = "")

fig <- pattern_figure(res, ligule_x = if (SP == "maize") 4,
  title = sprintf("%s leaf gradient: genes with opposite base-to-tip patterns", tools::toTitleCase(SP)),
  subtitle = sprintf(paste0(">= %g-fold base vs tip, extreme in basal third, tip within %g-fold of the other ",
                            "extreme, linear FDR < %g\nline = median, band = 10-90%%%s"),
                     MIN_FOLD, TIP_FOLD, FDR, if (SP == "maize") ", dotted = approx. leaf-2 ligule" else ""),
  caption = "Data: Wang et al. (2014) Nat Biotechnol 32:1158-1165, GEO GSE54274")
ggsave(sprintf("results/%s_gradient_patterns.png", SP), fig, width = 10, height = 4.8, dpi = 200, bg = "white")

cat(sprintf("%s: %d expressed genes (max section mean >= %g FPKM)\n", SP, nrow(res), MIN_FPKM))
print(table(pattern = res$pattern, shape = res$shape))
