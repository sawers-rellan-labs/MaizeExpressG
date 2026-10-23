# Data access helpers. Everything the app shows comes from data/leafgrad.rds,
# built by scripts/02_build_app_data.R.

LG <- readRDS("data/leafgrad.rds")
LG$samples <- as.data.frame(LG$samples)
LG$sections <- as.data.frame(LG$sections)
GENES <- LG$genes

SPECIES  <- c(Maize = "maize", Rice = "rice")
DATASETS <- c("GEO replicates - FPKM (5 maize / 4 rice reps)" = "geo",
              "Paper Supp. Table 13 - RPKM section means"     = "paper",
              "Paper UDM - unified maize/rice gradient"      = "udm")
XAXES    <- c("Distance from leaf base (cm)" = "pos_cm",
              "Section number (base = 1)"     = "section",
              "GEO section label"             = "sec_label",
              "Relative to leaf-2 ligule (cm, maize)" = "pos_ligule")
SCALES   <- c("log2(FPKM + 1)" = "log2", "FPKM" = "linear")
DEFAULT_MIN_READS <- 250000
PAL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
         "#882255", "#117733", "#332288", "#999933", "#AA4499", "#44AA99")

gene_cols <- function(n) rep_len(PAL, n)

x_label <- function(xaxis, dataset) {
  if (dataset == "udm") return("UDM unified gradient position (1-15)")
  names(XAXES)[XAXES == xaxis]
}

# Libraries (GEO) or section columns (paper/UDM) used for a species/dataset.
dataset_samples <- function(sp, dataset, min_reads = 0) {
  if (dataset == "geo") {
    s <- LG$samples[LG$samples$species == sp & LG$samples$mapped_reads >= min_reads, ]
  } else {
    m <- if (dataset == "paper") LG$paper[[sp]] else LG$udm[[sp]]
    k <- ncol(m)
    s <- if (dataset == "paper") LG$sections[LG$sections$species == sp, ] else
      data.frame(species = sp, section = seq_len(k), sec_label = seq_len(k),
                 pos_cm = seq_len(k), pos_ligule = NA_real_)
    s <- s[order(s$section), ]
    s$sample <- as.character(s$section)
    s$rep <- NA_integer_
    s$mapped_reads <- NA_real_
  }
  rownames(s) <- NULL
  s
}

dataset_matrix <- function(sp, dataset) {
  switch(dataset, geo = LG$expr[[sp]], paper = LG$paper[[sp]], udm = LG$udm[[sp]])
}

transform_values <- function(v, scale) if (scale == "log2") log2(pmax(v, 0) + 1) else v

# Long table of expression values: one row per gene x library (or section).
profile_long <- function(sp, dataset, ids, scale = "log2", min_reads = 0, xaxis = "pos_cm") {
  m <- dataset_matrix(sp, dataset)
  ids <- ids[ids %in% rownames(m)]
  if (!length(ids)) return(NULL)
  s <- dataset_samples(sp, dataset, min_reads)
  v <- m[ids, s$sample, drop = FALSE]
  xa <- if (dataset == "udm") "section" else xaxis
  if (xa == "pos_ligule" && sp != "maize") xa <- "pos_cm"
  data.frame(gene = rep(ids, times = nrow(s)),
             sample = rep(s$sample, each = length(ids)),
             rep = rep(s$rep, each = length(ids)),
             section = rep(s$section, each = length(ids)),
             x = rep(s[[xa]], each = length(ids)),
             raw = as.vector(v),
             value = transform_values(as.vector(v), scale),
             stringsAsFactors = FALSE)
}

summarise_sections <- function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  g <- interaction(d$gene, d$section, drop = TRUE)
  out <- do.call(rbind, lapply(split(d, g), function(z) data.frame(
    gene = z$gene[1], section = z$section[1], x = z$x[1], n = sum(!is.na(z$value)),
    mean = mean(z$value, na.rm = TRUE), sd = sd(z$value, na.rm = TRUE),
    stringsAsFactors = FALSE)))
  out$se <- out$sd / sqrt(out$n)
  out[order(out$gene, out$section), ]
}

# Per-gene z-score using the section means (applied to replicates too).
zscore_long <- function(d, means) {
  st <- do.call(rbind, lapply(split(means, means$gene), function(z)
    data.frame(gene = z$gene[1], mu = mean(z$mean), s = sd(z$mean))))
  st$s[!is.finite(st$s) | st$s == 0] <- 1
  i <- match(d$gene, st$gene); d$value <- (d$value - st$mu[i]) / st$s[i]
  j <- match(means$gene, st$gene)
  means$mean <- (means$mean - st$mu[j]) / st$s[j]
  means$sd <- means$sd / st$s[j]; means$se <- means$se / st$s[j]
  list(d = d, means = means)
}

# Section-mean matrix (genes x sections) for genome-wide operations.
section_mean_matrix <- function(sp, dataset, scale = "log2", min_reads = 0, ids = NULL) {
  m <- dataset_matrix(sp, dataset)
  if (!is.null(ids)) m <- m[intersect(ids, rownames(m)), , drop = FALSE]
  s <- dataset_samples(sp, dataset, min_reads)
  v <- transform_values(m[, s$sample, drop = FALSE], scale)
  secs <- sort(unique(s$section))
  A <- sapply(secs, function(k) (s$section == k) / sum(s$section == k))
  out <- v %*% A
  colnames(out) <- secs
  out
}

# ---- gene lookup -----------------------------------------------------------
gene_choices <- function(sp, idtype = "v5") {
  g <- GENES[GENES$species == sp, ]
  setNames(g$id, if (idtype == "v5") g$label_v5 else g$label)
}

ALIAS <- local({
  g <- GENES
  split_ids <- function(x) strsplit(ifelse(x == "", NA, x), ",\\s*")
  a4 <- split_ids(g$v4); a5 <- split_ids(g$v5)
  al <- rbind(
    data.frame(alias = toupper(g$id), id = g$id),
    data.frame(alias = toupper(unlist(a4)), id = rep(g$id, lengths(a4))),
    data.frame(alias = toupper(unlist(a5)), id = rep(g$id, lengths(a5))),
    data.frame(alias = toupper(g$symbol), id = g$id))
  al <- al[!is.na(al$alias) & al$alias != "", ]
  al[!duplicated(al), ]
})

# Resolve pasted identifiers (v3/v4/v5 maize IDs, symbols, MSU loci).
resolve_ids <- function(txt) {
  orig <- unique(trimws(unlist(strsplit(txt, "[,;\\s]+", perl = TRUE))))
  orig <- sub("_T[0-9]+$|\\.[0-9]+$|_P[0-9]+$", "", orig[orig != ""])
  tok <- toupper(orig)
  hit <- ALIAS[ALIAS$alias %in% tok, ]
  hit$alias <- orig[match(hit$alias, tok)]                     # report IDs as typed
  list(ids = unique(hit$id), missing = orig[!tok %in% ALIAS$alias], hits = hit)
}

gene_species <- function(ids) GENES$species[match(ids, GENES$id)]

gene_table <- function(ids) {
  g <- GENES[match(ids, GENES$id), c("id", "species", "symbol", "description", "chr", "start",
                                     "end", "v4", "v5", "v5_mapping", "orthologs")]
  g$chr <- ifelse(is.na(g$chr), "", paste0("v3 ", g$chr, ":", g$start, "-", g$end))
  g$start <- g$end <- NULL
  names(g) <- c("ID (v3)", "Species", "Symbol", "Description", "Location", "B73 v4", "B73 v5",
                "v3 -> v5", "Orthologs (Supp. Table 3)")
  g
}

# Name shown for a gene. In v5 mode maize genes show their v5 ID; genes that
# share a v5 model with another v3 gene also show the v3 ID.
display_id <- function(ids, idtype = "v5") {
  if (idtype != "v5") return(ids)
  i <- match(ids, GENES$id)
  v5 <- GENES$v5_display[i]
  out <- ifelse(is.na(v5) | v5 == "", ids, v5)
  amb <- GENES$v5_shared[i] %in% TRUE
  out[amb] <- paste0(out[amb], " [", ids[amb], "]")
  out
}

short_name <- function(ids, idtype = "v3") {
  sym <- GENES$symbol[match(ids, GENES$id)]
  d <- display_id(ids, idtype)
  ifelse(!is.na(sym) & sym != "", paste0(d, " (", sym, ")"), d)
}

# Insert v5 ID and mapping-type columns after a v3 ID column (maize tables only).
add_v5 <- function(df, idcol = "id", prefix = "") {
  i <- match(df[[idcol]], GENES$id)
  if (!any(GENES$species[i] == "maize", na.rm = TRUE)) return(df)
  new <- data.frame(GENES$v5[i], GENES$v5_mapping[i], stringsAsFactors = FALSE)
  names(new) <- paste0(prefix, c("v5_id", "v3_to_v5"))
  k <- match(idcol, names(df))
  cbind(df[, seq_len(k), drop = FALSE], new, df[, -seq_len(k), drop = FALSE])
}

download_csv <- function(data_fun, prefix) downloadHandler(
  filename = function() paste0(prefix, "_", format(Sys.Date()), ".csv"),
  content = function(file) write.csv(data_fun(), file, row.names = FALSE))
