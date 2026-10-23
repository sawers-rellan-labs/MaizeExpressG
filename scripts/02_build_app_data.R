# Build app/data/leafgrad.rds from the files in data-raw/.
# Run from the project root after 01_download_data.R:
#   Rscript scripts/02_build_app_data.R

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

raw <- function(f) file.path("data-raw", f)
read_gz <- function(f, ...) setDT(read.delim(gzfile(raw(f)), check.names = FALSE, ...))

# ---------------------------------------------------------------------------
# 1. Per-library FPKM from GEO (GSE54274)
# ---------------------------------------------------------------------------
read_geo <- function(f) {
  d <- read_gz(f)
  d <- d[, !grepl("^V[0-9]+$|^$", names(d)), with = FALSE]   # trailing tab column
  m <- as.matrix(d[, -1])
  rownames(m) <- d$Gene
  m[is.na(m)] <- 0
  m <- m[rowSums(m) > 0, , drop = FALSE]                       # drop never-detected genes
  m
}
expr <- list(maize = read_geo("GSE54274_maize_leaf_section_fpkm.txt.gz"),
             rice  = read_geo("GSE54274_rice_leaf_section_fpkm.txt.gz"))

lib <- fread(raw("SuppTable17_library_stats.csv"))

samples <- rbindlist(lapply(names(expr), function(sp) {
  cn <- colnames(expr[[sp]])
  data.table(sample    = cn,
             species   = sp,
             rep       = as.integer(sub(".*_rep([0-9]+)_.*", "\\1", cn)),
             sec_label = as.integer(sub(".*_sec([0-9]+)$", "\\1", cn)))
}))
# GEO/SuppTable17 maize labels run 1-4, 6-16 (no "sec5") for the fifteen 1-cm
# segments; section = rank of the label, which matches columns M1..M15 of
# Supplementary Table 13 (checked by correlation of section means).
samples[, section := match(sec_label, sort(unique(sec_label))), by = species]
seg_cm <- c(maize = 1, rice = 2)                               # segment length (GEO design)
samples[, pos_cm := (section - 0.5) * seg_cm[species]]
# Supp. Table 1 numbers maize sections -4..-1, 1..11 relative to the leaf-2
# ligule, i.e. the ligule sits at the boundary of sections 4 and 5 (4 cm).
samples[, pos_ligule := fifelse(species == "maize", pos_cm - 4, NA_real_)]
samples <- merge(samples, lib[, .(species, rep, sec_label, total_reads, mapped_reads)],
                 by = c("species", "rep", "sec_label"), all.x = TRUE, sort = FALSE)
setorder(samples, species, section, rep)
for (sp in names(expr)) expr[[sp]] <- expr[[sp]][, samples[species == sp, sample]]
stopifnot(!anyNA(samples$mapped_reads))

sections <- unique(samples[, .(species, section, sec_label, pos_cm, pos_ligule)])

# ---------------------------------------------------------------------------
# 2. Paper values: Supplementary Table 13 (RPKM section means + UDM)
# ---------------------------------------------------------------------------
st13 <- raw("SuppTable13_RPKM_gradients.xlsx")
sheet_matrix <- function(sheet, skip, id_col = 1) {
  d <- as.data.frame(read_excel(st13, sheet = sheet, skip = skip))
  m <- as.matrix(sapply(d[, -seq_len(id_col)], function(v) suppressWarnings(as.numeric(v))))
  rownames(m) <- d[[id_col]]
  m[!is.na(rownames(m)), , drop = FALSE]
}
paper <- list(maize = sheet_matrix("Maize", 2), rice = sheet_matrix("Rice", 0))
colnames(paper$maize) <- seq_len(ncol(paper$maize))
colnames(paper$rice)  <- seq_len(ncol(paper$rice))

u <- as.data.frame(read_excel(st13, sheet = "Unified"))
udm <- lapply(c(maize = "Maize", rice = "Rice"), function(sp) {
  x <- u[u$Species == sp, ]
  m <- as.matrix(x[, grep("^Gradient", names(x))])
  dimnames(m) <- list(x$ID, seq_len(ncol(m)))
  m
})

# ---------------------------------------------------------------------------
# 3. Orthologs (Supplementary Table 3: rice locus -> one or more maize genes)
# ---------------------------------------------------------------------------
o <- as.data.frame(read_excel(raw("SuppTable03_orthologs.xlsx"), col_names = FALSE, skip = 3,
                               col_types = "text", .name_repair = "minimal"))
o <- as.data.table(o)
setnames(o, 1, "rice")
orthologs <- melt(o, id.vars = "rice", value.name = "maize", na.rm = TRUE)[, .(rice, maize = trimws(maize))]
orthologs <- unique(orthologs[grepl("^LOC_Os", rice) & maize != ""])

# ---------------------------------------------------------------------------
# 4. Gene annotation
# ---------------------------------------------------------------------------
clean <- function(x) { x <- trimws(x); x[is.na(x)] <- ""; x }

v3 <- read_gz("Zea_mays.AGPv3.21.full_gene_data.txt.gz", quote = "", comment.char = "",
              strip.white = TRUE, colClasses = "character")
setnames(v3, trimws(names(v3)))
v3 <- unique(v3[, .(id = gene_model, chr = gene_model_chr, start = as.integer(gene_model_start),
                    end = as.integer(gene_model_end), symbol = clean(locus_symbol),
                    name = clean(locus_name), products = clean(gene_products))], by = "id")
v3[grepl("^Zm0000", symbol), symbol := ""]                      # v4 IDs used as placeholders
v3[, products := gsub('^\\{|\\}$|"', "", products)]
v3[, products := gsub(",", "; ", products)]

xref <- function(f, max_ids = 3) {
  x <- fread(raw(f), header = FALSE, col.names = c("v3", "ids"))
  x[, n := lengths(strsplit(ids, ","))]
  x[n <= max_ids, .(v3, ids = gsub(",", ", ", ids))]           # many-to-many blocks are not useful
}
x4 <- xref("B73v3_to_B73v4.tsv")

# v3 -> v5: MaizeGDB pan-gene cross-reference (covers the v3 filtered gene set).
# Keep every match and record its shape: 1:1, split (one v3 -> several v5),
# merged (several v3 -> one v5) or complex (several -> several).
x5_pairs <- fread(raw("B73v3_to_B73v5.tsv"), header = FALSE, col.names = c("v3", "ids"))[
  ids != "", .(v5 = strsplit(ids, ",")[[1]]), by = v3]
x5_pairs[, n_v3 := .N, by = v5]
x5 <- x5_pairs[, .(v5 = paste(v5, collapse = ", "),
                   v5_mapping = if (all(n_v3 == 1)) (if (.N == 1) "1:1" else "split") else
                     if (.N == 1) "merged" else "complex"), by = v3]

nlm_lines <- readLines(raw("Uniprot_ProtNLM_maize_v5.txt"), warn = FALSE)
nlm_lines <- nlm_lines[lengths(strsplit(nlm_lines, "\t")) == 4]          # file has one line broken in two
nlm <- fread(text = nlm_lines, sep = "\t", header = FALSE, quote = "",
             col.names = c("v5", "uniprot", "uniprot_name", "protnlm"))
nlm[, protnlm := gsub('^"|"$', "", protnlm)]
nlm <- unique(nlm, by = "v5")

rice_annot <- fread(raw("rice_MSU7_locus_brief_info.txt"), quote = "")
rice_annot <- rice_annot[is_representative == "Y",
  .(id = locus, chr, start, end = stop, description = annotation)]
rice_annot <- unique(rice_annot, by = "id")

c4 <- as.data.table(read_excel(raw("SuppTable07_C4_shuttle_genes.xlsx"), skip = 2))
c4 <- c4[!is.na(`Maize Gene ID`), .(id = `Maize Gene ID`, abbrev = Abbreviation, full_name = `Full Name`)]

maize_ids <- sort(unique(c(rownames(expr$maize), rownames(paper$maize), rownames(udm$maize))))
rice_ids  <- sort(unique(c(rownames(expr$rice),  rownames(paper$rice),  rownames(udm$rice))))

mg <- data.table(id = maize_ids, species = "maize")
mg <- merge(mg, v3, by = "id", all.x = TRUE, sort = FALSE)
mg <- merge(mg, x4[, .(id = v3, v4 = ids)], by = "id", all.x = TRUE, sort = FALSE)
mg <- merge(mg, x5[, .(id = v3, v5, v5_mapping)], by = "id", all.x = TRUE, sort = FALSE)
mg[is.na(v5_mapping), v5_mapping := "no v5 match"]
# Genes the cross-reference misses: position-based match (scripts/04_liftover_v3_to_v5.R)
lo_file <- file.path("data-derived", "v3_to_v5_liftover.tsv")
if (file.exists(lo_file)) {
  lo <- fread(lo_file, na.strings = "")
  mg <- merge(mg, lo[, .(id = v3, lo_v5 = v5_best, lo_class = liftover, lo_frac_v3 = frac_v3,
                         lo_frac_v5 = frac_v5)], by = "id", all.x = TRUE, sort = FALSE)
  use <- mg$v5_mapping == "no v5 match" & mg$lo_class %in% c("strong", "partial")
  mg[use, `:=`(v5 = lo_v5, v5_mapping = paste0("liftover (", lo_class, ")"))]
} else {
  mg[, `:=`(lo_v5 = NA_character_, lo_class = NA_character_, lo_frac_v3 = NA_real_, lo_frac_v5 = NA_real_)]
}
mg[, v5_first := sub(",.*", "", v5)]
mg <- merge(mg, nlm[, .(v5_first = v5, protnlm)], by = "v5_first", all.x = TRUE, sort = FALSE)
mg[, v5_first := NULL]
for (cl in c("symbol", "name", "products", "v4", "v5", "protnlm")) set(mg, which(is.na(mg[[cl]])), cl, "")
mg[, description := fifelse(name != "", name, fifelse(products != "", products, protnlm))]

rg <- merge(data.table(id = rice_ids, species = "rice"), rice_annot, by = "id", all.x = TRUE, sort = FALSE)
rg[is.na(description), description := ""]
rg[, `:=`(symbol = "", v4 = "", v5 = "", v5_mapping = "", lo_v5 = NA_character_, lo_class = NA_character_,
          lo_frac_v3 = NA_real_, lo_frac_v5 = NA_real_)]

orth_m <- orthologs[, .(orthologs = paste(sort(unique(rice)), collapse = ", ")), by = .(id = maize)]
orth_r <- orthologs[, .(orthologs = paste(sort(unique(maize)), collapse = ", ")), by = .(id = rice)]
mg <- merge(mg, orth_m, by = "id", all.x = TRUE, sort = FALSE)
rg <- merge(rg, orth_r, by = "id", all.x = TRUE, sort = FALSE)

# Maize genes with no description borrow their rice ortholog's MSU annotation
mz_from_rice <- merge(orthologs, rice_annot[, .(rice = id, description)], by = "rice")[
  description != "" & !grepl("expressed protein|hypothetical protein", description),
  .(rice_desc = description[1]), by = .(id = maize)]
mg <- merge(mg, mz_from_rice, by = "id", all.x = TRUE, sort = FALSE)
mg[description == "" & !is.na(rice_desc), description := paste0("[rice ortholog] ", rice_desc)]
mg[, rice_desc := NULL]

cols <- c("id", "species", "symbol", "description", "chr", "start", "end", "v4", "v5", "v5_mapping",
          "lo_v5", "lo_class", "lo_frac_v3", "lo_frac_v5", "orthologs")
genes <- rbind(mg[, ..cols], rg[, ..cols], fill = TRUE)
genes[is.na(orthologs), orthologs := ""]
genes[, c4 := id %in% c4$id]
genes[, in_geo   := id %in% c(rownames(expr$maize),  rownames(expr$rice))]
genes[, in_paper := id %in% c(rownames(paper$maize), rownames(paper$rice))]
genes[, in_udm   := id %in% c(rownames(udm$maize),   rownames(udm$rice))]
genes[, v5_n := fifelse(v5 == "", 0L, lengths(strsplit(v5, ", ")))]
genes[, v5_display := fifelse(v5_n == 0, "",
  paste0(sub(",.*", "", v5), fifelse(v5_n > 1, paste0(" +", v5_n - 1), "")))]
genes[, v5_n := NULL]
# v5 models claimed by more than one v3 gene (display then keeps the v3 ID)
sh <- genes[species == "maize" & v5 != "", .(v5 = strsplit(v5, ", ")[[1]]), by = id]
sh[, n := uniqueN(id), by = v5]
genes[, v5_shared := id %in% sh[n > 1, id]]
lab_tail <- genes[, paste0(fifelse(symbol != "", paste0(" | ", symbol), ""),
                           fifelse(description != "", paste0(" | ", substr(description, 1, 90)), ""))]
genes[, label := paste0(id, fifelse(v5_display != "", paste0(" | ", v5_display), ""), lab_tail)]
genes[, label_v5 := paste0(fifelse(v5_display != "", paste0(v5_display, " | ", id),
                                   fifelse(species == "maize", paste(id, "(no v5)"), id)), lab_tail)]
setkey(genes, id)

# ---------------------------------------------------------------------------
# 5. Metabolites / enzyme activities (Supplementary Tables 1 and 2)
# ---------------------------------------------------------------------------
read_metab <- function(f, n_sec) {
  d <- as.data.frame(read_excel(raw(f), col_names = FALSE, skip = 2, .name_repair = "minimal"))
  d <- d[!is.na(d[[1]]), ]
  num <- function(v) suppressWarnings(as.numeric(v))                   # "NA" strings -> NA
  mean <- as.matrix(sapply(d[, 3 + seq_len(n_sec)], num))
  sd   <- as.matrix(sapply(d[, 3 + n_sec + seq_len(n_sec)], num))
  trait <- make.unique(paste0(trimws(d[[1]]), " [", trimws(d[[3]]), "]"), sep = " #")
  dimnames(mean) <- dimnames(sd) <- list(trait, seq_len(n_sec))
  list(info = data.table(trait = trait, name = trimws(d[[1]]), unit = trimws(d[[2]]), platform = trimws(d[[3]])),
       mean = mean, sd = sd)
}
metabolites <- list(maize = read_metab("SuppTable01_maize_metabolites.xlsx", 15),
                    rice  = read_metab("SuppTable02_rice_metabolites.xlsx", 11))

# ---------------------------------------------------------------------------
leafgrad <- list(
  expr = expr, samples = as.data.frame(samples), sections = as.data.frame(sections),
  paper = paper, udm = udm, genes = genes, orthologs = orthologs,
  metabolites = metabolites, c4 = c4,
  info = list(built = format(Sys.time()), geo = "GSE54274", doi = "10.1038/nbt.3019")
)
dir.create("app/data", showWarnings = FALSE, recursive = TRUE)
saveRDS(leafgrad, "app/data/leafgrad.rds", compress = "xz")

cat(sprintf("maize: %d genes x %d libraries | rice: %d genes x %d libraries\n",
            nrow(expr$maize), ncol(expr$maize), nrow(expr$rice), ncol(expr$rice)))
cat(sprintf("paper ST13: maize %d, rice %d | UDM: maize %d, rice %d | ortholog pairs: %d\n",
            nrow(paper$maize), nrow(paper$rice), nrow(udm$maize), nrow(udm$rice), nrow(orthologs)))
cat(sprintf("annotated maize genes: %d/%d | rice: %d/%d\n",
            sum(mg$description != ""), nrow(mg), sum(rg$description != ""), nrow(rg)))
cat("wrote app/data/leafgrad.rds (", round(file.size("app/data/leafgrad.rds") / 2^20, 1), "MB )\n")
