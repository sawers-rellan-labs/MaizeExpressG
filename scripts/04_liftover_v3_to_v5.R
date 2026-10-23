# Position-based v3 -> v5 gene matching, for genes the MaizeGDB cross-reference misses.
#
#   1. exons of every B73 v3 (5b+) gene model are lifted AGPv3 -> AGPv4 -> v5 with the
#      Ensembl Plants assembly chain files;
#   2. lifted exons are intersected (same strand) with the exons of v5 gene models
#      (Zm00001eb.1 protein-coding + non-coding);
#   3. each v3 gene gets its best-overlapping v5 gene plus overlap statistics.
# The result is checked against the genes the cross-reference maps 1:1.
#
# Run from the project root after 01_download_data.R (needs gzip on the PATH):
#   Rscript scripts/04_liftover_v3_to_v5.R
# Output: data-derived/v3_to_v5_liftover.tsv (read by 02_build_app_data.R)

suppressPackageStartupMessages(library(data.table))
raw <- function(f) file.path("data-raw", f)
nochr <- function(x) sub("^chr", "", x)

# ---- gene models -> merged exons (0-based, half-open) ----------------------
# Keeps only exon lines and the transcript lines that link exons to genes.
read_exons <- function(f) {
  cmd <- sprintf("gzip -dc %s | awk -F'\\t' 'NF == 9 && $3 != \"repeat_region\" && $3 != \"CDS\" && $3 !~ /UTR/'",
                 shQuote(raw(f)))
  g <- fread(cmd = cmd, header = FALSE, sep = "\t", quote = "",
             col.names = c("chr", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
  g[, id := sub(".*(^|;)ID=([^;]+).*", "\\2", attr)]
  g[!grepl("(^|;)ID=", attr), id := NA]
  g[, parent := sub(".*(^|;)Parent=([^;]+).*", "\\2", attr)]
  g[!grepl("(^|;)Parent=", attr), parent := NA]
  tx <- g[type != "exon" & !is.na(id) & !is.na(parent), .(tx = id, gene = parent)]
  ex <- g[type == "exon", .(chr = nochr(chr), start = start - 1L, end, strand, tx = parent)]
  ex <- merge(ex, tx, by = "tx")[, tx := NULL]
  merge_intervals(unique(ex))
}

merge_intervals <- function(d) {
  setorder(d, gene, chr, strand, start)
  d[, cm := cummax(end), by = .(gene, chr, strand)]
  d[, blk := cumsum(c(TRUE, start[-1] > cm[-.N])), by = .(gene, chr, strand)]
  d[, .(start = min(start), end = max(end)), by = .(gene, chr, strand, blk)][, blk := NULL][]
}

# ---- chain files -> ungapped alignment blocks ------------------------------
read_chain <- function(f) {
  x <- readLines(gzfile(raw(f)))
  x <- x[x != ""]
  hdr <- startsWith(x, "chain")
  h <- tstrsplit(x[hdr], " ", fixed = TRUE)
  heads <- data.table(cid = seq_len(sum(hdr)), score = as.numeric(h[[2]]),
                      t_chr = nochr(h[[3]]), t_start = as.numeric(h[[6]]),
                      q_chr = nochr(h[[8]]), q_size = as.numeric(h[[9]]), q_strand = h[[10]],
                      q_start = as.numeric(h[[11]]))
  # the files list every chain twice ("1" and "chr1"); keep one copy
  heads <- heads[!duplicated(heads[, .(score, t_chr, t_start, q_chr, q_start)])]
  b <- tstrsplit(x[!hdr], "\t| ", perl = TRUE)
  blocks <- data.table(cid = cumsum(hdr)[!hdr], size = as.numeric(b[[1]]),
                       dt = as.numeric(b[[2]]), dq = as.numeric(b[[3]]))
  blocks[is.na(dt), `:=`(dt = 0, dq = 0)]
  blocks <- merge(blocks, heads, by = "cid")
  blocks[, `:=`(t_s = t_start + cumsum(c(0, head(size + dt, -1))),
                q_s = q_start + cumsum(c(0, head(size + dq, -1)))), by = cid]
  blocks[, .(cid, score, t_chr, t_s, t_e = t_s + size, q_chr, q_s, q_size, q_strand)]
}

# Lift intervals (chr, start, end, strand, gene) through one chain set.
lift <- function(iv, blocks) {
  x <- copy(iv)[, `:=`(fs = start, fe = end - 1)]
  b <- copy(blocks)[, `:=`(bs = t_s, be = t_e - 1)]
  setkey(b, t_chr, bs, be)
  o <- foverlaps(x, b, by.x = c("chr", "fs", "fe"), by.y = c("t_chr", "bs", "be"), nomatch = NULL)
  o[, `:=`(cs = pmax(start, t_s), ce = pmin(end, t_e))]
  o[, `:=`(qs = q_s + (cs - t_s), qe = q_s + (ce - t_s))]
  rev <- o$q_strand == "-"
  o[rev, `:=`(qs = q_size - qe, qe = q_size - qs)]
  o[, strand := fifelse(rev, chartr("+-", "-+", strand), strand)]
  merge_intervals(o[, .(gene, chr = q_chr, start = qs, end = qe, strand)])
}

message("reading gene models")
v3 <- read_exons("Zea_mays.AGPv3.21.gff3.gz")
v5 <- rbind(read_exons("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"),
            read_exons("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.nc.gff3.gz"))
message("reading chains")
c34 <- read_chain("AGPv3_to_B73_RefGen_v4.chain.gz")
c45 <- read_chain("B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain.gz")

message("lifting ", uniqueN(v3$gene), " v3 genes")
v3_bp <- v3[, .(exon_bp = sum(end - start), v3_chr = chr[1]), by = gene]
l5 <- lift(lift(v3, c34), c45)
lifted_bp <- l5[, .(lifted_bp = sum(end - start)), by = gene]

# ---- overlap with v5 exons (same strand) -----------------------------------
v5_bp <- v5[, .(v5_exon_bp = sum(end - start)), by = .(v5 = gene)]
a <- copy(l5)[, `:=`(fs = start, fe = end - 1)]
b <- v5[, .(v5 = gene, chr, strand, bs = start, be = end - 1, s5 = start, e5 = end)]
setkey(b, chr, strand, bs, be)
o <- foverlaps(a, b, by.x = c("chr", "strand", "fs", "fe"), by.y = c("chr", "strand", "bs", "be"),
               nomatch = NULL)
o[, ov := pmin(end, e5) - pmax(start, s5)]
pairs <- o[, .(overlap_bp = sum(ov)), by = .(gene, v5)]
pairs <- merge(merge(pairs, v3_bp[, .(gene, exon_bp)], by = "gene"), v5_bp, by = "v5")
pairs[, `:=`(frac_v3 = overlap_bp / exon_bp, frac_v5 = overlap_bp / v5_exon_bp)]
setorder(pairs, gene, -overlap_bp)

best <- pairs[, .(v5_best = v5[1], frac_v3 = frac_v3[1], frac_v5 = frac_v5[1],
                  v5_other = paste(v5[-1][frac_v3[-1] >= 0.2], collapse = ", ")), by = gene]
res <- merge(merge(v3_bp, lifted_bp, by = "gene", all.x = TRUE), best, by = "gene", all.x = TRUE)
res[is.na(lifted_bp), lifted_bp := 0]
res[, frac_lifted := lifted_bp / exon_bp]
# Evidence classes
res[, liftover := fcase(
  frac_lifted < 0.5,                                "not lifted",
  is.na(v5_best) | frac_v3 < 0.1,                   "lifted, no v5 gene",
  frac_v3 >= 0.5 & frac_v5 >= 0.2,                  "strong",
  frac_v3 >= 0.2,                                   "partial",
  default =                                         "weak")]

# ---- validate against the cross-reference ----------------------------------
x5 <- fread(raw("B73v3_to_B73v5.tsv"), header = FALSE, col.names = c("gene", "ids"))[
  ids != "", .(v5 = strsplit(ids, ",")[[1]]), by = gene]
x5[, n_v3 := .N, by = v5]
one <- x5[, if (.N == 1 && n_v3 == 1) .(v5_xref = v5), by = gene]
val <- merge(res, one, by = "gene")
cat("\nAgreement with the cross-reference on its 1:1 genes:\n")
print(val[, .(genes = .N, best_v5_agrees = round(mean(v5_best == v5_xref, na.rm = TRUE), 4),
              no_call = sum(is.na(v5_best))), by = liftover][order(liftover)])
res <- merge(res, x5[, .(xref_v5 = paste(v5, collapse = ", ")), by = gene], by = "gene", all.x = TRUE)

dir.create("data-derived", showWarnings = FALSE)
fwrite(res[, .(v3 = gene, v3_chr, exon_bp, frac_lifted = round(frac_lifted, 3), v5_best,
               frac_v3 = round(frac_v3, 3), frac_v5 = round(frac_v5, 3), v5_other, liftover, xref_v5)],
       "data-derived/v3_to_v5_liftover.tsv", sep = "\t", na = "")
cat("\nAll v3 genes by evidence class:\n"); print(res[, .N, by = liftover][order(-N)])
cat("\nGenes without a cross-reference match:\n"); print(res[is.na(xref_v5), .N, by = liftover][order(-N)])
