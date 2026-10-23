# Write the maize expression tables with B73 v5 gene IDs to exports/.
# Run from the project root after 02_build_app_data.R:
#   Rscript scripts/03_export_v5_tables.R
#
# Expression was quantified on the v3 gene models, so every table keeps one row
# per v3 gene and adds its v5 ID(s) and the mapping type (MaizeGDB cross-reference):
#   1:1          one v3 gene <-> one v5 gene
#   split        one v3 gene -> several v5 genes (all listed)
#   merged       several v3 genes -> one v5 gene (those rows share a v5 ID)
#   complex      several v3 genes <-> several v5 genes
#   liftover (strong / partial)
#                not in the cross-reference; matched by genome position instead
#                (scripts/04_liftover_v3_to_v5.R): lifted v3 exons overlap the v5
#                gene's exons by >= 50% (strong) or >= 20% (partial) of their length
#   no v5 match  neither (mostly v3 working-set models with no v5 gene model)
# For a table keyed by v5 ID, keep the rows with v3_to_v5 == "1:1".

lg <- readRDS("app/data/leafgrad.rds")
g <- lg$genes[lg$genes$species == "maize", ]
dir.create("exports", showWarnings = FALSE)

with_v5 <- function(m) {
  i <- match(rownames(m), g$id)
  data.frame(v5_id = g$v5[i], v3_id = rownames(m), v3_to_v5 = g$v5_mapping[i], symbol = g$symbol[i],
             signif(m, 6), check.names = FALSE, stringsAsFactors = FALSE)
}

mapping <- data.frame(v3_id = g$id, v5_id = g$v5, v3_to_v5 = g$v5_mapping, v4_id = g$v4,
                      symbol = g$symbol, description = g$description,
                      liftover_v5 = g$lo_v5, liftover_class = g$lo_class,
                      liftover_frac_v3_exons = g$lo_frac_v3, liftover_frac_v5_exons = g$lo_frac_v5,
                      in_GEO = g$in_geo, in_SuppTable13 = g$in_paper, in_UDM = g$in_udm)
write.csv(mapping, "exports/maize_v3_to_v5_mapping.csv", row.names = FALSE)

geo <- lg$expr$maize
s <- lg$samples[match(colnames(geo), lg$samples$sample), ]
colnames(geo) <- sprintf("%s|section%02d|%.1fcm", colnames(geo), s$section, s$pos_cm)
write.csv(with_v5(geo), "exports/maize_GEO_FPKM_by_library_v5.csv", row.names = FALSE)

paper <- lg$paper$maize; colnames(paper) <- paste0("M", colnames(paper))
write.csv(with_v5(paper), "exports/maize_SuppTable13_RPKM_v5.csv", row.names = FALSE)

udm <- lg$udm$maize; colnames(udm) <- paste0("UDM_", colnames(udm))
write.csv(with_v5(udm), "exports/maize_UDM_v5.csv", row.names = FALSE)

writeLines(c(
  "Source of the expression data in this folder - please cite:",
  "",
  paste("Wang L, Czedik-Eysenberg A, Mertz RA, Si Y, Tohge T, Nunes-Nesi A, Arrivault S, Dedow LK, Bryant DW,",
        "Zhou W, Xu J, Weissmann S, Studer A, Li P, Zhang C, LaRue T, Shao Y, Ding Z, Sun Q, Patel RV, Turgeon R,",
        "Zhu X, Provart NJ, Mockler TC, Fernie AR, Stitt M, Liu P, Brutnell TP (2014) Comparative analyses of C4",
        "and C3 photosynthesis in developing leaves of maize and rice. Nature Biotechnology 32(11):1158-1165.",
        "doi:10.1038/nbt.3019. PMID 25306245."),
  "",
  "maize_GEO_FPKM_by_library_v5.csv: GEO GSE54274 (GSE54274_maize_leaf_section_fpkm.txt.gz)",
  "maize_SuppTable13_RPKM_v5.csv, maize_UDM_v5.csv: Supplementary Table 13 of the paper",
  "",
  "B73 v5 IDs: MaizeGDB v3-to-v5 cross-reference (https://download.maizegdb.org/Pan-genes/B73_gene_xref/),",
  "plus position-based matches from Ensembl Plants assembly chain files and MaizeGDB Zm00001eb.1 gene models",
  "(scripts/04_liftover_v3_to_v5.R). Annotation: MaizeGDB AGPv3 gene data and UniProt ProtNLM names."),
  "exports/CITATION.txt")

tab <- function(ids) table(factor(g$v5_mapping[match(ids, g$id)],
                                  c("1:1", "split", "merged", "complex", "liftover (strong)",
                                    "liftover (partial)", "no v5 match")))
print(rbind(GEO = tab(rownames(lg$expr$maize)), SuppTable13 = tab(rownames(lg$paper$maize)),
            UDM = tab(rownames(lg$udm$maize))))
cat("wrote", list.files("exports", full.names = TRUE), sep = "\n  ")
