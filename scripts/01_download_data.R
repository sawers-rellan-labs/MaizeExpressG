# Download the expression data and annotation used by the app.
#
# Source study: Wang L, et al. (2014) Comparative analyses of C4 and C3
# photosynthesis in developing leaves of maize and rice. Nat Biotechnol 32:1158-1165.
# https://doi.org/10.1038/nbt.3019   GEO: GSE54274   SRA: SRP018823
#
# Run from the project root:  Rscript scripts/01_download_data.R

options(timeout = 600)
dest <- "data-raw"
dir.create(dest, showWarnings = FALSE)

esm <- function(n, ext) sprintf(
  "https://media.springernature.com/original/springer-static/esm/art%%3A10.1038%%2Fnbt.3019/MediaObjects/41587_2014_BFnbt3019_MOESM%d_ESM.%s", n, ext)

files <- c(
  # Per-library FPKM (cuffdiff v2.1.1; maize AGPv3, rice MSU7) - GEO GSE54274
  "GSE54274_maize_leaf_section_fpkm.txt.gz" = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54274/suppl/GSE54274_maize_leaf_section_fpkm.txt.gz",
  "GSE54274_rice_leaf_section_fpkm.txt.gz"  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54274/suppl/GSE54274_rice_leaf_section_fpkm.txt.gz",
  # Paper supplementary material
  "SuppText_Figures.pdf"                  = esm(6,  "pdf"),
  "SuppTable01_maize_metabolites.xlsx"    = esm(7,  "xlsx"),
  "SuppTable02_rice_metabolites.xlsx"     = esm(8,  "xlsx"),
  "SuppTable03_orthologs.xlsx"            = esm(9,  "xlsx"),
  "SuppTable07_C4_shuttle_genes.xlsx"     = esm(10, "xlsx"),
  "SuppTable13_RPKM_gradients.xlsx"       = esm(14, "xlsx"),
  "SuppTable18_UDM_anchor_genes.xlsx"     = esm(18, "xlsx"),
  # Gene annotation
  "Zea_mays.AGPv3.21.full_gene_data.txt.gz" = "https://download.maizegdb.org/B73_RefGen_v3/Zea_mays.AGPv3.21.full_gene_data.txt.gz",
  "B73v3_to_B73v4.tsv"                    = "https://download.maizegdb.org/Pan-genes/B73_gene_xref/B73v3_to_B73v4.tsv",
  "B73v3_to_B73v5.tsv"                    = "https://download.maizegdb.org/Pan-genes/B73_gene_xref/B73v3_to_B73v5.tsv",
  "Uniprot_ProtNLM_maize_v5.txt"          = "https://download.maizegdb.org/GeneFunction_and_Expression/ProtNLM/Uniprot_ProtNLM_maize_v5.txt",
  "rice_MSU7_locus_brief_info.txt"        = "https://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.locus_brief_info.7.0",
  # Position-based v3 -> v5 matching (scripts/04_liftover_v3_to_v5.R)
  "Zea_mays.AGPv3.21.gff3.gz"             = "https://download.maizegdb.org/B73_RefGen_v3/Zea_mays.AGPv3.21.gff3.gz",
  "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"    = "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz",
  "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.nc.gff3.gz" = "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.nc.gff3.gz",
  "AGPv3_to_B73_RefGen_v4.chain.gz"       = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/assembly_chain/zea_mays/AGPv3_to_B73_RefGen_v4.chain.gz",
  "B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain.gz" = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/assembly_chain/zea_mays/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain.gz"
)

for (f in names(files)) {
  out <- file.path(dest, f)
  if (file.exists(out) && file.size(out) > 0) { message("have  ", f); next }
  message("fetch ", f)
  download.file(files[[f]], out, mode = "wb", quiet = TRUE)
}

# SuppTable17_library_stats.csv (reads per library) was transcribed from
# Supplementary Table 17 in SuppText_Figures.pdf and ships with the project.
stopifnot(file.exists(file.path(dest, "SuppTable17_library_stats.csv")))
message("done")
