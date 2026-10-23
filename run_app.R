# Launch MaizeExpressG:  Rscript run_app.R   (from this folder)
# or, inside R/RStudio:                shiny::runApp("app")
for (p in c("shiny", "bslib", "plotly", "DT", "ggplot2", "data.table", "readxl"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
if (!file.exists("app/data/leafgrad.rds")) {
  source("scripts/01_download_data.R")
  source("scripts/04_liftover_v3_to_v5.R")
  source("scripts/02_build_app_data.R")
}
shiny::runApp("app", launch.browser = TRUE)
