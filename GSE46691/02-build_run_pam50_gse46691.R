# ---- PAM50 ----
# 
# Downloads and runs the PAM50 analysis on appropriately processed data. This
# script automatically downloads the PAM50 scripts, organizes input and output
# data files, and creates the script to run the PAM50 process. In general, this
# script can be modified to run the PAM50 analysis on any pre-processed data by
# changing the variables in the PAM50 Parameters section.


# ---- Libraries ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_install("ctc", force = FALSE)
pacman::p_install("heatmap.plus", force = FALSE)
source("functions/utils.R")

# ---- PAM50 Parameters ----
# PAM50 Params (see `PAM50/bioclassifier_example/bioclassifier_subtypeAssignment.R`)
calibrationParameters <- -1
hasClinical           <- FALSE
collapseMethod        <- "iqr"
# Data for this run
pam50_data_file       <- "out/gse46691_prepped_pam50.tsv"
data_id               <- "GSE46691"

# ---- Download and Extract PAM50 Scripts ----
if (!file.exists("PAM50.zip")) 
  download.file("https://genome.unc.edu/pubsup/breastGEO/PAM50.zip", "PAM50.zip")

pamdir <- ask_if_overwrite("PAM50")
unzip("PAM50.zip", exdir = pamdir)

# ---- Copy prepped data to PAM dir ----
dir.create(file.path(pamdir, data_id))
file.copy(pam50_data_file, file.path(pamdir, data_id))

# ---- Modify PAM50 example script to use our data ----
collapseMethod <- paste0('"', collapseMethod, '"')
pam50_R_code <- c(
  "## Modified from bioclassifier_example/bioclassifier_subtypeAssignment.R",
  "library(ctc)",
  "library(heatmap.plus)",
  "paramDir <- 'bioclassifier_R'",
  glue::glue('inputDir <- "{data_id}"'),
  glue::glue('inputFile <- "{basename(pam50_data_file)}"'),
  glue::glue('short <- "{data_id}"'),
  glue::glue('calibrationParameters <- {calibrationParameters}'),
  glue::glue('hasClinical <- {hasClinical}'),
  glue::glue('collapseMethod <- {collapseMethod}'),
  'source(file.path(paramDir, "subtypePrediction_functions.R"))',
  'source(file.path(paramDir, "subtypePrediction_distributed.R"))',
  glue::glue('save.image(file.path("{data_id}", "pam50_{data_id}_envir.RData"))')
)
pam50_R <- file.path(pamdir, glue::glue("{data_id}_PAM50.R"))
writeLines(pam50_R_code, pam50_R)

# ---- Run PAM50 or give instructions ----
if (interactive() && yesno::yesno2(crayon::yellow(cli::symbol$bullet), " Run PAM50 now?")) {
  run_pam50 <- function(path_pam50_R) {
    owd <- setwd(dirname(path_pam50_R))
    on.exit(setwd(owd))
    system(paste("Rscript", basename(path_pam50_R)), wait = TRUE)
  }
  pam50_rs <- run_pam50(pam50_R) # 0 = GOOD, otherwise bad
  if (pam50_rs) {
    cli::cat_bullet(
      "Errors occurred while running PAM50, review output in ", crayon::red(pamdir),
      bullet = "cross", bullet_col = "red"
    )
  } else {
    cli::cat_bullet(
      "PAM50 analysis complete, review output in ", crayon::green(file.path(getwd(), pamdir, data_id)),
      bullet = "tick", bullet_col = "green"
    )
  }
} else {
  cli::cat_bullet(
    "Run PAM50 script: ", crayon::cyan(filel.path(getwd(), pam50_R)),
    bullet = "pointer"
  )
}