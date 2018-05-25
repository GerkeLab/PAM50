## This version uses GEOquery but doesn't get all the files
# if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
# pacman::p_load(GEOquery)
# 
# dir.create(data_dir, showWarnings = FALSE)
# getGEOfile("GSE46691", AnnotGPL = TRUE, destdir = data_dir, amount = "full", getGPL = TRUE)
# getGEOfile("GPL5188", destdir = data_dir, amount = "full")
# getGEOSuppFiles(
#   GEO = "GSE46691", 
#   makeDirectory = FALSE, 
#   baseDir = data_dir, 
#   fetch_files = TRUE, 
#   filter_regex = "quantile_normalized")

cli::cat_rule("Gathering GSE46691 Data Files")
gse46991_urls <- c(
  "GSE46691_quantile_normalized.txt.gz" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46691&format=file&file=GSE46691%5Fquantile%5Fnormalized%2Etxt%2Egz",
  "GSE46691_family.soft.gz" = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46691/soft/GSE46691_family.soft.gz",
  "GSE46691_series_matrix.txt.gz" = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46691/matrix/GSE46691_series_matrix.txt.gz",
  "GPL5188.soft" = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL5188&form=text&view=full"
)
if (!exists("data_dir")) data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE)
cat(names(gse46991_urls), 
    "GSE46691_quantile_normalized.txt", 
    "GSE46691.soft",
    "GSE46691_quantile_normalized.txt.gz-ID_REF.txt",
    file = file.path(data_dir, ".gitignore"), sep = "\n")

download_file <- function(url, dest.path, dest.name) {
  if (file.exists(file.path(dest.path, dest.name))) {
    cli::cat_bullet("Skipping '", crayon::green(dest.name), "' because it's already downloaded in ", 
                    crayon::yellow(dest.path), bullet = "pointer")
  } else {
    cli::cat_bullet("Downloading '", crayon::blue(dest.name), "' into ", 
                    crayon::yellow(dest.path), bullet = "pointer")
    download.file(url, file.path(dest.path, dest.name))
    cli::cat_bullet(crayon::green(dest.name), " download complete", bullet = "tick")
  }
}

for (i in seq_along(gse46991_urls)) {
  download_file(gse46991_urls[i], data_dir, names(gse46991_urls)[i])
}

if (!"GSE46691_quantile_normalized.txt" %in% dir(data_dir)) {
  cli::cat_bullet("Unzipping ", crayon::green("GSE46691_series_matrix.txt.gz"),
                  bullet = "pointer")
  owd <- setwd(data_dir)
  system("gunzip -c GSE46691_quantile_normalized.txt.gz > GSE46691_quantile_normalized.txt")
  setwd(owd)
}

cli::cat_bullet("GSE46691 data files are available in ", 
                crayon::yellow(file.path(getwd(), data_dir)), bullet = "tick")