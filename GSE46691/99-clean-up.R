# Clean Up Scripts
if (!exists("data_dir")) data_dir <- "data"

# ---- Clean everything ----
gse46691_clean_all <- function() {
  cli::cat_line(
    "This will delete input and output data. You will need to re-download data after running this.",
    bullet = "square_small_filled", col = "red")
  if (!interactive() || yesno::yesno2("Do you want to continue?"))
  if (dir.exists(data_dir)) unlink(data_dir, recursive = TRUE)
  if (dir.exists("out")) unlink("out", recursive = TRUE)
  if (dir.exists("PAM50")) unlink("PAM50", recursive = TRUE)
  if (file.exists("PAM50.zip")) unlink("PAM50.zip")
}

# ---- Clean outputs ----
gse46691_clean_outputs <- function() {
  if (dir.exists("out")) unlink("out", recursive = TRUE)
  if (dir.exists("PAM50")) unlink("PAM50", recursive = TRUE)
}