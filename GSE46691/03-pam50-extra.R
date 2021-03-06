# ---- Preamble ----
if (!exists("data_id")) data_id <- "GSE46691"
if (!exists("pamdir")) pamdir <- "PAM50"
pamdir_results <- file.path(pamdir, data_id)
pam50_envir <- file.path(pamdir_results, "pam50_GSE46691_envir.RData")
pam50_scores_file <- file.path(pamdir_results, "GSE46691_pam50scores.txt")
if (!file.exists(pam50_envir)) {
  stop("Please run PAM50 scripts first.")
}

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse)

# ---- Read in PAM50 scores ----
pam50_scores <- suppressWarnings(read_tsv(pam50_scores_file)) %>% 
  rename(sample = X1)

pam50_scores_no_her2 <- pam50_scores
pam50_scores_no_her2$Call <- pam50_scores %>% 
  select(Basal:Normal) %>% 
  select(-Her2) %>% 
  apply(1, function(x) names(which.max(x)))

pam50_scores_no_her2_normal <- pam50_scores
pam50_scores_no_her2_normal$Call <- pam50_scores %>% 
  select(Basal:Normal) %>% 
  select(-Her2, -Normal) %>% 
  apply(1, function(x) names(which.max(x)))


# ---- Save Modified PAM50 Scores ----
saveRDS(pam50_scores, file.path("out", "gse46691_pam50_scores.rds"))
saveRDS(pam50_scores_no_her2, file.path("out", "gse46691_pam50_scores_no_her2.rds"))
saveRDS(pam50_scores_no_her2_normal, file.path("out", "gse46691_pam50_scores_no_her2_normal.rds"))

write_tsv(pam50_scores, file.path("out", "gse46691_pam50_scores.tsv"))
write_tsv(pam50_scores_no_her2, file.path("out", "gse46691_pam50_scores_no_her2.tsv"))
write_tsv(pam50_scores_no_her2_normal, file.path("out", "gse46691_pam50_scores_no_her2_normal.tsv"))

# done!
cli::cat_bullet(
  "Modified PAM50 scores available in ", crayon::green(file.path(getwd(), "out")),
  bullet = "tick", bullet_col = "green"
)