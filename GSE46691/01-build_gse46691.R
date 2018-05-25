# ---- Load Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(Biobase)
pacman::p_load(tidyverse)
pacman::p_load(GEOquery)
source("functions_gse46691.R")
source("pam50_genes.R")

# ---- Gather Data ----
data_dir <- "data"
source("gather_gse46691.R")

# ---- Load GPL and Series Matrix ----
gse_46691 <- build_gse_46691(file_exprs = NULL, data_dir = data_dir)

# ---- Parse and Tidy Gene Assignment ----
# Get gene_assignment from GPL5188 annotation
gse_46691_genes <- pData(gse_46691$featureData) %>% 
  select(ID, gene_assignment) %>% 
  as_tibble() %>% 
  tidy_gene_assignment(gene_assignment) %>% 
  # gene_assignment_1 and _2 are now a list_cols but we can splat them out
  # which also drops anything that didn't have an assignment, i.e. "---"
  tidyr::unnest() %>% 
  select(-gene_assignment)

# ---- Get probe to PAM50 gene map ----
gse_46691_pam50 <- gse_46691_genes %>% 
  filter(gene_assignment_1 %in% pam50_annotation$GeneName |
           gene_assignment_2 %in% pam50_annotation$GeneName) %>% 
  left_join(pam50_annotation[, c("pcrID", "GeneName")], by = c("gene_assignment_2" = "GeneName")) %>%
  rename(gene_name_pam50 = pcrID)

# ---- Read Quantile Normalized Data ----
# ...while filtering in IDs on the fly
gse_46691$exprs <- read_tsv_filtered(
  file.path(data_dir, "GSE46691_quantile_normalized.txt"),
  ID_REF %in% gse_46691_pam50$ID,
  chunk_size = 10000, 
  col_types = readr::cols(ID_REF = col_integer(), .default = col_double())
)

# ---- Save Checkpoint ----
dir.create("out")
saveRDS(gse_46691, file.path("out", "gse46691_checkpoint.rds"))
saveRDS(gse_46691$exprs, file.path("out", "gse46691_pam50_exprs.rds"))
# gse_46691 <- readRDS(file.path("out", "gse46691_checkpoint.rds"))
# gse_46691$exprs <- readRDS(file.path("out", "gse46691_pam50_exprs.rds"))


# ---- Process Data for PAM ----

# PAM50 Pre-processing Steps
gse_46691$exprs_prepped <- gse_46691$exprs %>% 
  tidyr::gather("sample", "value", -ID_REF) %>% 
  # PAM50 step 2: Log transform expression estimates
  # mutate(value = log(value)) %>% 
  # PAM50 step 3: Median center
  group_by(ID_REF) %>% median_center(value) %>% 
  ungroup()

# Merge probeset observations and gene information
gse_46691$qn <- gse_46691$exprs_prepped %>% 
  mutate(sample = str_replace(sample, "\\.CEL$", "")) %>% 
  left_join(gse_46691_pam50, by = c("ID_REF" = "ID"))

# PAM50: Output TSV File
gse_46691$qn_pam50 <- gse_46691$qn %>% 
  tidyr::spread("sample", "value") %>% 
  select(-contains("gene_assignment"), -ID_REF) %>% 
  select(gene = gene_name_pam50, everything())

write_tsv(gse_46691$qn_pam50, file.path("out", "gse46691_prepped_pam50.tsv"))


# ---- Additional formatting for our use ----

# Do our own summarization of values across probes
gse_46691$qn_summarized <- gse_46691$qn %>% 
  group_by(sample, gene_assignment_2) %>% 
  summarize(value = median(value))

# Clean up pheno data and then merge with gse
gse_46691$qn_ready <- pData(gse_46691$phenoData) %>% 
  as_tibble() %>% 
  biogroom::remove_common(description) %>% 
  biogroom:::clean_channel_vars() %>% 
  select(-matches("title|characteristics|supplementary")) %>% 
  mutate_all(as.character) %>% 
  readr::type_convert() %>% 
  left_join(gse_46691$qn_summarized, ., by = c("sample" = "description"))

saveRDS(gse_46691, file.path("out", "gse46691.rds"))
# gse_46691 <- readRDS(file.path("out", "gse46691.rds"))

