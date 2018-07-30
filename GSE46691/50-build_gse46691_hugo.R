# ---- OPTIONS ----
# 
# * Summarize Function *
# When HUGO Gene names are associated with multiple probesets, the probset
# expression values will be summarized using the following function:
SUMMARIZE_FUNCTION <- median

# * HUGO Dictionary URL*
# This URL points to most recent HUGO gene names information.
# See https://beta.genenames.org/download/custom/ for more information.
HUGO_DICT_URL <- "https://beta.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

# ---- Load Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(Biobase)
pacman::p_load(tidyverse)
pacman::p_load(GEOquery)

# Source functions for these scripts
source_files <- dir("functions", pattern = "\\.R$", full.names = TRUE)
purrr::walk(source_files, source, local = globalenv())

# ---- Gather Data ----
data_dir <- gather_gse46691("data")

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

gse_46691_genes2 <- tidyr::gather(gse_46691_genes, drop, gene_name, -ID) %>% 
  select(-drop)

# ---- Get Probeset Annotations ----
# http://www.affymetrix.com/support/technical/byproduct.affx?product=huexon-st
# Download: http://www.affymetrix.com/Auth/analysis/downloads/na36/wtexon/HuEx-1_0-st-v2.na36.hg19.probeset.csv.zip
# Save zip file in "data" and extract.
huex_annotation_path <- file.path("data", "HuEx-1_0-st-v2.na36.hg19.probeset.csv", "HuEx-1_0-st-v2.na36.hg19.probeset.csv")
huex_header_lines <- readLines(huex_annotation_path, n = 50)
huex_header_lines <- max(which(grepl("^#", huex_header_lines)))
huex_annotation <- readr::read_csv(
  huex_annotation_path, 
  col_types = cols_only(
    gene_assignment = col_character(), 
    probeset_id = col_integer()
  ), 
  skip = huex_header_lines
) %>% 
  tidy_gene_assignment(gene_assignment) %>% 
  tidyr::unnest() %>% 
  select(probeset_id, gene_assignment = gene_assignment_1)

# ---- Get Latest HUGO names ----
download_files(
  urls = c("hgnc_dict.tsv" = HUGO_DICT_URL),
  dest.dir = "data"
)

hgnc_dict <- readr::read_delim("data/hgnc_dict.tsv", 
    "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(`Approved symbol` = str_replace(`Approved symbol`, "~withdrawn", ""))

probe2hugo <- huex_annotation %>% 
  distinct() %>% 
  inner_join(hgnc_dict, by = c(gene_assignment = "RefSeq IDs"))

# # This part is not needed but kept for future reference. Uncomment if the
# # probset-hugo mapping contains un-approved symbols (check `Status` column)
# 
# probes_with_only_unapproved_names <- probe2hugo %>% 
#   group_by(ID, Status) %>%
#   count() %>%
#   tidyr::spread(Status, n, fill = 0) %>%
#   filter(Approved < 1, `Entry Withdrawn` > 0 | `Symbol Withdrawn` > 0)
# 
# # At this point, there are many HUGO names with Status = "Symbol withdrawn".
# # These entries point to a replacing symbol in the `Approved name` column
# # with the syntax "see XXX", e.g. "symbol withdrawn, see AGAP4" for "AGAP8".
# # So replace with the new name and re-merge with hgnc_dict.
# probe2hugo <- probe2hugo %>%
#   mutate(gene_assignment = if_else(
#     Status == "Symbol Withdrawn",
#     str_extract(`Approved name`, "(?<=see )(.+)$"),
#     gene_assignment
#   )) %>% 
#   distinct(probeset_id, gene_assignment) %>% 
#   inner_join(hgnc_dict, by = c(gene_assignment = "Approved symbol"))

gse_46691$exprs <- read_tsv_filtered(
  file.path(data_dir, "GSE46691_quantile_normalized.txt"),
  ID_REF %in% probe2hugo$probeset_id,
  chunk_size = 10000, 
  col_types = readr::cols(ID_REF = col_integer(), .default = col_double())
)

# ---- Save Checkpoint ----
dir.create("out")
saveRDS(gse_46691, file.path("out", "gse46691_hugo_checkpoint.rds"))
saveRDS(gse_46691$exprs, file.path("out", "gse46691_hugo_exprs.rds"))
saveRDS(probe2hugo, file.path("out", "probe2hugo.rds"))

## To restart from here without having to re-run the above:
# gse_46691 <- readRDS(file.path("out", "gse46691_hugo_checkpoint.rds"))
# probe2hugo <- readRDS(file.path("out", "probe2hugo.rds"))


# ---- Probe ID to HUGO ----
gse_46691_exprs_prepped <- gse_46691$exprs %>% 
  tidyr::gather("sample", "value", -ID_REF) %>% 
  mutate(sample = str_replace(sample, "\\.CEL$", "")) %>% 
  left_join(
    select(probe2hugo, probeset_id, hugo_name = `Approved symbol`), 
    ., 
    by = c(probeset_id = "ID_REF")
  ) %>% 
  group_by(sample, hugo_name) %>% 
  summarize(value = SUMMARIZE_FUNCTION(value)) %>% 
  tidyr::spread(hugo_name, value)

saveRDS(gse_46691_exprs_prepped, file.path("out", "gse466691_hugo_exprs_prepped.rds"))

write_tsv(gse_46691_exprs_prepped, file.path("out", "gse466691_hugo_exprs_prepped.tsv"))


# ---- Additional Patient Information ----

# Clean up pheno data and then merge with gse
gse_46691_pdata <- pData(gse_46691$phenoData) %>% 
  as_tibble() %>% 
  remove_common(description) %>% 
  clean_channel_vars() %>% 
  select(-matches("title|characteristics|supplementary")) %>% 
  mutate_all(as.character) %>% 
  readr::type_convert()

## To merge with expression dataset, run the following, which will prepend four
## columns before the exprs data: `geo_accession`, `sample`, `gleason score` and
## `metastatic event`.
#
# gse_46691_pdata %>% 
#   rename(sample = description) %>% 
#   left_join(gse_46691_exprs_prepped, by = "sample")
