#' Download or process GPL file and convert it to a featureData AnnotatedDataFrame
#' 
#' This code was extracted from [GEOquery::parseGSEMatrix](https://git.io/vhIYC).
#' 
#' @param gpl <GPL|chr|list> Either a GEOquery GPL object (a pre-loaded GPL), or
#'   a character string indicating the GPL that should be parsed. Also takes a
#'   list and applies `gpl2d` recursively to the list.
gpl2fd <- function(gpl) {
  if (inherits(gpl, "character")) gpl <- GEOquery:::parseGPL(gpl)
  if (inherits(gpl, "list")) return(purrr::map(gpl, gpl2fd))
  # At this point, throw an error if gpl is not a GPL object
  stopifnot(inherits(gpl, "GPL"))
  
  # This section from GEOquery
  gpl_meta <- Columns(gpl)
  gpl_dat  <- Table(gpl)
  
  # Match row names and column names between meta and data
  # - Add ID column of data (first col) to rownames of data (for AnnotatedDataFrame)
  # - Make sure colnames are unique for the data
  # - Match data columns to meta rows via meta rownames
  gpl_dat_ids <- if (ncol(gpl_dat)) as.character(gpl_dat[, 1]) else character()
  gpl_dat_ids <- ifelse(is.na(gpl_dat_ids), "NA", gpl_dat_ids)
  rownames(gpl_dat) <- make.unique(gpl_dat_ids)
  colnames(gpl_dat) <- make.unique(colnames(gpl_dat))
  rownames(gpl_meta) <- colnames(gpl_dat)
  
  # Return featureData as AnnotatedDataFrame
  new("AnnotatedDataFrame", data = gpl_dat, varMetadata = gpl_meta)
}

#' Build ExpressionSet or ExpressionSet-like list for GSE46691
#' 
#' Primarily used to build a list *like* the ExpressionSet for GSE46691.
#' Can also return an ExpressionSet, but you need lots of memory to load.
#' ExpressionSets get very mad when the rows and columns aren't all aligned
#' when you're building the ExpressionSet object. In this particular case, we're
#' loading the normalized probe data directly and would really like to be able
#' to load just the subset we care about. Loading the featureData is easy, but
#' subsetting is tough. For now, use the defaults and use the code below
#' to work with the data we want without the overhead of the ExpressionSet.
#' 
#' @param exprs_file <chr> Filename pointing to the normalized probe data.
#'   For GSE46691 this is `"GSE46691_quantile_normalized.txt.gz"`.
#' @param as_ExpressionSet <lgl> Try to build and return and ExpressionSet?
build_gse_46691 <- function(file_exprs = NULL, as_ExpressionSet = FALSE, data_dir = "data") {
  # Check params
  include_exprs <- !is.null(file_exprs)
  if (as_ExpressionSet && !include_exprs) {
    stop("Must load exprs to build ExpressionSet, retry by passing a file ",
         "name to `file_exprs`.")
  }
  
  # This part is easier if we're in the dir with the data
  if (!dir.exists(data_dir)) {
    stop("Directory ", data_dir, " does not exists. ", 
         'Use `getGEOfile("GSE46691", "', data_dir, '", amount = "full")`.')
  }
  owd <- setwd(data_dir)
  on.exit(setwd(owd))
  
  # Files we need
  file_gpl <- "GPL5188.soft"
  file_series_matrix <- "GSE46691_series_matrix.txt.gz"
  if (is.null(file_exprs)) file_exprs <- "GSE46691_quantile_normalized.txt.gz"
  stopifnot(file.exists(file_exprs))
  stopifnot(file.exists(file_gpl))
  stopifnot(file.exists(file_series_matrix))
  
  # I want to filter the features by the features included in the expression
  # matrix, but the expression matrix is HUGE.
  # So I used egrep to pull out the IDs in the first column of exprs file
  cli::cat_bullet("Extracting ", crayon::bold("feature IDs"), " (ID_REF) from ", 
                  crayon::green(file_exprs))
  file_exprs_ids <- paste0(file_exprs, "-ID_REF.txt")
  if (!file.exists(file_exprs_ids)) {
    grep_cmd <- if (grepl("gz$", file_exprs)) "zgrep -E" else "egrep"
    system(glue::glue('{grep_cmd} --only-matching "^\\d+\\b" {file_exprs} > {file_exprs_ids}'))
  }
  feature_ids <- scan(file_exprs_ids, character())

  
  cli::cat_bullet("Loading ", crayon::bold("feature data"), " from ", crayon::green(file_gpl))
  feature_data <- gpl2fd(file_gpl)
  feature_data <- feature_data[feature_ids, ] # limit to features in exprs
  
  if (as_ExpressionSet) {
    cli::cat_rule("Loading expression matrix and building expression set")
    new(
      "ExpressionSet",
      phenoData = phenoData(GEOquery:::parseGSEMatrix(file_series_matrix, getGPL=FALSE)$eset),
      annotation = "GPL5188",
      featureData = feature_data,
      exprs = as.matrix(read_tsv(
        file_exprs,
        col_types = cols(
          ID_REF = col_integer(),
          .default = col_double()
        )
      ))
    )
  } else {
    cli::cat_bullet("Loading series matrix", if (include_exprs) " and expression set")
    list(
      phenoData = phenoData(GEOquery:::parseGSEMatrix(file_series_matrix, getGPL=FALSE)$eset),
      featureData = feature_data,
      feature_ids = feature_ids,
      if (include_exprs) exprs = read_tsv(
        file_exprs, 
        col_types = cols(
          ID_REF = col_integer(), 
          .default = col_double()
        )
      )
    )
  }
}

#' Tidy up the gene assignment column
#' 
#' @param x <tbl|df> Feature Data (`fData(ExpressionSet)` or 
#'   `pData(featureData.AnnotateDataFrame)`)
#' @param col <bare> Unquoted column name containing the gene assignment info
#' @param sep_major <chr> Major delimiter between entries, default `"///"`
#' @param sep_minor <chr> Minor delimiter between entries, default `"//"`
tidy_gene_assignment <- function(x, col, sep_major = "///", sep_minor = "//") {
  col <- rlang::enquo(col)
  pull(x, !!col) %>% 
    str_match_all(glue::glue("\\s*(.*?)\\s*{sep_minor}\\s*(.+?)\\s*(?:{sep_major}|$)")) %>% 
    map(~ list(.[, 2], .[, 3])) %>% 
    transpose() %>% 
    set_names(., paste0(rlang::quo_text(col), "_", 1:length(.))) %>% 
    as_tibble() %>% 
    bind_cols(x, .)
}

#' Read TSV and Filter on the Fly
#' 
#' @param file <chr> Filename to read in
#' @param filter_expression <bare> The expression passed to `dplyr::filter()`
#' @inheritDotParams readr::read_tsv_chunked
read_tsv_filtered <- function(file, filter_expression, ...) {
  filter_expression <- rlang::enexpr(filter_expression)
  
  filter_by_id <- function(.data, pos) {
    dplyr::filter(.data, !!filter_expression)
  }
  
  readr::read_tsv_chunked(
    file, 
    readr::DataFrameCallback$new(filter_by_id),
    ...
  )
}

#' Median center a variable (respects groups)
#' 
#' @param x <grouped_df> A grouped tibble
#' @param value_var <bare> Unquote variable (col) containing value to center
median_center <- function(x, value_var) {
  stopifnot(inherits(x, "grouped_df"))
  value_var <- rlang::enquo(value_var)
  x_median <- dplyr::summarize(x, median = median(!!value_var, na.rm = TRUE))
  dplyr::left_join(x, x_median, by = dplyr::group_vars(x)) %>% 
    dplyr::mutate(!!rlang::quo_text(value_var) := !!value_var - median) %>% 
    dplyr::select(-median)
}