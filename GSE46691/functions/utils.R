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

#' Aask if we should overwrite a directory
#' 
#' Assumes that directory to be created will be in the current working directory
#' and cannot make recursive directories. 
ask_if_overwrite <- function(destdir, create_dir = FALSE, ...) {
  stopifnot(!grepl(destdir, .Platform$file.sep, fixed = TRUE))
  library(glue)
  library(dplyr)
  if (dir.exists(destdir)) {
    overwrite_ok <- interactive() && 
      yesno::yesno2(crayon::red(cli::symbol$cross), glue(" Directory {destdir} already exists. Overwrite?"))
    if (overwrite_ok) { 
      unlink(destdir, recursive = TRUE)
    } else {
      root_dir <- gsub(glue("({.Platform$file.sep}.+)$"), "", destdir)
      existing_dirs <- file.info(dir(root_dir, full.names = TRUE)) %>% 
        tibble::rownames_to_column("name") %>% 
        filter(isdir, grepl(destdir, name, fixed = TRUE)) %>% 
        pull(name)
      dirs <- make.unique(c(existing_dirs, destdir))
      destdir <- dirs[length(dirs)]
    }
  }
  if (create_dir) dir.create(destdir, ...)
  destdir
}