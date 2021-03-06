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

#' Ask if we should overwrite a directory
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

#' Download Files and Post Process
#' 
#' Download files from `urls` and alternatively post-process with the functions
#' provided in post_process. The `urls` should be a named vector of URLs where 
#' the name corresponds to the downloaded file name. The `post_processs`
#' argument accepts a list of functions, where the list entries are named to
#' match the downloaded file name. Each entry is itself a named list, where the
#' name provides the name of the file expected to be created as a result of
#' having applied the post processing. If the derivative name or the destination
#' file exist, no downloading occurrs.
#' 
#' @example 
#' download_files(
#'   urls = c("GSE46691_family.soft.gz" = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46691/soft/GSE46691_family.soft.gz"),
#'   post_process = list(
#'     "GSE46691_family.soft.gz" = list(
#'       "GSE46691_family.soft" = GEOquery::gunzip
#'     )
#'   ),
#'   dest.dir = tempdir()
#' )
download_files <- function(
  urls = NULL, 
  post_process = NULL, 
  dest.dir = getwd()
) {
  dir.create(dest.dir, showWarnings = FALSE)
  
  for (i in seq_along(urls)) {
    url <- urls[i]
    dest.file <- names(urls)[i]
    dest.alt  <- if (dest.file %in% names(post_process)) names(post_process[[dest.file]])
    skip_alt  <- !is.null(dest.alt) && file.exists(file.path(dest.dir, dest.alt))
    
    if (skip_alt) {
      cli::cat_bullet("Skipping '", crayon::green(dest.file), "' because ",
                      "'", crayon::green(dest.alt), "' is already available in ", 
                      crayon::yellow(dest.dir), bullet = "pointer")
    } else if (file.exists(file.path(dest.dir, dest.file))) {
      cli::cat_bullet("Skipping '", crayon::green(dest.file), "' because it's already downloaded in ", 
                      crayon::yellow(dest.dir), bullet = "pointer")
    } else {
      cli::cat_bullet("Downloading '", crayon::blue(dest.file), "' into ", 
                      crayon::yellow(dest.dir), bullet = "pointer")
      download.file(url, file.path(dest.dir, dest.file))
      cli::cat_bullet(crayon::green(dest.file), " download complete", bullet = "tick")
    }
    
    if (!is.null(dest.alt) && !dest.alt %in% dir(dest.dir)) {
      cli::cat_bullet("Post processing '", crayon::green(dest.file), "' to produce '",
                   crayon::yellow(dest.alt), "'")
      process_fun <- post_process[[dest.file]][[dest.alt]]
      owd <- setwd(dest.dir)
      process_fun(dest.file)
      setwd(owd)
    }
  }
  
  dest.dir
}