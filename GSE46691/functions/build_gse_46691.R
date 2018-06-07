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
#'   For GSE46691 this is `"GSE46691_quantile_normalized.txt.gz"` or 
#'   `"GSE46691_quantile_normalized.txt"`.
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
  if (is.null(file_exprs)) {
    file_exprs <- c("GSE46691_quantile_normalized.txt.gz",
                    "GSE46691_quantile_normalized.txt")
    file_exprs <- file_exprs[file.exists(file_exprs)][1]
  }
  stopifnot(is.na(file_exprs) || file.exists(file_exprs))
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

#' Gather GSE46691 Data
#' 
#' Download GSE46691 source data as needed to the provided `data_dir`.
gather_gse46691 <- function(data_dir = "data") {
  cli::cat_rule("Gathering GSE46691 Data Files")
  gse46691_urls <- c(
    "GSE46691_quantile_normalized.txt.gz" = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46691&format=file&file=GSE46691%5Fquantile%5Fnormalized%2Etxt%2Egz",
    "GSE46691_family.soft.gz" = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46691/soft/GSE46691_family.soft.gz",
    "GSE46691_series_matrix.txt.gz" = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46691/matrix/GSE46691_series_matrix.txt.gz",
    "GPL5188.soft" = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL5188&form=text&view=full"
  )
  
  dir.create(data_dir, showWarnings = FALSE)
  cat(names(gse46691_urls), 
      "GSE46691_quantile_normalized.txt", 
      "GSE46691.soft",
      "GSE46691_quantile_normalized.txt.gz-ID_REF.txt",
      file = file.path(data_dir, ".gitignore"), sep = "\n")
  
  download_files(
    urls = gse46691_urls,
    post_process = list(
      "GSE46691_quantile_normalized.txt.gz" = list(
        "GSE46691_quantile_normalized.txt" = GEOquery::gunzip
      )
    ),
    dest.dir = data_dir
  )
  
  cli::cat_bullet("GSE46691 data files are available in ", 
                  crayon::yellow(file.path(getwd(), data_dir)), bullet = "tick")
  data_dir
}
