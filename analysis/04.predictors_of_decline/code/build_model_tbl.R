build_model_tbl <- function(base_dir,
                            pattern = "\\.rds$",
                            recursive = TRUE) {
  
  files <- list.files(
    path = base_dir,
    pattern = pattern,
    recursive = recursive,
    full.names = TRUE
  )
  
  parse_model_file <- function(fp) {
    
    parts <- strsplit(fp, .Platform$file.sep)[[1]]
    scheme <- parts[length(parts) - 1]
    
    fname <- basename(fp)
    fname <- sub("\\.rds$", "", fname)
    
    pieces <- strsplit(fname, "_")[[1]]
    
    model_type <- pieces[1]
    variable <- ifelse(length(pieces) > 1, pieces[2], NA)
    
    tibble::tibble(
      file = fp,
      scheme = scheme,
      model_type = model_type,
      variable = variable
    )
  }
  
  purrr::map_dfr(files, parse_model_file)
}