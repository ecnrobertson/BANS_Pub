prep_genoscape_inputs <- function(Q_tibble, pca_df, all) {
  
  pca_df <- pca_df %>%
    dplyr::select(-V2)
  
  colnames(pca_df)[-1] <- paste0("PC", 1:(ncol(pca_df) - 1))
  names(pca_df)[1] <- "ind"
  
  Q_tibble$Sample <- pca_df$ind
  
  K <- ncol(Q_tibble) - 1
  
  ind_filtered <- all %>%
    dplyr::filter(BGP_ID %in% pca_df$ind) %>%
    dplyr::mutate(
      Region_site = paste(State, Lat, Long, sep = "_"),
      Region_pop = State
    )
  
  LatLong_tibble <- Q_tibble %>%
    dplyr::select(Sample) %>%
    dplyr::left_join(ind_filtered, by = c("Sample" = "BGP_ID")) %>%
    dplyr::select(Sample, Lat, Long)
  
  long_lat_tibble <- Q_tibble %>%
    dplyr::select(Sample) %>%
    dplyr::left_join(LatLong_tibble, by = "Sample") %>%
    dplyr::select(Long, Lat)
  
  long_lat_matrix <- as.matrix(long_lat_tibble)
  
  Q_matrix <- Q_tibble %>%
    dplyr::select(-Sample) %>%
    as.matrix()
  
  bad <- which(!is.finite(long_lat_matrix[,1]) | !is.finite(long_lat_matrix[,2]))
  
  return(list(
    Q_tibble = Q_tibble,
    pca_df = pca_df,
    ind_filtered = ind_filtered,
    LatLong_tibble = LatLong_tibble,
    long_lat_tibble = long_lat_tibble,
    long_lat_matrix = long_lat_matrix,
    Q_matrix = Q_matrix,
    K = K,
    bad = bad
  ))
}