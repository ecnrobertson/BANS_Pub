rasterize_genoscape_polygons <- function(genoscape_file,
                                         out_tif = NULL,
                                         cluster_col = "Cluster",
                                         plot_check = TRUE) {
  
  # read polygon file
  genoscape <- sf::st_read(genoscape_file, quiet = TRUE) %>%
    dplyr::mutate(
      !!cluster_col := as.factor(.data[[cluster_col]])
    )
  
  # quick polygon plot
  if (plot_check) {
    print(
      ggplot2::ggplot() +
        ggspatial::layer_spatial(genoscape, ggplot2::aes(fill = .data[[cluster_col]]))
    )
  }
  
  # rasterize selected cluster column
  genoscape_rast <- stars::st_rasterize(
    genoscape %>% dplyr::select(dplyr::all_of(cluster_col), geometry)
  )
  
  # optionally write to disk
  if (!is.null(out_tif)) {
    stars::write_stars(genoscape_rast, out_tif)
  }
  
  return(list(
    genoscape = genoscape,
    genoscape_rast = genoscape_rast
  ))
}