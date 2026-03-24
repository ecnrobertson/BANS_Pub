build_genoscape_surfaces <- function(Q_matrix,
                                     long_lat_matrix,
                                     breeding_na,
                                     cols,
                                     K,
                                     brick_names = NULL,
                                     out_tif = NULL,
                                     lamproj = "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") {
  
  keep <- is.finite(long_lat_matrix[,1]) & is.finite(long_lat_matrix[,2])
  
  Q_matrix_use <- Q_matrix[keep, , drop = FALSE]
  long_lat_matrix_use <- long_lat_matrix[keep, , drop = FALSE]
  
  genoscape_brick <- tess3r::tess3Q_map_rasters(
    x = Q_matrix_use,
    coord = long_lat_matrix_use,
    map.polygon = breeding_na,
    window = terra::ext(breeding_na)[1:4],
    resolution = c(1000,1000),
    col.palette = tess3r::CreatePalette(cols, K),
    method = "map.max",
    interpol = tess3r::FieldsKrigModel(10)
  )
  
  if (!is.null(brick_names)) {
    names(genoscape_brick) <- brick_names
  }
  
  if (!is.null(out_tif)) {
    terra::writeRaster(genoscape_brick, out_tif, overwrite = TRUE)
  }
  
  genoscape_rgba <- genoscapeRtools::qprob_rando_raster(
    TRB = genoscape_brick,
    cols = cols,
    alpha_scale = 3,
    abs_thresh = 0,
    alpha_exp = 1.55,
    alpha_chop_max = 230
  )
  
  terra::crs(genoscape_rgba) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  genoscape_spat <- terra::rast(genoscape_rgba)
  
  genoscape_projected <- terra::project(
    genoscape_spat,
    lamproj,
    method = "near"
  )
  
  return(list(
    genoscape_brick = genoscape_brick,
    genoscape_projected = genoscape_projected
  ))
}