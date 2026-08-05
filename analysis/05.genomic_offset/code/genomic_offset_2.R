#### Function to predict genomic offset from a RDA model
genomic_offset <- function(RDA, K, env_pres, env_fut, range = NULL, method = "loadings", scale_env, center_env){
  
  # Mask with the range if supplied
  if(!is.null(range)){
    env_pres <- raster::mask(env_pres, range)
    env_fut <- raster::mask(env_fut, range)
  }
  
  # Formatting and scaling environmental rasters for projection
  var_env_proj_pres <- as.data.frame(scale(raster::rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]])[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))
  var_env_proj_fut <- as.data.frame(scale(raster::rasterToPoints(env_fut[[row.names(RDA$CCA$biplot)]])[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))
  
  # Predicting pixels genetic component based on the loadings of the variables
  if (method == "loadings") {
    
    # Pre-allocate lists
    Proj_pres   <- vector("list", K)
    Proj_fut    <- vector("list", K)
    Proj_offset <- vector("list", K)
    names(Proj_pres)   <- paste0("RDA", 1:K)
    names(Proj_fut)    <- paste0("RDA", 1:K)
    names(Proj_offset) <- paste0("RDA", 1:K)
    
    # Which predictors are used in the RDA biplot loadings
    vars <- rownames(RDA$CCA$biplot)
    
    # Ensure matrices have the same variable order as the biplot rows
    Xpres <- as.matrix(var_env_proj_pres[, vars, drop = FALSE])
    Xfut  <- as.matrix(var_env_proj_fut[,  vars, drop = FALSE])
    
    # Indices of cells that are not NA in the template layer
    template_pres <- env_pres[[1]]
    template_fut  <- env_fut[[1]]
    idx_pres <- which(!is.na(raster::values(template_pres)))
    idx_fut  <- which(!is.na(raster::values(template_fut)))
    
    # Loop over axes
    for (i in 1:K) {
      
      # Loadings for axis i (vector length = number of variables)
      b <- as.matrix(RDA$CCA$biplot[vars, i, drop = FALSE])
      
      # Project environments onto axis i (fast; BLAS-backed)
      vals_pres <- as.vector(Xpres %*% b)
      vals_fut  <- as.vector(Xfut  %*% b)
      
      # Current climates raster
      ras_pres <- template_pres
      ras_pres[] <- NA_real_
      ras_pres[idx_pres] <- vals_pres
      names(ras_pres) <- paste0("RDA_pres_", i)
      Proj_pres[[i]] <- ras_pres
      
      # Future climates raster
      ras_fut <- template_fut
      ras_fut[] <- NA_real_
      ras_fut[idx_fut] <- vals_fut
      names(ras_fut) <- paste0("RDA_fut_", i)
      Proj_fut[[i]] <- ras_fut
      
      # Single axis genetic offset
      Proj_offset[[i]] <- abs(ras_pres - ras_fut)
      names(Proj_offset[[i]]) <- paste0("RDA_offset_", i)
    }
  }
  
  
  # Predicting pixels genetic component based on predict.RDA
  if(method == "predict"){ 
    # Prediction with the RDA model and both set of envionments 
    pred_pres <- predict(RDA, var_env_proj_pres[,-c(1,2)], type = "lc")
    pred_fut <- predict(RDA, var_env_proj_fut[,-c(1,2)], type = "lc")
    # List format
    Proj_offset <- list()    
    Proj_pres <- list()
    Proj_fut <- list()
    for(i in 1:K){
      # Current climates
      ras_pres <- raster::rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_pres[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      # Future climates
      ras_fut <- raster::rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_fut[,i])), crs = crs(env_pres))
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      Proj_fut[[i]] <- ras_fut
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Weights based on axis eigen values
  weights <- RDA$CCA$eig/sum(RDA$CCA$eig)
  
  # Weighing the current and future adaptive indices based on the eigen values of the associated axes
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) raster::rasterToPoints(Proj_pres[[x]])[,-c(1,2)]))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x]*weights[x])))
  Proj_offset_fut <- do.call(cbind, lapply(1:K, function(x) raster::rasterToPoints(Proj_fut[[x]])[,-c(1,2)]))
  Proj_offset_fut <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut[,x]*weights[x])))
  
  # Predict a global genetic offset, incorporating the K first axes weighted by their eigen values
  ras <- Proj_offset[[1]]
  ras[!is.na(ras)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut[x,]), method = "euclidean")))
  names(ras) <- "Global_offset"
  Proj_offset_global <- ras
  
  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset 
  return(list(Proj_pres = Proj_pres, Proj_fut = Proj_fut, Proj_offset = Proj_offset, Proj_offset_global = Proj_offset_global, weights = weights[1:K]))
}
