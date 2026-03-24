compare_model_aic <- function(data_path,
                              model_list,
                              response = "abd_trend",
                              esu_col = "ESU",
                              k_vals = c(ppt_june_delta = 15,
                                         soil_june_delta = 20,
                                         tmax_june_delta = 15,
                                         tmin_june_delta = 20),
                              method = "REML") {
  
  # read + clean data
  data <- read.csv(data_path)
  
  clean.data <- data %>%
    dplyr::filter(!is.na(.data[[esu_col]]))
  
  clean.data[[esu_col]] <- factor(clean.data[[esu_col]])
  
  # build species-level formula
  species_formula <- stats::as.formula(
    paste0(
      response, " ~ ",
      "s(ppt_june_delta, k=", k_vals["ppt_june_delta"], ") + ",
      "s(soil_june_delta, k=", k_vals["soil_june_delta"], ") + ",
      "s(tmax_june_delta, k=", k_vals["tmax_june_delta"], ") + ",
      "s(tmin_june_delta, k=", k_vals["tmin_june_delta"], ")"
    )
  )
  
  mod_species <- mgcv::gam(
    formula = species_formula,
    data = clean.data,
    method = method
  )
  
  # combine all models for AIC comparison
  all_models <- c(list(species = mod_species), model_list)
  
  aics <- do.call(AIC, all_models)
  aics$model <- rownames(aics)
  rownames(aics) <- NULL
  
  aics <- aics %>%
    dplyr::select(model, dplyr::everything()) %>%
    dplyr::mutate(
      delta = AIC - min(AIC),
      delta_vs_species = AIC - AIC[model == "species"]
    )
  
  return(list(
    data = data,
    clean.data = clean.data,
    mod_species = mod_species,
    aics = aics
  ))
}