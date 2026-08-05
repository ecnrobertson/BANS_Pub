plot_modS_modGS_predictions_all_schemes_adaptive_with_null <- function(
    data,
    schemes,
    cols,
    CU_names_ref,
    null_model,
    model = c("modS", "modGS"),
    focal_vars = c("tmax_june_delta"),
    response_label = "Predicted Abundance Trend",
    predictor_data_dir = "../../04.predictors_of_decline/scratch/spatial_predictors",
    model_dir = "../../04.predictors_of_decline/outputs/GAM_outputs/function_out",
    color_dir = "../colors",
    out_dir = "../figures_output",
    covars = c("ppt_june_delta", "soil_june_delta", "tmax_june_delta", "tmin_june_delta", "ws_june_delta"),
    n_points = 100,
    plot_dims = list(width = 7, height = 4),
    x_labels = NULL,
    save_plots = TRUE
) {
  
  model <- match.arg(model)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out <- vector("list", length(schemes))
  names(out) <- schemes
  
  for (scheme in schemes) {
    
    clean.data <- data %>%
      dplyr::filter(!is.na(AU)) %>%
      dplyr::mutate(ESU = factor(AU))
    
    scheme_out_dir <- file.path(
      out_dir,
      paste0(model, "_predictions_with_null_overlay"),
      scheme
    )
    
    dir.create(scheme_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    scheme_results <- vector("list", length(focal_vars))
    names(scheme_results) <- focal_vars
    
    for (focal_var in focal_vars) {
      
      other_covars <- setdiff(covars, focal_var)
      var_stub <- sub("_june_delta$", "", focal_var)
      
      model_file <- file.path(
        model_dir,
        scheme,
        paste0(model, "_", var_stub, ".rds")
      )
      
      mod_obj <- readRDS(model_file)
      
      # -----------------------------
      # AU-level prediction grid
      # -----------------------------
      
      pred_grid <- expand.grid(
        focal_seq = seq(
          min(clean.data[[focal_var]], na.rm = TRUE),
          max(clean.data[[focal_var]], na.rm = TRUE),
          length.out = n_points
        ),
        ESU = levels(clean.data$ESU)
      )
      
      names(pred_grid)[1] <- focal_var
      
      for (v in other_covars) {
        if (v %in% names(clean.data)) {
          pred_grid[[v]] <- median(clean.data[[v]], na.rm = TRUE)
        }
      }
      
      pred_grid$ESU <- factor(pred_grid$ESU, levels = levels(clean.data$ESU))
      
      p_full <- predict(
        mod_obj,
        newdata = pred_grid,
        se.fit = TRUE,
        type = "link"
      )
      
      pred_grid$fit_full <- p_full$fit
      pred_grid$se_full  <- p_full$se.fit
      
      pred_grid_named <- pred_grid %>%
        dplyr::left_join(CU_names_ref, by = c("ESU" = "AU"))
      
      # -----------------------------
      # Null / range-wide prediction grid
      # -----------------------------
      
      null_grid <- data.frame(
        focal_seq = seq(
          min(clean.data[[focal_var]], na.rm = TRUE),
          max(clean.data[[focal_var]], na.rm = TRUE),
          length.out = n_points
        )
      )
      
      names(null_grid)[1] <- focal_var
      
      for (v in other_covars) {
        if (v %in% names(clean.data)) {
          null_grid[[v]] <- median(clean.data[[v]], na.rm = TRUE)
        }
      }
      
      p_null <- predict(
        null_model,
        newdata = null_grid,
        se.fit = TRUE,
        type = "link"
      )
      
      null_grid$fit_null <- p_null$fit
      null_grid$se_null  <- p_null$se.fit
      
      # -----------------------------
      # Axis label
      # -----------------------------
      
      x_lab <- if (is.null(x_labels)) {
        focal_var
      } else if (focal_var %in% names(x_labels)) {
        x_labels[[focal_var]]
      } else {
        focal_var
      }
      
      # -----------------------------
      # Plot
      # -----------------------------
      
      p <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = pred_grid_named,
          ggplot2::aes(
            x = .data[[focal_var]],
            ymin = fit_full - 2 * se_full,
            ymax = fit_full + 2 * se_full,
            fill = ESU,
            group = ESU
          ),
          alpha = 0.15,
          colour = NA
        ) +
        ggplot2::geom_line(
          data = pred_grid_named,
          ggplot2::aes(
            x = .data[[focal_var]],
            y = fit_full,
            color = ESU,
            group = ESU
          )
        ) +
        ggplot2::geom_ribbon(
          data = null_grid,
          ggplot2::aes(
            x = .data[[focal_var]],
            ymin = fit_null - 2 * se_null,
            ymax = fit_null + 2 * se_null
          ),
          fill = "grey50",
          alpha = 0.15,
          colour = NA
        ) +
        ggplot2::geom_line(
          data = null_grid,
          ggplot2::aes(
            x = .data[[focal_var]],
            y = fit_null
          ),
          color = "black",
          linewidth = 1,
          linetype = "dashed"
        ) +
        ggplot2::theme_classic() +
        ggplot2::guides(fill = "none") +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::xlab(x_lab) +
        ggplot2::ylab(response_label)
      
      if (save_plots) {
        ggplot2::ggsave(
          filename = file.path(
            scheme_out_dir,
            paste0(
              model, "_",
              focal_var,
              "_predictions_with_null_overlay_",
              scheme,
              ".png"
            )
          ),
          plot = p,
          width = plot_dims$width,
          height = plot_dims$height
        )
      }
      
      scheme_results[[focal_var]] <- list(
        pred_grid = pred_grid,
        pred_grid_named = pred_grid_named,
        null_grid = null_grid,
        plot = p
      )
    }
    
    out[[scheme]] <- list(
      clean.data = clean.data,
      results_by_var = scheme_results
    )
  }
  
  out
}