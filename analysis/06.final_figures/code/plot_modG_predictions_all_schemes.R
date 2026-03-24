plot_modG_predictions_all_schemes <- function(
    schemes,
    focal_vars = c("tmax_june_delta"),
    response_label = "Predicted Abundance Trend",
    predictor_data_dir = "../../04.predictors_of_decline/scratch/spatial_predictors",
    model_dir = "../../04.predictors_of_decline/outputs/GAM_outputs/function_out",
    color_dir = "../colors",
    out_dir = "../figures_output/modG_predictions",
    model_file = "modG.rds",
    covars = c("ppt_june_delta", "soil_june_delta", "tmax_june_delta", "tmin_june_delta", "ws_june_delta"),
    n_points = 100,
    plot_dims = list(width = 6, height = 4),
    x_labels = NULL,
    save_plots = TRUE,
    include_pop_line = FALSE
) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out <- vector("list", length(schemes))
  names(out) <- schemes
  
  for (scheme in schemes) {
    
    CU_names_ref <- read.csv(
      file.path(color_dir, paste0("CU_", scheme, "_name_ref.csv"))
    )
    
    cols <- readRDS(
      file.path(color_dir, paste0("CU_grouping_colors_2_", scheme, "_fullnames.rds"))
    )
    
    data <- read.csv(
      file.path(
        predictor_data_dir,
        paste0("spatial_predictor_data_xy_from_ESU_", scheme, "_threshold45.csv")
      )
    )
    
    clean.data <- data %>%
      dplyr::filter(!is.na(ESU)) %>%
      dplyr::mutate(ESU = factor(ESU))
    
    modG <- readRDS(
      file.path(model_dir, scheme, model_file)
    )
    
    scheme_out_dir <- file.path(out_dir, scheme)
    dir.create(scheme_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    scheme_results <- vector("list", length(focal_vars))
    names(scheme_results) <- focal_vars
    
    for (focal_var in focal_vars) {
      
      other_covars <- setdiff(covars, focal_var)
      
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
        pred_grid[[v]] <- median(clean.data[[v]], na.rm = TRUE)
      }
      
      pred_grid$ESU <- factor(pred_grid$ESU, levels = levels(clean.data$ESU))
      
      p1 <- predict(modG, newdata = pred_grid, se.fit = TRUE, type = "link")
      pred_grid$fit <- p1$fit
      pred_grid$se <- p1$se.fit
      
      p2 <- predict(
        modG,
        newdata = pred_grid,
        se.fit = TRUE,
        type = "link",
        exclude = "s(ESU)"
      )
      pred_grid$fit_pop <- p2$fit
      pred_grid$se_pop <- p2$se.fit
      
      pred_grid_named <- pred_grid %>%
        dplyr::left_join(CU_names_ref, by = "ESU")
      
      x_lab <- if (is.null(x_labels)) {
        focal_var
      } else if (focal_var %in% names(x_labels)) {
        x_labels[[focal_var]]
      } else {
        focal_var
      }
      
      p <- ggplot2::ggplot(
        pred_grid_named,
        ggplot2::aes(
          x = .data[[focal_var]],
          group = Name,
          color = Name
        )
      ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            ymin = fit - 2 * se,
            ymax = fit + 2 * se,
            fill = Name
          ),
          alpha = 0.15,
          colour = NA
        ) +
        ggplot2::geom_line(ggplot2::aes(y = fit)) +
        ggplot2::theme_classic() +
        ggplot2::guides(fill = "none") +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::xlab(x_lab) +
        ggplot2::ylab(response_label)
      
      if (include_pop_line) {
        p <- p + ggplot2::geom_line(
          ggplot2::aes(y = fit_pop),
          linetype = "dashed",
          color= "black"
        )
      }
      
      if (save_plots) {
        ggplot2::ggsave(
          filename = file.path(
            scheme_out_dir,
            paste0("modG_", focal_var, "_predictions_", scheme, ".png")
          ),
          plot = p,
          width = plot_dims$width,
          height = plot_dims$height
        )
      }
      
      scheme_results[[focal_var]] <- list(
        pred_grid = pred_grid,
        pred_grid_named = pred_grid_named,
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