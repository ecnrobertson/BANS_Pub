plot_esu_data_vis <- function(clean.data,
                              out_prefix,
                              out_dir = ".",
                              esus_for_climate = NULL,
                              climate_var = "tmax_june_delta",
                              trend_var = "abd_trend",
                              esu_var = "ESU",
                              spatial_x = "x",
                              spatial_y = "y",
                              climate_vars = c("ppt_june_delta",
                                               "soil_june_delta",
                                               "tmax_june_delta",
                                               "tmin_june_delta",
                                               "ws_june_delta")) {
  
  # packages used:
  # dplyr, ggplot2, tidyr, broom
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  df <- clean.data
  
  if (is.null(esus_for_climate)) {
    esus_for_climate <- sort(unique(df[[esu_var]]))
  }
  
  df_climate <- df %>%
    dplyr::filter(.data[[esu_var]] %in% esus_for_climate)
  
  #-----------------------------
  # 1) Trend per ESU
  #-----------------------------
  p_trend_violin <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[esu_var]], y = .data[[trend_var]])
  ) +
    ggplot2::geom_violin(trim = FALSE, fill = "gray90") +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme_classic() +
    ggplot2::xlab(esu_var) +
    ggplot2::ylab(trend_var)
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, "_trend_by_ESU_violin.png")),
    plot = p_trend_violin
  )
  
  #-----------------------------
  # 2) Relationship with climate: linear
  #-----------------------------
  p_climate_linear <- ggplot2::ggplot(
    df_climate,
    ggplot2::aes(x = .data[[climate_var]], y = .data[[trend_var]])
  ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = TRUE) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", esu_var))) +
    ggplot2::theme_classic() +
    ggplot2::xlab(climate_var) +
    ggplot2::ylab(trend_var)
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, "_trend_by_", climate_var, "_linear.png")),
    plot = p_climate_linear
  )
  
  #-----------------------------
  # 3) Relationship with climate: smooth
  #-----------------------------
  p_climate_smooth <- ggplot2::ggplot(
    df_climate,
    ggplot2::aes(x = .data[[climate_var]], y = .data[[trend_var]])
  ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(se = TRUE) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", esu_var))) +
    ggplot2::theme_classic() +
    ggplot2::xlab(climate_var) +
    ggplot2::ylab(trend_var)
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, "_trend_by_", climate_var, "_smooth.png")),
    plot = p_climate_smooth
  )
  
  #-----------------------------
  # 4) Slopes by ESU
  #-----------------------------
  slope_formula <- stats::as.formula(paste(trend_var, "~", climate_var))
  
  slopes <- df %>%
    dplyr::group_by(.data[[esu_var]]) %>%
    dplyr::do(broom::tidy(stats::lm(slope_formula, data = .))) %>%
    dplyr::filter(term == climate_var)
  
  p_slopes <- ggplot2::ggplot(
    slopes,
    ggplot2::aes(x = .data[[esu_var]], y = estimate)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = estimate - std.error,
        ymax = estimate + std.error
      ),
      width = 0.2
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme_classic() +
    ggplot2::ylab(paste0(climate_var, " sensitivity (slope)")) +
    ggplot2::xlab(esu_var)
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, "_slope_by_ESU.png")),
    plot = p_slopes
  )
  
  #-----------------------------
  # 5) Overlap in climate variables
  #-----------------------------
  clim_long <- df %>%
    dplyr::select(dplyr::all_of(c(esu_var, climate_vars))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(climate_vars),
      names_to = "var",
      values_to = "value"
    ) %>%
    dplyr::filter(!is.na(value))
  
  p_climate_overlap <- ggplot2::ggplot(
    clim_long,
    ggplot2::aes(x = value, fill = factor(.data[[esu_var]]), alpha = 0.25)
  ) +
    ggplot2::geom_density(linewidth = 0.2) +
    ggplot2::facet_wrap(~ var, scales = "free", ncol = 2) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Value") +
    ggplot2::ylab("Density") +
    ggplot2::guides(alpha = "none")
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, "_climate_overlap.png")),
    plot = p_climate_overlap
  )
  
  #-----------------------------
  # 6) Which variables matter most for each ESU?
  #-----------------------------
  effect_formula <- stats::as.formula(
    paste(trend_var, "~", paste(climate_vars, collapse = " + "))
  )
  
  effects_df <- df %>%
    dplyr::group_by(.data[[esu_var]]) %>%
    dplyr::do(broom::tidy(stats::lm(effect_formula, data = .))) %>%
    dplyr::filter(term != "(Intercept)")
  
  p_variable_importance <- ggplot2::ggplot(
    effects_df,
    ggplot2::aes(x = term, y = .data[[esu_var]], fill = estimate)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(midpoint = 0) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Predictor") +
    ggplot2::ylab(esu_var)
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, "_variable_importance.png")),
    plot = p_variable_importance
  )
  
  #-----------------------------
  # 7) Trend over space
  #-----------------------------
  p_trend_space <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[spatial_x]],
      y = .data[[spatial_y]],
      color = .data[[trend_var]]
    )
  ) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", esu_var))) +
    ggplot2::scale_color_gradient2(midpoint = 0) +
    ggplot2::coord_equal() +
    ggplot2::theme_void()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, "_trend_over_space.png")),
    plot = p_trend_space
  )
  
  return(list(
    trend_violin = p_trend_violin,
    climate_linear = p_climate_linear,
    climate_smooth = p_climate_smooth,
    slopes_df = slopes,
    slopes_plot = p_slopes,
    climate_overlap = p_climate_overlap,
    effects_df = effects_df,
    variable_importance = p_variable_importance,
    trend_over_space = p_trend_space
  ))
}