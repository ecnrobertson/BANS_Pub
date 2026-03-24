plot_partial_effects_by_scheme <- function(
    schemes,
    model,
    clim_vars,
    color_dir = "../colors",
    predictor_data_dir = "../../04.predictors_of_decline/scratch/spatial_predictors",
    model_dir = "../../04.predictors_of_decline/outputs/GAM_outputs/function_out",
    out_dir = "../figures_output/partial_effects",
    width = 6,
    height = 2
) {
  
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
      dplyr::filter(!is.na(ESU))
    
    clean.data.names <- dplyr::left_join(clean.data, CU_names_ref, by = "ESU")
    
    scheme_out_dir <- file.path(out_dir, scheme)
    dir.create(scheme_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (var in clim_vars) {
      
      model_file <- file.path(
        model_dir,
        scheme,
        paste0(model, "_", var, ".rds")
      )
      
      mod_obj <- readRDS(model_file)
      
      sm_dat <- gratia::smooth_estimates(
        mod_obj,
        smooth = paste0("s(", var, "_june_delta,ESU)")
      )
      
      sm_dat.names <- dplyr::left_join(sm_dat, CU_names_ref, by = "ESU")
      
      p <- ggplot2::ggplot(
        sm_dat.names,
        ggplot2::aes(
          x = .data[[paste0(var, "_june_delta")]],
          y = .estimate,
          colour = Name
        )
      ) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~ Name, nrow = 1) +
        ggplot2::scale_colour_manual(values = cols) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          strip.text = ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          legend.position = "none",
          axis.title.y = ggplot2::element_text(size = 7),
          axis.text.x = ggplot2::element_text(
            size = 8,
            angle = 90,
            vjust = 0.5,
            hjust = 1
          )
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
        ggplot2::xlab(paste0("Change in ", var)) +
        ggplot2::ylab(expression("Partial effect on " ~ Delta ~ "Abundance"))
      
      ggplot2::ggsave(
        filename = file.path(
          scheme_out_dir,
          paste0("BANS_", var, "_", model, "_", scheme, "_partialeffect.png")
        ),
        plot = p,
        height = height,
        width = width
      )
    }
  }
}