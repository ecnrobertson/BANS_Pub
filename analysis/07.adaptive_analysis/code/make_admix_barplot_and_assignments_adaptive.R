make_admix_barplot_and_assignments <- function(K_target,
                                               cols,
                                               Q_tibble,
                                               threshold,
                                               pca_df,
                                               pops,
                                               group_label_order,
                                               geo_key,
                                               out_plot = NULL,
                                               out_assignments = NULL,
                                               plot_width = 10,
                                               plot_height = 2.2,
                                               plot_dpi = 300) {
  
  # make sure PCA IDs are formatted
  pca_df <- pca_df %>%
    dplyr::select(-V2)
  
  colnames(pca_df)[-1] <- paste0("PC", 1:(ncol(pca_df) - 1))
  names(pca_df)[1] <- "ind"
  
  # attach individual IDs to Q table
  Q_tibble$ind <- pca_df$ind
  
  # colors
  component_colors <- cols
  
  # cluster membership + labels
  pops <- pops %>%
    dplyr::transmute(ind = BGP_ID, Group)
  
  pops.cluster.letter <- dplyr::left_join(pops, geo_key, by = "Group")
  
  # attach grouping info
  Q_tibble2 <- Q_tibble %>%
    dplyr::left_join(pops.cluster.letter, by = "ind") %>%
    dplyr::filter(!is.na(Group_label))
  
  # identify ancestry columns to pivot
  AU_cols <- paste0("AU", seq_len(K_target))
  
  # long format for plotting
  q_long <- Q_tibble2 %>%
    dplyr::mutate(
      Group_label = factor(Group_label, levels = group_label_order)
    ) %>%
    dplyr::arrange(Group_label, ind) %>%
    dplyr::mutate(
      ind = factor(ind, levels = unique(ind))
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(AU_cols),
      names_to = "popGroup",
      values_to = "prob"
    ) %>%
    dplyr::mutate(
      group_label = factor(Group_label, levels = group_label_order),
      popGroup = factor(popGroup, levels = AU_cols)
    )
  
  # ADMIX barplot
  K_plot <- ggplot2::ggplot(
    q_long,
    ggplot2::aes(x = factor(ind), y = prob, fill = factor(popGroup))
  ) +
    ggplot2::geom_col(color = "gray", linewidth = 0.1) +
    ggplot2::facet_grid(~group_label, switch = "x", scales = "free", space = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "", y = "Ancestry", title = "") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 1)) +
    ggplot2::theme(
      panel.spacing.x = grid::unit(0.001, "lines"),
      axis.text.x = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = component_colors, guide = "none")
  
  # hard assignment table for PCA coloring / downstream use
  # admix_group <- q_long %>%
  #   dplyr::group_by(ind) %>%
  #   dplyr::slice_max(prob, n = 1, with_ties = FALSE) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::transmute(
  #     BGP_ID = ind,
  #     admix_group = popGroup
  #   )
  
  admix_group <- q_long %>%
    group_by(ind) %>%
    slice_max(prob, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      BGP_ID = ind,
      admix_group = if_else(prob >= threshold, as.character(popGroup), NA_character_)
    )
  
  # optional writing
  if (!is.null(out_plot)) {
    ggplot2::ggsave(
      filename = out_plot,
      plot = K_plot,
      width = plot_width,
      height = plot_height,
      units = "in",
      dpi = plot_dpi,
      bg = "white"
    )
  }
  
  if (!is.null(out_assignments)) {
    readr::write_tsv(admix_group, out_assignments)
  }
  
  return(list(
    Q_tibble = Q_tibble,
    Q_tibble2 = Q_tibble2,
    q_long = q_long,
    K_plot = K_plot,
    admix_group = admix_group
  ))
}