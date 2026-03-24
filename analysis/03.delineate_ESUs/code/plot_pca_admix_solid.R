plot_pca_admix_solid <- function(plink_pca,
                                 eigenval,
                                 cols,
                                 pops,
                                 admix_groups,
                                 out_plot = NULL,
                                 width = 8.5,
                                 height = 7) {
  
  # format PCA table
  pca_df <- plink_pca %>%
    dplyr::select(-V2)
  
  colnames(pca_df)[-1] <- paste0("PC", 1:(ncol(pca_df) - 1))
  names(pca_df)[1] <- "ind"
  
  # percent variance explained
  pve <- data.frame(
    PC = seq_along(eigenval),
    pve = eigenval / sum(eigenval) * 100
  )
  
  # standardize pops input
  pops2 <- pops %>%
    dplyr::transmute(BGP_ID, Group)
  
  # standardize admix group table
  admix_groups2 <- admix_groups
  if ("ind" %in% names(admix_groups2) && !"BGP_ID" %in% names(admix_groups2)) {
    admix_groups2 <- admix_groups2 %>% dplyr::rename(BGP_ID = ind)
  }
  
  pops.groups <- dplyr::left_join(pops2, admix_groups2, by = "BGP_ID") %>%
    dplyr::rename(ind = BGP_ID)
  
  # join PCA + metadata
  pca.cluster <- dplyr::left_join(pca_df, pops.groups, by = "ind")
  
  # plot
  p <- ggplot2::ggplot(
    pca.cluster,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      col = admix_group,
      text = paste0(
        "ID: ", ind, "<br>",
        "Group: ", Group, "<br>",
        "PC1: ", round(PC1, 4), "<br>",
        "PC2: ", round(PC2, 4)
      )
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(
      values = cols,
      name = "ESU Group"
    ) +
    ggplot2::theme_light() +
    ggplot2::xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
    ggplot2::ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
  
  if (!is.null(out_plot)) {
    ggplot2::ggsave(out_plot, plot = p, width = width, height = height)
  }
  
  return(list(
    pca_df = pca_df,
    pve = pve,
    pops.groups = pops.groups,
    pca.cluster = pca.cluster,
    pca_plot = p
  ))
}