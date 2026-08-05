plot_pca_location_gradient <- function(plink_pca,
                                       eigenval,
                                       base_cols,
                                       pops,
                                       admix_groups,
                                       geo_key,
                                       desired_order = NULL,
                                       out_plot = NULL,
                                       out_group_label_cols_rds = NULL,
                                       width = 8.5,
                                       height = 7) {
  
  # -----------------------------
  # Format PCA table
  # -----------------------------
  pca_df <- plink_pca %>%
    dplyr::select(-V2)
  
  colnames(pca_df)[-1] <- paste0("PC", 1:(ncol(pca_df) - 1))
  names(pca_df)[1] <- "ind"
  
  pve <- data.frame(
    PC = seq_along(eigenval),
    pve = eigenval / sum(eigenval) * 100
  )
  
  # -----------------------------
  # Standardize inputs
  # -----------------------------
  pops2 <- pops %>%
    dplyr::transmute(BGP_ID, Group)
  
  admix_groups2 <- admix_groups
  
  if ("ind" %in% names(admix_groups2) && !"BGP_ID" %in% names(admix_groups2)) {
    admix_groups2 <- admix_groups2 %>%
      dplyr::rename(BGP_ID = ind)
  }
  
  admix_groups2 <- admix_groups2 %>%
    dplyr::mutate(
      admix_group = tidyr::replace_na(admix_group, "No Assignment")
    )
  
  geo_key2 <- geo_key
  
  if ("Pop" %in% names(geo_key2) && !"Group" %in% names(geo_key2)) {
    geo_key2 <- geo_key2 %>%
      dplyr::rename(Group = Pop)
  }
  
  stopifnot(all(c("Group", "Group_label") %in% names(geo_key2)))
  
  geo_key2 <- geo_key2 %>%
    dplyr::mutate(
      Group_label = trimws(as.character(Group_label))
    )
  
  # -----------------------------
  # Join metadata
  # -----------------------------
  pops.groups.labels <- pops2 %>%
    dplyr::left_join(geo_key2, by = "Group") %>%
    dplyr::left_join(admix_groups2, by = "BGP_ID") %>%
    dplyr::mutate(
      admix_group = tidyr::replace_na(admix_group, "No Assignment"),
      Group_label = trimws(as.character(Group_label))
    )
  
  # -----------------------------
  # Determine owner ESU per location
  # -----------------------------
  base_cols <- base_cols
  base_cols["No Assignment"] <- "#BDBDBD"
  
  label_owner <- pops.groups.labels %>%
    dplyr::filter(!is.na(Group_label)) %>%
    dplyr::count(Group_label, admix_group, name = "n") %>%
    dplyr::group_by(Group_label) %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(owner_esu = admix_group) %>%
    dplyr::mutate(
      owner_esu = tidyr::replace_na(owner_esu, "No Assignment")
    )
  
  missing_cols <- setdiff(unique(label_owner$owner_esu), names(base_cols))
  if (length(missing_cols) > 0) {
    stop(
      "Missing colors for owner_esu values: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # -----------------------------
  # Build Group_label gradient colors
  # -----------------------------
  group_label_cols_df <- label_owner %>%
    dplyr::arrange(owner_esu, Group_label) %>%
    dplyr::group_by(owner_esu) %>%
    dplyr::group_modify(~{
      base <- base_cols[[as.character(.y$owner_esu)]]
      L <- .x$Group_label
      n <- length(L)
      
      shades <- grDevices::colorRampPalette(c("white", base))(n + 1)[-1]
      
      tibble::tibble(
        Group_label = L,
        col = shades
      )
    }) %>%
    dplyr::ungroup()
  
  group_label_cols <- stats::setNames(
    group_label_cols_df$col,
    group_label_cols_df$Group_label
  )
  
  stopifnot(!any(duplicated(names(group_label_cols))))
  
  if (!is.null(out_group_label_cols_rds)) {
    saveRDS(group_label_cols, file = out_group_label_cols_rds)
  }
  
  # -----------------------------
  # Add owner ESU to individuals
  # -----------------------------
  pops.groups.labels <- pops.groups.labels %>%
    dplyr::left_join(
      label_owner %>% dplyr::select(Group_label, owner_esu),
      by = "Group_label"
    )
  
  pops_for_join <- pops.groups.labels %>%
    dplyr::rename(ind = BGP_ID)
  
  # -----------------------------
  # Join PCA + metadata
  # -----------------------------
  pca.cluster <- dplyr::left_join(pca_df, pops_for_join, by = "ind") %>%
    dplyr::mutate(
      admix_group = tidyr::replace_na(admix_group, "No Assignment"),
      owner_esu = tidyr::replace_na(owner_esu, "No Assignment"),
      Group_label = trimws(as.character(Group_label))
    )
  
  # -----------------------------
  # Order legend / colors
  # -----------------------------
  if (!is.null(desired_order)) {
    desired_order <- as.character(desired_order)
    
    missing_from_palette <- setdiff(desired_order, names(group_label_cols))
    if (length(missing_from_palette) > 0) {
      warning(
        "These labels are in desired_order but not in group_label_cols: ",
        paste(missing_from_palette, collapse = ", ")
      )
    }
    
    labels_not_ordered <- setdiff(names(group_label_cols), desired_order)
    
    legend_order <- c(
      desired_order[desired_order %in% names(group_label_cols)],
      sort(labels_not_ordered)
    )
  } else {
    legend_order <- sort(names(group_label_cols))
  }
  
  group_label_cols <- group_label_cols[legend_order]
  
  pca.cluster <- pca.cluster %>%
    dplyr::mutate(
      Group_label = factor(Group_label, levels = legend_order)
    )
  
  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot2::ggplot(
    pca.cluster,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      col = Group_label,
      text = paste0(
        "ID: ", ind, "<br>",
        "Group: ", Group, "<br>",
        "Location: ", Group_label, "<br>",
        "Individual ESU: ", admix_group, "<br>",
        "Location owner ESU: ", owner_esu, "<br>",
        "PC1: ", round(PC1, 4), "<br>",
        "PC2: ", round(PC2, 4)
      )
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(
      values = group_label_cols,
      limits = legend_order,
      drop = FALSE,
      name = "Sampling location\nColored by ESU"
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
    ggplot2::ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
  
  if (!is.null(out_plot)) {
    ggplot2::ggsave(out_plot, plot = p, width = width, height = height)
  }
  
  return(list(
    pca_df = pca_df,
    pve = pve,
    pops.groups.labels = pops.groups.labels,
    label_owner = label_owner,
    group_label_cols = group_label_cols,
    pca.cluster = pca.cluster,
    pca_plot = p
  ))
}