plot_pca_location_gradient <- function(plink_pca,
                                       eigenval,
                                       cols,
                                       pops,
                                       admix_groups,
                                       letters,
                                       out_plot = NULL,
                                       out_letter_cols_rds = NULL,
                                       width = 8.5,
                                       height = 7) {
  
  #-----------------------------
  # standardize inputs
  #-----------------------------
  pca_df <- plink_pca %>%
    dplyr::select(-V2)
  
  colnames(pca_df)[-1] <- paste0("PC", 1:(ncol(pca_df) - 1))
  names(pca_df)[1] <- "ind"
  
  pve <- data.frame(
    PC = seq_along(eigenval),
    pve = eigenval / sum(eigenval) * 100
  )
  
  pops2 <- pops %>%
    dplyr::transmute(BGP_ID, Group)
  
  admix_groups2 <- admix_groups
  if ("ind" %in% names(admix_groups2) && !"BGP_ID" %in% names(admix_groups2)) {
    admix_groups2 <- admix_groups2 %>% dplyr::rename(BGP_ID = ind)
  }
  
  pops.groups <- dplyr::left_join(pops2, admix_groups2, by = "BGP_ID")
  
  # letters table may use Pop instead of Group
  letters2 <- letters
  if ("Pop" %in% names(letters2) && !"Group" %in% names(letters2)) {
    letters2 <- letters2 %>% dplyr::rename(Group = Pop)
  }
  
  stopifnot(all(c("Group", "cluster_letter") %in% names(letters2)))
  
  #-----------------------------
  # join letters
  #-----------------------------
  pops.groups.letters <- pops.groups %>%
    dplyr::left_join(letters2, by = "Group") %>%
    dplyr::mutate(cluster_letter = trimws(as.character(cluster_letter)))
  
  ESUs <- names(cols)
  
  # determine majority-owner ESU for each cluster letter
  letter_owner <- pops.groups.letters %>%
    dplyr::filter(
      !is.na(cluster_letter),
      !is.na(admix_group),
      admix_group %in% ESUs
    ) %>%
    dplyr::count(cluster_letter, admix_group, name = "n") %>%
    dplyr::group_by(cluster_letter) %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(owner_esu = admix_group)
  
  # build letter-level gradient palette within each owner ESU
  letter_cols_df <- letter_owner %>%
    dplyr::arrange(owner_esu, cluster_letter) %>%
    dplyr::group_by(owner_esu) %>%
    dplyr::group_modify(~{
      base <- cols[[.y$owner_esu]]
      L <- .x$cluster_letter
      n <- length(L)
      shades <- grDevices::colorRampPalette(c("white", base))(n + 1)[-1]
      tibble::tibble(cluster_letter = L, col = shades)
    }) %>%
    dplyr::ungroup()
  
  letter_cols <- stats::setNames(letter_cols_df$col, letter_cols_df$cluster_letter)
  
  stopifnot(!any(duplicated(names(letter_cols))))
  
  if (!is.null(out_letter_cols_rds)) {
    saveRDS(letter_cols, file = out_letter_cols_rds)
  }
  
  # add owner ESU back to individuals
  pops.groups.letters <- pops.groups.letters %>%
    dplyr::left_join(
      letter_owner %>% dplyr::select(cluster_letter, owner_esu),
      by = "cluster_letter"
    )
  
  # standardize individual ID name
  pops_for_join <- pops.groups.letters
  if ("BGP_ID" %in% names(pops_for_join) && !"ind" %in% names(pops_for_join)) {
    pops_for_join <- pops_for_join %>% dplyr::rename(ind = BGP_ID)
  }
  if ("sampleID" %in% names(pops_for_join) && !"ind" %in% names(pops_for_join)) {
    pops_for_join <- pops_for_join %>% dplyr::rename(ind = sampleID)
  }
  stopifnot("ind" %in% names(pops_for_join))
  
  # join PCA + metadata
  pca.cluster <- dplyr::left_join(pca_df, pops_for_join, by = "ind") %>%
    dplyr::mutate(cluster_letter = trimws(as.character(cluster_letter)))
  
  # lock palette mapping
  letter_cols <- letter_cols[order(names(letter_cols))]
  
  # plot
  p <- ggplot2::ggplot(
    pca.cluster,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      col = cluster_letter,
      text = paste0(
        "ID: ", ind, "<br>",
        "Group: ", Group, "<br>",
        "Individual ESU: ", admix_group, "<br>",
        "Letter owner ESU: ", owner_esu, "<br>",
        "PC1: ", round(PC1, 4), "<br>",
        "PC2: ", round(PC2, 4)
      )
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(
      values = letter_cols,
      limits = names(letter_cols),
      drop = FALSE,
      name = "Sampling location\nColored by ESU"
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
    pops.groups.letters = pops.groups.letters,
    letter_owner = letter_owner,
    letter_cols = letter_cols,
    pca.cluster = pca.cluster,
    pca_plot = p
  ))
}