## Replotting with an updated color/spatial ordering system
library(dplyr)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(readr)
library(data.table)

setwd("/Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/02.5.delineate_ESUs/")

################### spatial plotting
meta <- read.csv("../../data/BANS_all_sample_data_shifted.csv")
pop_list <- read.csv("pops.csv")

meta.sub <- left_join(meta, pop_list, by="BGP_ID")

cluster_geo <- meta.sub %>%
  group_by(Group) %>%
  summarise(
    cen_lon = mean(Long, na.rm = TRUE),
    cen_lat = mean(Lat, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(cen_lon, cen_lat) %>%   # west -> east, then south -> north
  mutate(cluster_letter = num_to_letters(row_number()))

# Save key
saveRDS(cluster_geo, "cluster_letter_key_geo_lonlat.rds")
write.csv(cluster_geo, "cluster_letter_key_geo_lonlat.csv", row.names = FALSE)

# Join back
meta.letter <- meta.sub %>%
  left_join(cluster_geo %>% select(Group, cluster_letter), by = "Group")

# Colors keyed by letters
K <- nrow(cluster_geo)
cols <- hcl.colors(21, palette = "Dark 3")
names(cols) <- LETTERS[1:21]

barplot(rep(1, 21), col = cols, border = NA)

cluster_colors_letters <- setNames(cols, cluster_geo$cluster_letter)
saveRDS(cluster_colors_letters, "cluster_colors_letters_geo_lonlat.rds")


# --- Plot (color + label are now LETTERS) ---

library(rnaturalearth)
library(raster)
library(terra)
map <- ne_states(country = c("United States of America", "Canada", "Mexico"), 
                 returnclass = 'sf')

p <- ggplot() +
  geom_sf(data = map, fill = "white", color = "black", linewidth = 0.1) +
  geom_point(
    data = cluster_geo,
    aes(x = cen_lon, y = cen_lat, color = cluster_letter),
    size = 1.5,
    position = position_jitter(width = 0.1, height = 0.1)
  ) +
  geom_label(
    data = cluster_geo,
    aes(x = cen_lon, y = cen_lat, label = cluster_letter),
    size = 2,
    vjust = -0.6
  ) +
  scale_color_manual(
    values = cluster_colors_letters,
    name = "Spatial clusters (PC1-ordered letters)",
    guide = "none"   # keep map clean; set to "legend" if you want a legend
  ) +
  coord_sf(xlim = c(-145, -55), ylim = c(15, 80)) +
  theme_minimal() +
  labs(
    title = "All clusters with ≥4 individuals (letter labels ordered by mean PC1)",
    x = "Longitude",
    y = "Latitude"
  )

p
ggsave("cluster_plot_letters_uninformed.png", p, width = 10, height = 6, dpi = 300)

################# TESTING WITH PCA ###################
plink_pca <- read.table("PCA_eigen/BANS.filtered.unrel.ld25-10-0.5.eigenvec", header = FALSE) %>% dplyr::select(-V2)
eigenval <- scan("PCA_eigen/BANS.filtered.unrel.ld25-10-0.5.eigenval")

pve05 <- data.frame(PC = 1:189, pve = eigenval/sum(eigenval)*100)
ggplot(pve05[1:20, ], aes(x = PC, y = pve)) +
  geom_col(width = 0.9) +
  labs(x = "PC", y = "PVE (%)") +
  theme_bw()

# set names
pca05 <- plink_pca
colnames(pca05)[-1] <- paste0("PC", 1:(ncol(pca05) - 1))
names(pca05)[1] <- "ind"
nrow(pca05)
#correct number of samples
pops <- meta.letter %>% filter(BGP_ID %in% pca05$ind) %>% dplyr::select(cluster_letter, BGP_ID)

# Join PCA data
pca.cluster <- left_join(pca05, pops %>% rename(ind=BGP_ID), by="ind")  # left_join PCA first to preserve PC columns

samp_unrel  <-  fread("sample_unrel_short.txt")
pca.cluster.norel <- pca.cluster %>% filter(ind %in% samp_unrel$IID)
nrow(pca.cluster.norel)

rel.nr <- ggplot(pca.cluster.norel, aes(
  x = PC1,
  y = PC2,
  col = cluster_letter,
  text = paste0(
    "ID: ", ind, "<br>",
    "Pop: ", " (", cluster_letter, ")<br>",
    "PC1: ", round(PC1, 4), "<br>",
    "PC2: ", round(PC2, 4)
  )
)) +
  geom_point(size = 3) +
  scale_color_manual(values = cluster_colors_letters) +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve05$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve05$pve[2], 3), "%)"))

ggplotly(rel.nr, tooltip = "text")

ggsave("plots/final_ESU/BANS_PCA_norel.png", width = 8.5, height = 7)
