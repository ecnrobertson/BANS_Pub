#rsync -avzP /Users/ericarobertson/Desktop/BANS_Pub/analysis/07.adaptive_analysis/scratch/admix_output/BANS_admix_grouping_adaptive_G3_2.tsv ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
#rsync -avzP /Users/ericarobertson/Desktop/BANS_Pub/analysis/06.final_figures/colors/AU_grouping_colors_K3.rds ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
library(dplyr)
library(ggplot2)
library(vegan)
library(tibble)
library(ggrastr)
library(readr)
library(tidyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
model_name <- args[1]

dir.create("plots/RDA", recursive = TRUE, showWarnings = FALSE)

model_files <- c(
  env_mem4 = "results/RDA_output/BANS.RDAresults.env_bio5.15.18.mem_4.RDS",
  env_only = "results/RDA_output/BANS.RDAresults.envonly_bio5.15.18.RDS",
  env_pc1  = "results/RDA_output/BANS.RDAresults.env_bio5.15.18.PC_1.RDS",
  env_mems = "results/RDA_output/BANS.RDAresults.env_bio5.15.18.mems_4.5.10.2.RDS",
  env_pc12 = "results/RDA_output/BANS.RDAresults.env_bio5.15.18.PC_1.2.RDS"
)

rda_obj <- readRDS(model_files[[model_name]])
print(model_files[[model_name]])

uncorselmem_ids <- readRDS("data/bans.mems.env.uncor_ids.rds")
pops <- read.csv("data/pops.csv")
geo_key <- readRDS("data/geo_informative_labels.RDS")
eig <- summary(rda_obj)$cont$importance
prop <- eig["Proportion Explained", grep("^RDA", colnames(eig))]

admix_groups <- data.table::fread("data/BANS_admix_grouping_adaptive_G3_2.tsv") %>%
  mutate(
    admix_group = tidyr::replace_na(admix_group, "No Assignment")
  )
cols <- readRDS("data/AU_grouping_colors_K3.rds")

cols["No Assignment"] <- "#BDBDBD"

print(cols)
print(table(admix_groups$admix_group))

site_scores <- as.data.frame(scores(rda_obj, choices = c(1, 2), display = "sites", scaling = 3)) %>%
  rownames_to_column("row") %>%
  rename(RDA1 = 2, RDA2 = 3) %>%
  mutate(BGP_ID = uncorselmem_ids$BGP_ID) %>%
  left_join(pops, by = "BGP_ID")

site_scores_admix <- site_scores %>% left_join(admix_groups, by = "BGP_ID")
  dplyr::mutate(
    admix_group = tidyr::replace_na(admix_group, "No Assignment")
  )
  
var_scores <- as.data.frame(scores(rda_obj, choices = c(1, 2), display = "bp", scaling = 3)) %>%
  rownames_to_column("variable") %>%
  rename(RDA1 = 2, RDA2 = 3)

rda_eigs <- vegan::eigenvals(rda_obj, model = "constrained")
rda_pct <- round(100 * rda_eigs / sum(rda_eigs), 1)
############ sample/site biplot ############

p_sites <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.8)) +
  geom_point(
    data = site_scores_admix,
    aes(x = RDA1, y = RDA2, color = admix_group),
    size = 3,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = cols,
    name = "AU"
  ) +
  geom_segment(
    data = var_scores,
    aes(x = 0, y = 0, xend = 10 * RDA1, yend = 10 * RDA2),
    arrow = arrow(length = unit(0.02, "npc")),
    linewidth = 0.5
  ) +
  geom_text(
    data = var_scores,
    aes(x = 3.25 * RDA1, y = 3.25 * RDA2, label = variable),
    fontface = "bold",
    size = 3
  ) +
  labs(
    title = paste("RDA sample/site biplot:", model_name),
    x = paste0("RDA1 (", rda_pct[1], "%)"),
    y = paste0("RDA2 (", rda_pct[2], "%)")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  paste0("plots/RDA/BANS.", model_name, ".sample_site_biplot_AU_colors.pdf"),
  p_sites,
  width = 8,
  height = 6,
  useDingbats = FALSE
)