# #!/bin/bash
# #SBATCH --job-name=RDA
# #SBATCH --output=RDA.%j.out
# #SBATCH --error=RDA.%j.err
# #SBATCH -t 8:00:00
# #SBATCH -p amilan
# #SBATCH --qos=normalm
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node 24
# #SBATCH --mem=90G
# #SBATCH --mail-type=ALL
# #SBATCH  --mail-user=ericacnr@colostate.edu
# 
# source ~/.bashrc
#
# ############ Load and Run RDA.R script ############
# conda activate R
# # Run your R script
# Rscript RDA_adaptive.R

############################# RDA.R ##################################
##### File Checklist #####
# rsync -avzP /Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/03.GEA/BANS_topRDA_vars1mil_narrow.csv ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
# rsync -avzP /Users/ericarobertson/Desktop/BANS_adaptive_units/data/BANS_all_sample_data.csv ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
# rsync -avzP /Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/02.5.delineate_ESUs/bg_colors.rds ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
# In data/
# BANS_sample_list_pop.tsv
# BANS.all.imputed.pruned_noheader.raw
# pruned_allindv.fam
# BANS_topRDA_vars1mil_narrow.csv
# BANS_all_sample_data.csv
# bg_colors.rds
# sample_unrel_short.txt
# env_narrow.rds
# env2_narrow.rds

############ Libraries ############ 
library(dplyr)       # select, rename, mutate, pipes (%>%)
library(readr)       # read_delim, read_csv
library(data.table)  # fread
library(sf)              # st_as_sf
library(sp)              # as(..., "Spatial")
library(geosphere)       # distm, distHaversine
library(adespatial)      # dbmem
library(rnaturalearth)   # ne_states
library(adegenet)   # read.PLINK, genlight, locNames, indNames
library(vegan)      # rda, ordiR2step, RsquareAdj, scores, anova.cca, eigenvals
library(ggplot2)    # ggplot
library(RColorBrewer) # brewer.pal

########### MAKING ENV #############
# indorder <- fread("BANS.unrel_ind.imputed4.1.fam") %>% rename(BGP_ID = V1) %>% dplyr::select(BGP_ID)
# env <- fread("BANS_topRDA_vars1mil_narrow.csv")
# env <- left_join(indorder, env, by = "BGP_ID")
# saveRDS(env, "env_narrow.rds")
# 
# env2 <- env %>% dplyr::select(-Long, -Lat)
# saveRDS(env2, "env2_narrow.rds")

############ PCA for pop structure variables ############ 
plink_pca <- read.table("../02.5.neutral_pop_struc/BANS.filtered.unrel.ld25-10-0.5.eigenvec", header = F) %>% select(-V1)
n_pcs <- ncol(plink_pca) - 2  # number of PC columns
colnames(plink_pca) <- c("FID", "IID", paste0("PC", 1:n_pcs))
plink_pca <- plink_pca %>% dplyr::select(-FID) %>% rename(BGP_ID = IID)
nrow(plink_pca)

print("read in PCA variables")

############ MEMS for spatial autocorrelation ############ 
unrel_ind <- fread("data/sample_unrel_short.txt", header=F) %>% rename(IID=V1)
indorder <- fread("BANS.unrel_ind.imputed4.1.fam") %>% rename(BGP_ID = V1) %>% dplyr::select(BGP_ID)
coords <- read_csv("data/BANS_all_sample_data.csv") %>%
  filter(BGP_ID %in% indorder$BGP_ID) %>% 
  dplyr::select(BGP_ID, Long, Lat) %>% filter(BGP_ID %in% unrel_ind$IID)
nrow(coords)

# reorder coordinates to match PCA data
print("reorder coordinates to match PCA data")
coords <- coords[(match(indorder$BGP_ID, coords$BGP_ID)),]

coords_sf <- st_as_sf(coords, coords = c("Long", "Lat"), crs = 4326)

gd <- geosphere::distm(as(coords_sf, "Spatial"), fun = geosphere::distHaversine)
dist <- as.dist(gd)

# generate moran's eigenvector maps (MEMs) based on coordinates
print("generate MEMS based on coordinates")
mem <- adespatial::dbmem(dist, MEM.autocor = "all")

#moran_test <- moran.randtest(mem)

# rda to assess how much MEMs explain genetic variation (PC1–PC2)
mem_rda <- vegan::rda(plink_pca[,c("PC1", "PC2")] ~ ., as.data.frame(mem[,1:10]))
summary(mem_rda)

# create a null (intercept-only) rda model as a baseline
mem_rda_full <- vegan::rda(plink_pca[,c("PC1", "PC2")] ~ 1, as.data.frame(mem[,1:10]))
summary(mem_rda_full)

# stepwise forward selection of mem variables based on adjusted R-square and significance
selmem <- vegan::ordiR2step(mem_rda_full, mem_rda, R2permutations=1000,Pin=0.01,R2scope=T)
summary(selmem)
# extract the names of the selected mem variables
selmem <- names(selmem$terminfo$ordered)

# subset the mem matrix to include only the mems that explain population structure
mem_popstr <- mem[,selmem]

saveRDS(mem_popstr, file = "data/BANS.mems.popstr_unrel.rds")

coords_mems_pc <- cbind(coords, mem, plink_pca[,c("PC1", "PC2")])

print("MEMS identified correctly")
############ RDA models prep ############ 

plink_data <- read.PLINK("BANS.unrel_ind.imputed4.1.raw", )
print(nrow(plink_data))

# Convert genlight to genotype matrix
geno_matrix <- as.matrix(plink_data)

# Keep SNP names
snp_names <- locNames(plink_data)  

# Recreate genlight object with SNP names
plink_subset <- new("genlight", geno_matrix, loc.names = snp_names, ind.names = indNames(plink_data))
print(nrow(plink_subset))

# Move to Dataframe format
geno_df <- as.data.frame(plink_subset)
rownames(geno_df) <- indNames(plink_data)
print(nrow(geno_df)==nrow(plink_subset))

# #pop <- read.table(file = "data/BANS_sample_list_pop.tsv", sep="\t", header = FALSE) %>%
# #  rename(BGP_ID = V1, Pop = V3) %>% dplyr::select(BGP_ID, Pop) %>% filter(BGP_ID %in% unrel_ind$IID)
# 
# env <- fread("data/BANS_topRDA_vars1mil_narrow.csv")
# env <- left_join(indorder, env, by = "BGP_ID")
# 
# print(env$BGP_ID==rownames(geno_df))
# bgp_ids <- data.frame(BGP_ID = rownames(geno_df))

env <- readRDS("data/env_narrow.rds")

#env <- cbind(pop, env)
#env2 <- env %>% dplyr::select(-Long, -Lat)
env2 <- readRDS("data/env2_narrow.rds")
nrow(env2)

# Convert BGP_ID to character (to avoid factor-related issues)
env2[, BGP_ID := as.character(BGP_ID)]

# Convert ClimGroup to a factor if it's categorical
env2[, Pop := as.factor(Pop)]

pred <- dplyr::select(env2, where(is.numeric))

pred_mems <- cbind(pred, as.data.frame(mem_popstr))

pred_pc_mems_env <- cbind(pred_mems, plink_pca %>% dplyr::select(PC1, PC2))


print("RDA data loaded")
############ Actual RDA Models ############ 

### 1. Fit RDAs ----
print("fit RDAs")
cond_term <- paste(selmem, collapse = " + ")

form <- as.formula(paste(
  "geno_df ~ Temperate.or.Subpolar.Shrubland + bio14 + bio09 + Wetland + clay + bio06 +",
  "Condition(", cond_term, ")"
))
rda_env_mems_uncor <- rda(form, data = pred_pc_mems_env, scale = TRUE)
saveRDS(rda_env_mems_uncor, file = "results/RDA_output/rda_env_mems_uncor1mil_narrow.rds")
print("rda_env_mems_uncor")

# Control for PCs
rda_env_pc1 <- rda(geno_df ~ Temperate.or.Subpolar.Shrubland + bio17 + Wetland +
                     clay + bio08 + Water +
                     Condition(PC1),
                   data = pred_pc_mems_env, scale = TRUE)
saveRDS(rda_env_pc1, file = "results/RDA_output/rda_PC1_1mil_narrow.rds")
rda_env_pc12 <- rda(geno_df ~ Temperate.or.Subpolar.Shrubland + bio17 + Wetland +
                      clay + bio08 + Water + 
                      Condition(PC1 + PC2),
                    data = pred_pc_mems_env, scale = TRUE)
saveRDS(rda_env_pc12, file = "results/RDA_output/rda_PC12_1mil_narrow.rds")

############ Getting Outlier SNPs ############
load.rda <- scores(rda_env_mems_uncor, choices=c(1:3), display="species")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)   #find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]             #locus names in these tails
}

cand1 <- outliers(load.rda[,1],3)
cand2 <- outliers(load.rda[,2],3)
cand3 <- outliers(load.rda[,3],3)

##Find total # candidates
ncand <- length(cand1) + length(cand2) + length(cand3)
print("number of candidate snps")
print(ncand)

df.cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
df.cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
df.cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(df.cand1) <- colnames(df.cand2) <- colnames(df.cand3) <- c("axis","snp","loading")

df.cand <- rbind(df.cand1, df.cand2, df.cand3)
df.cand$snp <- as.character(df.cand$snp)

df.cand$snp <- gsub("-", ".", df.cand$snp)

##Add environmental correlations to candidate snps
foo <- matrix(nrow=(ncand), ncol=ncol(pred))  #4 columns for 4 predictors
colnames(foo) <- colnames(pred)

for (i in 1:length(df.cand$snp)) {
  nam <- df.cand[i,2]
  snp.gen <- geno_df[[nam]]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(df.cand,foo)

length(cand$snp[duplicated(cand$snp)])
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) # 18202 duplicates on axis 2
table(foo[foo[,1]==3,2]) # 18787 duplicates on axis 3
cand <- cand[!duplicated(cand$snp),] # remove duplicate detection

cols <- as.numeric(ncol(cand))
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,(cols+1)] <- names(which.max(abs(bar[4:cols]))) # gives the variable
  cand[i,(cols+2)] <- max(abs(bar[4:cols]))              # gives the correlation
}

colnames(cand)[cols+1] <- "predictor"
colnames(cand)[cols+2] <- "correlation"

table(cand$predictor)

#rda_env_mems_uncor
table(cand$predictor)
write.table(cand, "env_only_RDA_cand.snps.txt", row.names=F,sep = "\t", quote=F)