############################# RDA_env_vars_cluster.R ##################################
##### File Checklist #####
# rsync -avzP /Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/03.GEA/Bans_env_pop.txt ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
# In data/
# Bans_env_pop.txt
# subsampled_noheader.raw

setwd("../BANS_Pub/analysis/05.genomic_offset/")
##### Libraries #####
library(vegan)
library(tidyverse)
library(data.table)

##### Running Foward Selection #####
unrel_ind <- fread("data/sample_unrel_short.txt", header=F) %>% rename(IID=V1)
env <- fread("data/Bans_envonly_pop.txt") %>% filter(BGP_ID %in% unrel_ind$IID)
print(nrow(env))

long_lat<-env %>% dplyr::select(Long,Lat)

#skipping the last column which doesn't have any values and biologically unimportant lc vars
pred <- env
# pred.vars <- colnames(pred[,6:24]) 
pred.vars <- c("bio01", "bio05", "bio06", "bio10", "bio12", "bio18", "bio04", "bio15")

print(pred.vars)
#scaling and centering them
pred.scale <- scale(pred[,..pred.vars], center = TRUE, scale = TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
pred.scale <- as.data.frame(pred.scale)
rownames(pred.scale) <- pred$BGP_ID

#with the 1mil subset
genotypes <- fread("subsampled_noheader.raw")
genotypes.filt <- genotypes %>% filter(genotypes$FID %in% pred$BGP_ID)
genotypes.clean <- genotypes.filt[,-(1:6)]
rownames(genotypes.clean) <- genotypes.filt$FID 

#null model
bans.mod0 <- rda(genotypes.clean  ~ 1, pred.scale) # Model with intercept only, R2 of 0

#all landcover
#all variables
# bans.rda.vars <- rda(genotypes.clean ~ bio01 + bio02 + bio03 + bio04 +
#                        bio05 + bio06 + bio07 + bio08 + bio09 + bio10 +
#                        bio13 + bio14 + bio15 + bio16 + bio17 + bio18 +
#                        bio19, pred.scale)
bans.rda.vars <- rda(genotypes.clean ~ bio01 + bio05 + bio06 +
                       bio10 + bio12 + bio18 + bio04 + bio15,
                     pred.scale)

print(RsquareAdj(bans.rda.vars))
saveRDS(bans.rda.vars,"results/RDA_output/BANS.RDAresults.envonly_b4_varsel1mil.RDS")

mod<-ordiR2step(bans.mod0,bans.rda.vars,direction="forward",R2permutations=1000,Pin=0.01,R2scope=T)  
print("ANOVA")
mod$anova

print("VIF")
vif.cca(mod)

#SAVE OUTPUT
saveRDS(mod,"results/RDA_output/BANS.varsel_envonly_results1mil.RDS") 