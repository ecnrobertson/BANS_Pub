############################# RDA_env_vars_cluster.R ##################################
##### File Checklist #####
# rsync -avzP /Users/ericarobertson/Desktop/BANS_adaptive_units/analysis/03.GEA/Bans_env_pop.txt ericacnr@colostate.edu@login.rc.colorado.edu:/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/data/
# In data/
# Bans_env_pop.txt
# subsampled_noheader.raw

##### Libraries #####
library(vegan)
library(tidyverse)
library(data.table)

##### Running Foward Selection #####
env <- fread("data/Bans_envonly_pop.txt")
print(nrow(env))

long_lat<-env %>% dplyr::select(Long,Lat)

pred <- env
pred.vars <- c("bio01", "bio04", "bio05", "bio06", "bio10", "bio12",  "bio15", "bio18")

print(pred.vars)
#scaling and centering them
pred.scale <- scale(pred[,..pred.vars], center = TRUE, scale = TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
pred.scale <- as.data.frame(pred.scale)
rownames(pred.scale) <- pred$BGP_ID
print(colnames(pred.scale))

#with the 1mil subset
genotypes <- fread("subsampled_unrelindv.1mil_noheader.raw")
genotypes.filt <- genotypes %>% filter(genotypes$FID %in% pred$BGP_ID)
genotypes.clean <- genotypes.filt[,-(1:6)]
rownames(genotypes.clean) <- genotypes.filt$FID 

#null model
bans.mod0 <- rda(genotypes.clean  ~ 1, pred.scale) # Model with intercept only, R2 of 0
print("Rsquare of null model:")
print(RsquareAdj(bans.mod0))

#all landcover
#all variables

bans.rda.vars <- rda(genotypes.clean ~ bio01 + bio05 +  + bio04 + bio06 +
                       bio10 + bio12 + bio15 + bio18,
                     pred.scale)

print("Rsquare of all variable model:")
print(RsquareAdj(bans.rda.vars))
saveRDS(bans.rda.vars,"results/RDA_output/BANS.RDAresults.envonly_bio1.4.5.6.10.12.15.18.RDS")

mod<-ordiR2step(bans.mod0,bans.rda.vars,direction="forward",R2permutations=1000,Pin=0.01,R2scope=T)  
print("ANOVA")
mod$anova

print("VIF")
vif.cca(mod)

#SAVE OUTPUT
saveRDS(mod,"results/RDA_output/BANS.RDAresults.envonly_bio1.4.5.6.10.12.15.18_varsel.RDS") 