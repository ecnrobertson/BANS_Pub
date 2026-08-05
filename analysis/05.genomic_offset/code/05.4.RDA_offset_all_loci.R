library(adegenet)
library(vegan)
library(data.table)
library(dplyr)
library(readr)
library(raster)
library(sf)
library(rnaturalearth)

options(warn = 1)  # print warnings as they occur
options(error = function() {
  traceback(2)
  quit(status = 1)
})

tmpbase <- "/scratch/alpine/ericacnr@colostate.edu/tmp"
jid <- Sys.getenv("SLURM_JOB_ID")
if (jid == "") jid <- paste0("pid_", Sys.getpid())
jobtmp <- file.path(tmpbase, jid)
dir.create(jobtmp, recursive=TRUE, showWarnings=FALSE)
raster::rasterOptions(tmpdir=jobtmp, progress="text", chunksize=1e8, maxmemory=2e10)

cat("raster tmpdir:", raster::rasterOptions()$tmpdir, "\n")
################## OUTPUTS ###################
date <- "2041"
RDS_output_result <- "results/RDA_proj_offset_all_loci_table_2041-2070_ssp126_bio5.15.18.rds"
future_path <- "data/ClimateNA_2041-2070/ssp126"

################## READING IN RDA ####################
rda_env_mems_uncor <- readRDS("/scratch/alpine/ericacnr@colostate.edu/BANS/03.GEA/results/RDA_output/RDA_output/BANS.RDAresults.envonly_bio5.15.18_varsel.RDS")

################## Setting up for genomic offset ####################
## Loading the climatic rasters
cont.paths <- list.files("data/ClimateNA_1981-2010", recursive = TRUE, full.names = TRUE)

cont.stack <- raster::stack(cont.paths)
names(cont.stack) <- c("bio05", "bio15", "bio18")

class(cont.stack)

# 2041-2070
futr.paths <- list.files(future_path, recursive = TRUE, full.names = T)
futr.stack <- raster::stack(futr.paths)
names(futr.stack) <- c("bio05", "bio15", "bio18")


class(futr.stack)

print("raster stacks projected")

## Loading the Breeding Range as a clipping mask
map <- ne_states(country = c("United States of America", "Canada"), returnclass = 'sf')
na_union <- st_union(map)
breeding.sf <- st_read("data/banswa_range_2023.gpkg") %>% filter(season=="breeding") %>% st_transform(crs(map))
breeding_na <- st_intersection(breeding.sf, na_union)
breeding_na <- st_make_valid(breeding_na)
range_sf <- st_union(breeding_na) |> st_make_valid() |> st_transform(4326)
range_sp <- as(range_sf, "Spatial")
plot(range_sp)

template <- cont.stack[[1]]

# rasterize: 1 inside breeding range, NA outside
range_raster <- rasterize(
  range_sp,
  template,
  field = 1,
  background = NA
)

## Source populations coordinates
unrel_ind <- fread("../03.GEA/data/sample_unrel_short.txt", header=F) %>% rename(IID=V1)
indorder <- fread("../03.GEA/BANS.unrel_ind.imputed4.1.fam") %>% rename(BGP_ID = V1) %>% dplyr::select(BGP_ID)
coords <- read_csv("data/BANS_all_sample_data_shifted.csv") %>%
  filter(BGP_ID %in% indorder$BGP_ID) %>%
  dplyr::select(BGP_ID, Long, Lat) %>% filter(BGP_ID %in% unrel_ind$IID)
nrow(coords)

print("reorder coordinates to match PCA data")
coords <- coords[(match(indorder$BGP_ID, coords$BGP_ID)),]

colnames(coords) <- c("BGP_ID", "Long", "Lat")

## Extracting environmental values for each source population
Env <- data.frame(raster::extract(cont.stack, coords[,2:3]))

## Standardization of the variables
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()

## Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')

## Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- c(coords$Population)
print("Centered and Scaled")
##################### Running the Genomic Offset #########################

## Running the function for 2070
source("src/genomic_offset.R")
#K should be the number of variables you're putting in, so I've got 3 enviro vars: K=3
res_RDA_proj <- genomic_offset(rda_env_mems_uncor, K = 3, env_pres = cont.stack, env_fut = futr.stack, range = range_raster, method = "loadings", scale_env =
                                 scale_env, center_env = center_env)
print("res_RDA_proj2041 done")
## Table global genetic offset predicted for 2050 and 2080
RDA_proj_offset <- data.frame(rasterToPoints(res_RDA_proj$Proj_offset_global), Date = rep(date, nrow(rasterToPoints(res_RDA_proj$Proj_offset_global))))
saveRDS(RDA_proj_offset, RDS_output_result)
print("saved RDA_proj_offset_all_loci_table.rds")