## HEADER --------------------------------------------------------------------
##
## Script name: Forest Canopy Gap Processing, Individual Tree Level
##      [03_snag-gaps_fcgap-calc.R]
##
## Purpose of script: 
##
## Author: Jessica M. Stitt
## Date created: 2021-05-22
## Email: jessmstitt@gmail.com
## 
## NOTES ---------------------------------------------------------------------
## 
##   
## 
## begin script --------------------------------------------------------------

##--------------------------------------------------------------------------##
## LOAD PACKAGES ----
##--------------------------------------------------------------------------##
## Forest canopy gaps
library(ForestGapR) #forest gap analyses**
## **If not yet installed, first install with devtools
# library(devtools)
# devtools::install_github("carlos-alberto-silva/ForestGapR")
##--------------------------------------------------------------------------##
## Spatial manipulation
library(raster)
library(rgdal) #for vector work
library(raster) #for metadata/attributes: vectors or rasters
library(rgeos)
library(sf)
library(sp)
##--------------------------------------------------------------------------##
## Data visualization & plotting
library(ggplot2)    #for plots & graphics
library(mapview)    #for maps
library(viridis)    #for color palettes
library(RColorBrewer)   #for color palettes
library(plotrix)    #for specialized plots
##--------------------------------------------------------------------------##
## Workflow & data organization
library(here)       #for filepath mgmt
library(beepr)      #for [unnecessary] sound notifications
library(tidyverse)  #for 'tidy' data
##--------------------------------------------------------------------------##
beep(10)
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
# ASSIGN FILE PATHWAYS ----
##--------------------------------------------------------------------------##
## Datasets
csv = here(path="01_datasets", "01_csv")    #field data files
gis = here(path="01_datasets", "03_gis")    #spatial data files
##--------------------------------------------------------------------------##
## Processing folders
### Raster processing folders
chm50m = here(path=gis,"02_asc_chm-r50m")   #chm raster R50m around plots
##--------------------------------------------------------------------------##
## Scripts & outputs
scripts = here(path="02_scripts")           #R scripts
tabs = here(path="03_results", "01_tables") #results tables
gapshp = here(path="03_results", "02_gis")  #gap spdf files
lasmx = here(path="03_results", "03_las-metrics") #lidar-derived metrics
figs = here(path="03_results", "04_figures")#graphic outputs
##--------------------------------------------------------------------------##
beep(10)
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## READ IN DATASETS ----
##--------------------------------------------------------------------------##
## Load relevant field data: IPNF 2017
## Data on 25m-radius survey plots, across sites
plots <- read_csv(file.path(csv, "snag-gaps_REFplots_ipnf2017.csv")) #n=53
## Data on standing dead trees found within survey plots, across sites
snags <- read_csv(file.path(csv, "snag-gaps_REFsnags_ipnf2017.csv")) #n=270
## Data on live trees approximated within survey plots, across sites
## **Based on individual tree detection method (script within project)
itdlive <- read_csv(file.path(tabs, "snag-gaps_ITDlive_ipnf2017.csv")) #n=2186

##--------------------------------------------------------------------------##
## Load buffers around plots
plotsClip_r25m = shapefile((paste0(gis, "/", "00_shp", "/",
                                   "plotsCLIP-r25m_ipnf2017.shp"))) 
##--------------------------------------------------------------------------##
## Assign labels, based on IDs
## Survey plot level (R25m)
plotNames <- as.character(plots$plotid)
nplots <- length(plotNames)
## Individual snags (from field surveys)
snagNames <- as.character(snags$snagid)
nsnags <- length(snagNames)
snagct <- 1:nsnags
## Individual live trees (from ITD, FWS 7x7)
treeNames <- as.character(itdlive$treeid)
nlivetrees <- length(treeNames)
livetreect <- 1:nlivetrees
##--------------------------------------------------------------------------##
beep(2) 
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
## Forest canopy gap analysis in R ----
##--------------------------------------------------------------------------##
## Set parameter values: canopy height & gap area thresholds for detection 
nthresholds <- seq(2,50, by = 2) #height (m) above ground 
size <- c(1,2500)                #m^2 for area cutoffs for canopy gap extent
rasterList <- list.files(paste0(chm50m, "/"), pattern = "*.asc") #chm rasters
# crs.ipnf <- CRS("+init=epsg:26911")  #coord ref sys
##--------------------------------------------------------------------------##
beep(10)
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
# Individual tree level: snags ----
##--------------------------------------------------------------------------##
## Crop CHM rasters for gap analysis of snags [STL]
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## SNAGS (4m x 4m) ----
##--------------------------------------------------------------------------##
## Create lists for snag chms & for canopy gaps by each chm, stacked by height
list_snagRasters_4x4m <- vector("list", nsnags)
stack_snagGaps_4x4m <- vector("list", nsnags)  

## Crop/mask indiv trees from larger raster & find gaps across ht thresholds
for(i in 1:length(snagNames)){
    chmName_tmp <- paste0(chm50m, "/", 
                          str_sub(snagNames[[i]], end = -6), ".asc")
    snagCHM_R50m <- raster(chmName_tmp)
    tmpCHM <- crop(snagCHM_R50m, extent(   #square crop of buffer around snag
        snags$xc[[i]] -2,
        snags$xc[[i]] +2,
        snags$yc[[i]] -2,
        snags$yc[[i]] +2))        
    # plot(tmpCHM, col = viridis(10),
    #      main = paste0("Plot ", i," out of ", length(snagNames), ", 4x4m"))
    list_snagRasters_4x4m[[i]] <- tmpCHM
    tmpstack <- stack()
    for (j in nthresholds){
        gaps_j <- getForestGaps(chm_layer = tmpCHM, threshold = j, size = size)
        names(gaps_j) <- sprintf("gaps_%02d", j)
        tmpstack <- stack(tmpstack, gaps_j)
        # plot(gaps_j, col="black", add=TRUE, 
        #      main="Forest Canopy Gaps", legend=FALSE)
    }
    crs(tmpstack) <- "+init=epsg:26911"
    stack_snagGaps_4x4m[[i]] <- tmpstack
    # writeRaster(tmpstack, format = "GTiff",
    #             filename=paste0(gapshp, "/", "snagGaps_4x4m_", 
    #                             snagNames[[i]], ".tif"))
    print((paste0("All 4x4m gaps found for snag: ", snagNames[[i]], 
                  " (", i, " out of ", length(snagNames), ")")))
}; beep(10)
##--------------------------------------------------------------------------##
## Rename gap layers & gap height layers
##--------------------------------------------------------------------------##
names(stack_snagGaps_4x4m) <- paste0(snagNames)
gaphtNames <- sprintf(
    "gaps_%02d", 
    nthresholds, 
    length(nthresholds), replace = TRUE)
snagGaps_4x4m <- stackSave(stack_snagGaps_4x4m,
                           filename=paste0(gapshp, "/",
                                           "snag-gaps_ipnf2017_snags-4x4m",
                                           ".stk"))
##--------------------------------------------------------------------------##
## Get gap statistics
##--------------------------------------------------------------------------##
gapStats_snags_4x4m <- data.frame()
for (i in 1:length(snagNames)){
    for (j in gaphtNames){
        tmpstats <- data.frame(
            gap_id=NA, gap_area=NA, chm_max=NA, chm_min=NA, 
            chm_mean=NA, chm_sd=NA, chm_gini=NA, chm_range=NA, 
            stringsAsFactors = F)
        tmpstats <- rbind(
            tmpstats, GapStats(gap_layer = stack_snagGaps_4x4m[[i]][[j]],
                               chm_layer = list_snagRasters_4x4m[[i]]))
        tmpstats$gap_id = paste(snagNames[[i]], " ", j)
        gapStats_snags_4x4m <- rbind(gapStats_snags_4x4m, tmpstats)
    }
    print((paste0("All statistics calculated for 4x4m gaps for snag: ", 
                  snagNames[[i]], " (",i," out of ",length(snagNames),")")))
}; beep(10)

## Save raw gap statistics to Results folder
write.csv(gapStats_snags_4x4m,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_snags-4x4m.csv"))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## SNAGS (6m x 6m) ----
##--------------------------------------------------------------------------##
## Create lists for snag chms & for canopy gaps by each chm, stacked by height
list_snagRasters_6x6m <- vector("list", nsnags)
stack_snagGaps_6x6m <- vector("list", nsnags)  

## Crop/mask indiv trees from larger raster & find gaps across ht thresholds
for(i in 1:length(snagNames)){
    chmName_tmp <- paste0(chm50m, "/", 
                          str_sub(snagNames[[i]], end = -6), ".asc")
    snagCHM_R50m <- raster(chmName_tmp)
    tmpCHM <- crop(snagCHM_R50m, extent(   #square crop of buffer around snag
        snags$xc[[i]] -3,
        snags$xc[[i]] +3,
        snags$yc[[i]] -3,
        snags$yc[[i]] +3))        
    plot(tmpCHM, col = viridis(10),
         main = paste0("Plot ", i," out of ", length(snagNames), ", 6x6m"))
    list_snagRasters_6x6m[[i]] <- tmpCHM
    tmpstack <- stack()
    for (j in nthresholds){
        gaps_j <- getForestGaps(chm_layer = tmpCHM, threshold = j, size = size)
        names(gaps_j) <- sprintf("gaps_%02d", j)
        tmpstack <- stack(tmpstack, gaps_j)
        # plot(gaps_j, col="black", add=TRUE, 
        #      main="Forest Canopy Gaps", legend=FALSE)
    }
    crs(tmpstack) <- "+init=epsg:26911"
    stack_snagGaps_6x6m[[i]] <- tmpstack
    # writeRaster(tmpstack, format = "GTiff",
    #             filename=paste0(gapshp, "/", "snagGaps_6x6m_", 
    #                             snagNames[[i]], ".tif"))
    print((paste0("All 6x6m gaps found for snag: ", snagNames[[i]], 
                  " (", i, " out of ", length(snagNames), ")")))
}; beep(10)
##--------------------------------------------------------------------------##
## Rename gap layers & gap height layers
##--------------------------------------------------------------------------##
names(stack_snagGaps_6x6m) <- paste0(snagNames)
gaphtNames <- sprintf(
    "gaps_%02d", 
    nthresholds, 
    length(nthresholds), replace = TRUE)
snagGaps_6x6m <- stackSave(stack_snagGaps_6x6m,
                           filename=paste0(gapshp, "/",
                                           "snag-gaps_ipnf2017_snags-6x6m",
                                           ".stk"))
##--------------------------------------------------------------------------##
## Get gap statistics
##--------------------------------------------------------------------------##
gapStats_snags_6x6m <- data.frame()
for (i in 1:length(snagNames)){
    for (j in gaphtNames){
        tmpstats <- data.frame(
            gap_id=NA, gap_area=NA, chm_max=NA, chm_min=NA, 
            chm_mean=NA, chm_sd=NA, chm_gini=NA, chm_range=NA, 
            stringsAsFactors = F)
        tmpstats <- rbind(
            tmpstats, GapStats(gap_layer = stack_snagGaps_6x6m[[i]][[j]],
                               chm_layer = list_snagRasters_6x6m[[i]]))
        tmpstats$gap_id = paste(snagNames[[i]], " ", j)
        gapStats_snags_6x6m <- rbind(gapStats_snags_6x6m, tmpstats)
    }
    print((paste0("All statistics calculated for 6x6m gaps for snag: ", 
                  snagNames[[i]], " (",i," out of ",length(snagNames),")")))
}; beep(10)

## Save raw gap statistics to Results folder
write.csv(gapStats_snags_6x6m,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_snags-6x6m.csv"))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## SNAGS (8m x 8m) ----
##--------------------------------------------------------------------------##
## Create lists for snag chms & for canopy gaps by each chm, stacked by height
list_snagRasters_8x8m <- vector("list", nsnags)
stack_snagGaps_8x8m <- vector("list", nsnags)  

## Crop/mask indiv trees from larger raster & find gaps across ht thresholds
for(i in 1:length(snagNames)){
    chmName_tmp <- paste0(chm50m, "/", 
                          str_sub(snagNames[[i]], end = -6), ".asc")
    snagCHM_R50m <- raster(chmName_tmp)
    tmpCHM <- crop(snagCHM_R50m, extent(   #square crop of buffer around snag
        snags$xc[[i]] -4,
        snags$xc[[i]] +4,
        snags$yc[[i]] -4,
        snags$yc[[i]] +4))        
    # plot(tmpCHM, col = viridis(10),
    #      main = paste0("Plot ", i," out of ", length(snagNames), ", 8x8m"))
    list_snagRasters_8x8m[[i]] <- tmpCHM
    tmpstack <- stack()
    for (j in nthresholds){
        gaps_j <- getForestGaps(chm_layer = tmpCHM, threshold = j, size = size)
        names(gaps_j) <- sprintf("gaps_%02d", j)
        tmpstack <- stack(tmpstack, gaps_j)
        # plot(gaps_j, col="black", add=TRUE, 
        #      main="Forest Canopy Gaps", legend=FALSE)
    }
    crs(tmpstack) <- "+init=epsg:26911"
    stack_snagGaps_8x8m[[i]] <- tmpstack
    # writeRaster(tmpstack, format = "GTiff",
    #             filename=paste0(gapshp, "/", "snagGaps_8x8m_", 
    #                             snagNames[[i]], ".tif"))
    print((paste0("All 8x8m gaps found for snag: ", snagNames[[i]], 
                  " (", i, " out of ", length(snagNames), ")")))
}; beep(10)
##--------------------------------------------------------------------------##
## Rename gap layers & gap height layers
##--------------------------------------------------------------------------##
names(stack_snagGaps_8x8m) <- paste0(snagNames)
gaphtNames <- sprintf(
    "gaps_%02d", 
    nthresholds, 
    length(nthresholds), replace = TRUE)
snagGaps_8x8m <- stackSave(stack_snagGaps_8x8m,
                           filename=paste0(gapshp, "/",
                                           "snag-gaps_ipnf2017_snags-8x8m",
                                           ".stk"))
##--------------------------------------------------------------------------##
## Get gap statistics
##--------------------------------------------------------------------------##
gapStats_snags_8x8m <- data.frame()
for (i in 1:length(snagNames)){
    for (j in gaphtNames){
        tmpstats <- data.frame(
            gap_id=NA, gap_area=NA, chm_max=NA, chm_min=NA, 
            chm_mean=NA, chm_sd=NA, chm_gini=NA, chm_range=NA, 
            stringsAsFactors = F)
        tmpstats <- rbind(
            tmpstats, GapStats(gap_layer = stack_snagGaps_8x8m[[i]][[j]],
                               chm_layer = list_snagRasters_8x8m[[i]]))
        tmpstats$gap_id = paste(snagNames[[i]], " ", j)
        gapStats_snags_8x8m <- rbind(gapStats_snags_8x8m, tmpstats)
    }
    print((paste0("All statistics calculated for 8x8m gaps for snag: ", 
                  snagNames[[i]], " (",i," out of ",length(snagNames),")")))
}; beep(10)

## Save raw gap statistics to Results folder
write.csv(gapStats_snags_8x8m,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_snags-8x8m.csv"))

##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
# Individual tree level: live trees ----
##--------------------------------------------------------------------------##
## Crop CHM rasters for gap analysis of live trees [STL]
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## ITD LIVE TREES (4m x 4m) ----
##--------------------------------------------------------------------------##
## Create lists for snag chms & for canopy gaps by each chm, stacked by height
list_liveRasters_4x4m <- vector("list", nlivetrees)
stack_liveGaps_4x4m <- vector("list", nlivetrees)  

## Crop/mask indiv trees from larger raster & find gaps across ht thresholds
for(i in 1:length(treeNames)){
    chmName_tmp <- paste0(chm50m, "/", 
                          str_sub(treeNames[[i]], end = -5), ".asc")
    liveCHM_R50m <- raster(chmName_tmp)
    tmpCHM <- crop(liveCHM_R50m, extent(   #square crop of buffer around snag
        itdlive$xc[[i]] -2,
        itdlive$xc[[i]] +2,
        itdlive$yc[[i]] -2,
        itdlive$yc[[i]] +2))        
    # plot(tmpCHM, col = viridis(10),
    #      main = paste0("Plot ", i," out of ", length(treeNames), ", 4x4m"))
    list_liveRasters_4x4m[[i]] <- tmpCHM
    tmpstack <- stack()
    for (j in nthresholds){
        gaps_j <- getForestGaps(chm_layer = tmpCHM, threshold = j, size = size)
        names(gaps_j) <- sprintf("gaps_%02d", j)
        tmpstack <- stack(tmpstack, gaps_j)
        # plot(gaps_j, col="black", add=TRUE, 
        #      main="Forest Canopy Gaps", legend=FALSE)
    }
    crs(tmpstack) <- "+init=epsg:26911"
    stack_liveGaps_4x4m[[i]] <- tmpstack
    # writeRaster(tmpstack, format = "GTiff",
    #             filename=paste0(gapshp, "/", "liveGaps_4x4m_", 
    #                             treeNames[[i]], ".tif"))
    print((paste0("All 4x4m gaps found for live tree: ", treeNames[[i]], 
                  " (", i, " out of ", length(treeNames), ")")))
}; beep(10)
##--------------------------------------------------------------------------##
## Rename gap layers & gap height layers
##--------------------------------------------------------------------------##
names(stack_liveGaps_4x4m) <- paste0(treeNames)
gaphtNames <- sprintf(
    "gaps_%02d", 
    nthresholds, 
    length(nthresholds), replace = TRUE)
liveGaps_4x4m <- stackSave(stack_liveGaps_4x4m,
                           filename=paste0(gapshp, "/",
                                           "snag-gaps_ipnf2017_live-4x4m",
                                           ".stk"))
##--------------------------------------------------------------------------##
## Get gap statistics
##--------------------------------------------------------------------------##
gapStats_live_4x4m <- data.frame()
for (i in 1:length(treeNames)){
    for (j in gaphtNames){
        tmpstats <- data.frame(
            gap_id=NA, gap_area=NA, chm_max=NA, chm_min=NA, 
            chm_mean=NA, chm_sd=NA, chm_gini=NA, chm_range=NA, 
            stringsAsFactors = F)
        tmpstats <- rbind(
            tmpstats, GapStats(gap_layer = stack_liveGaps_4x4m[[i]][[j]],
                               chm_layer = list_liveRasters_4x4m[[i]]))
        tmpstats$gap_id = paste(treeNames[[i]], " ", j)
        gapStats_live_4x4m <- rbind(gapStats_live_4x4m, tmpstats)
    }
    print((paste0("All statistics calculated for 4x4m gaps for live tree: ", 
                  treeNames[[i]], " (",i," out of ",length(treeNames),")")))
}; beep(10)

## Save raw gap statistics to Results folder
write.csv(gapStats_live_4x4m,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_live-4x4m.csv"))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## ITD LIVE TREES (6m x 6m) ----
##--------------------------------------------------------------------------##
## Create lists for snag chms & for canopy gaps by each chm, stacked by height
list_liveRasters_6x6m <- vector("list", nlivetrees)
stack_liveGaps_6x6m <- vector("list", nlivetrees)  

## Crop/mask indiv trees from larger raster & find gaps across ht thresholds
for(i in 1:length(treeNames)){
    chmName_tmp <- paste0(chm50m, "/", 
                          str_sub(treeNames[[i]], end = -5), ".asc")
    liveCHM_R50m <- raster(chmName_tmp)
    tmpCHM <- crop(liveCHM_R50m, extent(   #square crop of buffer around snag
        itdlive$xc[[i]] -3,
        itdlive$xc[[i]] +3,
        itdlive$yc[[i]] -3,
        itdlive$yc[[i]] +3))        
    plot(tmpCHM, col = viridis(10),
         main = paste0("Plot ", i," out of ", length(treeNames), ", 6x6m"))
    list_liveRasters_6x6m[[i]] <- tmpCHM
    tmpstack <- stack()
    for (j in nthresholds){
        gaps_j <- getForestGaps(chm_layer = tmpCHM, threshold = j, size = size)
        names(gaps_j) <- sprintf("gaps_%02d", j)
        tmpstack <- stack(tmpstack, gaps_j)
        # plot(gaps_j, col="black", add=TRUE, 
        #      main="Forest Canopy Gaps", legend=FALSE)
    }
    crs(tmpstack) <- "+init=epsg:26911"
    stack_liveGaps_6x6m[[i]] <- tmpstack
    # writeRaster(tmpstack, format = "GTiff",
    #             filename=paste0(gapshp, "/", "liveGaps_6x6m_", 
    #                             treeNames[[i]], ".tif"))
    print((paste0("All 6x6m gaps found for live tree: ", treeNames[[i]], 
                  " (", i, " out of ", length(treeNames), ")")))
}; beep(10)
##--------------------------------------------------------------------------##
## Rename gap layers & gap height layers
##--------------------------------------------------------------------------##
names(stack_liveGaps_6x6m) <- paste0(treeNames)
gaphtNames <- sprintf(
    "gaps_%02d", 
    nthresholds, 
    length(nthresholds), replace = TRUE)
liveGaps_6x6m <- stackSave(stack_liveGaps_6x6m,
                           filename=paste0(gapshp, "/",
                                           "snag-gaps_ipnf2017_live-6x6m",
                                           ".stk"))
##--------------------------------------------------------------------------##
## Get gap statistics
##--------------------------------------------------------------------------##
gapStats_live_6x6m <- data.frame()
for (i in 1:length(treeNames)){
    for (j in gaphtNames){
        tmpstats <- data.frame(
            gap_id=NA, gap_area=NA, chm_max=NA, chm_min=NA, 
            chm_mean=NA, chm_sd=NA, chm_gini=NA, chm_range=NA, 
            stringsAsFactors = F)
        tmpstats <- rbind(
            tmpstats, GapStats(gap_layer = stack_liveGaps_6x6m[[i]][[j]],
                               chm_layer = list_liveRasters_6x6m[[i]]))
        tmpstats$gap_id = paste(treeNames[[i]], " ", j)
        gapStats_live_6x6m <- rbind(gapStats_live_6x6m, tmpstats)
    }
    print((paste0("All statistics calculated for 6x6m gaps for live tree: ", 
                  treeNames[[i]], " (",i," out of ",length(treeNames),")")))
}; beep(10)

## Save raw gap statistics to Results folder
write.csv(gapStats_live_6x6m,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_live-6x6m.csv"))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## ITD LIVE TREES (8m x 8m) ----
##--------------------------------------------------------------------------##
## Create lists for snag chms & for canopy gaps by each chm, stacked by height
list_liveRasters_8x8m <- vector("list", nlivetrees)
stack_liveGaps_8x8m <- vector("list", nlivetrees)  

## Crop/mask indiv trees from larger raster & find gaps across ht thresholds
for(i in 1:length(treeNames)){
    chmName_tmp <- paste0(chm50m, "/", 
                          str_sub(treeNames[[i]], end = -5), ".asc")
    liveCHM_R50m <- raster(chmName_tmp)
    tmpCHM <- crop(liveCHM_R50m, extent(   #square crop of buffer around snag
        itdlive$xc[[i]] -4,
        itdlive$xc[[i]] +4,
        itdlive$yc[[i]] -4,
        itdlive$yc[[i]] +4))        
    # plot(tmpCHM, col = viridis(10),
    #      main = paste0("Plot ", i," out of ", length(treeNames), ", 8x8m"))
    list_liveRasters_8x8m[[i]] <- tmpCHM
    tmpstack <- stack()
    for (j in nthresholds){
        gaps_j <- getForestGaps(chm_layer = tmpCHM, threshold = j, size = size)
        names(gaps_j) <- sprintf("gaps_%02d", j)
        tmpstack <- stack(tmpstack, gaps_j)
        # plot(gaps_j, col="black", add=TRUE, 
        #      main="Forest Canopy Gaps", legend=FALSE)
    }
    crs(tmpstack) <- "+init=epsg:26911"
    stack_liveGaps_8x8m[[i]] <- tmpstack
    # writeRaster(tmpstack, format = "GTiff",
    #             filename=paste0(gapshp, "/", "liveGaps_8x8m_", 
    #                             treeNames[[i]], ".tif"))
    print((paste0("All 8x8m gaps found for live tree: ", treeNames[[i]], 
                  " (", i, " out of ", length(treeNames), ")")))
}; beep(10)
##--------------------------------------------------------------------------##
## Rename gap layers & gap height layers
##--------------------------------------------------------------------------##
names(stack_liveGaps_8x8m) <- paste0(treeNames)
gaphtNames <- sprintf(
    "gaps_%02d", 
    nthresholds, 
    length(nthresholds), replace = TRUE)
liveGaps_8x8m <- stackSave(stack_liveGaps_8x8m,
                           filename=paste0(gapshp, "/",
                                           "snag-gaps_ipnf2017_live-8x8m",
                                           ".stk"))
##--------------------------------------------------------------------------##
## Get gap statistics
##--------------------------------------------------------------------------##
gapStats_live_8x8m <- data.frame()
for (i in 1:length(treeNames)){
    for (j in gaphtNames){
        tmpstats <- data.frame(
            gap_id=NA, gap_area=NA, chm_max=NA, chm_min=NA, 
            chm_mean=NA, chm_sd=NA, chm_gini=NA, chm_range=NA, 
            stringsAsFactors = F)
        tmpstats <- rbind(
            tmpstats, GapStats(gap_layer = stack_liveGaps_8x8m[[i]][[j]],
                               chm_layer = list_liveRasters_8x8m[[i]]))
        tmpstats$gap_id = paste(treeNames[[i]], " ", j)
        gapStats_live_8x8m <- rbind(gapStats_live_8x8m, tmpstats)
    }
    print((paste0("All statistics calculated for 8x8m gaps for live tree: ", 
                  treeNames[[i]], " (",i," out of ",length(treeNames),")")))
}; beep(10)

## Save raw gap statistics to Results folder
write.csv(gapStats_live_8x8m,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_live-8x8m.csv"))

##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## Process gap statistics by footprint size ----
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## ALL SNAGS & LIVE TREES: 4m x 4m ----
##--------------------------------------------------------------------------##
## Calculate snag statistics (fp4x4m)
fp4x4_snagStats_proc <- as_tibble(gapStats_snags_4x4m) %>%
    separate("gap_id", c("plotid", "snagid", "gaps", "gap_ht")) %>%
    select(-c(gaps)) %>%
    unite("treeid", "plotid":"snagid", 
          sep = "_", remove=TRUE) %>% 
    add_column(tree_type = "snag", .after = "treeid"); fp4x4_snagStats_proc
fp4x4_snagStats_clean <- na.omit(fp4x4_snagStats_proc)

##--------------------------------------------------------------------------##
## Calculate live tree statistics (fp4x4m)
fp4x4_liveStats_proc <- as_tibble(gapStats_live_4x4m) %>%
    separate("gap_id", c("plotid", "treeid", "gaps", "gap_ht")) %>%
    select(-c(gaps)) %>%
    unite("treeid", "plotid":"treeid", 
          sep = "_", remove=TRUE) %>% 
    add_column(tree_type = "live", .after = "treeid"); fp4x4_liveStats_proc
fp4x4_liveStats_clean <- na.omit(fp4x4_liveStats_proc)

##--------------------------------------------------------------------------##
## Combine all tree statistics (fp4x4m)
fp4x4_gapStats_clean <- bind_rows(fp4x4_snagStats_clean, fp4x4_liveStats_clean)
write.csv(fp4x4_gapStats_clean,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_all-trees_4x4m.csv"))

##--------------------------------------------------------------------------##
## Aggregate & summarize all gaps per slice
fp4x4_gapStats_summ <- fp4x4_gapStats_clean %>%
    group_by(treeid, gap_ht) %>%
    summarise(gap_area_tot = sum(gap_area), 
              chm_min=min(chm_min), chm_max=max(chm_max), 
              tree_type=unique(tree_type))
write.csv(fp4x4_gapStats_summ,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats-summary_4x4m.csv"))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## ALL SNAGS & LIVE TREES: 6m x 6m ----
##--------------------------------------------------------------------------##
## Calculate snag statistics (fp6x6m)
fp6x6_snagStats_proc <- as_tibble(gapStats_snags_6x6m) %>%
    separate("gap_id", c("plotid", "snagid", "gaps", "gap_ht")) %>%
    select(-c(gaps)) %>%
    unite("treeid", "plotid":"snagid", 
          sep = "_", remove=TRUE) %>% 
    add_column(tree_type = "snag", .after = "treeid"); fp6x6_snagStats_proc
fp6x6_snagStats_clean <- na.omit(fp6x6_snagStats_proc)

##--------------------------------------------------------------------------##
## Calculate live tree statistics (fp6x6m)
fp6x6_liveStats_proc <- as_tibble(gapStats_live_6x6m) %>%
    separate("gap_id", c("plotid", "treeid", "gaps", "gap_ht")) %>%
    select(-c("gaps")) %>%
    unite("treeid", "plotid":"treeid", 
          sep = "_", remove=TRUE) %>% 
    add_column(tree_type = "live", .after = "treeid"); fp6x6_liveStats_proc
fp6x6_liveStats_clean <- na.omit(fp6x6_liveStats_proc)

##--------------------------------------------------------------------------##
## Combine all tree statistics (fp6x6m)
fp6x6_gapStats_clean <- bind_rows(fp6x6_snagStats_clean, fp6x6_liveStats_clean)
write.csv(fp6x6_gapStats_clean,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_all-trees_6x6m.csv"))

##--------------------------------------------------------------------------##
## Aggregate & summarize all gaps per slice
fp6x6_gapStats_summ <- fp6x6_gapStats_clean %>%
    group_by(treeid, gap_ht) %>%
    summarise(gap_area_tot = sum(gap_area), 
              chm_min=min(chm_min), chm_max=max(chm_max), 
              tree_type=unique(tree_type))
write.csv(fp6x6_gapStats_summ,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats-summary_6x6m.csv"))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## ALL SNAGS & LIVE TREES: 8m x 8m ----
##--------------------------------------------------------------------------##
## Calculate snag statistics (fp8x8m)
fp8x8_snagStats_proc <- as_tibble(gapStats_snags_8x8m) %>%
    separate("gap_id", c("plotid", "snagid", "gaps", "gap_ht")) %>%
    select(-c(gaps)) %>%
    unite("treeid", "plotid":"snagid", 
          sep = "_", remove=TRUE) %>% 
    add_column(tree_type = "snag", .after = "treeid"); fp8x8_snagStats_proc
fp8x8_snagStats_clean <- na.omit(fp8x8_snagStats_proc)

##--------------------------------------------------------------------------##
## Calculate live tree statistics (fp8x8m)
fp8x8_liveStats_proc <- as_tibble(gapStats_live_8x8m) %>%
    separate("gap_id", c("plotid", "treeid", "gaps", "gap_ht")) %>%
    select(-c(gaps)) %>%
    unite("treeid", "plotid":"treeid", 
          sep = "_", remove=TRUE) %>% 
    add_column(tree_type = "live", .after = "treeid"); fp8x8_liveStats_proc
fp8x8_liveStats_clean <- na.omit(fp8x8_liveStats_proc)

##--------------------------------------------------------------------------##
## Combine all tree statistics (fp8x8m)
fp8x8_gapStats_clean <- bind_rows(fp8x8_snagStats_clean, fp8x8_liveStats_clean)
write.csv(fp8x8_gapStats_clean,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats_all-trees_8x8m.csv"))

##--------------------------------------------------------------------------##
## Aggregate & summarize all gaps per slice
fp8x8_gapStats_summ <- fp8x8_gapStats_clean %>%
    group_by(treeid, gap_ht) %>%
    summarise(gap_area_tot = sum(gap_area), 
              chm_min=min(chm_min), chm_max=max(chm_max), 
              tree_type=unique(tree_type))
write.csv(fp8x8_gapStats_summ,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats-summary_8x8m.csv"))

##--------------------------------------------------------------------------##
##

##--------------------------------------------------------------------------##
## Combine gap statistics for all buffer sizes ----
##--------------------------------------------------------------------------##
## Process stats from each buffer size, add buffer ID column 
gapStats_proc_fp4x4 <- add_column(fp4x4_gapStats_summ, fp_m2 = "16m2", 
                                  .before = "tree_type" ) %>%
    mutate(site = case_when(
        startsWith(treeid, "BR") ~ "cda",
        startsWith(treeid, "HS") ~ "stj"
    ))
    
gapStats_proc_fp6x6 <- add_column(fp6x6_gapStats_summ, fp_m2 = "36m2", 
                                  .before = "tree_type") %>%
    mutate(site = case_when(
        startsWith(treeid, "BR") ~ "cda",
        startsWith(treeid, "HS") ~ "stj"
    ))
gapStats_proc_fp8x8 <- add_column(fp8x8_gapStats_summ, fp_m2 = "64m2", 
                                  .before = "tree_type") %>%
    mutate(site = case_when(
        startsWith(treeid, "BR") ~ "cda",
        startsWith(treeid, "HS") ~ "stj"
    ))

##--------------------------------------------------------------------------##
## Combine all FPS
##--------------------------------------------------------------------------##
gapStats_clean <- bind_rows(gapStats_proc_fp4x4,
                            gapStats_proc_fp6x6, 
                            gapStats_proc_fp8x8) %>%
    dplyr::relocate(site, .after = "gap_ht") %>%
    dplyr::relocate(fp_m2, .after = "site")  %>%
    dplyr::relocate(tree_type, .after = "fp_m2") %>%
    dplyr::mutate(fp_num = recode(fp_m2,
                           16 = `16m2`,
                           36 = `36m2`,
                           64 = `64m2`)) %>%
    dplyr::mutate(prop_gap = gap_area_tot / fp_num); gapStats_clean

# write.csv(gapStats_clean,
#           file.path(tabs, "snag-gaps_ipnf2017_gapstats_all-fp.csv"))
# gaps <- gapStats_clean
gaps <- read_csv(file.path(tabs, "snag-gaps_ipnf2017_gapstats_all-fp.csv")) %>%
    dplyr::select(-X1)
##--------------------------------------------------------------------------##
## Summarize gap stats for all FPS across live and dead trees ----
##--------------------------------------------------------------------------##
gapStats_summ <- gaps %>% 
    group_by(tree_type, fp_m2, gap_ht) %>%
    summarise(gap_area_tot = sum(gap_area_tot), 
              tree_type = unique(tree_type),
              site = unique(site),
              fp_m2 = unique(fp_m2)); gapStats_summ

write.csv(gapStats_summ,
          file.path(tabs, "snag-gaps_ipnf2017_gapstats-aggregate.csv"))


##--------------------------------------------------------------------------##
## Analyze differences withingap stats across live and dead trees ----
##--------------------------------------------------------------------------##
gapStats_analysis <- gaps %>%
    group_by(tree_type, site, fp_m2, gap_ht) %>%
    mutate(chm_range = chm_max - chm_min) %>% 
    summarise(across(c(prop_gap, gap_area_tot, chm_range),
        list(avg = mean, se = std.error))); gapStats_analysis

write.csv(gapStats_analysis,
         file.path(tabs, "snag-gaps_ipnf2017_gapstats_avgHT-se_all-fp.csv"))

gapStats_graphics <- gapStats_analysis[,c(1:6)]  %>%
    rename("prop_mu" = prop_gap_avg)  %>%
    rename("prop_se" = prop_gap_se); gapStats_graphics

## Compute summary statistics by groups: mean and sd
gaps %>%
    group_by(fp_m2, tree_type, site) %>%
    get_summary_stats(prop_gap, type = "common") 

##--------------------------------------------------------------------------##
beep(11)

##--------------------------------------------------------------------------##