## HEADER --------------------------------------------------------------------
##
## Script name: Individual Tree Detection, Accuracy Assessment (Quality Check)
##    [02_snag-gaps_itd-qaqc.R]
## 
## Purpose of script: Calibrate individual tree detection (ITD) methods to 
##  estimate number and location of live trees in all 25m-radius survey plots
##    Inputs: CHMs of field survey plots (R25m) + subplot sample data (R11.3m)
##    Outputs: accuracy assessments of ITD methods + live tree coordinates
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
## Lidar analysis in R
library(lidR)       #for lidar data manipulation**
library(rLiDAR)     #for lidar data manipulation & ITD**
## **If not yet installed, first install with devtools
# library(devtools)
# devtools::install_github("Jean-Romain/lidR")
# devtools::install_github("carlos-alberto-silva/rLiDAR")
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
laspc = here(path="01_datasets", "02_las")  #lidar point cloud files
gis = here(path="01_datasets", "03_gis")    #spatial data files
##--------------------------------------------------------------------------##
## Processing folders
### Raster processing folders
pfc = here(path=gis,"01_tif_pitfree-chm")   #pit-free chm R50m around plots
chm50m = here(path=gis,"02_asc_chm-r50m")   #chm raster R50m around plots
chmPlot = here(gis,"03_asc_chm-plot")       #chm raster from chm to R25m
rasterRMRS = here(gis,"04_asc_chm-rmrs")    #chm raster from chm to R11.3m
rasterTree = here(gis,"05_asc_chm-tree")    #chm raster from chm to ITL: 
rasterTree_4x4m_snag = here(path=rasterTree, "fp_4x4m", "snag") #4x4m chm snag
rasterTree_6x6m_snag = here(path=rasterTree, "fp_6x6m", "snag") #6x6m chm snag
rasterTree_8x8m_snag = here(path=rasterTree, "fp_8x8m", "snag") #8x8m chm snag
rasterTree_4x4m_live = here(path=rasterTree, "fp_4x4m", "live") #4x4m chm live
rasterTree_6x6m_live = here(path=rasterTree, "fp_6x6m", "live") #6x6m chm live
rasterTree_8x8m_live = here(path=rasterTree, "fp_8x8m", "live") #8x8m chm live
##--------------------------------------------------------------------------##
## Scripts & outputs
scripts = here(path="02_scripts")           #R scripts
tabs = here(path="03_results", "01_tables") #results tables
figs = here(path="03_results", "02_figures")#graphic outputs
lasmx = here(path="03_results", "03_las-metrics") #lidar-derived metrics
##--------------------------------------------------------------------------##
beep(10)
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## READ IN DATASETS ----
##--------------------------------------------------------------------------##
## Load relevant field data: IPNF 2017
## Data on 25m-radius survey plots, across sites
REFplots <- read.csv(file.path(csv, "snag-gaps_REFplots_ipnf2017.csv")) #n=53
plots <- as_tibble(REFplots)
## Data on USFS RMRS 11.3m-radius surveys, within a subset of survey plots
rmrs <- read.csv(file.path(csv, "snag-gaps_subplots-r11m_ipnf2017.csv"))#n=32
qaqc <- read_csv(file.path(csv, "snag-gaps_subplots-r11m_ipnf2017.csv"))#tibble
##--------------------------------------------------------------------------##
## Load buffers around plots (R50m, R25m circles)
# plotsClip_r50m = shapefile((paste0(gis, "/", "00_shp", "/",
#                                    "plotsCLIP-r50m_ipnf2017.shp"))) 
plotsClip_r25m = shapefile((paste0(gis, "/", "00_shp", "/",
                                   "plotsCLIP-r25m_ipnf2017.shp"))) 
plotsClip_r11m = shapefile((paste0(gis, "/", "00_shp", "/",
                                   "subplotsCLIP-r11m_ipnf2017.shp")))
##--------------------------------------------------------------------------##
## Assign names
## Set labels, based on field survey plot IDs
plotNames <- as.character(plots$plotid)
subplots <- subset(plots, subplot=="y")
subplotNames <- as.character(subplots$plotid)
## Identify number of plots
nplots <- length(plotNames)
nsubplots <- length(subplotNames)
##--------------------------------------------------------------------------##
beep(2) 
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
# INDIVIDUAL TREE DETECTION ----
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## Run ITD on live trees to generate list for gap analysis ----
##--------------------------------------------------------------------------##
## List of CHM rasters (R50m)
refRasters <- as.list(paste0(chm50m,"/", plotNames, ".asc"))

## Setup 
crs.ipnf <- CRS("+init=epsg:26911")  #coord ref sys

##--------------------------------------------------------------------------##
## 3x3m Scale ----
##--------------------------------------------------------------------------##
## Parameters for ITD function {rLidar}
fws <- 3   #fixed window size = dimensions (m); default = 3 (3x3)
minht <- 2 #min height above ground for the detection break (m)

## Dataframe to add ITD trees to (reset to empty before use)
treesxyList <- data.frame() 
## Run ITD method for all trees within 50m grid {rLidar}
for(i in 1:length(refRasters)){
    ## Load raster CHM for given plot (R50m)
    plotCHM <- raster(refRasters[[i]])
    ## Find individual trees within given plot
    treesPlot <- FindTreesCHM(plotCHM, fws, minht)
    treesPlot$treeid <- 0 # does this reset treeID?
    names(treesPlot)[names(treesPlot)=="x"] <- "xc"
    names(treesPlot)[names(treesPlot)=="y"] <- "yc"
    treesPlot <- treesPlot[,c(4,1,2,3)]
    for(j in 1:length(treesPlot$treeid)){
        ## Alter treeid name to include plotNames[[i]], "_", i
        treesPlot$treeid[[j]]<-paste0(plotNames[[i]], sprintf("_%03d",j))
    }
    ## Add all trees within given plot to treesxyList
    treesxyList <- rbind(treesxyList, treesPlot)
    print((paste0("Completed ITD for all trees in plot: ", plotNames[[i]], 
                  " (", i, " out of ", length(plotNames), ")")))
}; beep(11)

tidytrees_fws3x3m <- as.tibble(treesxyList) %>% 
    separate(treeid, c("plotid"), sep = "_", remove=FALSE) %>%
    relocate(plotid, .before = "treeid")
tidytrees_fws3x3m

# write.csv(tidytrees_fws3x3m,
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R50m_FWS3x3.csv"))
tidytrees_r50m_fws3x3 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R50m_FWS3x3.csv")) 

## Convert R50m ITD tree list to spatial points
ITD_r50m_fws3x3 <- SpatialPointsDataFrame(tidytrees_r50m_fws3x3[,3:4], 
                                          tidytrees_r50m_fws3x3)
# coordinates(ITD_r50m_fws3x3)
proj4string(ITD_r50m_fws3x3) <- crs.ipnf  # define projection system of data
# is.projected(ITD_r50m_fws3x3)
summary(ITD_r50m_fws3x3)
# plot(ITD_r50m_fws3x3)

## Crop from R50m to R11m for all plots
ITD_r11m_fws3x3 <- raster::crop(ITD_r50m_fws3x3, plotsClip_r11m)

# write.csv(ITD_r11m_fws3x3@data,
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS3x3.csv"))
tidytrees_r11m_fws3x3 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS3x3.csv")) %>%
    select(-X1)
##--------------------------------------------------------------------------##
## 5x5m Scale ----
##--------------------------------------------------------------------------##
## Parameters for ITD function {rLidar}
fws <- 5   #fixed window size = dimensions (m); default = 3 (3x3)
minht <- 2 #min height above ground for the detection break (m)

## Dataframe to add ITD trees to
treesxyList <- data.frame() 
## Run ITD method for all trees within 50m grid {rLidar}
for(i in 1:length(refRasters)){
    ## Load raster CHM for given plot (R50m)
    plotCHM <- raster(refRasters[[i]])
    ## Find individual trees within given plot
    treesPlot <- FindTreesCHM(plotCHM, fws, minht)
    treesPlot$treeid <- 0 # does this reset treeID?
    names(treesPlot)[names(treesPlot)=="x"] <- "xc"
    names(treesPlot)[names(treesPlot)=="y"] <- "yc"
    treesPlot <- treesPlot[,c(4,1,2,3)]
    for(j in 1:length(treesPlot$treeid)){
        ## Alter treeid name to include plotNames[[i]], "_", i
        treesPlot$treeid[[j]]<-paste0(plotNames[[i]], sprintf("_%03d",j))
    }
    ## Add all trees within given plot to treesxyList
    treesxyList <- rbind(treesxyList, treesPlot)
    print((paste0("Completed ITD for all trees in plot: ", plotNames[[i]], 
                  " (", i, " out of ", length(plotNames), ")")))
}; beep(11)

tidytrees_fws5x5m <- as.tibble(treesxyList) %>% 
    separate(treeid, c("plotid"), sep = "_", remove=FALSE) %>%
    relocate(plotid, .before = "treeid")
tidytrees_fws5x5m

# write.csv(tidytrees_fws5x5m,
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R50m_FWS5x5.csv"))
tidytrees_r50m_fws5x5 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R50m_FWS5x5.csv")) 

## Convert R50m ITD tree list to spatial points
ITD_r50m_fws5x5 <- SpatialPointsDataFrame(tidytrees_r50m_fws5x5[,3:4], 
                                          tidytrees_r50m_fws5x5)
coordinates(ITD_r50m_fws5x5)
proj4string(ITD_r50m_fws5x5) <- crs.ipnf  # define projection system of data
is.projected(ITD_r50m_fws5x5)
summary(ITD_r50m_fws5x5)
# plot(ITD_r50m_fws5x5)

## Crop from R50m to R11m for all plots
ITD_r11m_fws5x5 <- raster::crop(ITD_r50m_fws5x5, plotsClip_r11m)
# write.csv(ITD_r11m_fws5x5@data,
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS5x5.csv"))
tidytrees_r11m_fws5x5 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS5x5.csv")) %>%
    select(-X1)

##--------------------------------------------------------------------------##
## 7x7m Scale ----
##--------------------------------------------------------------------------##
## Parameters for ITD function {rLidar}
fws <- 7   #fixed window size = dimensions (m); default = 3 (3x3)
minht <- 2 #min height above ground for the detection break (m)
maxcrown=10.0 #maximum individual tree crown radius expected
exclusion=0.5 #% height value below maxht of pixels to exclude; 
##             default = 0.3 (30%)
## List of CHM rasters (R50m)
refRasters <- as.list(paste0(chm50m,"/", plotNames, ".asc"))
## Dataframe to add ITD trees to
treesxyList <- data.frame() 
## Run ITD method for all trees within 50m grid {rLidar}
for(i in 1:length(refRasters)){
    ## Load raster CHM for given plot (R50m)
    plotCHM <- raster(refRasters[[i]])
    ## Find individual trees within given plot
    treesPlot <- FindTreesCHM(plotCHM, fws, minht)
    treesPlot$treeid <- 0 # does this reset treeID?
    names(treesPlot)[names(treesPlot)=="x"] <- "xc"
    names(treesPlot)[names(treesPlot)=="y"] <- "yc"
    treesPlot <- treesPlot[,c(4,1,2,3)]
    for(j in 1:length(treesPlot$treeid)){
        ## Alter treeid name to include plotNames[[i]], "_", i
        treesPlot$treeid[[j]]<-paste0(plotNames[[i]], sprintf("_%03d",j))
    }
    ## Delineate boundaries for each tree crown
    # itdcrowns <- ForestCAS(plotCHM, loc, maxcrown, exclusion)
    # boundaryTrees <- itdcrowns[[1]]
    # crownList <- itdcrowns[[2]]
    # crownList$crad<-sqrt(crownList$ca/pi)
    # write.csv(crownList, paste0(tabs, "/", plotNames[[i]], 
    # "_ITDcrowns7x7.csv"))
    ## Add all trees within given plot to treesxyList
    treesxyList <- rbind(treesxyList, treesPlot)
    print((paste0("Completed ITD for all trees in plot: ", plotNames[[i]], 
                  " (", i, " out of ", length(plotNames), ")")))
}; beep(11)

tidytrees_fws7x7m <- as.tibble(treesxyList) %>% 
    separate(treeid, c("plotid"), sep = "_", remove=FALSE) %>%
    relocate(plotid, .before = "treeid")
tidytrees_fws7x7m

# write.csv(tidytrees_fws7x7m,
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS7x7.csv"))
tidytrees_r50m_fws7x7 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R50m_FWS7x7.csv")) 

## Convert R50m ITD tree list to spatial points
ITD_r50m_fws7x7 <- SpatialPointsDataFrame(tidytrees_r50m_fws7x7[,3:4], 
                                          tidytrees_r50m_fws7x7)
coordinates(ITD_r50m_fws7x7)
proj4string(ITD_r50m_fws7x7) <- crs.ipnf  # define projection system of data
is.projected(ITD_r50m_fws7x7)
summary(ITD_r50m_fws7x7)
# plot(ITD_r50m_fws7x7)

## Crop from R50m to R11m for all plots
ITD_r11m_fws7x7 <- raster::crop(ITD_r50m_fws7x7, plotsClip_r11m)

# write.csv(ITD_r11m_fws7x7@data, 
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS7x7.csv"))
tidytrees_r11m_fws7x7 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS7x7.csv")) %>%
    select(-X1)
##--------------------------------------------------------------------------##
## 9x9m Scale ----
##--------------------------------------------------------------------------##
## Parameters for ITD function {rLidar}
fws <- 9   #fixed window size = dimensions (m); default = 3 (3x3)
minht <- 2 #min height above ground for the detection break (m)
## List of CHM rasters (R50m)
refRasters <- as.list(paste0(chm50m,"/", plotNames, ".asc"))
## Dataframe to add ITD trees to
treesxyList <- data.frame() 
## Run ITD method for all trees within 50m grid {rLidar}
for(i in 1:length(refRasters)){
    ## Load raster CHM for given plot (R50m)
    plotCHM <- raster(refRasters[[i]])
    ## Find individual trees within given plot
    treesPlot <- FindTreesCHM(plotCHM, fws, minht)
    treesPlot$treeid <- 0 # does this reset treeID?
    names(treesPlot)[names(treesPlot)=="x"] <- "xc"
    names(treesPlot)[names(treesPlot)=="y"] <- "yc"
    treesPlot <- treesPlot[,c(4,1,2,3)]
    for(j in 1:length(treesPlot$treeid)){
        ## Alter treeid name to include plotNames[[i]], "_", i
        treesPlot$treeid[[j]]<-paste0(plotNames[[i]], sprintf("_%03d",j))
    }
    ## Add all trees within given plot to treesxyList
    treesxyList <- rbind(treesxyList, treesPlot)
    print((paste0("Completed ITD for all trees in plot: ", plotNames[[i]], 
                  " (", i, " out of ", length(plotNames), ")")))
}; beep(11)

tidytrees_fws9x9m <- as.tibble(treesxyList) %>% 
    separate(treeid, c("plotid"), sep = "_", remove=FALSE) %>%
    relocate(plotid, .before = "treeid")
tidytrees_fws9x9m

# write.csv(tidytrees_fws9x9m,
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R50m_FWS9x9.csv"))
tidytrees_r50m_fws9x9 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R50m_FWS9x9.csv")) 

## Convert R50m ITD tree list to spatial points
ITD_r50m_fws9x9 <- SpatialPointsDataFrame(tidytrees_r50m_fws9x9[,3:4], 
                                          tidytrees_r50m_fws9x9)
coordinates(ITD_r50m_fws9x9)
proj4string(ITD_r50m_fws9x9) <- crs.ipnf  # define projection system of data
is.projected(ITD_r50m_fws9x9)
summary(ITD_r50m_fws9x9)
# plot(ITD_r50m_fws9x9)

## Crop from R50m to R11m for all plots
ITD_r11m_fws9x9 <- raster::crop(ITD_r50m_fws9x9, plotsClip_r11m)
# write.csv(ITD_r11m_fws9x9@data, 
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS9x9.csv"))

tidytrees_r11m_fws9x9 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS9x9.csv")) %>%
    select(-X1)

##--------------------------------------------------------------------------##
# Process all ITD datasets to test accuracy against census data ----
##--------------------------------------------------------------------------##
## Clean up census (REF) data from USFS RMRS subplot surveys (32 of 53 plots)
qaqcREF <- qaqc %>% select(plotid, 
                           subplot.live_ALL,
                           subplot.live_20mHT,
                           subplot.live_30mHT) %>%
    left_join(plots, by="plotid") %>% 
    select(-c(6:11)) %>%
    relocate(site, .after = "plotid") %>%
    rename(
        REF_subALL = subplot.live_ALL,
        REF_sub20mHT = subplot.live_20mHT,
        REF_sub30mHT = subplot.live_30mHT
    ); qaqcREF

##--------------------------------------------------------------------------##
## Aggregate & tally all ITD-generated trees per survey plot
qaqcITD_3x3 <- tidytrees_r11m_fws3x3 %>%
    group_by(plotid) %>% tally() %>%
    rename("ITD_fwd3x3" = n)
qaqcITD_5x5 <- tidytrees_r11m_fws5x5 %>%
    group_by(plotid) %>% tally() %>%
    rename("ITD_fwd5x5" = n)
qaqcITD_7x7 <- tidytrees_r11m_fws7x7 %>%
    group_by(plotid) %>% tally() %>%
    rename("ITD_fwd7x7" = n)
qaqcITD_9x9 <- tidytrees_r11m_fws9x9 %>%
    group_by(plotid) %>% tally() %>%
    rename("ITD_fwd9x9" = n)

## Combine tallies across FWSizes
qaqcITD <- plyr::join_all(list(
    qaqcITD_3x3,qaqcITD_5x5,qaqcITD_7x7,qaqcITD_9x9), 
         by='plotid', type='left')
##--------------------------------------------------------------------------##
## Calculate differences between REF & ITD tree counts per plot (ERROR)
qaqcACC <- left_join(qaqcREF, qaqcITD, by='plotid') %>%
    mutate(
        ERROR_allv3 = ITD_fwd3x3 - REF_subALL,
        ERROR_allv5 = ITD_fwd5x5 - REF_subALL,
        ERROR_allv7 = ITD_fwd7x7 - REF_subALL,
        ERROR_allv9 = ITD_fwd9x9 - REF_subALL,
        ERROR_20mv3 = ITD_fwd3x3 - REF_sub20mHT,
        ERROR_20mv5 = ITD_fwd5x5 - REF_sub20mHT,
        ERROR_20mv7 = ITD_fwd7x7 - REF_sub20mHT,
        ERROR_20mv9 = ITD_fwd9x9 - REF_sub20mHT,
        ERROR_30mv3 = ITD_fwd3x3 - REF_sub30mHT,
        ERROR_30mv5 = ITD_fwd5x5 - REF_sub30mHT,
        ERROR_30mv7 = ITD_fwd7x7 - REF_sub30mHT,
        ERROR_30mv9 = ITD_fwd9x9 - REF_sub30mHT,
    ); qaqcACC

## Generate summary statistics on accuracy results across FWS & height filters
itdREPORT <- qaqcACC[,c(1:2,10:21)] %>%
    gather(key, value, ERROR_allv3:ERROR_30mv9 )%>%
    separate(key, c("census_filter", "fws"), sep="v") %>% 
    group_by(site, fws, census_filter) %>% 
    summarise(mean = mean(value), se = std.error(value)) %>%
    pivot_wider(names_from = census_filter, values_from = c(mean, se)) %>%
    rename("all_x" = mean_ERROR_all,
           "all_se" = se_ERROR_all,
           "20m_x" = mean_ERROR_20m,
           "20m_se" = se_ERROR_20m,
           "30m_x" = mean_ERROR_30m,
           "30m_se" = se_ERROR_30m) %>%
    relocate(site, fws, 
             "all_x", "all_se", "20m_x", "20m_se", "30m_x", "30m_se")


colMeans(itdREPORT[,c(3,5,7)])
itdREP_cda <- colMeans(itdREPORT[1:4,c(3,5,7)])
itdREP_stj <- colMeans(itdREPORT[5:8,c(3,5,7)])


itdERR <- qaqcACC[,c(1:3,10:21)] 
itdSUBSET <- itdERR %>% 
    select(c(plotid, REF_subALL, ERROR_30mv5, ERROR_30mv7, ERROR_20mv5)) %>%
    summarise(across(starts_with("ERROR"), 
                     list(mean = mean, se = std.error))) %>%
    rename(mean_30mv5 = ERROR_30mv5_mean,
           se_30mv5 = ERROR_30mv5_se,
           mean_30mv7 = ERROR_30mv7_mean,
           se_30mv7 = ERROR_30mv7_se,
           mean_20mv5 = ERROR_20mv5_mean,
           se_20mv5 = ERROR_20mv5_se); itdSUBSET

itdSUBSET_bySite <- itdERR %>% 
    group_by(site) %>%
    select(c(plotid, REF_subALL, ERROR_30mv5, ERROR_30mv7, ERROR_20mv5)) %>%
    summarise(across(starts_with("ERROR"), 
                     list(mean = mean, se = std.error))) %>%
    rename(mean_30mv5 = ERROR_30mv5_mean,
           se_30mv5 = ERROR_30mv5_se,
           mean_30mv7 = ERROR_30mv7_mean,
           se_30mv7 = ERROR_30mv7_se,
           mean_20mv5 = ERROR_20mv5_mean,
           se_20mv5 = ERROR_20mv5_se); itdSUBSET_bySite

##--------------------------------------------------------------------------##
## Select most accurate data and crop from R50m to R25m for all plots
##--------------------------------------------------------------------------##
## Crop from R50m to R11m for all plots
ITD_r25m_fws7x7 <- raster::crop(ITD_r50m_fws7x7, plotsClip_r25m)
# write.csv(ITD_r11m_fws7x7@data,
#           file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS7x7.csv"))

tidytrees_r11m_fws7x7 <- read_csv(
    file.path(tabs, "snag-gaps_ITDlive-acc_ipnf2017_R11m_FWS7x7.csv")) %>%
    select(-X1)
