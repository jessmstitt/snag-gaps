## HEADER --------------------------------------------------------------------
##
## Script name: Lidar pre-processing steps
##    [01_snag-gaps_als-proc.R]
## 
## Purpose of script: Process airborne lidar data acquisitions 
##    Inputs: raw point clouds from IPNF (2016) + field survey data (2017)
##    Outputs: height-normalized clips of survey plots + rasterized CHMs +
##             lidar-derived metrics + DEMs for topography
## 
## Author: Jess M. Stitt
## Date created: 2021-05-22
## Email: jessmstitt@gmail.com
## 
## NOTES ---------------------------------------------------------------------
## 
##   
## 
## begin script --------------------------------------------------------------

##--------------------------------------------------------------------------##
# LOAD PACKAGES ----
##--------------------------------------------------------------------------##
## Lidar analysis in R
library(lidR)       #for lidar data manipulation**
## **If not yet installed, first install with devtools
# library(devtools)
# devtools::install_github("Jean-Romain/lidR")
##--------------------------------------------------------------------------##
## Spatial manipulation
library(raster)
library(rgdal) #for vector work
library(sf)
library(sp)
##--------------------------------------------------------------------------##
## Data visualization & plotting
library(ggplot2)    #for plots & graphics
library(viridis)    #for color palettes
##--------------------------------------------------------------------------##
## Workflow & data organization
library(here)       #for filepath mgmt
library(tidyverse)  #for 'tidy' data
library(beepr)      #for [unnecessary] sound notifications
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
### Raw lidar datasets [not stored within project files**]
## **Access to full lidar point clouds can be found at: 
## <https://doi.org/10.3334/ORNLDAAC/1766> 
## 
### COEUR D'ALENE Study Site (cda) lidar
# rawLAS_cda = file.path(path/to/file/BuncoRoushRossi2016.las)
### SAINT JOE Study Site (stj) lidar
# rawLAS_stj = file.path(path/to/file/StJoe.las)
##--------------------------------------------------------------------------##
## Processing folders
### Lidar processing folders
las50m = here(path=laspc,"00_las_clip-r50m")#raw point cloud R50m
ground = here(path=laspc,"01_las_ground")   #ground filter R50m
dem = here(path=laspc,"02_las_dem")         #digital elevation model R50m
znorm = here(path=laspc,"03_las_znorm")     #height-normalized R50m
znf = here(path=laspc,"04_las_znfil")       #znorm filtered (>0, C==1L) R50m
clipPlot = here(path=laspc,"05_las_clip-plot") #lasclipCircle from znf to R25m
clipRMRS = here(path=laspc,"06_las_clip-rmrs") #lasclipCircle from znf to R11.3m
##--------------------------------------------------------------------------##
### Raster processing folders
pfc = here(path=gis,"01_tif_pitfree-chm")   #pit-free chm R50m around plots
chm50m = here(path=gis,"02_asc_chm-r50m")      #chm raster R50m around plots
chmPlot = here(gis,"03_asc_chm-plot")       #chm raster from chm to R25m
rasterRMRS = here(gis,"04_asc_chm-rmrs")      #chm raster from chm to R11.3m
##--------------------------------------------------------------------------##
## Scripts & outputs folders
scripts = here(path="02_scripts")           #R scripts
lasmx = here(path="03_results", "03_las-metrics") #lidar-derived metrics
##--------------------------------------------------------------------------##
beep(1)
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
# READ IN DATASETS ----
##--------------------------------------------------------------------------##
## Load relevant field data: IPNF 2017
## Data on 25m-radius survey plots, across sites
plots <- read.csv(file.path(csv, "snag-gaps_REFplots_ipnf2017.csv")) #n=53
## Data on USFS RMRS 11.3m-radius surveys, within a subset of survey plots
rmrs <- read.csv(file.path(csv, "snag-gaps_subplots-r11m_ipnf2017.csv"))#n=32
##--------------------------------------------------------------------------##
## Load buffers around plots (R50m, R25m circles)
plotsClip_r50m = shapefile((paste0(gis, "/", "00_shp", "/",
                                   "plotsCLIP-r50m_ipnf2017.shp"))) 
plotsClip_r25m = shapefile((paste0(gis, "/", "00_shp", "/",
                                   "plotsCLIP-r25m_ipnf2017.shp"))) 
## Was easier to clip all plots at R11.3m (for QAQC steps) vs only sampled
plotsClip_r11m = shapefile((paste0(gis, "/", "00_shp", "/",
                                   "plotsCLIP-r11m_ipnf2017.shp")))
# plotsClip_r11m = shapefile((paste0(gis, "/", "00_shp", "/",
#                                    "subplotsCLIP-r11m_ipnf2017.shp")))
##--------------------------------------------------------------------------##
## Assign names, based on field survey plot IDs; and number of plots
plotNames <- as.character(plots$plotid)
nplots <- length(plotNames)
##--------------------------------------------------------------------------##
beep(2) 
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
# LIDAR PROCESSING ----
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## 00. Clip raw las pointclouds for REF plots (R50m) ----
##--------------------------------------------------------------------------##
## Isolate regions of interest (ROIs): subset REF plots+buffer (R50m) from tiles
ctgR50m <- readLAScatalog(rawLAS) #catalog files of raw las from cda
opt_output_files(ctgR50m) <- (paste0(las50m,"/{plotid}")) #set catalog options
plotLAS_R50m <- lasclip(ctgR50m, plotsClip_r50m) #clip just the REF plots R50m
ctgROI <- readLAScatalog(las50m)    #catalog all files in folder
lidR:::catalog_laxindex(ctgROI)     #index all files & create lax files
las_check(ctgROI)

##--------------------------------------------------------------------------##
## 01. Ground points only for REF plots (R50m) ----
##--------------------------------------------------------------------------##
opt_output_files(ctgROI) <- (paste0(ground,"/{ORIGINALFILENAME}"))
ctgGRD <- classify_ground(ctgROI, algorithm = csf())
las_check(ctgGRD)
lidR::filter_duplicates(ctgGRD)
lidR:::catalog_laxindex(ctgGRD)
las_check(ctgGRD)
##--------------------------------------------------------------------------##
### Check quality with sample plot
classcheck <- lidR::readLAS(paste0(ground, "/BR210B.las"))
las_check(classcheck)
plot(classcheck, color = "Classification", size = 3, bg = "white")
## **Function for 2D cross-section view of point cloud
##  from: https://jean-romain.github.io/lidRbook/io.html#plot-crossection
plot_crossection <- function(las,
                             p1 = c(min(las@data$X), mean(las@data$Y)),
                             p2 = c(max(las@data$X), mean(las@data$Y)),
                             width = 4, colour_by = NULL)
{
    colour_by <- enquo(colour_by)
    data_clip <- clip_transect(las, p1, p2, width)
    p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 0.5) +
         coord_equal() + theme_minimal()

    if (!is.null(colour_by))
        p <- p + aes(color = !!colour_by) + labs(color = "")

    return(p)
}
plot_crossection(classcheck, colour_by = factor(Classification))  +
    theme_bw()
##--------------------------------------------------------------------------##
### Check quality with sample plot
groundcheck <- filter_ground(classcheck)
plot(groundcheck, size = 3, bg = "white")
plot_crossection(groundcheck, colour_by = factor(Classification))  +
    theme_bw()

##--------------------------------------------------------------------------##
## 02. DEM for REF plots (R50m) ----
##--------------------------------------------------------------------------##
opt_output_files(ctgROI) <- (paste0(dem,"/{ORIGINALFILENAME}"))
demR50m <- grid_terrain(ctgROI, algorithm = tin())
##--------------------------------------------------------------------------##
### Check quality with sample plot
demcheck <- paste0(dem, "/BR210B.tif") #cda; HS221B for stj
r1 <- raster(demcheck); plot(r1)

##--------------------------------------------------------------------------##
## 03. Create z-normalized las pointclouds for REF plots (R50m) ----
##--------------------------------------------------------------------------##
opt_output_files(ctgGRD) <- (paste0(znorm,"/{ORIGINALFILENAME}_zn"))
htnormalized <- normalize_height(ctgGRD, tin())
ctgZN <- readLAScatalog(znorm)
lidR:::catalog_laxindex(ctgZN)
las_check(ctgZN)
##--------------------------------------------------------------------------##
### Check quality with sample plot
zncheck <- lidR::readLAS(paste0(znorm, "/BR210B_zn.las"))
plot(zncheck, size = 2, bg = "white")
plot_crossection(zncheck, colour_by = factor(Classification))

##--------------------------------------------------------------------------##
## 04. Filter "noise" from las pointclouds for REF plots (R50m) ----
##--------------------------------------------------------------------------##
## Loop through to filter las files from z-normalized folder
lasList_ref <- list.files(paste0(znorm, "/"), pattern="*.las")
refList <- list()
for(i in lasList_ref){
    las <- lidR::readLAS(paste0(znorm,"/",i))   #read each las file
    f = filter_poi(las, Z > 0, Z <= 70,
                   Classification==1L)          #filter pointcloud
    refList[[i]]<-f                             #add to list of las+zn+filter
}
plotsZNF <- as.list(refList)
plotsZNF <- paste0(znf,"/",plotNames,"_znf.las")#replace filenames
for(i in 1:length(plotsZNF)){
    writeLAS(refList[[i]],file=plotsZNF[[i]])   #write LAS files from list
}; beep(5)
names(refList) <- plotNames                     #update refList with plot IDs

ctgZNF <- readLAScatalog(znf) 
lidR:::catalog_laxindex(ctgZNF)
las_check(ctgZNF)
##--------------------------------------------------------------------------##
### Check quality with sample plot
znfcheck <- lidR::readLAS(paste0(znf, "/BR210B_znf.las"))
plot(znfcheck, size = 2, bg = "white")
plot_crossection(znfcheck, colour_by = factor(Classification))
##--------------------------------------------------------------------------##
lasList_znf <- list.files(paste0(znf, "/"), pattern="*.las")
znfList <- list()
for(i in lasList_znf){
    lasznf <- lidR::readLAS(paste0(znf,"/",i))     #read each las file
    znfList[[i]]<-lasznf                         #add to list of las+zn+filter
}; beep(5)
plotsZNF <- as.list(znfList)
lidR::plot(plotsZNF[[1]], color="Z", axis=T)

##--------------------------------------------------------------------------##
## 05. Clip pointclouds to only REF survey plot extent (R25m) -----
##--------------------------------------------------------------------------##
ctgPlotClip <- readLAScatalog(znf)
lidR:::catalog_laxindex(ctgPlotClip)
opt_independent_files(ctgPlotClip) <- TRUE
opt_chunk_buffer(ctgPlotClip) <- 0
opt_chunk_size(ctgPlotClip)   <- 0
opt_output_files(ctgPlotClip) <- (paste0(clipPlot,"/{plotid}"))
plotLAS_R25m = clip_roi(ctgPlotClip, plotsClip_r25m)
ctgPlots <- readLAScatalog(clipPlot)  #catalog all files in folder
lidR:::catalog_laxindex(ctgPlots)     #index all files & create lax files
lascheck(ctgPlots)                    #inspect the files in catalog
##--------------------------------------------------------------------------##
### Check quality with sample plot
lasList_R25m <- list.files(paste0(clipPlot, "/"), pattern="*.las")
R25mList <- list()
plotcheckList <- list()
for(i in lasList_R25m){
    las25m <- lidR::readLAS(paste0(clipPlot,"/",i))     #read each las file
    R25mList[[i]]<-las25m                         #add to list of las+zn+filter
    plotcheckList[[i]] <-
        plot_crossection(las25m, colour_by = factor(Classification))
}
plotsR25m <- as.list(R25mList)
lidR::plot(plotsR25m[[1]], color="Z", axis=T)
names(plotcheckList) <- plotNames 
plotcheckList[[1]]

##--------------------------------------------------------------------------##
## 06. Clip subset of survey plots to USFS RMRS subplots (R11.3m) -----
##--------------------------------------------------------------------------##
ctgRMRSClip <- readLAScatalog(znf)
lidR:::catalog_laxindex(ctgRMRSClip)
opt_chunk_buffer(ctgRMRSClip) <- 0
opt_chunk_size(ctgRMRSClip)   <- 0
opt_output_files(ctgRMRSClip) <- (paste0(clipRMRS,"/{plotid}_sp-r11m"))
plotLAS_R11m = clip_roi(ctgRMRSClip, plotsClip_r11m)
ctgRMRS <- readLAScatalog(clipRMRS)  #catalog all files in folder
lidR:::catalog_laxindex(ctgRMRS)     #index all files & create lax files
lascheck(ctgRMRS)                    #inspect the files in catalog
##--------------------------------------------------------------------------##
### Check quality with sample plot
lasList_R11m <- list.files(paste0(clipRMRS, "/"), pattern="*.las")
R11mList <- list()
QAQCcheckList <- list()
for(i in lasList_R11m){
    las11m <- lidR::readLAS(paste0(clipRMRS,"/",i))     #read each las file
    R11mList[[i]]<-las11m                         #add to list of las+zn+filter
    QAQCcheckList[[i]] <-
        plot_crossection(las11m, colour_by = factor(Classification))
}
plotsR11m <- as.list(R11mList)
lidR::plot(plotsR11m[[1]], color="Z", axis=T)
subplotNames <- subset(plots, subplot=="y")
names(QAQCcheckList) <- subplotNames 
QAQCcheckList[[1]]
##--------------------------------------------------------------------------##
beep(1)
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
# CHM PROCESSING ----
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## 01. Create pit-free canopy height models for grids (R50m) <TIF> ----
##--------------------------------------------------------------------------##
## Set up: catalog details
ctgCHM <- readLAScatalog(znf)
lidR:::catalog_laxindex(ctgCHM)
opt_output_files(ctgCHM) <- (paste0(pfc, "/{ORIGINALFILENAME}_pfc"))
## Establish algorithm parameters & run CHM function
algo = pitfree(thresholds = c(0,2,5,10,20,50), subcircle = 0.2)
canopy = grid_canopy(ctgCHM, 0.5, algo)
lidR::plot(canopy, col = viridis(10))
summary(ctgCHM)

##--------------------------------------------------------------------------##
## 02. Create rasters of CHMs for grids (R50m) <ASCII> ----
##--------------------------------------------------------------------------##
## Pull pit-free CHMs from folder and change format
pfcList <- list.files(paste0(pfc, "/"), pattern = "*.tif")
pfcList <- paste0(pfc, "/", pfcList)
list_plotRasters <- list()
## Loop through CHMs to turn into rasters
for(i in pfcList) {
    f <- i
    plotCHM_R50m <- raster(f)
    plot(plotCHM_R50m, col = viridis(10))
    writeRaster(plotCHM_R50m, paste0(chm50m, "/", 
                                     str_sub(i, start = -18, end = -13),
                          ".asc"), format = "ascii") #get plotIDs
    list_plotRasters[[i]] <- plotCHM_R50m
}; beep(2)
names(list_plotRasters) <- plotNames      #update refList with plot IDs

##--------------------------------------------------------------------------##
## 03. Clip rasters to only REF survey plot extent (R25m) ----
##--------------------------------------------------------------------------##
## Loop through CHMs to turn into rasters
for(i in pfcList) {
    f <- i
    plotCHM_R50m <- raster(f)
    r <- mask(plotCHM_R50m, plotsClip_r25m)
    plot(r, col = viridis(10))
    writeRaster(r, paste0(chmPlot, "/", 
                          str_sub(i, start = -18, end = -13),
                          ".asc"), format = "ascii") #get plotIDs
    list_plotRasters[[i]] <- r
}; beep(2)
names(list_plotRasters) <- plotNames      #update refList with plot IDs
beep(2) 

##--------------------------------------------------------------------------##
## 04. Clip subset of CHMs to USFS RMRS subplots (R11.3m) -----
##--------------------------------------------------------------------------##
## Loop through CHMs to turn into rasters
lasList_R11m <- list.files(paste0(clipRMRS, "/"), pattern="*.las")
R11mList <- list()
QAQCcheckList <- list()

rmrsList <- list.files(paste0(pfc, "/"), pattern = "*.tif")
pfcList <- list.files(paste0(pfc, "/"), pattern = "*.tif")
pfcList <- paste0(pfc, "/", pfcList)
list_plotRasters <- list()

list_subplotRasters <- list()
for(i in pfcList) {
    f <- i
    plotCHM_R50m <- raster(f)
    r <- mask(plotCHM_R50m, plotsClip_r11m)
    plot(r, col = viridis(10))
    writeRaster(r, paste0(rasterRMRS, "/", 
                          str_sub(i, start = -18, end = -13),
                          ".asc"), format = "ascii") #get plotIDs
    list_subplotRasters[[i]] <- r
}; beep(2)
names(list_subplotRasters) <- plotNames      #update refList with plot IDs
beep(4) 
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
# LIDAR METRICS ----
##--------------------------------------------------------------------------##
## Calculate lidar-derived metrics at survey plot level (R25m)
plotLAS_list <- list.files(paste0(clipPlot, "/"), pattern="*.las")

plotMx <- lapply(plotLAS_list, function(file) {
    las <- lidR::readLAS(paste0(clipPlot,"/",file))
    cloud_metrics(las, .stdmetrics)
})
plotMx <- data.table::rbindlist(plotMx)
plotsMx <- as.tibble(cbind(plotid=c(plotNames), plotMx))
write.csv(plotsMx, file.path(lasmx,
                             "snag-gaps_REFplots_cloud-metrics_ipnf2017.csv"))
##--------------------------------------------------------------------------##
beep(4) 
##--------------------------------------------------------------------------##

## end script ----------------------------------------------------------------