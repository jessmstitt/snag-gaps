# snag-gaps README
 
Project: Forest canopy gaps around individual snags vs live trees

Purpose: Overview of the structure and files for the snag-gaps
 project, built in tandem with a manuscript entitled
      "Characterizing individual tree-level snags using airborne 
      lidar-derived forest canopy gaps within closed-canopy conifer forests"
      submitted to Methods in Ecology & Evolution (MEE), March 2021

Corresponding author: Jess M. Stitt
Email: jessmstitt [at] gmail.com

This repository contains all files related to a project investigating forest canopy gap dynamics around individual standing dead trees (snags) vs. live trees in a closed-canopy conifer forest. 

# *Project objectives*
 - Determine if forest canopy gaps can be used to improve snag detection / prediction in forests
 - Focus on two aspects: 
    - height in canopy where gaps occur & gap size for snags vs. live trees
    - extent of footprint size around all individual trees 

Canopy gaps are those areas in the overstory where trees are absent - often times they are created by treefalls, but may also represent still-standing dead trees, with the gap primarily composed of the area previously taken up by the live canopy of the now-dead snag (yet still containing the standing dead trunk and possibly some branches).

# *Project structure*
The files for this project are organized into folders, grouped into four (4) subfolders:
 - 01_datasets
 - 02_scripts
 - 03_results
 - 04_reports
    
## 01_datasets
Within this folder are the raw data used for this project, including:
 - 01_csv = spreadsheets (.csv) for plots, snags, and subplots
 - 02_las = airborne lidar data (**not included in this repository**)
 - 03_gis = all files for spatial information, further subdivided:
    - 00_shp = shapefiles (.shp) for the study sites and survey plots
    - 01:07 folders = raster processing steps; only step included here is 
    - 02_asc_chm-r50m = 50m-radius clips of pit-free canopy height models (chm) around each survey plot (n=53), at 0.5m resolution

## 02_scripts
Within this folder are the four (4) R scripts (named following a template of ORDER_PROJECT_FUNCTION):
 - 01_snag-gaps_als-proc.R = processing raw lidar data into chm for each plot
 - 02_snag-gaps_itd-qaqc.R = individual tree detection (itd) accuracy tests
 - 03_snag-gaps_fcgap-calc.R = individual tree-level gap detection & analyses
 - 04_snag-gaps_graphics.R =  generation of figures based on fc gap results

## 03_results
Within this folder are all results from processing done within the project:
 - 01_tables = csv files including the itd-generated live tree locations, gap statistics, and full itd accuracy test datasets
 - 02_gis = any shapefiles of gap boundaries (**not currently included**)
 - 03_las-metrics = lidar-derived metrics at the plot & tree level
 - 04_figures = graphics generated based on the fc gap results

## 04_reports
Within this folder are any larger "reports" including R Markdown files, maps created within QGIS, and additional figures

## **R Code: going from forest-wide canopy height model to canopy gap statistics for individual trees**
Much of the code contained within this project derives from the [ForestGapR github page](https://github.com/carlos-alberto-silva/ForestGapR) and other sources online, such as [the lidR package book](https://jean-romain.github.io/lidRbook/index.html) and [a Quantitative Ecology tutorial](http://quantitativeecology.org/using-rlidar-and-fusion-to-delineate-individual-trees-through-canopy-height-model-segmentation/). 

## Data Source & Access
Data were collected in large part by USFS RMRS, which provided the airborne lidar data of the Idaho Panhandle National Forests and conducted the subplot surveys across both sites.
