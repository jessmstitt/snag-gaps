## HEADER --------------------------------------------------------------------
##
## Script name: Generating plots and graphics for snag gap statistics
##      [04_snag-gaps_graphics.R]
##
## Purpose of script:
##
## Author: Jessica M. Stitt
##
## Date Created: 2021-05-22
##
## Copyright (c) Jessica M. Stitt, 2021
## Email: jessmstitt@gmail.com
##
## NOTES ---------------------------------------------------------------------
##
##   
##
## ---------------------------------------------------------------------------
##--------------------------------------------------------------------------##
## LOAD PACKAGES ----
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
library(rstatix)
library(ggpubr)
library(gridExtra)
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
## Read in datasets
allGaps_clean <- read.csv(
    file.path(tabs, "snag-gaps_ipnf2017_full-dataset-clean.csv"))

gapStats_analysis <- read.csv(
    file.path(tabs, "snag-gaps_gapstats-analysis-bytree_ipnf2017.csv"))
gapFigs <- as_tibble(gapStats_analysis) %>%
    convert_as_factor(gapht_m)
gapFigs <- gapFigs[,c(2:5,8:9)]
str(gapFigs)
summary(gapFigs)

gapStats_diffs <- read.csv(
    file.path(tabs, "snag-gaps_gapstats-analysis-propdiff_ipnf2017.csv"))
gapDiff <- as_tibble(gapStats_diffs) %>%
    convert_as_factor(gapht_m)
gapDiff <- gapDiff[,2:5]
str(gapDiff)
summary(gapDiff)

gaps <- allGaps_clean %>%
    convert_as_factor(fp_m2, gapht_m)
str(gaps)
summary(gaps)

summstats <- group_by(gaps, tree_type, site) %>%
    summarise(
        count.tree = n()/75,
        mean.gap.prop = mean(prop_gap, na.rm = TRUE),
        se = std.error(prop_gap, na.rm = TRUE)
    )

gaps %>%
    group_by(tree_type) %>%
    get_summary_stats(prop_gap, type="full") #median & interquartile range


summgaps <- 
    gaps %>%
    convert_as_factor(gapht_m, fp_m2) %>%
    group_by(tree_type, gapht_m, site, fp_m2) %>%
    summarise(
        gapmu = mean(prop_gap, na.rm = TRUE),
        se = std.error(prop_gap, na.rm = TRUE)) 

setrees <- 
    gaps %>%
    convert_as_factor(gapht_m, fp_m2) %>%
    group_by(tree_type, gapht_m, site, fp_m2) %>%
    summarise(
        se = std.error(prop_gap, na.rm = TRUE)) %>%
    spread(tree_type, se)
setrees <- setrees[,4:5]

difftrees <- 
    gaps %>%
    convert_as_factor(gapht_m, fp_m2) %>%
    group_by(tree_type, gapht_m, site, fp_m2) %>%
    summarise(
        gapmu = mean(prop_gap, na.rm = TRUE)) %>%
    spread(tree_type, gapmu) %>%
    mutate(
        diff = snag - live,
        pct.diff = snag - live, pct.diff = scales::percent(pct.diff))

diffsbysite <- 
    summgaps %>%
    group_by(gapht_m, site, fp_m2) %>%
    dplyr::select(-c(live, snag, pct.diff)) %>%
    spread(site, diff) %>%
    mutate(
        site.diff = cda - stj)

diffsbysite %>%
    group_by(fp_m2) %>%
    summarise(avgsitediff = mean(site.diff))

avgdiffspersite <- 
    summgaps %>%
    group_by(gapht_m, site, fp_m2) %>%
    # select(-c(live, snag, pct.diff)) %>%
    group_by(site) %>%
    summarise(avgdiff = mean(diff))

avgdiffsperfp <- 
    difftrees %>%
    # group_by(gapht_m, site, fp_m2) %>%
    # select(-c(live, snag, pct.diff)) %>%
    group_by(fp_m2) %>%
    summarise(avgdiff = mean(diff))

sitediffXht <- gaps %>%
    convert_as_factor(gapht_m, fp_m2) %>%
    group_by(tree_type, gapht_m, site, fp_m2) %>%
    summarise(
        gapmu = mean(prop_gap, na.rm = TRUE)) %>%
    spread(site, gapmu) %>%
    mutate(sitediff = cda - stj) %>%
    group_by(tree_type, gapht_m) %>%
    summarise(avgsitediff = mean(sitediff)) %>%
    spread(tree_type, avgsitediff)

##--------------------------------------------------------------------------##
beep(2) 
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
plotting <- gapFigs
# New facet label names for site variable
site.labels <- c("CdA site", "StJ site")
names(site.labels) <- c("CdA", "StJ")

# New facet label names for fp variable
fp.labels <- c("16m2", "36m2", "64m2")
names(fp.labels) <- c("16", "36", "64")

## Visualize canopy gap area across heights: between sites + tree type + fps
bysite_vert <-
    ggplot(plotting, aes(y = as.numeric(gapht_m)*2,
                         x = prop_mu,
                         group = interaction(type,fp_m2),
                         color = type)) + 
    facet_wrap(~site, ncol = 2) +
    geom_hline(yintercept = c(10,20,30,40,50), show.legend=F,
               color = "grey", size=0.4, lty="dotted") +
    geom_vline(xintercept = c(0,0.25,0.5,0.75,1), show.legend=F,
               color = "grey", size=0.4, lty="dotted") +
    # geom_vline(xintercept = c(1), show.legend=F,
    #           color = "blue", size=1.8, alpha=0.2) +
    scale_color_manual(values = c("#1B9E77", "#D95F02"),
                       name = "Tree Type",
                       labels = c("Live Tree (CdA n=1187 | StJ n=999)",
                                  "Snag (CdA n=180 | StJ n=90)")
    ) +
    geom_line(lwd = 0.6, alpha = 0.6) + 
    geom_point(size = 2, aes(shape=fp_m2)) + 
    scale_shape_discrete(name = expression(bold("Footprint Size "~(m^2), "")),
                         labels = c("Small (16)","Medium (36)", "Large (64)")) +
    geom_errorbar(aes(xmin = prop_mu - prop_se,
                      xmax = prop_mu + prop_se), width = 0.6,
                  position = position_dodge(0.05)) +
    xlim(0,1) + ylim(2,50) +
    labs(
        # title = "Average proportion of canopy gap area (+/- se) 
        # across forest canopy heights for snags vs. live trees",
        # subtitle = expression(italic(
        #   "Grouped by tree type & footprint size, split by site; 
        #   single tree-level (STL)")),
        y = "Forest canopy height (m)",
        x = "Mean (±s.e.) canopy gap proportion per unit area") +
    # theme_pubclean()
    theme_pubr(border = T, legend = "bottom") +
    labs_pubr(base_size = 14) + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.justification=c(0.5,0), legend.position=c(0.84,0.02))
bysite_vert

# ggsave(plot = bysite_vert, "snag-gaps_graphic_typeXsite-FIG.pdf",
#        device="pdf", width=36, height=24, units="cm",
#        path = figs)
##--------------------------------------------------------------------------##
# DIFFERENCE PLOTS
diff_bysite_vert <-
    ggplot(gapDiff, aes(x = as.numeric(gapht_m)*2,
                        y = diffmu_propcga,
                        group = fp_m2,
                        color = fp_m2)) + 
    facet_wrap(~site, ncol=2) +
    geom_vline(xintercept = c(0,10,20,30,40,50), show.legend=F,
               color = "grey", size=0.4, lty="dotted") +
    geom_hline(yintercept = c(0,0.1,0.2,0.3,0.4, 0.5), show.legend=F,
               color = "grey", size=0.4, lty="dotted") +
    geom_line(lwd = 0.8, alpha = 0.8) + 
    geom_point(size = 2, aes(shape=fp_m2)) + 
    # scale_shape_discrete(name = expression(bold("Footprint Size "~(m^2), "")),
    #                      labels = c("Small (16)","Medium (36)", 
    #                      "Large (64)")) +
    scale_color_discrete(name = expression(bold("Footprint Size "~(m^2), "")),
                         labels = c("Small (16)","Medium (36)", "Large (64)")) +
    scale_shape_discrete(name = expression(bold("Footprint Size "~(m^2), "")),
                         labels = c("Small (16)","Medium (36)", "Large (64)")) +
    ylim(0,0.5) + xlim(2,50) +
    labs(
        # title = "Differences in mean gap area between snags & live trees 
        # across canopy heights",
        #      subtitle = expression(italic("Grouped by study site, 
        #      single tree-level (STL)")),
        x = "Forest canopy height (m)",
        y = "Difference in mean canopy gap proportion (snag - live)") +
    theme_pubr(border=T) +
    labs_pubr(base_size = 14)+ 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.justification=c(0.5,0), legend.position=c(0.84,0.04)) +
    coord_flip()
diff_bysite_vert

# ggsave(plot = diff_bysite_vert, "snag-gaps_graphic_diffXsite-FIG.pdf",
#        device="pdf", width=36, height=24, units="cm",
#        path = figs)
##--------------------------------------------------------------------------##
