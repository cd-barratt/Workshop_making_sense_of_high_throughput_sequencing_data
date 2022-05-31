library(grImport)
library(Rcpp) 
library(RcppEigen)
library(geosphere)
library(raster)
library(rgeos)
library(sp)
library(rgdal)
library(rworldmap)
library(rworldxtra)
library(RColorBrewer)
#library(devtools)
#install_github("dipetkov/eems/plotting/rEEMSplots")
library(rEEMSplots)


##########################################
########### Plot png files ##############
#Look at results across chains

Sys.setenv(R_GSCMD = '/usr/local/Cellar/ghostscript/9.27_1/bin/gs')
        
              
setwd=("/Users/chris/Desktop/Workshop_2/Outputs/EEMS/")
mcmcpaths = c("/Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/Outputs/EEMS/Lflav_nDemes700-chain1/","/Users/chris/Desktop/Workshop_2/Outputs/EEMS/Lflav_nDemes700-chain2/")
name.figures.to.save = "/Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/Outputs/EEMS/Lflav-nDemes700-PLOTS"

eems.plots(mcmcpath = mcmcpaths, plotpath = paste(name.figures.to.save,"-default",sep=""), longlat = FALSE,
           add.demes = TRUE, add.map = TRUE, projection.in = "+proj=longlat +datum=WGS84", projection.out = "+proj=merc +datum=WGS84", out.png = TRUE)

##########################################
########### Plot pdf files ##############
#Look at results across chains

eems.plots(mcmcpath = mcmcpaths, plotpath = paste(name.figures.to.save,"-default",sep=""), longlat = FALSE,
           add.demes = TRUE, add.map = TRUE, projection.in = "+proj=longlat +datum=WGS84", projection.out = "+proj=merc +datum=WGS84", out.png = FALSE)

##########################################
########### Plot png files with grids ##############
name.figures.to.save = "/Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/Outputs/EEMS/Lflav-nDemes700-PLOTS_grid"

eems.plots(mcmcpath = mcmcpaths, plotpath = paste(name.figures.to.save,"-default",sep=""), longlat = FALSE,
           add.demes = TRUE, add.map = TRUE, projection.in = "+proj=longlat +datum=WGS84", projection.out = "+proj=merc +datum=WGS84", out.png = TRUE, add.grid = TRUE, lwd.grid = 0.5)



######################################
#set color scheme

color_pal2= c("White","#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8" ,"#807DBA", "#6A51A3" ,"#54278F")
color_pal2

name.figures.to.save = "/Users/chris/Desktop/Workshop_2/work/Lflavomaculatus/Outputs/EEMS/Lflav-nDemes700-C2PLOTS"

eems.plots(mcmcpath = mcmcpaths, plotpath = paste(name.figures.to.save,"-default",sep=""), longlat = FALSE,
           add.demes = TRUE, add.map = TRUE, projection.in = "+proj=longlat +datum=WGS84", projection.out = "+proj=merc +datum=WGS84", out.png = TRUE, eems.colors = color_pal2)


