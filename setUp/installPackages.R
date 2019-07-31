#!/usr/bin/Rscript

#-------------------------------------------------------
#
# sourcable files
#
#-------------------------------------------------------
PredictingProteinInteractionsFolder = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/"


addSource <- function(table, name, fName){
  sources2 = data.frame(matrix(0,ncol = 2, nrow = 1))
  colnames(sources2) = c("name", "path")
  
  sources2[1, ] = c(name,paste(PredictingProteinInteractionsFolder,fName, sep = ""))
  
  table = rbind(table, sources2)
  return(table)
}

sources = data.frame(matrix(0,ncol = 2, nrow = 0))
colnames(sources) = c("name", "path")

sources = addSource(sources, "helperFunctions", "/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R")
sources = addSource(sources, "UltraQuickRepeatedSubSampling", "/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R")
sources = addSource(sources, "TriangulateIsoSurface", "/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R")
sources = addSource(sources, "QuickRepeatedSubSampling", "/MetricGeometry/QuickRepeatedSubSampling/QuickRepeatedSubSampling.R")


write.table(sources, paste(PredictingProteinInteractionsFolder, "setUp/SourcableFiles.txt", sep = ""), row.names = FALSE)
#-------------------------------------------------------
#
# install packages
#
#-------------------------------------------------------


is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}

if(!is.installed("plot3D")){install.packages("plot3D")}
if(!is.installed("funr")){install.packages("funr")}
if(!is.installed("emdist")){install.packages("emdist")}
if(!is.installed("ggplot2")){install.packages("ggplot2")}
if(!is.installed("ggdendro")){install.packages("ggdendro")}
if(!is.installed("plot3D")){install.packages("plot3D")}
if(!is.installed("tcltk")){install.packages("tcltk")}
if(!is.installed("getopt")){install.packages("getopt")}
if(!is.installed("rfUtilities")){install.packages("rfUtilities")}
if(!is.installed("rbenchmark")){install.packages("rbenchmark")}


if(!is.installed("raster")){install.packages("raster")}
if(!is.installed("foreach")){install.packages("foreach")}
if(!is.installed("doMC")){install.packages("doMC")}
if(!is.installed("FNN")){install.packages("FNN")}
if(!is.installed("RSNNS")){install.packages("RSNNS")}


# misc3d
# geometry
# xtable
# spatstat
# igraph
# gplots
# lubridate


library(rgl)

library(readobj)




# R.Version()
# install.packages("installr")
# library(installr)
# updateR()



# for yupiter ...
# https://www.datacamp.com/community/blog/jupyter-notebook-r?utm_source=adwords_ppc&utm_campaignid=898687156&utm_adgroupid=48947256715&utm_device=c&utm_keyword=&utm_matchtype=b&utm_network=g&utm_adpostion=1t1&utm_creative=255798340456&utm_targetid=dsa-473406574235&utm_loc_interest_ms=&utm_loc_physical_ms=1004446&gclid=CjwKCAjw-ITqBRB7EiwAZ1c5U2ue8Ma0thAn3YtgCI-B15RYUtSIZSxCqGx4dPPT8BoKx7gYbnHrchoCL4AQAvD_BwE



