#!/usr/bin/Rscript
# options <- commandArgs(trailingOnly = TRUE)


#-------------------------------------------------------
#
# sourcable files
#
#-------------------------------------------------------
# PredictingProteinInteractionsFolder = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/"
# PredictingProteinInteractionsFolder = options[1]

path = funr::get_script_path()
vec = strsplit(path, "/setUp")

PredictingProteinInteractionsFolder =  paste(vec[[1]][1],"/", sep ="")
print(PredictingProteinInteractionsFolder)


addSource <- function(table, name, fName){
  sources2 = data.frame(matrix(0,ncol = 2, nrow = 1))
  colnames(sources2) = c("name", "path")
  
  sources2[1, ] = c(name,paste(PredictingProteinInteractionsFolder,fName, sep = ""))
  
  table = rbind(table, sources2)
  return(table)
}

print("creating sources-table ...")
sources = data.frame(matrix(0,ncol = 2, nrow = 0))
colnames(sources) = c("name", "path")

sources = addSource(sources, "helperFunctions", "/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R")
sources = addSource(sources, "UltraQuickRepeatedSubSampling", "/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R")
sources = addSource(sources, "TriangulateIsoSurface", "/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R")
sources = addSource(sources, "QuickRepeatedSubSampling", "/MetricGeometry/QuickRepeatedSubSampling/QuickRepeatedSubSampling.R")
sources = addSource(sources, "BoostedKNN", "/Classification/NNClassification/optimizeDifferentModels/BoostedKNN.R")

sources = addSource(sources, "extrinsicDistances", "/Classification/NNClassification/additionalScripts/extrinsicDistances.R")
sources = addSource(sources, "isoFaces", "/Classification/NNClassification/additionalScripts/isoFaces.R")


write.table(sources, paste(PredictingProteinInteractionsFolder, "setUp/SourcableFiles.txt", sep = ""), row.names = FALSE)

#-------------------------------------------------------
#
# install packages
#
#-------------------------------------------------------

# print("installing packages ...")
# 
# is.installed <- function(mypkg){
#   is.element(mypkg, installed.packages()[,1])
# }
# 
# if(!is.installed("plot3D")){install.packages("plot3D")}
# if(!is.installed("funr")){install.packages("funr")}
# if(!is.installed("emdist")){install.packages("emdist")}
# if(!is.installed("ggplot2")){install.packages("ggplot2")}
# if(!is.installed("ggdendro")){install.packages("ggdendro")}
# if(!is.installed("plot3D")){install.packages("plot3D")}
# if(!is.installed("tcltk")){install.packages("tcltk")}
# if(!is.installed("getopt")){install.packages("getopt")}
# if(!is.installed("rfUtilities")){install.packages("rfUtilities")}
# if(!is.installed("rbenchmark")){install.packages("rbenchmark")}
# 
# 
# if(!is.installed("raster")){install.packages("raster")}
# if(!is.installed("foreach")){install.packages("foreach")}
# if(!is.installed("doMC")){install.packages("doMC")}
# if(!is.installed("FNN")){install.packages("FNN")}
# if(!is.installed("RSNNS")){install.packages("RSNNS")}

# install.packages("rdist")
# misc3d
# geometry
# xtable
# spatstat
# igraph
# gplots
# lubridate
# install.packages("clue")

# install.packages("pracma")

# library(rgl)
# 
# library(readobj)


# install.packages("keras")
# library(keras)
# 
# # install.packages("lime")
# library(lime)
# 
# # install.packages("tidyquant")
# library(tidyquant)
# 
# # install.packages("rsample")
# # library(rsample)
# 
# # install.packages("recipes")
# library(recipes)
# 
# # install.packages("yardstick")
# library(yardstick)
# 
# # install.packages("corrr")
# # library(corrr)
# 
# library(permute)
# library(xtable)
# install.packages("DescTools")



# R.Version()
# install.packages("installr")
# library(installr)
# updateR()



# for yupiter ...
# https://www.datacamp.com/community/blog/jupyter-notebook-r?utm_source=adwords_ppc&utm_campaignid=898687156&utm_adgroupid=48947256715&utm_device=c&utm_keyword=&utm_matchtype=b&utm_network=g&utm_adpostion=1t1&utm_creative=255798340456&utm_targetid=dsa-473406574235&utm_loc_interest_ms=&utm_loc_physical_ms=1004446&gclid=CjwKCAjw-ITqBRB7EiwAZ1c5U2ue8Ma0thAn3YtgCI-B15RYUtSIZSxCqGx4dPPT8BoKx7gYbnHrchoCL4AQAvD_BwE



