
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


write.table(sources, paste(PredictingProteinInteractionsFolder, "setUp/SourcableFiles.txt", sep = ""), row.names = FALSE)
#-------------------------------------------------------
#
# install packages
#
#-------------------------------------------------------

install.packages("rfUtilities")
library(rfUtilities)

library(rgl)

library(readobj)

R.Version()


install.packages("installr")
library(installr)


updateR()






