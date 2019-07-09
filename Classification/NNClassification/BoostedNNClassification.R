#!/usr/bin/Rscript
#----------------------------------------------------------------------------------
# 7.7.2019
# Willy Bruhn
#
#----------------------------------------------------------------------------------
# Input:
# model_f_name              ... /home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/evAlg/bestModelLOO.RData
#                                 name of the file with the model, a ".RData"-file
#
#                                 This file has to be produced with the function BoostedKNN.saveModelToFile()
#                                 The known labels are also stored in this file.
#
# path_to_test_proteins       ... /home/willy/Schreibtisch/PPItest100/Output/
#                                 path to the proteins that predictions should be made for
#                                 Needs the output of MutComp, and produced pts-files.
#
# predictions_folder          ... /home/willy/Schreibtisch/PPItest100/KNNPredictions/
#                                 path to were the predictions should be stored.
#
#----------------------------------------------------------------------------------
# install.packages("getopt")
library(getopt)

options(warn=-1)

#----------------------------------------------------------------------------------
# Input
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'model_f_name'  , 'f', 2, "character",
  'path_to_test_proteins'   , 't', 2, "character",
  'predictions_folder'   , 'p', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
if ( is.null(opt$model_f_name    ) ) { opt$model_f_name    = "/home/willy/Schreibtisch/PPIToy/ModelTrain/bestModel/bestModel.RData"     }
if ( is.null(opt$path_to_test_proteins    ) ) { opt$path_to_test_proteins    = "/home/willy/Schreibtisch/PPIToy/Test/Proteins/Output/"    }
if ( is.null(opt$predictions_folder    ) ) { opt$predictions_folder    = "/home/willy/Schreibtisch/PPIToy/Predictions/"    }
if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }

# print some progress messages to stderr, if requested.
if ( opt$verbose ) { write("writing...",stderr()) }

print("--------------------------------------------------------------------------")
print(opt)
print("--------------------------------------------------------------------------")
#----------------------------------------------------------------------------------
# source("/home/willy/RedoxChallenges/MasterThesis/memoliModels/scripts/TriangulateIsoSurface.R")

# /home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts
s1 = paste(funr::get_script_path(),"/additionalScripts/TriangulateIsoSurface.R", sep ="")
source(s1)

s2 = paste(funr::get_script_path(),"/additionalScripts/TriangulateIsoSurface.R", sep ="")
# source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/isoFaces.R")
source(s2)

# source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/extrinsicDistances.R")

s3 = paste(funr::get_script_path(),"/additionalScripts/extrinsicDistances.R", sep ="")
source(s3)

print("----------------------")
print(funr::get_script_path())

s4 = paste(funr::get_script_path(),"/optimizeDifferentModels/BoostedKNN.R", sep ="")
source(s4)

print("done sourcing ...")

#----------------------------------------------------------------------------------
# load the model
#----------------------------------------------------------------------------------
print("loading model ...")

# parameters
# model_f_name = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/evAlg/bestModelLOO.RData"
# path_to_test_proteins = "/home/willy/Schreibtisch/PPItest100/Output/"
# predictions_folder = "/home/willy/Schreibtisch/PPItest100/KNNPredictions/"

model_f_name = opt$model_f_name
path_to_test_proteins = opt$path_to_test_proteins
predictions_folder = opt$predictions_folder

# RepeatedSamplingPath = "/home/willy/PredictingProteinInteractions/MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/"
RepeatedSamplingPath = paste(funr::get_script_path(),"/../../MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/", sep = "")
RepeatedSamplingExe = "./main"

model = loadModelFromFile(model_f_name)

# check if the model is valid
if(checkIfModelValid(model) == 1) q(status=1)

print(model)


copyDirectoriesInDirectory <- function(from,to){
  # copies only the contents of the directory
  # not the directory itself
  #---------------------------------------
  
  dirs = list.dirs(path = from, full.names = TRUE,recursive = FALSE)
  for(i in 1:length(dirs)){
    file.copy(dirs[i],to,recursive = TRUE,overwrite = FALSE)
    
    print(paste("copying",dirs[i]))
  }
}

#----------------------------------------------------------------------------------
# copy the needed proteins into the destinationfolder
#----------------------------------------------------------------------------------
proteinsFolder = paste(predictions_folder, "/proteins/", sep = "")
if(!dir.exists(predictions_folder))dir.create(predictions_folder)
if(!dir.exists(proteinsFolder)) {  dir.create(proteinsFolder)}
copyDirectoriesInDirectory(path_to_test_proteins, proteinsFolder)
copyDirectoriesInDirectory(model$pathToProteinFiles, proteinsFolder)
print("... done copying")

distances_folder = paste(predictions_folder, "/distances/", sep = "")
if(!dir.exists(distances_folder)) {  dir.create(distances_folder)}

# create the names-file
print("Creating file with the names ...")
names_file = paste(proteinsFolder,"/names.txt", sep ="")
proteinnames = list.dirs(path = proteinsFolder, recursive = FALSE, full.names = FALSE)
write(file = names_file,proteinnames)

# store the distance-matrices in a list
print("Calculating all distances. This might take a while ...")
model_distances = list()
for(i in 1:length(model$distance_files)){
  parameters = getRepSampParametersFromFileName(fName = model$distance_files[[i]])
  
  # getRepeatedSampling(RepeatedSamplingPath = RepeatedSamplingPath,
  #                     RepeatedSamplingExe = RepeatedSamplingExe,
  #                     path = proteinsFolder,
  #                     outPath = distances_folder,
  #                     proteinsToCompareFile_target = names_file,
  #                     proteinsToCompareFile = names_file,
  #                     measure = parameters$measure,
  #                     number_of_selected_points = parameters$n,
  #                     rounds = parameters$m,
  #                     c1 = parameters$c1,
  #                     c2 = parameters$c2,
  #                     c3 = parameters$c3)
  
  
  model_distances[[i]] = getRepeatedSampling(RepeatedSamplingPath = RepeatedSamplingPath,
                                              RepeatedSamplingExe = RepeatedSamplingExe,
                                              path = proteinsFolder,
                                              outPath = distances_folder,
                                              proteinsToCompareFile_target = names_file,
                                              proteinsToCompareFile = names_file,
                                              measure = parameters$measure,
                                              number_of_selected_points = parameters$n,
                                              rounds = parameters$m,
                                              c1 = parameters$c1,
                                              c2 = parameters$c2,
                                              c3 = parameters$c3)
}
print("... done calculating distances.")

# get the names of the proteins
protNames = colnames(model_distances[[1]])

# create the vector with the labels
# NA are the proteins that are to be predicted
y = rep(NA,length(protNames))
for(i in 1:length(protNames)){
  if(protNames[i] %in% model$labels$name){
    y[i] = as.character(model$labels$label[which(model$labels$name  == protNames[i])])
  }
}

print("Calculating predictions ...")
predictions = KnnBoosted(y = y, distances_list = model_distances, k_list = model$k_list, distWeight_list = model$distWeight_list, classFreq_list = model$classFreq_list,weight_list = model$weight_list)
print("... done")

print(paste("writing predictions to ", predictions_folder,"predictions.txt", sep =""))
write.table(predictions, file = paste(predictions_folder,"predictions.txt", sep = ""), row.names = FALSE)

print(predictions)

#----------------------------------------------------------------------------------
# signal success and exit.
q(status=0)
