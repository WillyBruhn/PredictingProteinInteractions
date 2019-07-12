#!/usr/bin/Rscript
#
# Willy Bruhn 25.6.19
#
# This script puts together the functionality of all the other scripts.
#
# 1. Call MutComp to produce dx/obj-files
# 2. Calculate pts-files
# Then either
#           3.1 Train a model
#               - calculate all distances with createAllDistances.R
#               - find optimal boosted model evolutionaryAlgOptParam.R
# or
#           3.2 Make predictions on new data
#               - run BoostedNNClassification.R
#
# 4. Build a clustering
#
#   ------------------------------------------------------------------
#   ------------------------------------------------------------------
#   Input:
#
#     pdb_folder          ... pdb-files to be processed
#
#     MutCompOutPut       ... output of MutComp (dx,pts, ...)
#
#     doMutComp           ... compute MutComp
#   
#     doClustering        ... compute a clustering
#
#     mode(Train/Predict/SingleDistance) ... either Train a model or make predictions
#
#     ------------------------------------------------------------------
#     mode == Train 
#     -------------------------------------------------------------------
#     createAllDistances.R
#     -------------------------------------------------------------------
#     MutCompOutPut         ... path to all proteins as produced by MutComp
#
#     distances_train       ... folder in which all distance-matrices will be stored
#
#     numberOfPoints        ... number of points to select (see parameters of RepeatedSampling)
#
#     rounds                ... number of repetitions
#
#     -------------------------------------------------------------------
#     evolutionary algorithm to find optimal boosted model
#     -------------------------------------------------------------------
#     distances_train       ... folder with all available distances
# 
#     MutCompOutPut         ... folder with the proteins, that means
#                               as output of MutComp (dx,pts-files neccessary)
#   
#     labels_train          ... a file with "names" and "labels" sepcifying the functions of the proteins
#                               the model will be build only on the proteins that are mentioned in this file
#                               That means the names must occur in the column-names and row-names of the distances
#                               in the folder "Distances"
#   
#     bestModelDirectory    ... directory to store the best model
#   
#     bestModelFileName     ... name of the best model
#   
#     popSize_train         ... size of the population of the evolutionary algorithm
#   
#     generations_train     ... iterations to run the ev alg
#   
#     royalNumber_train     ... number of individuals that are kept without mutation for the next generation
#   
#     randoms_train         ... number of individuals that are randomly generated in each generation
#   
#   
#     ------------------------------------------------------------------
#     mode == Predict 
#     ------------------------------------------------------------------
#     bestModelDirectory          ... directory to store the best model
#   
#     bestModelFileName           ... name of the best model
#
#                                     name of the file with the model, a ".RData"-file
#
#                                     This file has to be produced with the function BoostedKNN.saveModelToFile()
#                                     The known labels are also stored in this file.
#
#     MutCompOutPut            ... /home/willy/Schreibtisch/PPItest100/Output/
#                                     path to the proteins that predictions should be made for
#                                     Needs the output of MutComp, and produced pts-files.
#
#     predictions_folder          ... /home/willy/Schreibtisch/PPItest100/KNNPredictions/
#                                     path to were the predictions should be stored.
#
#----------------------------------------------------------------------------------
#
#
# 
#-------------------------------------------------------------------------

is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}

if(!is.installed("plot3D")){install.packages("plot3D")}
library(plot3D)
if(!is.installed("funr")){install.packages("funr")}
library(funr)


if(!is.installed("emdist")){install.packages("emdist")}
if(!is.installed("ggplot2")){install.packages("ggplot2")}
if(!is.installed("ggdendro")){install.packages("ggdendro")}
if(!is.installed("plot3D")){install.packages("plot3D")}
if(!is.installed("emdist")){install.packages("emdist")}
if(!is.installed("tcltk")){install.packages("tcltk")}


library("cluster")
library("tcltk")
library("ggplot2")
library("ggdendro")
library("plot3D")
library("emdist")

if(!is.installed("getopt")){install.packages("getopt")}
library("getopt")

rm(is.installed)


#----------------------------------------------------------------------------------
library(getopt)

options(warn=-1)
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'pdb_folder'  , 'a', 2, "character",
  'MutCompOutPut'   , 'b', 2, "character",
  'doMutComp'   , 'c', 2, "logical",
  'doClustering'   , 'd', 2, "logical",
  'mode'   , 'e', 2, "character",

  'distances_train'   , 'f', 2, "character",
  'numberOfPoints'   , 'g', 2, "numeric",
  'rounds'   , 'Q', 2, "numeric",
  
  'labels_train'   , 'i', 2, "character",
  'bestModelDirectory'   , 'j', 2, "character",
  'bestModelFileName'   , 'k', 2, "character",
  'popSize_train'   , 'l', 2, "numeric",
  'generations_train'   , 'm', 2, "numeric",
  'royalNumber_train'   , 'n', 2, "numeric",
  'randoms_train'   , 'o', 2, "numeric",
  
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
if ( is.null(opt$pdb_folder    ) ) { opt$pdb_folder    = paste(funr::get_script_path(),"/../QuickStart/ModelTrain/pdb/",sep = "")     }
if ( is.null(opt$MutCompOutPut    ) ) { opt$MutCompOutPut   = paste(funr::get_script_path(),"/../QuickStart/ModelTrain/Proteins/",sep = "")   }
if ( is.null(opt$doMutComp    ) ) { opt$doMutComp    = TRUE    }
if ( is.null(opt$doClustering    ) ) { opt$doClustering    = FALSE    }
if ( is.null(opt$mode    ) ) { opt$mode    = "Train"    }

if ( is.null(opt$distances_train    ) ) { opt$distances_train    =  paste(funr::get_script_path(),"/../QuickStart/ModelTrain/RepeatedSubSampling/",sep = "")    }
if ( is.null(opt$numberOfPoints   ) ) { opt$numberOfPoints    = 2   }
if ( is.null(opt$rounds    ) ) { opt$rounds    = 5  }

if ( is.null(opt$labels_train    ) ) { opt$labels_train    =    paste(funr::get_script_path(),"/../QuickStart/ModelTrain/labels.txt",sep = "") }
if ( is.null(opt$bestModelDirectory    ) ) { opt$bestModelDirectory    = paste(funr::get_script_path(),"/../QuickStart/ModelTrain/bestModel/",sep = "")    }
if ( is.null(opt$bestModelFileName    ) ) { opt$bestModelFileName    = "bestModel.RData"    }
if ( is.null(opt$popSize_train    ) ) { opt$popSize_train    = 20    }
if ( is.null(opt$generations_train    ) ) { opt$generations_train    = 10    }
if ( is.null(opt$royalNumber_train    ) ) { opt$royalNumber_train    = 1    }
if ( is.null(opt$randoms_train    ) ) { opt$randoms_train    = 2    }

if ( is.null(opt$predictions_folder    ) ) { opt$predictions_folder    = paste(funr::get_script_path(),"/../QuickStart/Predictions/",sep = "")   }

if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }


print(opt)


s1 = paste(funr::get_script_path(), "/NNClassification/optimizeDifferentModels/BoostedKNN.R", sep = "")
print(s1)
source(s1)

#-------------------------------------------------------------------------
CreatALLdx <- function(ListOfProtNames,PathToProtData)
{
  print("Creating .pts-files from dx-files. This might take a while ...")
  
  for(i in 1:NROW(ListOfProtNames)){
    # needs only name without extension
    fileName <- ListOfProtNames[i]
    
    # needs full path
    inPath=paste(PathToProtData,"/",ListOfProtNames[i],sep="")
    
    # needs full path
    outPath=inPath
    
    if(!file.exists(paste(outPath,"/",fileName,"_pot_positive.pts",sep="")) | 
       !file.exists(paste(outPath,"/",fileName,"_pot_negative.pts",sep="")))
    {
      print(paste("Creating pts-files for ", fileName, " (",i,"/",NROW(ListOfProtNames),")",sep =""))
      
      dxData <- read.csv(file=paste(inPath,"/",fileName,"_pot",".dx",sep=""),
                         sep=' ', skip= 11, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
      dxData <- head(dxData,-10)
      
      v1 = as.numeric(dxData$V1)
      v2 = as.numeric(dxData$V2)
      v3 = as.numeric(dxData$V3)
      
      size = 129
      x <- c(1:size)
      y <- c(1:size)
      z <- c(1:size)
      
      merged <- as.vector(rbind(v1,v2,v3)) 
      V <- array(merged, c(size,size,size))
      
      # print("creating isosurface ...")
      iso <- createisosurf(x, y, z, V, level = 1.0)
      iso2 <- createisosurf(x, y, z, V, level = -1.0)
      
      iso <- iso[!duplicated(iso), ]
      iso2 <- iso2[!duplicated(iso2), ]
      # print(paste("writing to file ", outPath, "/", fileName, ".pts ...", sep = ""))
      write.table(iso, file = paste(outPath,"/",fileName,"_pot_positive.pts",sep=""),
                  row.names = F,na = "",sep = ";",dec = ".")
      write.table(iso2, file = paste(outPath,"/",fileName,"_pot_negative.pts",sep=""),
                  row.names = F,na = "",sep = ";",dec = ".")
    }
  }
  
  print("... done creating .pts-files from dx-files.")
}


createFolderHierarchy <- function(folderName, pdbs){
  # Creates a folder-hierarchy that is needed in the next steps
  # Input -pdb
  #       -pqr
  #
  # Output
  # 
  # the pdbs in the folder "pdbs" are copied to Input/pdb/
  #
  
  
  if(!dir.exists(paste(folderName))) dir.create(folderName)
  
  Input = paste(folderName, "/Input/", sep = "")
  if(!dir.exists(Input)) dir.create(Input)
  if(!dir.exists(paste(Input, "/pdb/", sep = "")))  dir.create(paste(Input, "/pdb/", sep = ""))
  if(!dir.exists(paste(Input, "/pqr/", sep = "")))  dir.create(paste(Input, "/pqr/", sep = ""))
  
  
  Output = paste(folderName, "/Output/", sep = "")
  if(!dir.exists(Output)) dir.create(Output)
  
  list.of.files <- list.files(pdbs, ".pdb", full.names = TRUE)
  file.copy(list.of.files, paste(Input, "/pdb/", sep = ""))
  
}
#-------------------------------------------------------------------------
# 0. Parameter-processing to this script
# pdbFolder = "/home/willy/Schreibtisch/PPIToy/ModelTrain/pdb/"
# PPIoutputFolder = "/home/willy/Schreibtisch/PPIToy/ModelTrain/Proteins/"
# 
# doMutComp = TRUE
# doRepeatedSampling = FALSE
# doClustering = FALSE
# doClassification = FALSE
# 
# labels_train = "/home/willy/Schreibtisch/PPIToy/ModelTrain/labels.txt"

# pdbFolder = "/home/willy/Schreibtisch/PPIToy/Test/pdb/"
# PPIoutputFolder = "/home/willy/Schreibtisch/PPIToy/Test/Proteins/"
# 
# doMutComp = TRUE
# doRepeatedSampling = FALSE
# doClustering = FALSE
# doClassification = FALSE
# 
# labels_train = "/home/willy/Schreibtisch/PPIToy/ModelTrain/labels.txt"


pdbFolder = opt$pdb_folder
PPIoutputFolder = opt$MutCompOutPut

pathToProteins = paste(PPIoutputFolder, "/Output/",sep="")

# print(pathToProteins)
# 
# return()

doMutComp = opt$doMutComp
doRepeatedSampling = FALSE
doClustering = opt$doClustering
mode = opt$mode

labels_train = opt$labels_train

# the dummy-file containing the needed parameters for MutComp
# funr::get_script_path() does not work within Rstudio
# MutCompSettings = "/home/willy/PredictingProteinInteractions/Classification/predictingProteinInteractionsSettings/MutCompParametersDummy.txt"
MutCompSettings = paste(funr::get_script_path(),"/predictingProteinInteractionsSettings/MutCompParametersDummy.txt",sep = "")

createFolderHierarchy(PPIoutputFolder, pdbFolder)
file.copy(labels_train, paste(PPIoutputFolder, "/Output/labels.txt", sep = ""))

# df -h --total
#-------------------------------------------------------------------------

copyMutCompParameterDummy <- function(MutCompParametersPath,MutCompParametersFile, MutCompDummy){
  t2 = c()
  
  if(!file.exists(paste(MutCompParametersPath,MutCompParametersFile, sep = ""))) {
    t = read.csv2(MutCompDummy, header = FALSE)
    t2 = data.frame(matrix(rep(0,nrow(t)), nrow = 1))
    colnames(t2) = t$V1
    t2[1,] = t$V2
  } else {
    t = read.csv2(paste(MutCompParametersPath,MutCompParametersFile, sep = ""), header = FALSE)
    t2 = data.frame(matrix(rep(0,nrow(t)), nrow = 1))
    colnames(t2) = t$V1
    t2[1,] = t$V2
  }
  
  return(t2)
}

s_d <- function(x, k = 4) trimws(format(round(x, k), nsmall=k))

# s_d(8)

#-------------------------------------------------------------------------
# 1. Call MutComp to produce dx/obj-files
MutCompPath = paste(funr::get_script_path(),"/../PreProcessingProteins/MutComp/",sep = "")
MutCompExe = "./process.sh"
MutCompParametersPath = PPIoutputFolder
MutCompParametersFile = "parameters.txt"


MutCompParameters = copyMutCompParameterDummy(MutCompParametersPath,MutCompParametersFile, MutCompSettings)

# change the parameters accordingly
MutCompParameters$parametersPath = PPIoutputFolder
MutCompParameters$MutCompInputFolder = paste(PPIoutputFolder, "/Input/", sep = "")
MutCompParameters$MutCompOutputFolder = paste(PPIoutputFolder, "/Output/", sep = "")
MutCompParameters$MutCompLog = paste(PPIoutputFolder, "/Output/mutComp.log", sep = "")
MutCompParameters$MutCompError = paste(PPIoutputFolder, "/Output/mutComp.err", sep = "")

M2 = t(MutCompParameters[1,])
# MutCompSettings = paste("/home/willy/PredictingProteinInteractions/Classification","/predictingProteinInteractionsSettings/MutCompParametersDummy.txt",sep = "")

# write the new parameters to the file
write.table(M2, file = paste(MutCompParametersPath,MutCompParametersFile, sep =""), sep = ";", row.names = TRUE, quote = FALSE, col.names = FALSE)

# execute MutComp
if(doMutComp) system2(paste(MutCompPath,MutCompExe, sep = ""), args = paste(MutCompParametersPath,MutCompParametersFile, sep = ""))

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# 2. Call a programm that calculates all pairwise distances (RepeatedSampling, 2step-downsampling, ICP, ...)
# In case of CompareProteins:
CompareProteinsPath = paste(funr::get_script_path(),"/../MetricGeometry/ComparingProteins/",sep = "")
CompareProteinsExe = "./CompareIsosurfaces.R"
# MutCompParametersPath = "/home/willy/Schreibtisch/ExampleHierarchy/"
# MutCompParametersFile = "parameters.txt"

# system2(paste(MutCompPath,MutCompExe, sep = ""), args = paste(MutCompParametersPath,MutCompParametersFile, sep = ""))
#-------------------------------------------------------------------------
# In case of RepeatedSampling:
# RepeatedSamplingPath = "/home/willy/PredictingProteinInteractions/MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/"
# 
# RepeatedSamplingPath = paste(funr::get_script_path(),"/../MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/", sep = "")
# RepeatedSamplingExe = "./main"
# RepeatedSamplingArguments = "/home/willy/Schreibtisch/ExampleHierarchy/"
# 
# 
# RepeatedSamplingArguments = list("path" = paste(PPIoutputFolder, "/Output/",sep=""),
#      "outPath" = paste(PPIoutputFolder, "/RepSubOutput/",sep=""),
#      "proteinsToCompareFile_target" = "/home/willy/PredictingProteinInteractions/data/additionalPDBS_1/Output/names.txt",
#      "proteinsToCompareFile" = "/home/willy/PredictingProteinInteractions/data/additionalPDBS_1/Output/names.txt",
#      "measure" = s_d(1),
#      "number_of_selected_points" = "100",
#      "rounds" = "5000",
#      "c1" = s_d(1),
#      "c2" = s_d(0),
#      "c3" = s_d(0),
#      "emd_list_id" = "test",
#      "allParameterCombinations" = "0",
#      "NNtoActCent" = "0")
# 
# if(!dir.exists(paste(PPIoutputFolder, "/RepSubOutput/",sep=""))) dir.create(paste(PPIoutputFolder, "/RepSubOutput/",sep=""))
# 
# # create a file with all the directory-names
proteinNames = list.dirs(path = paste(PPIoutputFolder, "/Output/",sep=""), recursive = FALSE,full.names = FALSE)
# proteinNamesFile = paste(RepeatedSamplingArguments$path,"names.txt",sep="")
# 
# write(proteinNames, file = proteinNamesFile, ncolumns = 1, append = FALSE, sep = " ")
# RepeatedSamplingArguments$proteinsToCompareFile_target = proteinNamesFile
# RepeatedSamplingArguments$proteinsToCompareFile = proteinNamesFile
# 
# 
# # s << outpath << "/EMD_" << n << "_" << rounds << "_" << distribution
# # << "_" << c1 << "_" << c2 << "_" << c3 << "_id_" << id << "_NNact_" << NNToActCent << ".csv";
# 
# 
# 
# RepeatedSamplingArguments_distanceMatrix = paste("EMD_", RepeatedSamplingArguments$number_of_selected_points,
#                                                  "_", RepeatedSamplingArguments$rounds, "_", RepeatedSamplingArguments$measure,
#                                                  "_", RepeatedSamplingArguments$c1, "_", RepeatedSamplingArguments$c2, "_", RepeatedSamplingArguments$c3,
#                                                  "_id_", RepeatedSamplingArguments$emd_list_id, "_NNact_", RepeatedSamplingArguments$NNtoActCent, ".csv", sep ="")

CreatALLdx(ListOfProtNames = proteinNames, PathToProtData = paste(PPIoutputFolder, "/Output/",sep=""))

#-------------------------------------------------------------------------
# 3.2 Training
#-------------------------------------------------------------------------

makeCustomCall <-function(pathToExe,call, parameterList, print = TRUE){
  
  arguments = rep(0,2*length(parameterList))
  
  for(i in 1:length(parameterList)){
    arguments[1+ (i-1)*2] = names(parameterList[i])
    arguments[2+ (i-1)*2] = as.character(parameterList[i])
  }
  
  print(paste(pathToExe,call,sep = ""))
  print(arguments)
  
  system2(paste(pathToExe,call,sep = ""), args = arguments)
}

# create the names-file
print("Creating file with the names ...")
names_file = paste(pathToProteins,"/names.txt", sep ="")
proteinnames = list.dirs(path = pathToProteins, recursive = FALSE, full.names = FALSE)
print(proteinnames)
write(file = names_file,proteinnames)

if(mode == "Train"){
  print("Creating all distances with createAllDistances.R ...")
  pathToCreateAllDistances = paste(funr::get_script_path(),"/NNClassification/optimizeDifferentModels/",sep = "")
  createAllDistancesExe = "./createAllDistances.R"

  createAllDistancesParameters = list(
    "--ProteinsPath" = pathToProteins,
    "--outPutRepSubSamp" = opt$distances_train,
    "--n" = opt$numberOfPoints,
    "--m" = opt$rounds)
  
  # makeCustomCall(pathToCreateAllDistances,createAllDistancesExe,createAllDistancesParameters)
  
  #-------------------------------------------------------------------------
  print("Training with evolutionaryAlgOptParam.R ...")
  pathToEvAlg = paste(funr::get_script_path(),"/NNClassification/optimizeDifferentModels/",sep = "")
  evalgExe = "./evolutionaryAlgOptParam.R"
  
  
  evAlgParameters = list(
    "--Proteins" = pathToProteins,
    "--distances_folder" = opt$distances_train,
    "--labels" = labels_train,
    "--bestModelDirectory" = opt$bestModelDirectory,
    "--bestModelFileName" = opt$bestModelFileName,
    "--popSize" = opt$popSize,
    "--generations" = opt$generations,
    "--royalNumber" = opt$royalNumber,
    "--randoms" = opt$randoms)
  
  makeCustomCall(pathToEvAlg,evalgExe,evAlgParameters)
}
#-------------------------------------------------------------------------
# 3.2 Classification
#-------------------------------------------------------------------------

if(mode == "Predict"){
  # pathToBoostedNNClassification = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/"
  pathToBoostedNNClassification = paste(funr::get_script_path(),"/NNClassification/",sep = "")
  BoostedNNClassificationCall = "./BoostedNNClassification.R"
  
  BoostedNNClassificationParameters = list(
       "--model_f_name" = paste(opt$bestModelDirectory, "/", opt$bestModelFileName, sep = ""),
       "--path_to_test_proteins" = paste(opt$MutCompOutPut,"/Output/",sep=""),
       "--predictions_folder" = opt$predictions_folder)
  
  # BoostedNNarguments = rep(0,2*length(BoostedNNClassificationParameters))
  # 
  # for(i in 1:length(BoostedNNClassificationParameters)){
  #   BoostedNNarguments[1+ (i-1)*2] = names(BoostedNNClassificationParameters[i])
  #   BoostedNNarguments[2+ (i-1)*2] = as.character(BoostedNNClassificationParameters[i])
  # }
  # 
  # system2(paste(pathToBoostedNNClassification,BoostedNNClassificationCall,sep = ""), args = BoostedNNarguments)
  
  makeCustomCall(pathToBoostedNNClassification,BoostedNNClassificationCall,BoostedNNClassificationParameters)
}  



#-------------------------------------------------------------------------
# 4. Clustering

if(mode == "SingleDistance"){
  #-------------------------------------------------------------------------
  # comput one distance-matrix with RepSubSample
  #-------------------------------------------------------------------------
  
  RepeatedSamplingPath = paste(funr::get_script_path(),"/../MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/",sep = "")
  RepeatedSamplingExe = "./main"
  
  positive = getRepeatedSampling(RepeatedSamplingPath,RepeatedSamplingExe,
                                 path = pathToProteins,
                                 outPath = opt$distances_train, 
                                 proteinsToCompareFile_target = names_file, 
                                 proteinsToCompareFile = names_file,
                                 measure = 1,
                                 number_of_selected_points = opt$numberOfPoints,
                                 rounds = opt$rounds,
                                 c1 = 1, c2 = 0, c3 = 0)
  
  
  # negative = getRepeatedSampling(RepeatedSamplingPath,RepeatedSamplingExe,
  #                                path = pathToProteins,
  #                                outPath = opt$distances_train,
  #                                proteinsToCompareFile_target = names_file,
  #                                proteinsToCompareFile = names_file,
  #                                measure = 1,
  #                                number_of_selected_points = opt$numberOfPoints,
  #                                rounds = opt$rounds,
  #                                c1 = 0, c2 = 1, c3 = 0)
  
  
  AllvsAll.Cluster <- function(outPath, distance_matrix, fname)
  {  
    mydendrogramplot <- function(clust,xlim=NULL,ylim=NULL, title=NULL)
    {
      
      dendrogram <- as.dendrogram(clust)
      dendro.data <- dendro_data(dendrogram)
      
      p <- ggplot() +
        geom_segment(data = dendro.data$segments,
                     aes_string(x = "x", y = "y", xend = "xend", yend = "yend"))+
        theme_dendro()+
        scale_x_continuous(breaks = seq_along(dendro.data$labels$label), 
                           labels = dendro.data$labels$label) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
        ggtitle(title)
      
      if(is.null(xlim) &is.null(ylim))
      {
        p <- p +  coord_cartesian(xlim = xlim, ylim = ylim)
        
      }
      p
    }
    # data = distance_matrix
    # 
    # ProtList <- unique(c(as.character(data[,1]),as.character(data[,2])))
    # print(ProtList)
    # 
    # print(data$emd_distance)
    # 
    # matr.Neg <- matrix(0,nrow = NROW(ProtList),ncol = NROW(ProtList), dimnames = list(ProtList,ProtList))
    # 
    # 
    # 
    # k <- 1
    # for(i in 1:(NROW(ProtList)-1))
    # {
    #   for(j in (i+1):NROW(ProtList))
    #   {
    #     if(k <= NROW(data)){
    #       matr.Neg[i,j] <- data[k,3]
    #       matr.Neg[j,i] <- data[k,3]
    #     }
    #     k <- k+1
    #   }
    # }
    
    agnes.average.Neg <- agnes(x = distance_matrix, diss = T,method = "average",keep.diss = F,keep.data = F)
    mydendrogramplot(agnes.average.Neg,title = "UPGMA")
    ggsave(filename = paste(outPath,"/Dendrogram_", fname, ".pdf",sep=""),height=7, width = 14)
  }
  
  positive_name = paste("positive_n_", opt$numberOfPoints, "_m_", opt$rounds, sep = "")
  # negative_name = paste("negative_n_", opt$numberOfPoints, "_m_", opt$rounds, sep = "")
  
  if(doClustering) AllvsAll.Cluster(outPath = pathToProteins, distance_matrix = positive, positive_name)
  # if(doClustering) AllvsAll.Cluster(outPath = pathToProteins, distance_matrix = negative, negative_name)
  
  
  #-------------------------------------------------------------------------
  # Compute some Summary-statistics
  #-------------------------------------------------------------------------
  write.table(file = paste(pathToProteins, "/summary_",positive_name,".txt",sep=""), summaryFromDistanceMatrix(positive), row.names = FALSE)
  # write.table(file = paste(pathToProteins, "/summary_",negative_name,".txt",sep=""), summaryFromDistanceMatrix(negative), row.names = FALSE)
  
}












