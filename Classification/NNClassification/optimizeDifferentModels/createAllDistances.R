#!/usr/bin/Rscript
#---------------------------------------------------------------------------------
# script to produce all the distance-matrices with all different parameters
# Willy Bruhn 7.7.2019
# 
#   ProteinsPath        ... path to all proteins as produced by MutComp
#
#   outPutRepSubSamp    ... folder in which all distance-matrices will be stored
#
#   n                   ... number of points to select (see parameters of RepeatedSampling)
#
#   m                   ... number of repeatitions
#
#---------------------------------------------------------------------------------
library(getopt)

options(warn=-1)

#----------------------------------------------------------------------------------
# Input
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'ProteinsPath'  , 'p', 2, "character",
  'outPutRepSubSamp'   , 'r', 2, "character",
  'n'   , 'n', 2, "integer",
  'm'   , 'm', 2, "integer"
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
if ( is.null(opt$ProteinsPath    ) ) { opt$ProteinsPath    = "/home/willy/PredictingProteinInteractions/data/106Redoxins/Output/"     }
if ( is.null(opt$outPutRepSubSamp    ) ) { opt$outPutRepSubSamp    = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/"    }
if ( is.null(opt$predictions_folder    ) ) { opt$predictions_folder    = "/home/willy/Schreibtisch/PPItest100/KNNPredictions/"    }
if ( is.null(opt$n    ) ) { opt$n    = 100   }
if ( is.null(opt$m    ) ) { opt$m    = 500  }
if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }

#----------------------------------------------------------------------------------
thisLocation = funr::get_script_path()
print(thisLocation)

split = strsplit(thisLocation,split = "/")
oneUp = split[[1]][1:(length(split[[1]])-1)]

oneUpPath = ""
for(i in 2:length(oneUp)){
  oneUpPath = paste(oneUpPath, oneUp[i],sep = "/")
}
oneUpPath = paste(oneUpPath, "/", sep ="")

print(oneUpPath)



s1 = paste(oneUpPath,"additionalScripts/TriangulateIsoSurface.R", sep ="")
# s1 = paste(thisLocation, "../additionalScripts/TriangulateIsoSurface.R", sep ="")
print(paste("sourcing ", s1, sep =""))

print("-----------------------------------------------------------")
source(s1)
# /home/sysgen/Documents/LWB/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/
# source("/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R")

s2 = paste(oneUpPath,"/additionalScripts/isoFaces.R", sep ="")
print(paste("sourcing ", s2, sep =""))
source(s2)
# source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/isoFaces.R")

# source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/extrinsicDistances.R")

s3 = paste(oneUpPath,"/additionalScripts/extrinsicDistances.R", sep ="")
print(paste("sourcing ", s3, sep =""))
source(s3)

# source("/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/BoostedKNN.R")

s4 = paste(oneUpPath,"/optimizeDifferentModels/BoostedKNN.R", sep ="")
print(paste("sourcing ", s4, sep =""))
source(s4)

print("done sourcing")

#--------------------------------------------------------------------------------------------------


plotKnnBootstrap <- function(F1,ACC,TPR,TNR,fname = "NONE"){
  if(fname != "NONE"){
    pdf(file = fname)
  }
  Kmax = length(F1)
  
  plot(y = F1, x = c(1:Kmax), type = "l", col = "blue", ylim = c(0,1), ylab = "", xlab = "k nearest neighbors", main = "bootstrap")
  points(y = ACC, x = c(1:Kmax), type = "l", col = "red")
  points(y = TPR, x = c(1:Kmax), type = "l", col = "green")
  points(y = TNR, x = c(1:Kmax), type = "l", col = "darkgreen")
  abline(v = which.max(F1))
  legend(x = 15, y = 0.4, legend=c("F_1", "ACC", "TPR", "TNR"),
         col=c("blue", "red", "green", "darkgreen"), lty=rep(1,4), cex=0.8)
  
  if(fname != "NONE"){
    dev.off()
  }
}

optimizeAndPlotKnnBootstrap <- function(d,labels,Kmax = 25, n = 1000, fname = "NONE", method = "Memoli", distWeight = TRUE, classFreq = TRUE){
  f1_scores = rep(0,Kmax)
  accs = rep(0,Kmax)
  TPR = rep(0,Kmax)
  TNR = rep(0,Kmax)
  
  TP = rep(0,Kmax)
  TN = rep(0,Kmax)
  FN = rep(0,Kmax)
  FP = rep(0,Kmax)
  
  for(k in c(1:Kmax)){
    print(k)
    
    conf = c()
    if(method == "Memoli"){
      conf =  memoliNNclassificationErrorEstimate(d, labels, n = n, kNN = k,normalized = FALSE)
    } else {
      conf = CrossValidation(X = d, y = as.vector(labels), bsTimes = n, bsSize = 1,equalSizes = TRUE , k= k, distWeight = distWeight, classFreq = classFreq)
    }

    print(conf)
    
    F1 = F1_score_confusion(conf)
    
    f1_scores[k] = F1$F1
    accs[k] = F1$ACC
    
    TPR[k] = F1$TPR
    TNR[k] = F1$TNR
    
    TP[k] = F1$TP
    TN[k] = F1$TN
    FP[k] = F1$FP
    FN[k] = F1$FN
  }
  
  # if(fname != "NONE"){
  #   pdf(file = fname)
  # }
  # plot(y = f1_scores, x = c(1:Kmax), type = "l", col = "blue", ylim = c(0,1), ylab = "", xlab = "k nearest neighbors", main = "bootstrap")
  # points(y = accs, x = c(1:Kmax), type = "l", col = "red")
  # points(y = TPR, x = c(1:Kmax), type = "l", col = "green")
  # points(y = TNR, x = c(1:Kmax), type = "l", col = "darkgreen")
  # abline(v = which.max(f1_scores))
  # legend(x = 15, y = 0.4, legend=c("F_1", "ACC", "TPR", "TNR"),
  #        col=c("blue", "red", "green", "darkgreen"), lty=rep(1,4), cex=0.8)
  # 
  # if(fname != "NONE"){
  #   dev.off()
  # }
  
  plotKnnBootstrap(F1 = f1_scores, ACC = accs, TPR = TPR, TNR = TNR,fname = fname)
  
  return(list("F1"=f1_scores, "ACC"=accs, "TPR" = TPR, "TNR" = TNR, "TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN))
}


#--------------------------------------------------------------------------------------------------
OutputPath = opt$ProteinsPath 

OutputFolderRepSamp = opt$outPutRepSubSamp
c1_vals = c(0,1,0.5,0.1)
c2_vals = c(0,1,0.5,0.1)
c3_vals = c(0,1,0.5,0.1)
n = opt$n
m = opt$m

n_vals = c(2,5,10,50)


print(paste(n, m, sep = ""))

# n_vals = c(2,10,500)

# RepeatedSamplingPath = "/home/willy/PredictingProteinInteractions/MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/"

RepeatedSamplingPath = paste(funr::get_script_path(),"/../../../MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/", sep = "")
RepeatedSamplingExe = "./main"


if(!dir.exists(OutputFolderRepSamp)) dir.create(OutputFolderRepSamp)

print(OutputFolderRepSamp)

count = 0
for(c1 in c1_vals){
  for(c2 in c2_vals){
    for(c3 in c3_vals){
      for(n in n_vals){
      
      if(c1 == 0 && c2 == 0 && c3 == 0) break;
      
        # if(count >= 2) break;
  
        RepeatedSamplingArguments = list("path" = OutputPath,
                                         "outPath" = OutputFolderRepSamp,
                                         "proteinsToCompareFile_target" = paste(OutputPath,"names.txt",sep=""),
                                         "proteinsToCompareFile" = paste(OutputPath,"names.txt",sep=""),
                                         "measure" = s_d(1),
                                         "number_of_selected_points" = n,
                                         "rounds" = m,
                                         "c1" = s_d(c1),
                                         "c2" = s_d(c2),
                                         "c3" = s_d(c3),
                                         "emd_list_id" = "opt",
                                         "allParameterCombinations" = "0",
                                         "NNtoActCent" = "0")
        
        
        RepeatedSamplingArguments_distanceMatrix = paste("EMD_", RepeatedSamplingArguments$number_of_selected_points,
                                                         "_", RepeatedSamplingArguments$rounds, "_", RepeatedSamplingArguments$measure,
                                                         "_", RepeatedSamplingArguments$c1, "_", RepeatedSamplingArguments$c2, "_", RepeatedSamplingArguments$c3,
                                                         "_id_", RepeatedSamplingArguments$emd_list_id, "_NNact_", RepeatedSamplingArguments$NNtoActCent, sep ="")
        
        RepeatedSamplingArguments_distanceMatrix_csv = paste(RepeatedSamplingArguments_distanceMatrix, ".csv", sep = "")
        
        count = count + 1
        print(paste("Calculating ",n, c1, c2, c3, count, "/", length(c1_vals)*length(c2_vals)*length(c3_vals)*length(n_vals),sep = " "))
        if(!file.exists(paste(OutputFolderRepSamp,"/",RepeatedSamplingArguments_distanceMatrix_csv, sep = ""))) system2(paste(RepeatedSamplingPath,RepeatedSamplingExe, sep = ""), args = unlist(RepeatedSamplingArguments))
   
      }
    }
  }
}

print("-----------------------------------------------------------")
print("Done calculating RepSub-distances ...")
print("-----------------------------------------------------------")



