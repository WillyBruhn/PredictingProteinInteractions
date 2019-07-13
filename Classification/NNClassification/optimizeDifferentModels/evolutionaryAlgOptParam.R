#!/usr/bin/Rscript

print(funr::get_script_path())
s1 = paste(funr::get_script_path(), "/BoostedKNN.R", sep = "")
print(s1)
source(s1)

#------------------------------------------------------------------------------------------------------------------------
# evolutionary algorithm to find optimal boosted model
#
#   Distances             ... folder with all available distances
#
#   Proteins              ... folder with the proteins, that means
#                             as output of MutComp (dx,pts-files neccessary)
#
#   labels                ... a file with "names" and "labels" sepcifying the functions of the proteins
#                             the model will be build only on the proteins that are mentioned in this file
#                             That means the names must occur in the column-names and row-names of the distances
#                             in the folder "Distances"
#
#   bestModelDirectory    ... directory to store the best model
#
#   bestModelFileName     ... name of the best model
#
#   popSize               ... size of the population of the evolutionary algorithm
#
#   generations           ... iterations to run the ev alg
#
#   royalNumber           ... number of individuals that are kept without mutation for the next generation 
#
#   randoms               ... number of individuals that are randomly generated in each generation
# 
#------------------------------------------------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------------------------------------------------

library(getopt)

options(warn=-1)

#----------------------------------------------------------------------------------
# Input
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'Proteins'  , 'f', 2, "character",
  'labels'   , 'q', 2, "character",
  'bestModelDirectory'   , 't', 2, "character",
  'bestModelFileName'   , 'b', 2, "character",
  'popSize'   , 'S', 2, "numeric",
  'generations'   , 'g', 2, "numeric",
  'royalNumber'   , 'r', 2, "numeric",
  'randoms'   , 'R', 2, "numeric",
  'distances_folder'  , 'p', 2, "character"
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
if ( is.null(opt$distances_folder    ) ) { opt$distances_folder    =     "/home/willy/Schreibtisch/PPIToy/ModelTrain/RepeatedSubSampling/"}
if ( is.null(opt$Proteins    ) ) { opt$Proteins    =     "/home/willy/Schreibtisch/PPIToy/ModelTrain/Proteins/Output/"}
if ( is.null(opt$bestModelDirectory    ) ) { opt$bestModelDirectory    = "/home/willy/Schreibtisch/PPIToy/ModelTrain/bestModel/"    }
if ( is.null(opt$labels    ) ) { opt$labels    = "/home/willy/Schreibtisch/PPIToy/ModelTrain/Proteins/Output/labels.txt"    }
if ( is.null(opt$bestModelFileName    ) ) { opt$bestModelFileName    = "bestModel.RData"    }

if ( is.null(opt$popSize    ) ) { opt$popSize    = 20    }
if ( is.null(opt$generations    ) ) { opt$generations    = 100    }
if ( is.null(opt$royalNumber    ) ) { opt$royalNumber    = 1    }
if ( is.null(opt$randoms    ) ) { opt$randoms    = 2    }

if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }

#----------------------------------------------------------------------------------

listOfFoldersWithDistances = c(opt$distances_folder)


if(!dir.exists(opt$bestModelDirectory)) dir.create(opt$bestModelDirectory)

bestModelDirectory = opt$bestModelDirectory
bestModelFileName = opt$bestModelFileName

# this directory has to contain all the proteins used for the model
# as an output of MutComp (all dx, pts, have to be present)
pathToMutCompOutput = opt$Proteins

protNames = list.dirs(path = pathToMutCompOutput, full.names = FALSE, recursive = FALSE)
# print(protNames)

labels = readLabels(file = opt$labels)

popSize = opt$popSize
generations = opt$generations
royalNumber = opt$royalNumber
randoms = opt$randoms


#------------------------------------------------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------------------------------------------------

readInAllDistances <- function(path_list, type = 3){
  # path_list ... list of folders from which all csv-files are read in
  # type ... format of the distance-matrix
  #-------------------------------------------------------------------
  print("Reading in all distances ...")
  
  distances_list = list()
  ind = 0
  
  for(i in 1:length(path_list)){
    files_csv = list.files(path_list[i], pattern = ".csv", full.names = TRUE, recursive = FALSE)
    
    for(j in 1:length(files_csv)){
      
      ind = ind + 1
      # print(files_csv[j])
      
      if(type == 3){
        distances_list[[ind]] = readDistanceMatrix3(files_csv[j])
      }    
    }
  }
  
  return(list("all_distances" = distances_list, "file_names" = files_csv))
}

# read in all distances
all_distances_obj = readInAllDistances(listOfFoldersWithDistances)

all_distances = all_distances_obj$all_distances
all_distances_file_names = all_distances_obj$file_names


generateRandomBoostedKnn <- function(numOfDistances, k_range = c(1:35), maxNumberOfModels = 7, maxNumHard = FALSE){
  # generate a random boosted KNN-model
  #------------------------------------------------
  
  numberOfModels = sample(c(1:maxNumberOfModels),size = 1)
  
  if(maxNumHard == TRUE) numberOfModels = maxNumberOfModels
  
  distance_indices = sample(c(1:numOfDistances), size = numberOfModels, replace = TRUE)
  k_list = sample(k_range, size = numberOfModels, replace = TRUE)
  distWeight_list = sample(c(TRUE,FALSE), size = numberOfModels, replace = TRUE)
  classFreq_list = sample(c(TRUE,FALSE), size = numberOfModels, replace = TRUE)
  weight_list = sample(c(1:100), size = numberOfModels, replace = TRUE)
  F1 = -1  # to indicate that this model has not been evaluated yet
  
  return(list("distance_indices" = distance_indices, "k_list" = k_list, "distWeight_list" = distWeight_list, "classFreq_list" = classFreq_list, "weight_list" = weight_list, "F1" = F1))
}


# saveModelToFile(individual = individual, dir = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/evAlg/", file = "bestModel.RData")
# l = loadModelFromFile(paste("/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/evAlg/","bestModel.RData", sep = ""))

loadModelAndBoostedLOO <- function(dir, file, all_distances){
  
  ind = loadModelFromFile(paste(dir,"/",file, sep =""))
  conf = LOOErrorEstimateBoosted(y = y,
                                 distances_list = all_distances[ind$distance_indices],
                                 k_list = ind$k_list,
                                 distWeight_list = ind$distWeight_list,
                                 classFreq_list = ind$classFreq_list,
                                 weight_list = ind$weight_list)
  
  eval = F1_score_confusion(conf)
  
  return(eval)
}

# loadModelAndBoostedLOO(dir = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/evAlg/", file = "bestModel.RData", all_distances)

randomSearch <- function(all_distances, Trials = 100, method = "LOO", bsTimes = 1000, bsSize = 1){
  # Random different Trials to search for optimal parameters
  #--------------------------------------------------
  
  AllVals = c()
  maximalF1Score = 0
  best = c()
  bestModel = c()
  
  bestConf = c()
  
  conf = c()
  
  for(i in 1:Trials){
    randomBoostedKnn = generateRandomBoostedKnn(all_distances, maxNumberOfModels = 7)
    
    if(method == "LOO"){
      conf = LOOErrorEstimateBoosted(y,distances_list = all_distances[randomBoostedKnn$distance_indices], k_list = randomBoostedKnn$k_list,distWeight_list = randomBoostedKnn$distWeight_list, classFreq_list = randomBoostedKnn$classFreq_list,weight_list = randomBoostedKnn$weight_list)
    } else if (method == "CV"){
      conf = CrossValidationBoosted(y,bsTimes = bsTimes,bsSize = bsSize,equalSizes = FALSE,distances_list = all_distances[randomBoostedKnn$distance_indices], k_list = randomBoostedKnn$k_list,distWeight_list = randomBoostedKnn$distWeight_list, classFreq_list = randomBoostedKnn$classFreq_list,weight_list = randomBoostedKnn$weight_list)
      
    }
    
    f1_val = F1_score_confusion(conf)
    if(!is.na(f1_val$F1) && f1_val$F1 > maximalF1Score){
      maximalF1Score = f1_val$F1
      best = f1_val
      
      bestConf = conf
      
      
      bestModel = randomBoostedKnn
    }
    print(paste("Trial:", i, sep = ""))
    print(bestConf)
    print(bestModel)
    print("------------------------------------")
  }
  
  return(list("evaluation" = best, "model" = bestModel))
}

# randomSearch(all_distances = all_distances,Trials = 100)

crossOver <- function(individual1, individual2){
  #-----------------------------------------------------
  # An individual consisists of:
  # ind$distance_indices
  # ind$k_list
  # ind$distWeight_list
  # ind$classFreq_list
  # ind$weight_list
  #
  # The models of the two individuals are combined
  #-----------------------------------------------------
  
  numberOfModels1 = length(individual1$distance_indices)
  numberOfModels2 = length(individual2$distance_indices)
  
  numberOfModelsFinal = round((numberOfModels1+numberOfModels2)/2 + rnorm(n = 1))
  
  if(numberOfModelsFinal <= 0) numberOfModelsFinal = 1
  
  finalModelsIndices = sample(c(1:(numberOfModels1+numberOfModels2)), size = numberOfModelsFinal, replace = TRUE)
  
  
  return(numberOfModelsFinal)
}

mergeIndividuals <- function(individual1, individual2){
  ind_ret = list("distance_indices" = c(individual1$distance_indices,individual2$distance_indices),
                 "k_list" = c(individual1$k_list,individual2$k_list),
                 "distWeight_list" = c(individual1$distWeight_list,individual2$distWeight_list),
                 "classFreq_list" = c(individual1$classFreq_list, individual2$classFreq_list),
                 "weight_list" = c(individual1$weight_list, individual2$weight_list))
  
  return(ind_ret)
}

reduceIndividual <- function(individual1, indices){
  ind_ret = list("distance_indices" = individual1$distance_indices[indices],
                 "k_list" = individual1$k_list[indices],
                 "distWeight_list" = individual1$distWeight_list[indices],
                 "classFreq_list" = individual1$classFreq_list[indices],
                 "weight_list" = individual1$weight_list[indices])
  
  return(ind_ret)
}

circularNext <- function(a,b,val){
  # a,b are the boarders
  # val is the current value
  #---------------------------
  
  val = val  + sign(rnorm(n = 1))
  
  if(val > b) return(a)
  
  if(val < a) return(b)
  return(val)
}

mutation <- function(individual, distance_num, k_range, mutationProbs = c(1,5,50,25,25,50)){
  #-----------------------------------------------------
  # An individual consisists of:
  # length(ind$distance_indices)
  # ind$distance_indices
  # ind$k_list
  # ind$distWeight_list
  # ind$classFreq_list
  # ind$weight_list
  #
  # Only one property is changed
  # mutationProbs ... specifies the probability that this property is changed
  #----------------------------------------------------------------------------
  
  # property to change
  property_ind = sample(c(1:6), size = 1, prob = mutationProbs)
  
  # change number of models
  if(property_ind == 1){
    # print("change number of models")
    
    num = round(rnorm(n = 1) + length(individual$distance_indices))
    if(num <= 0) num = 1
    
    # print(num)
    
    if(num > length(individual$distance_indices)){
      numOfNewModels = num - length(individual$distance_indices)
      
      additionalModels = generateRandomBoostedKnn(numOfDistances = distance_num, k_range = k_range,maxNumberOfModels = numOfNewModels,maxNumHard = TRUE)
      
      return(mergeIndividuals(individual,additionalModels))
      
    } else {
      inds_to_keep = sample(c(1:length(individual$distance_indices)), size = num, replace = FALSE)
      
      return(reduceIndividual(individual, inds_to_keep))
    }
  } else if(property_ind == 2){
    # change distance_matrix
    # print("change distance matrix")
    
    ind_to_change = sample(c(1:length(individual$distance_indices)), size = 1)
    additionalModels = generateRandomBoostedKnn(numOfDistances = distance_num, k_range = k_range,maxNumberOfModels = 1,maxNumHard = TRUE)
    
    inds_to_keep = c(1:length(individual$distance_indices))
    inds_to_keep[-ind_to_change]
    ret = reduceIndividual(individual, inds_to_keep)
    
    return(mergeIndividuals(ret,additionalModels))
  } else if(property_ind == 3){
    # change k value
    # print("change k value")
    
    ind_to_change = sample(c(1:length(individual$distance_indices)), size = 1)
    
    individual$k_list[ind_to_change] = circularNext(k_range[1], k_range[2], individual$k_list[ind_to_change])
    
    return(individual)
  } else if(property_ind == 4){
    # change distWeight value
    # print("change distWeight")
    
    ind_to_change = sample(c(1:length(individual$distance_indices)), size = 1)
    
    individual$distWeight_list[ind_to_change] = !individual$distWeight_list[ind_to_change]
    return(individual)
  }else if(property_ind == 5){
    # change classFreq_list value
    # print("change classFreq_list")
    
    ind_to_change = sample(c(1:length(individual$distance_indices)), size = 1)
    
    individual$classFreq_list[ind_to_change] = !individual$classFreq_list[ind_to_change]
    return(individual)
  } else if(property_ind == 6){
    # change weight_list value
    # print("change weight_list value")
    
    ind_to_change = sample(c(1:length(individual$distance_indices)), size = 1)
    
    individual$weight_list[ind_to_change] = circularNext(1, 100, individual$weight_list[ind_to_change])
    
    return(individual)
  }
}

# mutation(individual = population[[2]],distance_num = 30,k_range = c(1,30))

evolutionaryAlgorithm <- function(labels,protNames, pathToProteinFiles, all_distances, all_distances_file_names, popSize = 20, generations = 100, royalNumber = 2, randoms = 2, loadBestModel = TRUE, bestModelDirectory, bestModelFileName, numToLeaveOut = 1){
  
  # create the vector with the labels
  # NA are the proteins that are to be predicted
  y = rep(NA,length(protNames))
  for(i in 1:length(protNames)){
    if(protNames[i] %in% labels$name){
      y[i] = as.character(labels$label[which(labels$name  == protNames[i])])
    }
  }
  
  # reducedLabels = labels[,]
  indices = rep(0,length(protNames))
  for(i in 1:length(protNames)){
    if(protNames[i] %in% labels$name){
      indices[i] = which(labels$name == protNames[i])
    }
  }
  
  labels = labels[indices,]
  
  #---------------------------------------------------------------------------------
  # randomly initialize the population
  #---------------------------------------------------------------------------------
  print("Initializing population ...")
  population = list()
  for(i in 1:popSize){
    population[[i]] = generateRandomBoostedKnn(length(all_distances), k_range = c(1,length(y)),maxNumberOfModels = 7)
  }
  
  if(loadBestModel == TRUE){
    if(file.exists(paste(bestModelDirectory,"/", bestModelFileName,sep =""))) population[[1]] = loadModelFromFile(paste(bestModelDirectory,"/", bestModelFileName,sep =""))

  }
  
  print("Starting evolution ...")
  for(iteration in 1:generations){
    #---------------------------------------------------------------------------------
    # Evaluate each model
    #---------------------------------------------------------------------------------
    evaluations = rep(0,length(population))
    for(i in 1:popSize){
      # print(i)
      if(numToLeaveOut == 1){
        conf = LOOErrorEstimateBoosted(y = y,
                                       distances_list = all_distances[population[[i]]$distance_indices],
                                       k_list = population[[i]]$k_list,
                                       distWeight_list = population[[i]]$distWeight_list,
                                       classFreq_list = population[[i]]$classFreq_list,
                                       weight_list = population[[i]]$weight_list)
        

      } else if(numToLeaveOut > 1){
        conf = CrossValidationBoostedModelWraper(y = y,bsTimes = 100,bsSize = numToLeaveOut,equalSizes = FALSE,model = population[[i]],all_distances = all_distances)
      }
      
      evaluations[i] =  F1_score_confusion(conf)$F1
      if(is.na(evaluations[i])) evaluations[i] = 0
      population[[i]]$F1 = evaluations[i]
    }
    
    print("-----------------------------------------------------------------------------")
    print(paste("Gen ", iteration, ": Best F1-score: ", max(evaluations) , sep = ""))
    # print(population[[which.max(evaluations)]])
  
    if(numToLeaveOut == 1){
      saveModelToFile(individual = population[[which.max(evaluations)]], all_distances_file_names, dir = bestModelDirectory, file = bestModelFileName, pathToProteinFiles, labels)
    } else {
      # saveModelToFile(individual = population[[which.max(evaluations)]], all_distances_file_names, dir = bestModelDirectory, file = paste("bestModelL", numToLeaveOut,"O.RData",sep =""), pathToProteinFiles, labels)
    }

    
    print("-----------------------------------------------------------------------------")
    
    #---------------------------------------------------------------------------------
    # Royal selection
    #---------------------------------------------------------------------------------
    newPop = list()
    royals = which.maxn(evaluations[i], n = royalNumber)
    for(j in 1:royalNumber){
      newPop[[j]] = population[[royals[j]]]
    }
    
    
    #---------------------------------------------------------------------------------
    # Fitness proportionate selection
    #---------------------------------------------------------------------------------
    totalFitness = sum(evaluations)
    
    # probability for each individual to be selected
    selectionProbability = evaluations/totalFitness
    
    # print(selectionProbability)
    
    # randomly select two individuals for cross-over, proportionally to selectionProbability
    # twoIndividuals = sample(c(1:popSize), size = 2, prob = selectionProbability, replace = FALSE)
    
    #---------------------------------------------------------------------------------
    # Cross-over, not implemented yet
    #---------------------------------------------------------------------------------
    
    
    #---------------------------------------------------------------------------------
    # Mutation
    #---------------------------------------------------------------------------------
    for(i in (royalNumber+1):(popSize-randoms)){
      pop_ind = sample(c(1:popSize), size = 1, prob = selectionProbability, replace = TRUE)
      
      newPop[[i]] = mutation(individual = population[[pop_ind]],distance_num = length(all_distances),k_range = c(1,length(y)))
    }
    
    for(i in (popSize-randoms+1):(popSize)){
      newPop[[i]] = generateRandomBoostedKnn(length(all_distances), k_range = c(1,length(y)), maxNumberOfModels = 7)
    }
    
    population = newPop
  }
  
  return(population[[1]])
}


#------------------------------------------------------------------------------------------------------------------------
# evolutionary algorithm to find optimal boosted model
#
#------------------------------------------------------------------------------------------------------------------------




evolutionaryAlgorithm(labels,protNames,pathToProteinFiles = pathToMutCompOutput,all_distances = all_distances,all_distances_file_names = all_distances_file_names, popSize = popSize,generations = generations,royalNumber = royalNumber,randoms = randoms, 
                      loadBestModel = TRUE, bestModelDirectory = bestModelDirectory, bestModelFileName = bestModelFileName)


# loadModelAndBoostedLOO(bestModelDirectory, bestModelFileName, all_distances)
# 
# 
# 
# 
# m = loadModelFromFile(paste(bestModelDirectory,bestModelFileName, sep = ""))
# m
# 
# conf = CrossValidationBoostedModelWraper(y = y,bsTimes = 1000,bsSize = 2,equalSizes = FALSE,model = m,all_distances = all_distances)
# conf
# 
# F1_score_confusion(conf)



