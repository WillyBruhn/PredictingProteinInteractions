#---------------------------------------------------------------------------------
# this file is meant to be loaded by other scripts
# provides the functions needed to perform boosted k-nearest-neighbors optimization
# Willy Bruhn 7.7.2019
#---------------------------------------------------------------------------------
library(doBy)

readDistanceMatrix1 <- function(file = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/d_matrix.csv"){
  d = read.csv(file = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/d_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
  d = data.matrix(frame = d)
  
  return(d)
}

readDistanceMatrix2 <- function(file = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/FLB_repeated_sampling/FelixMasterThesis/ListOfEMD_positive_100.csv"){
  # Read the distance from a file formated with three columns
  # ProtA ProtB EMDdistances
  #
  # returns a quadratic matrix
  
  q = read.csv(file)
  names = unique(q$NameforEMDA)
  
  d = matrix(0, length(names), length(names))
  colnames(d) = names
  rownames(d) = names
  for(i in 1:nrow(q)){
    
    d[which(rownames(d) == q$NameforEMDA[i]),which(colnames(d) == q$NameforEMDB[i])] = q$EMD[i]
    d[which(colnames(d) == q$NameforEMDB[i]),which(rownames(d) == q$NameforEMDA[i])] = q$EMD[i]
  }
  return(d)
}

readDistanceMatrix3 <- function(file = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Example/EMD_1000_100_10.0000_1_1.0000_1.0000_1.0000_id_test.csv
"){
  q = read.csv(file, sep = ";")
  names = unique(q$protein_name)
  
  # print(names)
  
  d = matrix(0, length(names), length(names))
  colnames(d) = names
  rownames(d) = names
  for(i in 1:nrow(q)){
    
    d[which(rownames(d) == q$protein_name_target[i]),which(colnames(d) == q$protein_name[i])] = q$emd_distance[i]
    d[which(colnames(d) == q$protein_name[i]),which(rownames(d) == q$protein_name_target[i])] = q$emd_distance[i]
    
  }
  return(d)
}

readLabels <- function(file = "/home/willy/PredictingProteinInteractions/data/106Redoxins/Labels/labels.txt"){
  # reads the labels (functional, not_functional) of the proteins.
  
  t = read.table(file = file, header = TRUE)
  return(t)
}

s_d <- function(x, k = 4) trimws(format(round(x, k), nsmall=k))

F1_score_confusion <- function(confusion, positive.class="functional") {
  #                functional_pred not_functional_pred
  # functional               0.923               0.077
  # not_functional           0.133               0.867
  #
  # TP FN
  # FP TN
  
  # check for correct format
  positive.class_pred = paste(positive.class,"_pred", sep = "")
  
  # if(which(colnames(confusion) == positive.class_pred) != 1 || which(rownames(confusion) == positive.class) != 1){
  #   print("Error: confusion-matrix does not have the correct format. Aborting ...")
  #   return(NULL)
  # }
  
  pos_pred_col = which(colnames(confusion) == positive.class_pred)
  pos_ind = which(rownames(confusion) == positive.class)
  
  # transform confusion so that in the top-left corner the True-positives are placed
  # TP FN
  # FP TN
  
  confusion2 = confusion
  confusion[1,] = confusion2[pos_ind,]
  confusion[2,] = confusion2[-pos_ind,]
  
  confusion2 = confusion
  confusion[,1] = confusion2[,pos_pred_col]
  confusion[,2] = confusion2[,-pos_pred_col]
  
  TP = confusion[1,1]
  TN = confusion[2,2]
  FN = confusion[1,2]
  FP = confusion[2,1]
  
  
  PPV = TP/(TP + FP)
  TPR = TP/(TP + FN)
  
  TNR = TN/(TN + FP)
  
  F1 = 2*(PPV*TPR)/(PPV+TPR)
  
  ACC = (TP + TN) /(sum(confusion))
  
  # print(ACC)
  
  return(list("F1" = F1, "ACC" = ACC, "TPR" = TPR, "TNR" = TNR, "TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN))
}

readTestNames <- function(file = "/home/willy/PredictingProteinInteractions/data/106Redoxins/Labels/testNames.txt"){
  # reads the labels (functional, not_functional) of the proteins.
  
  t = read.table(file = file, header = TRUE, colClasses = c("character"))
  return(t)
}
#--------------------------------------------------------------------------------------------------
# Further optimization with boosted KNN
#
#--------------------------------------------------------------------------------------------------

KnnPredictions <- function(X,y,k = 1, method = "majority", distWeight = TRUE, classFreq = TRUE){
  # X ... matrix of all pairwise distances
  # y ... vector of all known labels
  # a "NA" indicates that this value should be predicted
  #
  # method:
  #         "majority" ... a majority voting is performed
  #         "probability" ... a class probability is calculated
  #
  # invDist: TRUE, then each point is weighted with its inverse distance
  #                 that means closer points have a higher voting-power
  #
  # classFreq: True, then each point is multiplied with the inverse of the 
  #                   freuqency of the points with the same label in the 
  #                   trainings-set
  #------------------------------------------------------
  
  test_indices = which(is.na(y) == TRUE)
  if(length(test_indices) == 0){
    print("Error: Nothing to predict. Specify with NA.")
    return(NULL)
  }
  
  
  # print(y)
  # we train row-wise and predict column-wise
  X_train = X[,-test_indices]
  y_train = y[!is.na(y)]
  
  # print(y_train)
  
  y_levels = unique(y_train)
  if(length(y_levels) < 2) {
    print("Error: Need at least 2 different class-labels!")
    return(NULL);
  }
  
  # class-frequencies in the training-set
  freq = rep(0,length(y_levels))
  for(j in 1:length(y_levels)){
    freq[j] = length(which(y_train == y_levels[j]))/length(y_train)
  }
  
  # print(freq)
  
  # the predictions are returned in a data-frame
  predictions = data.frame(matrix(0, ncol = length(y_levels)+1, nrow = length(test_indices)))
  colnames(predictions) = c("ind",y_levels)
  
  # print(test_indices)
  # print(X_train)
  
  for(i in 1:length(test_indices)){
    predictions[i,1] = colnames(X)[test_indices[i]]
    
    # get all distances to all instances
    distances = X_train[test_indices[i],]
    
    # k+1 because the dist to itself is included
    min_indices = which.minn(distances,n = k)
    
    # print(min_indices)
    
    y_cand = y_train[min_indices]
    
    # print(y_cand)
    
    occ = rep(0,length(y_levels))
    for(j in 1:length(y_levels)){
      occ[j] = length(which(y_cand == y_levels[j]))
    }
    
    # print(y_levels)
    # print(occ)
    
    weights = occ
    
    # each point is weighted with its inverse distance
    # that means closer points have a higher voting-power
    dist_weights = rep(1,k)
    if(distWeight == TRUE){
      # dist_weights = 1/distances[min_indices]
      
      min_d = min(distances[min_indices])
      max_d = max(distances[min_indices])
      dif = max_d - min_d
      if(dif != 0){
        for(j in 1:k){
          dist_weights[j] = (max_d - distances[min_indices[j]])/dif
        }
      }
      
    }
    
    # print(distances[min_indices])
    # print(dist_weights)
    
    invFreq_weights = rep(1,length(y_levels))
    if(classFreq == TRUE){
      invFreq_weights = 1/freq
    }
    
    for(j in 1:length(y_levels)){
      weights[j] = sum(dist_weights[which(y_cand == y_levels[j])]*invFreq_weights[j])
    }
    
    if(method == "majority"){
      maj = which.max(weights)
      predictions[i,maj+1] = 1
    }
    
    if(method == "probability"){
      su = sum(weights)
      for(j in 1:length(weights)){
        predictions[i,j+1] = weights[j]/su
      }
      
    }
  }
  
  return(predictions)
}

# d0 = readDistanceMatrix2("/home/willy/RedoxChallenges/MasterThesis/MACD/Output106_FLB_100/ListOfEMD_negative_100.csv")
# d1 = readDistanceMatrix3("/home/willy/Schreibtisch/106Test/RepSubOutput/EMD_100_500_1.0000_1.0000_1.0000_1.0000_id_test_NNact_0.csv")
# d2 = readDistanceMatrix3("/home/willy/Schreibtisch/106Test/RepSubOutput/EMD_100_5000_1.0000_1.0000_0.0000_0.0000_id_test_NNact_0.csv")
# d3 = readDistanceMatrix3("/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/EMD_100_500_1.0000_0.0000_0.1000_0.0000_id_opt_NNact_0.csv")

KnnBoosted <- function(y, distances_list, k_list, distWeight_list, classFreq_list, weight_list = rep(1,length(k_list))){
  # y ... labels, entries with "NA" are predicted with the rest used for training
  # distances_list ... list of distance-matrices
  # k_list ... list of values for k
  # distWeight_list ... TRUE/FALSE specifies for each model if the distance to the test-instance
  #                       should be weighted
  # classFreq_list ... TRUE/FALSE specifies for each model if the frequencies of the classes should 
  #                     be used for prediction
  # weight_list ... how much each models vote is weighted in prediction
  #-------------------------------------------------------------------------------------------------
  
  if(length(distances_list) != length(k_list) || length(k_list) != length(distWeight_list) || length(k_list) != length(classFreq_list) || length(k_list) != length(weight_list)){
    print("Error: Not all model-parameters set!")
    return(NULL)
  }
  
  if(sum(weight_list) == 0) weight_list = rep(1, length(weight_list))
  
  # get prediction from first model
  preds = KnnPredictions(distances_list[[1]], y, k = k_list[1], method = "probability", distWeight = distWeight_list[1], classFreq = classFreq_list[1])
  preds[,2:3] = preds[,2:3]*weight_list[1]
  
  # get predictions from the other models
  if(length(distances_list) > 1){
    for(i in 2:length(distances_list)){
      new_preds = KnnPredictions(distances_list[[i]], y, k = k_list[i], method = "probability", distWeight = distWeight_list[i], classFreq = classFreq_list[i])
      new_preds[,2:3] = new_preds[,2:3]*weight_list[i]
      preds = rbind(preds,new_preds)
    }
  }  
  
  # print(preds)
  # average over the predictions
  p1 = aggregate(preds[,2], by=list(preds$ind), sum)
  p2 = aggregate(preds[,3], by=list(preds$ind), sum)
  
  ret = data.frame(p1,p2[,2])
  colnames(ret) = colnames(preds)
  
  # normalize predictions
  for(i in 1:nrow(ret)){
    s = ret[i,2] + ret[i,3]
    ret[i,2] = ret[i,2]/s
    ret[i,3] = ret[i,3]/s
  }
  
  return(ret)
}

# y = as.vector(labels$label)
# y[c(7,6,35,18,20)] = NA
# KnnBoosted(y,list("d1" = d1, "d2" = d2), c(10,10), rep(TRUE,2), rep(TRUE,2), weight_list = c(0,0))

LOOErrorEstimateBoosted <- function(y, distances_list, k_list, distWeight_list, classFreq_list, weight_list, confMatrix = TRUE){
  # calculate the Leave one out error
  # confMatrix ... FALSE, then caclulate the error 
  #-------------------------------------------------------------
  
  y = as.vector(y)
  
  conf = data.frame(matrix(0, nrow = 2, ncol = 2))
  
  positive = "functional"
  negative = "not_functional"
  colnames(conf) = c(positive,negative)
  rownames(conf) = c(paste(positive,"_pred",sep = ""),paste(negative,"_pred",sep = ""))
  
  err = 0
  
  for(i in 1:length(y)){
    y_i = y
    y_i[i] = NA
    pred = KnnBoosted(y_i,distances_list, k_list,distWeight_list,classFreq_list, weight_list)
    
    # print(pred)
    
    if(confMatrix == FALSE){
      pred = pred[,-1]
      if(colnames(pred)[which.max(pred[1,])] != y[i]) err = err + 1
    } else {
      predicted_labels = convertToLabel(pred)
      
      # print(pred)
      # print("prediction")
      # print(predicted_labels)
      # print("groundTruth")
      # print(y[i])
      # print("------------------")
      
      conf = getConfusion(predicted_labels,y[i]) + conf
    }
    
    
    
  }
  
  if(confMatrix == FALSE){
    err = err/length(y)
    
    return(err)
  } else {
    return(data.frame(conf))
  }
}

CrossValidationBoosted <- function(y,bsTimes,bsSize,equalSizes = FALSE, distances_list,k_list,distWeight_list, classFreq_list,weight_list){
  
  y = as.vector(y)
  
  conf = data.frame(matrix(0, nrow = 2, ncol = 2))
  
  positive = "functional"
  negative = "not_functional"
  colnames(conf) = c(positive,negative)
  rownames(conf) = c(paste(positive,"_pred",sep = ""),paste(negative,"_pred",sep = ""))
  
  for(i in 1:bsTimes){
    v = sample(1:length(y) ,bsSize,replace = FALSE)
    
    if(equalSizes == TRUE){
      v_pos = sample(which(y == positive) ,bsSize,replace = FALSE)
      v_neg = sample(which(y == negative) ,bsSize,replace = FALSE)
      
      v = sort(c(v_pos,v_neg))
    }
    
    y_train = y
    y_train[v] = NA
    
    # print(v)
    # print(v)
    # pred = KnnPredictions(X, y_train, k = k, method = "probability", distWeight = distWeight, classFreq = classFreq)
    pred = KnnBoosted(y_train,distances_list, k_list,distWeight_list,classFreq_list,weight_list)
    
    predicted_labels = convertToLabel(pred)
    
    conf = getConfusion(predicted_labels,y[v]) + conf
    # print(conf)
  }
  
  return(data.frame(conf))
}


CrossValidationBoostedModelWraper <- function(y, bsTimes,bsSize,equalSizes = FALSE, model, all_distances){
  return(CrossValidationBoosted(y = y,bsTimes = bsTimes,bsSize = bsSize, equalSizes = equalSizes, distances_list = all_distances[model$distance_indices], model$k_list, distWeight_list = model$distWeight_list, model$classFreq_list, model$weight_list))
}

# CrossValidationBoosted(y,1000,1,equalSizes = FALSE,list("d1" = d2), c(17), rep(FALSE,1), rep(FALSE,1))
# CrossValidationBoosted(y,1000,1,equalSizes = FALSE,list("d1" = d2, "d2" = d2), c(25,25), c(TRUE,FALSE), c(TRUE,FALSE))
# CrossValidationBoosted(y,1000,1,equalSizes = FALSE,list("d1" = d2, "d2" = d2), c(17,17), c(TRUE,FALSE), c(TRUE,FALSE))
# 
# d0 = readDistanceMatrix2("/home/willy/RedoxChallenges/MasterThesis/MACD/Output106_FLB_100/ListOfEMD_negative_100.csv")
# d5  = readDistanceMatrix2("/home/willy/RedoxChallenges/MasterThesis/MACD/Output106_FLB_100/ListOfEMD_positive_100.csv")
# 
# y2 = y[-106]
# CrossValidationBoosted(y2,1000,1,equalSizes = FALSE,list("d1" = d0, "d2" = d5), c(5,5), c(TRUE,FALSE), c(TRUE,FALSE))
# CrossValidationBoosted(y2,1000,1,equalSizes = FALSE,list("d1" = d5), c(17), c(TRUE), c(TRUE))
# CrossValidationBoosted(y,1000,1,equalSizes = FALSE,list("d1" = d2), c(17), c(TRUE), c(TRUE))


LOOErrorEstimate <- function(X,y, k, invDist, classFreq){
  # calculate the Leave one out error
  
  y = as.vector(y)
  
  err = 0
  for(i in 1:nrow(X)){
    y_i = y
    y_i[i] = NA
    pred = KnnPredictions(X,y_i, k,method = "majority",invDist,classFreq)[,-1]
    
    # print(pred)
    
    if(colnames(pred)[which(pred[1,] == 1)] != y[i]) err = err + 1
    
  }
  
  err = err/nrow(X)
  
  return(err)
}

convertToLabel <- function(pred){
  lab = rep(0,nrow(pred))
  
  for(i in 1:nrow(pred)){
    lab[i] = colnames(pred)[which.max(pred[i,-1]) + 1]
  }
  
  return(lab)
}

getConfusion <- function(pred, groundTruth, positive = "functional", negative = "not_functional"){
  if(length(pred) != length(groundTruth)){
    print("Error: Differing lengths")
    return(NULL)
  }
  
  conf = data.frame(matrix(0, nro = 2, ncol = 2))
  rownames(conf) = c(positive,negative)
  colnames(conf) = c(paste(positive,"_pred",sep = ""),paste(negative,"_pred",sep = ""))
  
  for(i in 1:length(pred)){
    # TP
    if(pred[i] == groundTruth[i] && groundTruth[i] == positive) conf[1,1] = conf[1,1] +1
    
    # TN
    if(pred[i] == groundTruth[i] && groundTruth[i] == negative) conf[2,2] = conf[2,2] +1
    
    # FP
    if(pred[i] != groundTruth[i] && groundTruth[i] == positive) conf[1,2] = conf[1,2] +1
    
    # FN
    if(pred[i] != groundTruth[i] && groundTruth[i] == negative) conf[2,1] = conf[2,1] +1
  }
  
  return(conf)
}

CrossValidation <- function(X,y,bsTimes,bsSize,equalSizes = FALSE, k,distWeight, classFreq){
  
  y = as.vector(y)
  
  conf = data.frame(matrix(0, nrow = 2, ncol = 2))
  
  positive = "functional"
  negative = "not_functional"
  colnames(conf) = c(positive,negative)
  rownames(conf) = c(paste(positive,"_pred",sep = ""),paste(negative,"_pred",sep = ""))
  
  for(i in 1:bsTimes){
    v = sample(1:length(y) ,bsSize,replace = FALSE)
    
    if(equalSizes == TRUE){
      v_pos = sample(which(y == positive) ,bsSize,replace = FALSE)
      v_neg = sample(which(y == negative) ,bsSize,replace = FALSE)
      
      v = sort(c(v_pos,v_neg))
    }
    
    y_train = y
    y_train[v] = NA
    
    # print(v)
    # print(v)
    pred = KnnPredictions(X, y_train, k = k, method = "probability", distWeight = distWeight, classFreq = classFreq)
    
    predicted_labels = convertToLabel(pred)
    
    conf = getConfusion(predicted_labels,y[v]) + conf
    # print(conf)
  }
  
  return(data.frame(conf))
}


saveModelToFile <- function(individual, all_distances_file_names ,dir,file, pathToProteinFiles, labels){
  if(!dir.exists(dir)) dir.create(dir)
  
  individual = list("F1" = individual$F1,
                    "pathToProteinFiles" = pathToProteinFiles,
                    "labels" = labels,
                    "distance_files" = all_distances_file_names[individual$distance_indices], 
                    "distance_indices" = individual$distance_indices,
                    "k_list" = individual$k_list,
                    "distWeight_list" = individual$distWeight_list,
                    "classFreq_list" = individual$classFreq_list,
                    "weight_list" = individual$weight_list)
  
  saveRDS(individual, file=paste(dir, "/", file, sep = ""),ascii = TRUE)
}

loadModelFromFile <- function(file){
  return(readRDS(file))
}

checkIfModelValid <- function(model){
  print("Checking if model is valid ...")
  
  if(!dir.exists(model$pathToProteinFiles)){
    print(paste("Error: pathToProteinFiles does not exist. ", model$pathToProteinFiles, " Exiting ..."))
    return(1)
  }
  
  print(paste("reading ", model$distance_files[1],sep =""))
  d = readDistanceMatrix3(model$distance_files[1])
  if(!(setequal(colnames(d),model$labels$name))){
    print(paste("Error: Label names do not fit to names in distance-matrix. Exiting ...." ))
    return(1)
  }
  return(0)
}

getRepSampParametersFromFileName <- function(fName){
  
  s1 = strsplit(x = fName,split = "EMD_")[[1]][2]
  s1
  s2 = strsplit(x = s1,split = "_id")[[1]][1]
  lastSplit = strsplit(x = s2,split = "_")[[1]]
  
  n = lastSplit[1]
  m = lastSplit[2]
  measure = lastSplit[3]
  c1 = lastSplit[4]
  c2 = lastSplit[5]
  c3 = lastSplit[6]
  
  return(list("n" = n, "m" = m, "measure" = measure, "c1" = c1, "c2" = c2, "c3" = c3))
}

# saveModelToFile(population[[1]], all_distances_file_names = all_distances_file_names, 
#                 dir = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/evAlg/",
#                 file = "bestModel.RData")
# 
# l = loadModelFromFile(file ="/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/RepSubSamp/evAlg/bestModel.RData")


getRepeatedSampling <- function(RepeatedSamplingPath,RepeatedSamplingExe, path,outPath, proteinsToCompareFile_target, proteinsToCompareFile, measure, number_of_selected_points, rounds, c1, c2, c3){
  RepeatedSamplingArguments = list("path" = path,
                                   "outPath" = outPath,
                                   "proteinsToCompareFile_target" = proteinsToCompareFile_target,
                                   "proteinsToCompareFile" = proteinsToCompareFile,
                                   "measure" = s_d(as.numeric(measure)),
                                   "number_of_selected_points" = as.numeric(number_of_selected_points),
                                   "rounds" = as.numeric(rounds),
                                   "c1" = s_d(as.numeric(c1)),
                                   "c2" = s_d(as.numeric(c2)),
                                   "c3" = s_d(as.numeric(c3)),
                                   "emd_list_id" = "opt",
                                   "allParameterCombinations" = "0",
                                   "NNtoActCent" = "0")
  
  # print(RepeatedSamplingArguments)
  
  # return()
  RepeatedSamplingArguments_distanceMatrix = paste("EMD_", RepeatedSamplingArguments$number_of_selected_points,
                                                   "_", RepeatedSamplingArguments$rounds, "_", RepeatedSamplingArguments$measure,
                                                   "_", RepeatedSamplingArguments$c1, "_", RepeatedSamplingArguments$c2, "_", RepeatedSamplingArguments$c3,
                                                   "_id_", RepeatedSamplingArguments$emd_list_id, "_NNact_", RepeatedSamplingArguments$NNtoActCent, sep ="")
  
  RepeatedSamplingArguments_distanceMatrix_csv = paste(RepeatedSamplingArguments_distanceMatrix, ".csv", sep = "")
  
  if(!file.exists(paste(outPath,"/",RepeatedSamplingArguments_distanceMatrix_csv, sep = ""))) system2(paste(RepeatedSamplingPath,RepeatedSamplingExe, sep = ""), args = unlist(RepeatedSamplingArguments))
  
  return(readDistanceMatrix3(paste(outPath,"/",RepeatedSamplingArguments_distanceMatrix_csv, sep = "")))
}
