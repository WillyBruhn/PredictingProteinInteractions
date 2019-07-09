#!/usr/bin/Rscript
#----------------------------------------------------------------------------------
# 15.4.2019
# Willy Bruhn
#
#----------------------------------------------------------------------------------
# Input:
# distance_f_name("")           ... name of the file with the distances
#
# distance_format(1)
#                           = 1 ... a simple distance-matrix
#                           = 2 ... a 3 column-data-frame. First column ProtA
#                                   second column ProtB, third column distance.
#
# labels_train                   ... name of the file with the labels
#                         
# evaluate(TRUE)        
#                           = 1 ... with crossvalidation the model is tested
#                           = 2 ... the model is tested on the testset
#                           = 0 ... the model makes predictions on new data
#                                   to which not all classlabels are given.
#
# evaluate_n                    ... number of times to repeat the classification.
#                                   Higher number leads to more accurate estimation.
#
# kNN                           ... number of nearest neighbors to use for the majority-
#                                   voting.
#
# test_split                    ... if(evaluate == 2) then this file specifies
#                                   the names to be tested with a model built from 
#                                   the remaining names.
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
  'distance_f_name'  , 'd', 2, "character",
  'distance_format'   , 'f', 2, "integer",
  'evaluate'   , 'e', 2, "integer",
  'evaluate_n'   , 'n', 2, "integer",
  'kNN'   , 'k', 2, "integer",
  'test_split'   , 't', 2, "character",
  'labels_train'   , 'a', 2, "character"
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
if ( is.null(opt$distance_f_name    ) ) { opt$distance_f_name    = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/FLB_repeated_sampling/FelixMasterThesis/ListOfEMD_positive_100.csv"     }
if ( is.null(opt$distance_format    ) ) { opt$distance_format    = 2     }
if ( is.null(opt$evaluate ) ) { opt$evaluate = 0 }
if ( is.null(opt$evaluate_n ) ) { opt$evaluate_n = 100 }
if ( is.null(opt$kNN ) ) { opt$kNN = 10 }
if ( is.null(opt$labels_train    ) ) { opt$labels_train    = ""}
if ( is.null(opt$test_split ) ) { opt$test_split = "" }
if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }

# print some progress messages to stderr, if requested.
if ( opt$verbose ) { write("writing...",stderr()) }

print(opt)
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

print("done sourcing ...")

# labels = getFunctionalProteins()
# label = getProteinLabels(d_only_pos)
# lf = data.frame(colnames(d_train), label)
# 
# names(lf) = c("name", "label")
# write.table(lf, file = "/home/willy/PredictingProteinInteractions/data/106Redoxins/Labels/labels.txt", row.names = FALSE)
# 
# onlyNames = data.frame(lf$name)
# names(onlyNames) = c("names")
# write.table(onlyNames, file = "/home/willy/PredictingProteinInteractions/data/106Redoxins/Labels/AllNames.txt", row.names = FALSE)

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

readTestNames <- function(file = "/home/willy/PredictingProteinInteractions/data/106Redoxins/Labels/testNames.txt"){
  # reads the labels (functional, not_functional) of the proteins.
  
  t = read.table(file = file, header = TRUE, colClasses = c("character"))
  return(t)
}


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
  
  
  
  PPV = confusion[1,1]/(confusion[1,1] + confusion[2,1])
  TPR = confusion[1,1]/(confusion[1,1] + confusion[1,2])
  
  TNR = confusion[2,2]/(confusion[2,2] + confusion[2,1])
  
  F1 = 2*(PPV*TPR)/(PPV+TPR)
  
  ACC = (confusion[1,1] + confusion[2,2]) /(sum(confusion))
  
  # print(ACC)
  
  return(list("F1" = F1, "ACC" = ACC, "TPR" = TPR, "TNR" = TNR))
}

getPredictions <- function(d_train,test_names, labels_train,kNN){
  test_indices = rep(0,length(test_names))
  for(i in 1:length(test_indices)){
    test_indices[i] = which(colnames(d_train) == test_names[i])
  }
  
  lab = rep("NA", ncol(d_train))
  
  for(i in 1:ncol(d_train)){
    ind = which(labels_train$name == names(d_train[,1])[i])
    if(length(ind) > 0 )  lab[i] = as.character(labels_train$label[ind])
  }
  
  predictions = memoliNNclassification(d_train,lab,test_indices,n=kNN)
  df = data.frame(test_names,predictions)
  
  return(df)
}


LeaveOneOutCV <- function(d,classlabels, kNN = 1, positive.class = "functional", negative.class = "not_functional"){
  labels_un = sort(unique(classlabels))
  random_sample = matrix(0,nrow = n, ncol = length(labels_un))
  colnames(random_sample) = labels_un
  for(i in 1:length(labels_un)){
    random_camel = which(classlabels == labels_un[i])
    # print(random_camel)
    random_sample[,i] = sample(random_camel,n, replace = TRUE)
  }
  
  conf = matrix(0,ncol = length(labels_un), nrow = length(labels_un))
  colnames(conf) = paste(labels_un,"_pred",sep="")
  rownames(conf) = labels_un
  
  # print(conf)
  
  for(i in 1:nrow(d)){
    pred = memoliNNclassification(d,classlabels,i,n=kNN)
    
    # print(pred)
    
    # TP
    if(pred == classlabels[i] && classlabels[i] == positive.class){
      conf[1,1] = conf[1,1] + 1
    }
    
    # TN
    if(pred == classlabels[i] && classlabels[i] == negative.class){
      conf[2,2] = conf[2,2] + 1
    }
    
    # FN
    if(pred != classlabels[i] && classlabels[i] == positive.class){
      conf[1,2] = conf[1,2] + 1
    }
    
    # FP
    if(pred != classlabels[i] && classlabels[i] == negative.class){
      conf[2,1] = conf[2,1] + 1
    }
  }
  
  if(sum(conf) != nrow(d)) return(NULL)
  
  return(conf)
}


optimizeAndPlotKnn <- function(d,labels,Kmax = 25){
  f1_scores = rep(0,Kmax)
  accs = rep(0,Kmax)
  TPR = rep(0,Kmax)
  TNR = rep(0,Kmax)
  
  for(k in c(1:Kmax)){
    conf =  LeaveOneOutCV(d,labels, kNN = k)
    
    F1 = F1_score_confusion(conf)
    
    f1_scores[k] = F1$F1
    accs[k] = F1$ACC
    
    TPR[k] = F1$TPR
    TNR[k] = F1$TNR
  }
  
  plot(y = f1_scores, x = c(1:Kmax), type = "l", col = "blue", ylim = c(0,1), ylab = "", xlab = "k nearest neighbors")
  points(y = accs, x = c(1:Kmax), type = "l", col = "red")
  points(y = TPR, x = c(1:Kmax), type = "l", col = "green")
  points(y = TNR, x = c(1:Kmax), type = "l", col = "darkgreen")
  abline(v = which.max(f1_scores))
  legend(x = 15, y = 0.4, legend=c("F_1", "ACC", "TPR", "TNR"),
         col=c("blue", "red", "green", "darkgreen"), lty=rep(1,4), cex=0.8)
}

# the matrix with which the model will be trained
# d_train = readDistanceMatrix2(file = opt$distance_f_name)
d_train = c()
if(opt$distance_format == 2) {
  d_train = readDistanceMatrix2(file = opt$distance_f_name)
} else if(opt$distance_format == 3) {
  d_train = readDistanceMatrix3(file = opt$distance_f_name)
} else {
  d_train = readDistanceMatrix1(file = opt$distance_f_name)
}

labels_train = readLabels(file = opt$labels_train)

flags =  labels_train$name %in% colnames(d_train)
if(length(which(flags == "TRUE")) != nrow(labels_train)) print("Warning: no label specified for some proteins!")

print(paste(length(which(flags == "TRUE")),nrow(labels_train), sep = " "))
print(colnames(d_train))
print("-----------------")
print(colnames(labels_train$name))


labels_train = labels_train[flags,]


if(opt$evaluate == 1 || opt$evaluate == 2){
  
  # if no test-file is specified then the prediciton-accuracy is estimated with cross-validation
  if(opt$evaluate == 1){
    conf = memoliNNclassificationErrorEstimate(d_train, labels_train$label, n = opt$evaluate_n, kNN = opt$kNN,normalized = TRUE)
    
    print(conf)
  } else {
    test_names = readTestNames(file = opt$test_split)
    
    # print(test_names)
    
    df = getPredictions(d_train,test_names = test_names[,1], labels_train = labels_train, kNN = opt$kNN)
    
    trueLabels = rep("f",nrow(test_names))
    
    for(i in 1:length(trueLabels)){
      trueLabels[i] = as.character(labels_train$label[which(labels_train$name == test_names[i,1])])
    }

    df = data.frame(df, trueLabels)
    names(df) = c("names", "predictions", "trueLabel")
    
    # print(df)
    
    conf_matrix = data.frame(matrix(0,ncol=2,nrow=2))
    colnames(conf_matrix) = c("functional_pred", "not_functional_pred")
    rownames(conf_matrix) = c("functional", "not_functional")
    
    conf_matrix
    
    for(i in 1:nrow(df)){
      if(df$trueLabel[i] == "functional"){
        if(df$predictions[i] == df$trueLabel[i]) conf_matrix[1,1] = conf_matrix[1,1] + 1
        else conf_matrix[1,2] = conf_matrix[1,2] + 1
      } else {
        if(df$predictions[i] == df$trueLabel[i]) conf_matrix[2,2] = conf_matrix[2,2] + 1
        else conf_matrix[2,1] = conf_matrix[2,1] + 1
      }
    }
    
    print(conf_matrix)
    
    


  }
  
} else {
  # the names not mentioned in labels_train are used for testing
  # print(colnames(d_train))
  
  test_names = setdiff(colnames(d_train), labels_train$name)
  train_names = setdiff(colnames(d_train), test_names)
  
  # test_indices = rep(0,length(test_names))
  # for(i in 1:length(test_indices)){
  #   test_indices[i] = which(colnames(d_train) == test_names[i])
  # }
  # 
  # lab = rep("NA", ncol(d_train))
  # 
  # for(i in 1:ncol(d_train)){
  #   ind = which(labels_train$name == names(d_train[,1])[i])
  #   if(length(ind) > 0 )  lab[i] = as.character(labels_train$label[ind])
  # }
  # 
  # predictions = memoliNNclassification(d_train,lab,test_indices,n=opt$kNN)
  # df = data.frame(test_names,predictions)
  

  print("Trained model with the following proteins:")
  print(train_names)
  
  print("Predictions:")
  # print(test_names)
  
  if(length(test_names) == 0) {
    print(paste("All names are used for training. No names left for testing. Change the file ", opt$labels_train, sep =""))
    q(status=1)
    
  }
  
  df = getPredictions(d_train,test_names = test_names, labels_train = labels_train, kNN = opt$kNN)
  
  print(df)
}





#----------------------------------------------------------------------------------
# signal success and exit.
q(status=0)
