#!/usr/bin/Rscript
#---------------------------------------------------------------------------------------------
# Willy Bruhn, 18.7.19
#
# pathToProjection  ... path to the folder where the projection-file is found
#
# n                 ... n, number of points selected in each sample
#
# m                 ... m, number of samples that were drawn from each object
#
# q                 ... number of quantiles to approximate the distribution
#
#
# projectionName    ... fileName, e.g. "proj"
#                       /home/willy/PredictingProteinInteractions/data/animals//proj_n_500_m_200_q_200.csv
#
#                       The projection has the following structure
#                       name  n q_1, ..., q_q
#                     
#                       name          ... name of the subclass of the object
#                       n             ... number of points that were used to generate the distribution
#                       q_1, ..., q_q ... the quantiles to approximate the distribution
#
# trainSplitRatio   ... one number specifying how much data is going into the trainging
#                       (default = 0.7)
#
#  
# ModelName         ... name of the file in which the trained model is saved
#
# clustNum          ... number of clusters to run the training on
#
#----------------------------------------------------------------------------------------------

library(getopt)

options(warn=-1)

#----------------------------------------------------------------------------------
# Input
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'pathToProjection'  , 'p', 2, "character",
  'projectionName'   , 'r', 2, "character",
  'n'   , 'n', 2, "integer",
  'm'   , 'm', 2, "integer",
  'q'   , 'q', 2, "integer",
  'trainSplitRatio'   , 't', 2, "integer",
  'clustNum'   , 'c', 2, "integer",
  'ModelName'  , 'o', 2, "character",
  'mode'  , 'M', 2, "character"
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
if ( is.null(opt$pathToProjection    ) ) { opt$pathToProjection    = "/home/willy/PredictingProteinInteractions/data/animals/"     }
if ( is.null(opt$projectionName    ) ) { opt$projectionName    = "proj"    }
if ( is.null(opt$n    ) ) { opt$n    = 500   }
if ( is.null(opt$m    ) ) { opt$m    = 500  }
if ( is.null(opt$q    ) ) { opt$q    = 10  }
if ( is.null(opt$trainSplitRatio    ) ) { opt$trainSplitRatio    = 0.7  }
if ( is.null(opt$clustNum    ) ) { opt$clustNum    = 5  }
if ( is.null(opt$ModelName    ) ) { opt$ModelName    = paste(opt$pathToProjection,"/KnnArea.RData", sep = "")  }
if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }
if ( is.null(opt$mode    ) ) { opt$mode    = "nn"  }
#----------------------------------------------------------------------------------
thisLocation = funr::get_script_path()

#----------------------------------------------------------------------------------------------

library(randomForest)
library(mlbench)
library(caret)
library(doParallel)
library(doBy)

library(caret)

s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
source(s1)

# generateFileName <- function(n,m,q,fName,path){
#   return(paste(path, "/", fName,"_n_", n, "_m_", m, "_q_", q, ".csv", sep = ""))
# }
# 
# writeProjectionToFile <- function(proj,n,m,q,path = "/home/willy/PredictingProteinInteractions/data/animals/", fName = "proj"){
#   fName_final = generateFileName(n=n,m=m,q=q,fName = fName,path = path)
#   print(paste("writing projection to ",fName_final, sep =""))
#   
#   write.csv(proj, file = fName_final, row.names = FALSE)
# }
# 
# readProjectionFromFile <- function(n,m,q,path,fName){
#   fName_final = generateFileName(n=n,m=m,q=q,fName = fName,path = path)
#   print(paste("loading projection from ", fName_final, sep =""))
#   
#   return(read.csv(file = fName_final))
# }

trainKnnOnEachDistribution <- function(df_proj_train, k_params = c(1,10,30), numClusters = 5){
  # for each point the probability is learned that it belongs to a certain subclass
  # later the classes are predicted based on these subclasses
  #-------------------------------------------------------------------
  print("Training model k-nearest-neighbor ...")
  
  # 1st col = name
  # 2nd col = number of points (n)
  x = df_proj_train[,3:ncol(df_proj_train)]
  y = as.character(df_proj_train[,1])
  
  
  data_train = cbind(x,y)
  # data_train[1:5,]
  
  control <- trainControl(method="repeatedcv", number=5, repeats=10)
  seed <- 71
  metric <- "Accuracy"
  set.seed(seed)
  
  
  mtry <- data.frame(matrix(0, ncol = 1, nrow = 0))
  colnames(mtry) = c("k")
  mtry[1:length(k_params),] = k_params
  
  
  cl <- makePSOCKcluster(numClusters)
  registerDoParallel(cl)
  knn <- train(y~., data=data_train, method="knn", metric=metric, tuneGrid=mtry, trControl=control)
  stopCluster(cl)
  
  return(knn)
}

trainNeuralNetOnEachDistribution <- function(df_proj_train, splitPattern, numClusters = 5, size =c(10), thresh = 0.75, k = 5, number = 10, decay = c(0.1)){
  # for each point the probability is learned that it belongs to a certain subclass
  # later the classes are predicted based on these subclasses
  #-------------------------------------------------------------------
  print("Training model neural-network ...")
  
  # 1st col = name
  # 2nd col = number of points (n)
  x = df_proj_train[,3:ncol(df_proj_train)]
  y = as.character(df_proj_train[,1])
  
  y = getClassNamesFromSubClasses(y,splitPattern = splitPattern)
  
  data_train = cbind(x,y)
  
  seed <- 71
  metric <- "Accuracy"
  set.seed(seed)
  
  numFolds <- trainControl(method = 'cv', number = number, classProbs = TRUE, verboseIter = TRUE, preProcOptions = list(thresh = thresh, ICAcomp = 3, k = k))
  
  cl <- makePSOCKcluster(numClusters)
  registerDoParallel(cl)
  nn <- train(y~., data=data_train, method="nnet", metric=metric, trControl=numFolds, tuneGrid=expand.grid(size=size, decay=decay))
  stopCluster(cl)
  
  return(nn)
}


transformToInputNodes <- function(matRow){
  
  name  = as.character(unique(matRow[,1])[1])
  # print(name)
  
  matRow = matRow[do.call(order,matRow),]

  matRow = matRow[,-c(1:2)]
  # print(matRow)
  
  m = nrow(matRow)
  q = ncol(matRow)
  inputNodes = data.frame(matrix(0,ncol = m*q+1, nrow = 1 ))
  
  # print(inputNodes)
  
  for(i in 1:m){
    
    start= (i-1)*q+1+1
    end = start+q-1
    
    # print(paste(i,q,start,end))
    
    # print(matRow[i,])
    # print(inputNodes)
    inputNodes[1,start:end] = matRow[i,]
  }
  
  inputNodes[1,1] = name
  # colnames(inputNodes)[1] = "name"

  return(inputNodes)
}


# dd[order(-dd[,4], dd[,1]), ]
# 
# order(df_proj[1:m,])
# 
# df_proj[1:m,]
# 
# transformToInputNodes(df_proj[1:10,])
# 
# 
# 
# require(stats)
# 
# (ii <- order(x <- c(1,1,3:1,1:4,3), y <- c(9,9:1), z <- c(2,1:9)))
# ## 6  5  2  1  7  4 10  8  3  9
# rbind(x, y, z)[,ii] # shows the reordering (ties via 2nd & 3rd arg)
# 
# 
# 
# dd <- transform(data.frame(x, y, z),z = factor(z, labels = LETTERS[9:1]))
# ## Either as above {for factor 'z' : using internal coding}:
# dd[ order(x, -y, z), ]
# ## or along 1st column, ties along 2nd, ... *arbitrary* no.{columns}:
# dd[ do.call(order, dd), ]
# sub = df_proj[1:10,]



# t = transformToInputNodes(df_proj[1:10,])
# t

transformAllToInputNodes <- function(df_proj,m){
  q_temp = ncol(df_proj)-2
  inputNodes = data.frame(matrix(0,ncol = m*q_temp+1, nrow = nrow(df_proj)/m ))
  
  for(i in 1:nrow(inputNodes)){
    print(i/nrow(inputNodes))
    start = (i-1)*m+1
    end = start+m-1
    
    inputNodes[i,] = transformToInputNodes(df_proj[start:end,])
  }
  
  colnames(inputNodes) = c("name", make.names(seq(1:(q_temp*m))))
  
  return(inputNodes)
}



# transformAllToInputNodes(df_proj[1:100,], 10)


trainNeuralNetOnEachDistributionQuantiles <- function(df_proj_train, m, splitPattern, numClusters = 5, size =c(10), thresh = 0.75, k = 2, number = 2, decay = c(0.1)){
  # for each point the probability is learned that it belongs to a certain subclass
  # later the classes are predicted based on these subclasses
  #-------------------------------------------------------------------
  print("Training model neural-network ...")
  
  print("Preprocessing ...")
  nodes = transformAllToInputNodes(df_proj_train, m)
  
  x = nodes[,2:ncol(nodes)]
  y = as.character(nodes[,1])
  
  y = getClassNamesFromSubClasses(y,splitPattern = splitPattern)
  
  print(y)
  
  print(x[1:5,])
  
  data_train = cbind(x,y)
  
  write.csv(x = data_train, file = "/home/willy/Schreibtisch/test.csv", row.names = FALSE)
  
  print("... done preprocessing")
  
  seed <- 71
  metric <- "Accuracy"
  # set.seed(seed)
  
  
  numFolds <- trainControl(method = 'cv', number = number, classProbs = TRUE, verboseIter = TRUE, preProcOptions = list(thresh = thresh, ICAcomp = 3, k = k))
  
  cl <- makePSOCKcluster(numClusters)
  registerDoParallel(cl)
  nn <- caret::train(y~., data=data_train, method="nnet", metric=metric, trControl=numFolds, tuneGrid=expand.grid(size=size, decay=decay))
  stopCluster(cl)
  
  return(nn)
}

evaluateAccuracyQuantileNN <- function(model = model,df_proj_test, m = m){
  print("Preprocessing ...")
  nodes = transformAllToInputNodes(df_proj_test, m)
  
  x = nodes[,2:ncol(nodes)]
  y = as.character(nodes[,1])
  
  y = getClassNamesFromSubClasses(y,splitPattern = splitPattern)
  
  print(y)
  
  print(x[1:5,])
  
  data_train = cbind(x,y)
  
  print("... done preprocessing")
  
  accuracy(predict(object = model, newdata = x),y)
}

getClassNamesFromSubClasses <- function(subClasses, splitPattern = "-"){
  # gets the className from a subclass
  # e.g.
  # Lion-01 --> Lion
  #-------------------------------------------------------
  
  classNames = rep("", length(subClasses))
  for(i in 1:length(subClasses)){
    classNames[i] = strsplit(as.character(subClasses[i]),split = splitPattern)[[1]][1]
  }
  
  return(classNames)
}

splitTrainTest <- function(df_proj, trainSplitRatio){
  numOfObjects = length(unique(df_proj[,1]))
  
  print(numOfObjects)
  
  trainNames = sample(unique(df_proj[,1]), size = numOfObjects*trainSplitRatio,replace = FALSE)
  
  trainIndices = which(df_proj[,1] %in% trainNames)
  
  df_train = df_proj[trainIndices,]
  df_test = df_proj[-trainIndices,]
  
  list("train" = df_train, "test" = df_test)
}

calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

calculate_class_probabilities <- function(subClassProbs, classNamesAll) {
  
  classProbs = data.frame(matrix(0,ncol = length(classNamesAll), nrow =1))
  colnames(classProbs) = classNamesAll
  
  # print(classProbs)
  
  classNames = getClassNamesFromSubClasses(names(subClassProbs))


  # classProbs = rep(0,length(classNames))
  for(i in 1:length(classNamesAll)){

    indices = which(classNames == classNamesAll[i])
    
    su = sum(subClassProbs[indices])

    classProbs[1, i] = su
  }
  
  return(classProbs)
}

evaluateAccuracy <- function(model, x_test, y_test, m, probs = TRUE){
  print("Evaluating accuracy ...")
  
  y_pred_subClass = c()
  if(probs == TRUE){
    y_pred_subClass = predict(model, newdata = x_test, type = "prob")
  } else {
    y_pred_subClass = predict(model, newdata = x_test)
  }
  
  # for each row we have on object and in each column we have another 
  # sample from the object. For each sample a sub-class-probability is calculated
  # These sub-class-probabilities are converted to class-probabilities and
  # then averaged to get the prediction for the class
  y_pred_class = matrix(getClassNamesFromSubClasses(subClasses = y_pred_subClass), ncol = m, byrow = TRUE)
  
  y_pred_final =rep(0,length(y_test)/m)
  
  # these are the actual classes that need to be predicted correctly
  # of each subclass m instances are present
  y_test_classes = getClassNamesFromSubClasses(subClasses = y_test)[seq(1,length(y_test),m)]
  
  y_levels = unique(y_test_classes)
  

  if(probs == TRUE){
    classProbs_test_objs = data.frame(matrix(0,ncol = length(y_levels), nrow = length(y_test_classes)))
    colnames(classProbs_test_objs)  = y_levels
    
    classProbs = data.frame(matrix(0,ncol = length(y_levels), nrow = m))
    colnames(classProbs) = y_levels
    # weighted-version
    for(i in 1:length(y_pred_final)){
      
      for(j in 1:m){
        
        # print(y_pred_subClass[(i-1)*m+j,])
        # print(cl)
  
        cl = calculate_class_probabilities(y_pred_subClass[(i-1)*m+j,], y_levels)
        classProbs[j,] = cl[1,]
      }
      
      # print(classProbs)
      classProbs_test_objs[i,] = colSums(classProbs)/sum(classProbs)
      
      # print(classProbs_test_objs[i,] )
      # return()
    }
  
    for(i in 1:length(y_pred_final)){
      y_pred_final[i] = y_levels[which.maxn(classProbs_test_objs[i,], n = 1)]
      
      # print(classProbs_test_objs[i,])
      # print(y_pred_final[i])
    }
  
  } else {
    # mode-version
    for(i in 1:length(y_pred_final)){
      y_pred_final[i] = calculate_mode(y_pred_class[i,])
    }
  }
  
  dfOut = cbind(y_pred_final,y_test_classes)
  
  accuracy = length(which((y_pred_final == y_test_classes) == TRUE))/length(y_test_classes)
  print(paste("accuracy:",accuracy))
  
  
  list("all" = dfOut,
       "wrongPred" = dfOut[which(dfOut[,1] != dfOut[,2]),],
       "accuracy" = accuracy,
       "wrongPredSubClasses" = NULL)
}


#------------
# manual tests
#-----------


# parameters
n = opt$n
m = opt$m
q = opt$q

ModelName = opt$ModelName

# so much of the data goes into training
trainSplitRatio = opt$trainSplitRatio

# path to the projection-file
path = opt$pathToProjection

#------------
# manual tests
n = 30
m = 100
q = 10
opt$mode = "nnQuantiles"
path = "/home/willy/PredictingProteinInteractions/data/ModelNet10/projections/"
trainSplitRatio = 0.8
numClusters = 5

ModelName = paste(path,"/NNArea.RData", sep = "") 
splitPattern = "_"
#-----------
#------------
# manual tests 2
n = 100
m = 20
q = 2
opt$mode = "nnQuantiles"
path = "/home/willy/PredictingProteinInteractions/data/animals/"
trainSplitRatio = 0.8
numClusters = 5

ModelName = paste(path,"/NNArea.RData", sep = "") 
splitPattern = "-"
#-----------



# read in the projection of all distributions of all objects
df_proj = readProjectionFromFile(n = n, m = m, q = q, path = path, fName = "proj")

#-----------------------------------------------------------------------------
# # produce models with ...
# models_all = getAllModel_F_approximations(model_vec, 72, n = n,m = m,q = q)
# df = getManhattanProjection(models_all)
# writeProjectionToFile(proj = df,n = n,m = m,q = q,path = path,fName = "proj")
#-----------------------------------------------------------------------------
numOfObjects = length(unique(df_proj[,1]))
print(paste("found ", numOfObjects, " sub-Classes", sep =""))

TrainAndTest = c()
if(opt$mode == "knn"){
  TrainAndTest = splitTrainTest(df_proj = df_proj,trainSplitRatio = trainSplitRatio)
  
  # Train model
  knnSubClasses = trainKnnOnEachDistribution(df_proj_train = TrainAndTest$train)
  
  # save the model
  saveRDS(knnSubClasses, file = ModelName)
} else if (opt$mode == "nn"){
  TrainAndTest = splitTrainTest(df_proj = df_proj,trainSplitRatio = trainSplitRatio)
  
  # Train model
  nn = trainNeuralNetOnEachDistribution(df_proj_train = TrainAndTest$train, splitPattern = splitPattern)
  
  print(paste("saving model to ", ModelName, " ...", sep =""))
  # save the model
  saveRDS(nn, file = ModelName)
}else if (opt$mode == "nnQuantiles"){

  TrainAndTest = splitTrainTest(df_proj = df_proj,trainSplitRatio = trainSplitRatio)
  
  # Train model
  nn = trainNeuralNetOnEachDistributionQuantiles(df_proj_train = TrainAndTest$train, m = m, splitPattern = splitPattern)
  
  print(paste("saving model to ", ModelName, " ...", sep =""))
  # save the model
  saveRDS(nn, file = ModelName)
}



model = readRDS(ModelName)



# evaluate model on test-set
# 1st column = subClassName
# 2nd column = n
x_test = TrainAndTest$test[,3:ncol(TrainAndTest$test)]
y_test = as.factor(TrainAndTest$test[,1])

y_test = getClassNamesFromSubClasses(y_test, splitPattern = splitPattern)

# evaluateAccuracy(knn, x_test, y_test, m)
evaluateAccuracy(model = model,x_test = x_test, y_test = y_test, m = m, probs = FALSE)



evaluateAccuracyQuantileNN(model = model,df_proj_test = TrainAndTest$test, m = m)
  