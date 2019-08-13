# s3 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
# source(s3)

wsPath = as.character(paste(funr::get_script_path(), "/../../setUp/SourceLoader.R", sep = ""))

source(wsPath)
sourceFiles(c("helperFunctions"))


# install.packages("keras")
library(keras)

# install_keras(tensorflow = "gpu")
# pip install --upgrade pip

# install.packages("lime")
# library(lime)

# install.packages("tidyquant")
# library(tidyquant)

# install.packages("rsample")
# library(rsample)

# install.packages("recipes")
# library(recipes)

# install.packages("yardstick")
# library(yardstick)

# install.packages("corrr")
# library(corrr)

# install.packages("permute")
# library(permute)

# install.packages("xtable")
library(xtable)

# library(rgl)

library(permute)

transformToNodes <- function(quantiles, sSize){
  # all the ones that are merged together
  name = unique(as.character(quantiles[,1]))
  quantiles = quantiles[,-1]
  
  quantiles = quantiles[ do.call(order, quantiles), ]
  qPlusTwo = ncol(quantiles)
  
  out = data.frame(matrix(0, nrow = 1, ncol = nrow(quantiles)*qPlusTwo+1))
  
  for(i in 1:nrow(quantiles)){
    start = (i-1)*qPlusTwo+1+1
    end = start+qPlusTwo-1
    
    out[1,start:end] = as.vector(quantiles[i,])
  }
  
  out[1,1] = name
  
  colnames(out) = c("name", as.character(seq(1:(qPlusTwo*sSize))))
  return(out)
}

plotQuantilesProteins<- function(quantiles, functionals){
  t = strsplit(classes[10], split = "_")[[1]]
  paste(t[1:(length(t)-1)], collapse = '_')
  
  classes = quantiles[,1] %in% functionals
  classLevels = unique(classes)
  numOfClasses = length(classLevels)
  
  print(classLevels)
  
  print(paste("number of classes: ", numOfClasses))
  
  colMap = c("red", "blue", "green", "black")
  
  
  for(i in 1:numOfClasses){
    inds = which(classes == classLevels[i])
    points3d(quantiles[inds,-1], col = colMap[i])  
  }
}


getNNInputFromQuantiles <- function(quantiles,m, sampleSize=m, sampleTimes=1){
  # m number of distributions
  # sampleSize how many of the m distributions to select
  # sampleTimes how many times to select
  
  quantilesAsNNInput = data.frame(matrix(0,ncol = (ncol(quantiles)-1)*sampleSize+1, nrow = sampleTimes*nrow(quantiles)/m))
  for(i in 1:(nrow(quantiles)/m)){
    print(i*sampleTimes/(nrow(quantilesAsNNInput)))
    
    start = (i-1)*m+1
    end = start+m-1
    
    for(j in 1:sampleTimes){
      inds = sample(c(start:end),size = sampleSize,replace = FALSE)
      quantilesAsNNInput[(i-1)*(sampleTimes)+j,] = transformToNodes(quantiles[inds,],sampleSize)
    }
  }
  
  return(quantilesAsNNInput)
}

getTrainAndTest <- function(quantiles, sampleSize, sampleTimes, m, fName, reDo = FALSE){
  if(!file.exists(fName) || reDo == TRUE){
    subClassNames = unique(quantiles[,1])
    numObjects = length(subClassNames)
    
    numClasses = length(unique(getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")))
    
    subClassNames_test_ind = sample(c(1:numObjects), size = numObjects*0.1, replace = FALSE)
    subClassNames_test =subClassNames[subClassNames_test_ind]
    subClassNames_train =subClassNames[-subClassNames_test_ind]
    
    # return(list("train" = subClassNames_train, "test" = subClassNames_test))
    
    # NNInput_train = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_train),],
    #                                         m = m,
    #                                         sampleSize = sampleSize,
    #                                         sampleTimes = sampleTimes)
    
    NNInput_train = quantiles[which(quantiles[,1] %in% subClassNames_train),]
    
    # return(NNInput_train)
    
    # install.packages("permute")
    library(permute)
    NNInput_train = NNInput_train[shuffle(1:nrow(NNInput_train)), ]
    
    # NNInput_test = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_test),],
    #                                        m,sampleSize = sampleSize,sampleTimes = sampleTimes)
    
    NNInput_test = quantiles[which(quantiles[,1] %in% subClassNames_test),]
    
    y_train = getClassNamesFromSubClasses(NNInput_train[,1],splitPattern = "_")
    
    classLevels = unique(y_train)
    y_train = as.numeric(as.factor(y_train))-1
    x_train_prev = as.matrix(NNInput_train[,-1])
    y_train_prev <- to_categorical(y_train, numClasses)
    
    q = 1
    quants = q +2
    # sampSize = 1
    
    x_train = list()
    y_train = y_train_prev[1:(nrow(x_train_prev)/sampleSize),]
    for(i in 1:(nrow(x_train_prev)/sampleSize)){
      start_ind = (i-1)*sampleSize+1
      end_ind = start_ind+sampleSize-1
      
      # print(paste(start_ind,end_ind, i))
      x_train[[i]] = array_reshape(x_train_prev[start_ind:end_ind,],dim = c(sampleSize,quants,2), order = "F")
      
      # ?array_reshape
      
      y_train[i,] = y_train_prev[start_ind,]
      
      print("y_train_prev")
      print(y_train_prev[start_ind,])
      print("y_train")
      print(y_train[i,])
    }
    
    
    
    y_test = getClassNamesFromSubClasses(NNInput_test[,1],splitPattern = "_")
    y_test = as.numeric(as.factor(y_test))-1
    x_test_prev = as.matrix(NNInput_test[,-1])
    y_test_prev <- to_categorical(y_test, numClasses)
    
    x_test = list()
    y_test = y_test_prev[1:(nrow(x_test_prev)/sampleSize),]
    for(i in 1:(nrow(x_test_prev)/sampleSize)){
      start_ind = (i-1)*sampleSize+1
      end_ind = start_ind+sampleSize-1
      
      print(paste(start_ind,end_ind, i))
      x_test[[i]] = array_reshape(x_test_prev[start_ind:end_ind,],dim = c(sampleSize,quants,2), order = "F")
      
      y_test[i,] = y_test_prev[start_ind,]
    }
    
    
    l = list("x_train" = x_train, "y_train" = y_train, "x_test" = x_test, "y_test" = y_test, "numClasses" = numClasses)
    
    saveRDS(list(list("x_train" = x_train, "y_train" = y_train, "x_test" = x_test, "y_test" = y_test, "numClasses" = numClasses)),file = fName)
  } else {
    l = readRDS(fName)
  }
  
  return(l)
}


sampleMultipleTimes  <- function(quantiles,sampleSize, sampleTimes, numPermutations = 1){
  
  # q = ncol(quantiles)-2-1
  
  # name in first column
  quants = ncol(quantiles)-1
  
  subClassNames_train = unique(quantiles[,1])
  
  quantTrain = data.frame(matrix(0, ncol = sampleSize*quants+1, nrow = numPermutations*sampleTimes*length(subClassNames_train)))
  
  lastInd = 0
  
  colnames(quantTrain)[1] = "name"
  
  # preallocate mem
  sampled_indices = rep(1, sampleSize)
  unordered = quantiles[sampled_indices,c(1:quants)+1]
  
  # for each subclassName (that means object) draw sampleSize big sets for sampleTimes times
  for(i in 1:length(subClassNames_train)){
    subClassName_tmp = as.character(subClassNames_train[i])
    
    print(paste(subClassName_tmp, i/length(subClassNames_train)))
    
    # available points (distributions)
    train_indices = which(quantiles[,1] == subClassName_tmp)
    for(j in 1:sampleTimes){
      sampled_indices = sample(train_indices, size = sampleSize, replace = FALSE)
      
      # print(as.matrix(quantiles[sampled_indices,2:(q+3)]))
      unordered = quantiles[sampled_indices,(c(1:quants)+1)]
      if(numPermutations == 1){
        ordered = as.matrix(unordered[order(rowSums(unordered)),])
        quantTrain[j+(i-1)*sampleTimes,2:ncol(quantTrain)] = array_reshape(ordered, dim = c(1,sampleSize*quants), order = "C")
        quantTrain[j+(i-1)*sampleTimes,1] = subClassName_tmp
      } else {
        
        iInd = (i-1)*sampleTimes*numPermutations
        jInd = (j-1)*numPermutations
        for(k in 1:numPermutations){
          ind = k + ((j-1)*numPermutations + ((i-1)*sampleTimes*numPermutations))
          
          ordered = as.matrix(unordered[shuffle(nrow(unordered)),])
          quantTrain[ind,2:ncol(quantTrain)] = array_reshape(ordered, dim = c(1,sampleSize*quants), order = "C")
          quantTrain[ind,1] = subClassName_tmp
          
          if((ind)-lastInd !=1) {
            print("ERROR")
            print((ind)-lastInd )
          }
          lastInd = ind
        }
      }
    }
  }
  return(quantTrain)
}
# 
# 4*3*10
# quantos = sampleMultipleTimes(quantiles = quantiles[1:1000,],sampleSize = 3,sampleTimes = 4, numPermutations = 3)[,1]
# quantos
# 
# un = unique(quantos)
# v = rep(0,length(un))
# for(i in 1:length(un)){
#   v[i] = length(which(quantos == un[i]))
# }
# v


library(foreach)
sampleMultipleTimesParallel  <- function(quantiles,sampleSize, sampleTimes, numPermutations = 1, m = 100){
  
  # q = ncol(quantiles)-2-1
  
  # name in first column
  quants = ncol(quantiles)-1
  
  subClassNames_train = unique(quantiles[,1])
  
  quantTrain = data.frame(matrix(0, ncol = sampleSize*quants+1, nrow = numPermutations*sampleTimes*length(subClassNames_train)))
  
  lastInd = 0
  
  colnames(quantTrain)[1] = "name"
  
  relevantQuantiles = as.matrix(quantiles[,(c(1:quants)+1)])
  
  # for each subclassName (that means object) draw sampleSize big sets for sampleTimes times
  for(i in 1:length(subClassNames_train)){
    subClassName_tmp = as.character(subClassNames_train[i])
    
    print(paste(subClassName_tmp, i/length(subClassNames_train)))
    
    train_indices = ((i-1)*m+1):(i*m)
    for(j in 1:sampleTimes){
      sampled_indices = sample(train_indices, size = sampleSize, replace = FALSE)
      
      ordered = relevantQuantiles[sampled_indices,]
      quantTrain[j+(i-1)*sampleTimes,2:ncol(quantTrain)] = as.vector(t(ordered))
      quantTrain[j+(i-1)*sampleTimes,1] = subClassName_tmp
      
    }
  }
  return(quantTrain)
}

sampleMultipleTimesParallel2  <- function(quantiles,sampleSize, sampleTimes, numPermutations = 1, m = 100){
  # name in first column
  quants = ncol(quantiles)-1
  
  subClassNames_train = unique(quantiles[,1])
  
  quantTrain = data.frame(matrix(0, ncol = sampleSize*quants+1, nrow = numPermutations*sampleTimes*length(subClassNames_train)))
  colnames(quantTrain)[1] = "name"
  
  relevantQuantiles = as.matrix(quantiles[,(c(1:quants)+1)])
  
  # # for each subclassName (that means object) draw sampleSize big sets for sampleTimes times
  # for(i in 1:length(subClassNames_train)){
  #   subClassName_tmp = as.character(subClassNames_train[i])
  #   
  #   print(paste(subClassName_tmp, i/length(subClassNames_train)))
  #   
  #   train_indices = ((i-1)*m+1):(i*m)
  #   for(j in 1:sampleTimes){
  #     sampled_indices = sample(train_indices, size = sampleSize, replace = FALSE)
  #     
  #     ordered = relevantQuantiles[sampled_indices,]
  #     quantTrain[j+(i-1)*sampleTimes,2:ncol(quantTrain)] = as.vector(t(ordered))
  #     quantTrain[j+(i-1)*sampleTimes,1] = subClassName_tmp
  #     
  #   }
  # }
  
  for(k in 1:nrow(quantTrain)){
    class_ind = ceiling(k/sampleTimes)
    
    subClassName_tmp = as.character(subClassNames_train[class_ind])
    print(paste(subClassName_tmp, k/nrow(quantTrain)))
    
    train_indices = ((class_ind-1)*m+1):(class_ind*m)
    
    sampled_indices = sample(train_indices, size = sampleSize, replace = FALSE)
    
    ordered = relevantQuantiles[sampled_indices,]
    quantTrain[k,2:ncol(quantTrain)] = as.vector(t(ordered))
    quantTrain[k,1] = subClassName_tmp
  }
  
  return(quantTrain)
}


sampleMultipleTimesParallel3  <- function(quantiles,sampleSize, sampleTimes, numPermutations = 1, m = 100, sort = FALSE){
  # name in first column
  quants = ncol(quantiles)-1
  
  subClassNames_train = unique(quantiles[,1])
  
  relevantQuantiles = as.matrix(quantiles[,(c(1:quants)+1)])
  
  quantTrainMat = t(apply(matrix(c(1:(numPermutations*sampleTimes*length(subClassNames_train))), ncol = 1),1,
                          FUN = function(k){
                            class_ind= ceiling(k/sampleTimes)
                            
                            subClassName_tmp = as.character(subClassNames_train[class_ind])
                            
                            train_indices = ((class_ind-1)*m+1):(class_ind*m)
                            
                            sampled_indices = sample(train_indices, size = sampleSize, replace = FALSE)
                            
                            ordered = relevantQuantiles[sampled_indices,]
                            
                            # if(sort == TRUE){
                            #   print(ordered)
                            #   ordered = ordered[ do.call(order, as.list(ordered)), ]
                            # 
                            #   print("ordered")
                            #   print(ordered)
                            # }
                            
                            c(subClassName_tmp,as.vector(t(ordered)))
                          }))
  
  quantTrain = data.frame(quantTrainMat, stringsAsFactors=FALSE)
  colnames(quantTrain)[1] = "name"
  
  return(quantTrain)
}

# 
# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDirStandard/All_ind_Distances_nE_1000_nD_100_n_40_m_100_q_1.csv"
# quantiles = read.csv(file =fName, header = TRUE)
# 
# sampleMultipleTimesParallel3(quantiles = quantiles[1:200,-2] , sampleSize = 2,sampleTimes = 5,numPermutations = 1, m = 100)
# 
# 
# library(rbenchmark)
# num = 2000
# sampleSize = 10
# sampleTimes = 50
# numPermutations = 10
# benchmark("Parallel3" = sampleMultipleTimesParallel3(quantiles = quantiles[1:num,-2], sampleSize = sampleSize,sampleTimes = sampleTimes,numPermutations = 1),
#           "Parallel2" = sampleMultipleTimesParallel2(quantiles = quantiles[1:num,-2], sampleSize = sampleSize,sampleTimes = sampleTimes,numPermutations = 1),
#           "NotParallel" = sampleMultipleTimes(quantiles = quantiles[1:num,-2], sampleSize = sampleSize,sampleTimes = sampleTimes,numPermutations = 1),
#           replications = 1,
#           columns = c("test", "replications", "elapsed",
#                       "relative", "user.self", "sys.self"))
# 
# library(permute)
# 
# perm = diag(10)[shuffle(10),]
# 
# m = matrix(seq(1:10), ncol = 10)
# 
# m %*% perm
# 
# 
# ?aperm
# 
# library(utils)
# 
# num = 1000
# sampleSize = 10
# sampleTimes = 50
# numPermutations = 10
# 
# Rprof(tmp <- tempfile())
# m = sampleMultipleTimesParallel3(quantiles = quantiles[1:10000,], sampleSize = 10,sampleTimes = 100,numPermutations = 1, m = 100)
# Rprof()
# summaryRprof(tmp)
# 
# nrow(quantiles)
# 
# 
# sampleMultipleTimesParallel(quantiles = quantiles[1:2,], sampleSize = 2,sampleTimes = 1,numPermutations = 1, m = 2)
# sampleMultipleTimes(quantiles = quantiles[1:2,], sampleSize = 2,sampleTimes = 1,numPermutations = 1)
# 
# quantiles[1:2,]


sampleMultipleTimesMoments  <- function(quantiles,sampleSize, sampleTimes){
  moments = c(1)
  
  # q = ncol(quantiles)-2-1
  
  # name in first column
  quants = ncol(quantiles)-1
  
  subClassNames_train = unique(quantiles[,1])
  
  quantTrain = data.frame(matrix(0, ncol = quants*length(moments)+1, nrow = sampleTimes*length(subClassNames_train)))
  
  colnames(quantTrain)[1] = "name"
  # for each subclassName (that means object) draw sampleSize big sets for sampleTimes times
  for(i in 1:length(subClassNames_train)){
    subClassName_tmp = as.character(subClassNames_train[i])
    
    print(paste(subClassName_tmp, i/length(subClassNames_train)))
    
    # available points (distributions)
    train_indices = which(quantiles[,1] == subClassName_tmp)
    for(j in 1:sampleTimes){
      sampled_indices = sample(train_indices, size = sampleSize, replace = FALSE)
      quantTrain[j+(i-1)*sampleTimes,1] = subClassName_tmp
      
      # print(as.matrix(quantiles[sampled_indices,2:(q+3)]))
      unordered = quantiles[sampled_indices,(c(1:quants)+1)]
      
      
      mom = colMeans(unordered)
      
      quantTrain[j+(i-1)*sampleTimes,2:ncol(quantTrain)] = array_reshape(mom, dim = c(1,length(moments)*quants), order = "C")
    }
  }
  
  return(quantTrain[shuffle(1:nrow(quantTrain)),])
}

# m = matrix(seq(1:25), ncol = 5)
# colMeans(m)
# 
# qps = sampleMultipleTimesMoments(quantiles[,1:4], sampleSize = 3, sampleTimes = 1)
# 
# points3d(qps[,-1])
# 
# 
# # 
# # sampleMultipleTimes(quantiles = quantiles[,1:4],sampleSize = 2,sampleTimes = 2)
# 
# shuffle(seq(1:10))


getTrainAndTestOnlySurf <- function(quantiles, sampleSize, sampleTimes, numPermutations = 1, fName = "NOPE", reDo = FALSE,
                                    euklid = FALSE,
                                    TrainTestSplit = 0.1){
  # if(!file.exists(fName) || reDo == TRUE){
  
  x_test = NULL
  y_test = NULL
  
  if(euklid == FALSE){
    quantiles = quantiles[,1:(1+(ncol(quantiles)-1)/2)]
  }
  
  subClassNames = unique(quantiles[,1])
  numObjects = length(subClassNames)
  
  numClasses = length(unique(getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")))
  
  subClassNames_test = c()
  subClassNames_train = subClassNames
  if(TrainTestSplit != 0) {
    subClassNames_test_ind = sample(c(1:numObjects), size = numObjects*TrainTestSplit, replace = FALSE)
    subClassNames_test =subClassNames[subClassNames_test_ind]
    subClassNames_train =subClassNames[-subClassNames_test_ind]
  }
  
  print("starting sampling ...")
  sampledQuantiles = sampleMultipleTimes(quantiles = quantiles,sampleSize = sampleSize, sampleTimes = sampleTimes, numPermutations = numPermutations)
  
  un = unique(sampledQuantiles[,1])
  v = rep(0,length(un))
  for(i in 1:length(un)){
    v[i] = length(which(sampledQuantiles[,1] == un[i]))
  }
  
  print(v)
  if(var(v) != 0) return(NULL)
  
  print("generating train and test ...")
  
  library(permute)
  # sampledQuantiles = sampledQuantiles[shuffle(1:nrow(sampledQuantiles)), ]
  y_all = getClassNamesFromSubClasses(sampledQuantiles[,1],splitPattern = "_")
  
  train_indices = which(sampledQuantiles[,1] %in% subClassNames_train)
  test_indices = which(sampledQuantiles[,1] %in% subClassNames_test)
  
  y_train = y_all[train_indices]
  x_train = as.matrix(sampledQuantiles[train_indices,-1])
  
  classLevels = unique(y_all)
  y_train = as.numeric(as.factor(y_train))-1
  y_train <- to_categorical(y_train, numClasses)
  
  y_train_original_names = sampledQuantiles[train_indices,1]
  
  y_test_original_names = "0"
  if(TrainTestSplit > 0){
    y_test= y_all[test_indices]
    x_test = as.matrix(sampledQuantiles[test_indices,-1])
    
    y_test = as.numeric(as.factor(y_test))-1
    y_test <- to_categorical(y_test, numClasses)
    
    y_test_original_names = sampledQuantiles[test_indices,1]
    
    un = unique(y_test_original_names)
    v = rep(0,length(un))
    for(i in 1:length(un)){
      v[i] = length(which(y_test_original_names == un[i]))
    }
    if(length(v) > 0 && var(v) != 0) return(NULL)
  }
  
  if(is.null(x_test)) {
    x_test = "onlyTrain"
    y_test = "onlyTrain"
  }
  
  return(list("x_train" = x_train, "y_train" = y_train, "x_test" = x_test, "y_test" = y_test, "numClasses" = numClasses, "y_test_original_names" = y_test_original_names ,
              "y_train_original_names" = y_train_original_names))
}


getSamplesSurf2 <- function( quantiles,
                             sampleSize, 
                             sampleTimes,
                             numPermutations = 1,
                             fName = "NOPE",
                             reDo = FALSE,
                             euklid = FALSE,
                             numClasses,
                             m,
                             splitPattern = "_",
                             sort = FALSE){
  
  # the rest of quantiles is euclidean distances
  if(euklid == FALSE){
    quantiles = quantiles[,1:(1+(ncol(quantiles)-1)/2)]
  }
  
  print("starting sampling ...")
  sampledQuantiles = sampleMultipleTimesParallel3(quantiles = quantiles,
                                                  sampleSize = sampleSize,
                                                  sampleTimes = sampleTimes, 
                                                  numPermutations = numPermutations,
                                                  m = m,
                                                  sort = sort)
  
  un = unique(sampledQuantiles[,1])
  v = rep(0,length(un))
  for(i in 1:length(un)){
    v[i] = length(which(sampledQuantiles[,1] == un[i]))
  }
  
  print(v)
  if(var(v) != 0) return(NULL)
  
  print("generating x and y ...")
  
  y_all_class_names = getClassNamesFromSubClasses(sampledQuantiles[,1],splitPattern = splitPattern)
  
  print(y_all_class_names)
  
  y = y_all_class_names
  X = as.matrix(sampledQuantiles[,-1])
  
  classLevels = unique(y)
  y = as.numeric(as.factor(y))-1
  y <- to_categorical(y, numClasses)
  
  y_original_names = sampledQuantiles[,1]
  
  return(list("X" = X, "y" = y, "y_original_names" = y_original_names))
}

# getTrainAndTestOnlySurf(quantiles[1:2000,],sampleSize = 3,sampleTimes = 3, numPermutations = 2)
# 
# quantiles[1,]


convModel1 <- function(TrainTest, sampleSize, sampleTimes, m, q, epochs = 30){
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_reshape(target_shape = c(sampleSize, q+2,1),input_shape = c(ncol(x_train))) %>% 
    layer_conv_2d(filter = 30, kernel_size = c(3,1), padding = 'same', input_shape = c(sampleSize, q+2,1)) %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.25) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(4,1), padding = 'same', input_shape = c(sampleSize, q+2,1)) %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.25) %>%
    
    layer_flatten() %>%
    layer_dense(100) %>%
    layer_dense(100) %>%
    layer_activation("relu") %>%
    layer_dropout(0.5) %>%
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.5) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  # ?layer_global_average_pooling_1d
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = 10, 
    validation_split = 0.2
  )
  
  
  
  # $loss
  # [1] 0.7399024
  # 
  # $acc
  # [1] 0.6863636
  print(model %>% evaluate(x_test, y_test))
  return(model)
}

convModel2 <- function(TrainTest, sampleSize, sampleTimes, m, q, epochs = 30){
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_reshape(target_shape = c(sampleSize, q+2,1),input_shape = c(ncol(x_train))) %>% 
    
    layer_conv_2d(filter = 30, kernel_size = c(2,1), padding = 'same', input_shape = c(sampleSize, q+2,1)) %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(5,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(6,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(6,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(6,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(6,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(6,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(6,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(6,3), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(6,1), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(4,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_flatten() %>%
    layer_dense(30) %>%
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  # ?layer_global_average_pooling_1d
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = 16, 
    validation_split = 0.2
  )
  
  
  
  # $loss
  # [1] 0.7399024
  # 
  # $acc
  # [1] 0.6863636
  print(model %>% evaluate(x_test, y_test))
  return(model)
}


convModel3 <- function(TrainTest, sampleSize, sampleTimes, m, q, epochs = 30){
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_reshape(target_shape = c(sampleSize, q+2,1),input_shape = c(ncol(x_train))) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(2,1), padding = 'same', input_shape = c(sampleSize, q+2,2)) %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(2,1), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(2,1), padding = 'same') %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_flatten() %>%
    layer_dense(100) %>%
    layer_dense(50) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  # ?layer_global_average_pooling_1d
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = 32, 
    validation_split = 0.2
  )
  
  
  
  # $loss
  # [1] 0.7399024
  # 
  # $acc
  # [1] 0.6863636
  print(model %>% evaluate(x_test, y_test))
  return(model)
}


convModel4 <- function(TrainTest, sampleSize, sampleTimes, q, epochs = 30, batch_size = 64){
  
  
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_reshape(target_shape = c(sampleSize, ncol(x_train)/sampleSize,1),input_shape = c(ncol(x_train))) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same', input_shape = c(sampleSize, ncol(x_train)/sampleSize,1)) %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 10, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 10, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 10, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 10, kernel_size = c(4,4), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(3,3),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    
    layer_flatten() %>%
    layer_dense(400) %>%
    layer_dense(200) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    
    layer_dense(100) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(100) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(units = numClasses, activation = 'softmax')
  
  # ?layer_global_average_pooling_1d
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    validation_split = 0.2
  )
  
  
  
  # $loss
  # [1] 0.7399024
  # 
  # $acc
  # [1] 0.6863636
  print(model %>% evaluate(x_test, y_test))
  
  
  
  
  # pred <- data.frame(y = predict(model, as.matrix(x_test)))
  
  testNamesOrig = unique(Tr$y_test_original_names)
  
  # correct = 0
  # 
  # for(i in 1:length(testNamesOrig)){
  #   origNameInds = which(Tr$y_test_original_names == testNamesOrig[i])
  #   prediction = which.max(colSums(pred[origNameInds,])/length(origNameInds))
  #   actual = which.max(colSums(y_test[origNameInds,])/length(origNameInds))
  #   if(prediction == actual) correct = correct +1
  # }
  # 
  # print(correct/length(testNamesOrig))
  
  
  
  return(model)
}

model5 <- function(TrainTest, epochs = 30, batch_size = 64, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model5 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 500, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 400, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    validation_split = 0.2
  )
  
  print(model %>% evaluate(x_test, y_test))
  
  return(model)
}



model6 <- function(TrainTest, epochs = 30, batch_size = 64, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model6 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 400, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 300, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    validation_split = 0.2
  )
  
  if(!is.null(x_test)) print(model %>% evaluate(x_test, y_test))
  
  return(model)
}


model7 <- function(TrainTest, epochs = 30, batch_size = 64, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model7 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 200, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 50, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    validation_split = 0.2
  )
  
  print(model %>% evaluate(x_test, y_test))
  
  return(model)
}

library(DescTools)

reverseToCategorical <- function(oneHot, names){
  names_out = rep("", nrow(oneHot))
  for(i in 1:nrow(oneHot)){
    names_out[i] = names[which.max(oneHot[i,])]
  }
  
  return(names_out)
}



