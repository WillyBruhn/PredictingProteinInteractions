s3 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
source(s3)


library(keras)

# install.packages("lime")
library(lime)

# install.packages("tidyquant")
library(tidyquant)

# install.packages("rsample")
# library(rsample)

# install.packages("recipes")
library(recipes)

# install.packages("yardstick")
library(yardstick)

# install.packages("corrr")
# library(corrr)

library(permute)



# #-------------------------------------------------------------------
# # MNist
# #-------------------------------------------------------------------
# mnist <- dataset_mnist()
# x_train <- mnist$train$x
# y_train <- mnist$train$y
# x_test <- mnist$test$x
# y_test <- mnist$test$y
# 
# 
# 
# # reshape
# x_train <- array_reshape(x_train, c(nrow(x_train), 784))
# x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# # rescale
# x_train <- x_train / 255
# x_test <- x_test / 255
# 
# 
# y_train <- to_categorical(y_train, 10)
# y_test <- to_categorical(y_test, 10)
# 
# 
# 
# model <- keras_model_sequential() 
# model %>% 
#   layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>% 
#   layer_dropout(rate = 0.4) %>% 
#   layer_dense(units = 128, activation = 'relu') %>%
#   layer_dropout(rate = 0.3) %>%
#   layer_dense(units = 10, activation = 'softmax')
# 
# 
# model %>% compile(
#   loss = 'categorical_crossentropy',
#   optimizer = optimizer_rmsprop(),
#   metrics = c('accuracy')
# )
# 
# 
# history <- model %>% fit(
#   x_train, y_train, 
#   epochs = 30, batch_size = 128, 
#   validation_split = 0.2
# )
# 
# plot(history)
# 
# 
# model %>% evaluate(x_test, y_test)
# 
# model %>% predict_classes(x_test)
# #-------------------------------------------------------------------
# # MNist
# #-------------------------------------------------------------------


#----------------------------------------------------------
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

#------------------------------------------------------------------
library(rgl)

n = 40
m = 100
q = 1

sampleSize = 20
sampleTimes = 200

# path = "/home/willy/PredictingProteinInteractions/data/animals/"
# fName = "quantiles"
# quantiles = readQuantilesFromFile(path = path, fName = fName, n = n, m= m, q = q)

quantiles = read.csv(file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_200_nS_50_n_48_m_100_q_1.csv")


# df = quantiles
# #
# # df = readProjectionFromFile(n = n,m = m,q = q,path = "/home/willy/PredictingProteinInteractions/data/animals/",fName = "proj")
# classes = unique(getClassNamesFromSubClasses(df[,1]))
# colors = as.numeric(as.factor(classes))
# colorMap = c("red", "blue", "green", "yellow", "black", "pink", "orange")
# for(i in 1:length(colorMap)){
#   print(paste(classes[i], colorMap[i]))
# 
#   if(!(classes[i] == "head")){
#     inds = which(getClassNamesFromSubClasses(df[,1]) == classes[i])
#     # print(inds)
#     points3d(x = df[inds,2], y = df[inds,3], z = df[inds,4],colorMap[i], add = TRUE,size = 10)
# 
#     geox = sum(df[inds,2])/length(inds)
#     geoy = sum(df[inds,3])/length(inds)
#     geoz = sum(df[inds,4])/length(inds)
# 
#     textCoord = c(geox-0.1,geoy-0.1,geoz)
#     if(classes[i] == "horse") {
#       textCoord = c(geox+0.1,geoy+0.1,geoz)
#       text3d(x = textCoord[1]+0.05, y = textCoord[2]+0.05, z = textCoord[3],texts = classes[i],cex = 2)
#       
#       arrow3d(p0 = textCoord, c(geox+0.05,geoy+0.05,geoz), col = "black")
#     }else{
#       text3d(x = textCoord[1]-0.05, y = textCoord[2]-0.05, z = textCoord[3],texts = classes[i],cex = 2)
#       
#       arrow3d(p0 = textCoord, c(geox-0.05,geoy-0.05,geoz), col = "black")
#     }
#     
# 
# 
#   }
# }
# 
# rgl.snapshot("/home/willy/PredictingProteinInteractions/Results/Images/animals3dProjectionExample.png")




subClassNames = unique(quantiles[,1])
numObjects = length(subClassNames)

numClasses = length(unique(getClassNamesFromSubClasses(quantiles[,1],splitPattern = "-")))

subClassNames_test_ind = sample(c(1:numObjects), size = numObjects*0.3, replace = FALSE)
subClassNames_test =subClassNames[subClassNames_test_ind]
subClassNames_train =subClassNames[-subClassNames_test_ind]


NNInput_train = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_train),],
                                        m,sampleSize = sampleSize,sampleTimes = sampleTimes)

# install.packages("permute")
library(permute)
NNInput_train = NNInput_train[shuffle(1:nrow(NNInput_train)), ]

NNInput_test = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_test),],
                                        m,sampleSize = sampleSize,sampleTimes = sampleTimes)

y_train = getClassNamesFromSubClasses(NNInput_train[,1],splitPattern = "-")

classLevels = unique(y_train)
y_train = as.numeric(as.factor(y_train))-1
x_train = as.matrix(NNInput_train[,-1])
y_train <- to_categorical(y_train, numClasses)


# unique(NNInput_test[,1])

y_test = getClassNamesFromSubClasses(NNInput_test[,1],splitPattern = "-")
y_test = as.numeric(as.factor(y_test))-1
x_test = as.matrix(NNInput_test[,-1])
y_test <- to_categorical(y_test, numClasses)


# y_test <- to_categorical(y_test, 10)

#---------------------------------------------------------
model <- keras_model_sequential()
model %>% 
  layer_dense(units = 300, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
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

# history <- model %>% fit(
#   x_train, y_train, 
#   epochs = 30, batch_size = 10, 
#   validation_split = 0.01
# )

history <- model %>% fit(
  x_train, y_train, 
  epochs = 30, batch_size = 10, 
  validation_split = 0.2
)

model %>% evaluate(x_test, y_test)

# 4200/4200 [==============================] - 0s 27us/sample - loss: 0.0453 - acc: 0.9936
# $loss
# [1] 0.04529072
# 
# $acc
# [1] 0.9935714

#-------------------------------------------------------------------------------------------------------------
# ModelNet10
# example bathtub vs toilet
#

sampleSize = 20
m = 100
sampleTimes = 100


# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_90_m_3_q_1.csv"
#

# this worked !
# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_12_m_3_q_1.csv"

fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_12_m_3_q_10.csv"
quantiles = read.csv(file =fName, header = TRUE, row.names = 1)

sub = which(getClassNamesFromSubClasses(quantiles[,1], splitPattern = "_") %in% c("bathtub", "toilet"))
quantiles = quantiles[sub,]

subClassNames = unique(quantiles[,1])
numObjects = length(subClassNames)

numClasses = length(unique(getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")))

subClassNames_test_ind = sample(c(1:numObjects), size = numObjects*0.1, replace = FALSE)
subClassNames_test =subClassNames[subClassNames_test_ind]
subClassNames_train =subClassNames[-subClassNames_test_ind]


NNInput_train = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_train),],
                                        m = m,
                                        sampleSize = sampleSize,
                                        sampleTimes = sampleTimes)

# install.packages("permute")
library(permute)
NNInput_train = NNInput_train[shuffle(1:nrow(NNInput_train)), ]

NNInput_test = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_test),],
                                       m,sampleSize = sampleSize,sampleTimes = sampleTimes)

y_train = getClassNamesFromSubClasses(NNInput_train[,1],splitPattern = "_")

classLevels = unique(y_train)
y_train = as.numeric(as.factor(y_train))-1
x_train = as.matrix(NNInput_train[,-1])
y_train <- to_categorical(y_train, numClasses)



y_test = getClassNamesFromSubClasses(NNInput_test[,1],splitPattern = "_")
y_test = as.numeric(as.factor(y_test))-1
x_test = as.matrix(NNInput_test[,-1])
y_test <- to_categorical(y_test, numClasses)


# y_test <- to_categorical(y_test, 10)

#---------------------------------------------------------
model <- keras_model_sequential()
model %>% 
  layer_dense(units = 100, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
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

# history <- model %>% fit(
#   x_train, y_train, 
#   epochs = 30, batch_size = 10, 
#   validation_split = 0.01
# )

history <- model %>% fit(
  x_train, y_train, 
  epochs = 50, batch_size = 10, 
  validation_split = 0.2
)

model %>% evaluate(x_test, y_test)
# $loss
# [1] 0.3235773
# 
# $acc
# [1] 0.865

#-------------------------------------------------------------------------------------------------------------
# ModelNet10
# example bathtub vs toilet
#

sampleSize = 20
m = 100
sampleTimes = 10


# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_90_m_3_q_1.csv"
#
fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_12_m_3_q_10.csv"
quantiles = read.csv(file =fName, header = TRUE, row.names = 1)

sub = which(getClassNamesFromSubClasses(quantiles[,1], splitPattern = "_") %in% c("bathtub", "toilet", "sofa"))
quantiles = quantiles[sub,]

subClassNames = unique(quantiles[,1])
numObjects = length(subClassNames)

numClasses = length(unique(getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")))

subClassNames_test_ind = sample(c(1:numObjects), size = numObjects*0.1, replace = FALSE)
subClassNames_test =subClassNames[subClassNames_test_ind]
subClassNames_train =subClassNames[-subClassNames_test_ind]


NNInput_train = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_train),],
                                        m = m,
                                        sampleSize = sampleSize,
                                        sampleTimes = sampleTimes)

# install.packages("permute")
library(permute)
NNInput_train = NNInput_train[shuffle(1:nrow(NNInput_train)), ]

NNInput_test = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_test),],
                                       m,sampleSize = sampleSize,sampleTimes = sampleTimes)

y_train = getClassNamesFromSubClasses(NNInput_train[,1],splitPattern = "_")

classLevels = unique(y_train)
y_train = as.numeric(as.factor(y_train))-1
x_train = as.matrix(NNInput_train[,-1])
y_train <- to_categorical(y_train, numClasses)



y_test = getClassNamesFromSubClasses(NNInput_test[,1],splitPattern = "_")
y_test = as.numeric(as.factor(y_test))-1
x_test = as.matrix(NNInput_test[,-1])
y_test <- to_categorical(y_test, numClasses)


# y_test <- to_categorical(y_test, 10)

#---------------------------------------------------------
model <- keras_model_sequential()
model %>% 
  layer_dense(units = 300, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 100, activation = 'relu') %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 100, activation = 'relu') %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10, activation = 'relu') %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = numClasses, activation = 'softmax')

# ?layer_global_average_pooling_1d

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

# history <- model %>% fit(
#   x_train, y_train, 
#   epochs = 30, batch_size = 10, 
#   validation_split = 0.01
# )

history <- model %>% fit(
  x_train, y_train, 
  epochs = 50, batch_size = 10, 
  validation_split = 0.2
)

model %>% evaluate(x_test, y_test)

# $loss
# [1] 0.4404805
# 
# $acc
# [1] 0.8

#---------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
# ModelNet10
#

sampleSize = 20
m = 100
sampleTimes = 10
q = 1


# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_90_m_3_q_1.csv"
#
fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDirStandard/All_ind_0_nE_1000_nD_100_n_80_m_100_q_1.csv"
quantiles = read.csv(file =fName, header = TRUE, row.names = 1)

sub = which(getClassNamesFromSubClasses(quantiles[,1], splitPattern = "_") %in% c("bathtub", "toilet"))
quantiles = quantiles[sub,]

subClassNames = unique(quantiles[,1])
numObjects = length(subClassNames)

numClasses = length(unique(getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")))

subClassNames_test_ind = sample(c(1:numObjects), size = numObjects*0.1, replace = FALSE)
subClassNames_test =subClassNames[subClassNames_test_ind]
subClassNames_train =subClassNames[-subClassNames_test_ind]


quantiles[1:200,]


NNInput_train = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_train),],
                                        m = m,
                                        sampleSize = sampleSize,
                                        sampleTimes = sampleTimes)

# install.packages("permute")
library(permute)
NNInput_train = NNInput_train[shuffle(1:nrow(NNInput_train)), ]

NNInput_train[1,]

NNInput_test = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_test),],
                                       m,sampleSize = sampleSize,sampleTimes = sampleTimes)

y_train = getClassNamesFromSubClasses(NNInput_train[,1],splitPattern = "_")

classLevels = unique(y_train)
y_train = as.numeric(as.factor(y_train))-1
x_train = as.matrix(NNInput_train[,-1])
y_train <- to_categorical(y_train, numClasses)



y_test = getClassNamesFromSubClasses(NNInput_test[,1],splitPattern = "_")
y_test = as.numeric(as.factor(y_test))-1
x_test = as.matrix(NNInput_test[,-1])
y_test <- to_categorical(y_test, numClasses)


# y_test <- to_categorical(y_test, 10)


?layer_reshape
#---------------------------------------------------------
model <- keras_model_sequential()
model %>% 
  layer_reshape(target_shape = c(sampleSize, q+2,1),input_shape = c(ncol(x_train))) %>% 
  layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same', input_shape = c(sampleSize, q+2,1)) %>% 
  layer_activation("relu") %>%
  layer_max_pooling_2d(pool_size = c(2,2),padding = 'same') %>%
  layer_dropout(0.25) %>%
  
  layer_conv_2d(filter = 10, kernel_size = c(5,1), padding = 'same', input_shape = c(sampleSize, q+2,1)) %>% 
  layer_activation("relu") %>%
  layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_dropout(0.25) %>%
  
  layer_flatten() %>%
  layer_dense(512) %>%
  layer_activation("relu") %>%
  layer_dropout(0.5) %>%
  layer_dense(units = numClasses, activation = 'softmax')

# ?layer_global_average_pooling_1d

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

# history <- model %>% fit(
#   x_train, y_train, 
#   epochs = 30, batch_size = 10, 
#   validation_split = 0.01
# )

history <- model %>% fit(
  x_train, y_train, 
  epochs = 50, batch_size = 10, 
  validation_split = 0.2
)

model %>% evaluate(x_test, y_test)
# $loss
# [1] 0.3933986
# 
# $acc
# [1] 0.875

#---------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
# ModelNet10
#

sampleSize = 20
m = 200
sampleTimes = 10
q = 1


# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_90_m_3_q_1.csv"
#
fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDirStandard/All_ind_0_nE_1000_nD_100_n_10_m_200_q_1.csv"
quantiles = read.csv(file =fName, header = TRUE, row.names = 1)

nrow(quantiles)

# plotQuantiles(quantiles)

# TrainTest2 = TrainTest

TrainTest = getTrainAndTest(quantiles, sampleSize = sampleSize, sampleTimes = sampleTimes, m = m,
                            fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/tmpTrainTest2.Rdata")


# convModel1(TrainTest, sampleSize = sampleSize, sampleTimes = sampleTimes, m = m, q = q, epochs = 30)
# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDirStandard/All_ind_0_nE_1000_nD_100_n_90_m_100_q_1.csv"
# "toilet", "sofa"
#
# $loss
# [1] 0.3118481
# 
# $acc
# [1] 0.8571429

# ncol(TrainTest[[1]]$x_train)
600 /sampleSize


convModel2(TrainTest, sampleSize = sampleSize, sampleTimes = sampleTimes, m = m, q = q, epochs = 30)

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

# m = matrix(seq(1:25), nrow = 5)
# array_reshape(m, dim = c(1,25), order = "F")
# 
# as.matrix(quantiles[1:3,2:4])
# 
# array_reshape(as.matrix(quantiles[1:3,2:4]), dim = c(1,9), order = "C")


# order(rowSums(quantiles[,-1]))



sampleMultipleTimes  <- function(quantiles,sampleSize, sampleTimes){
  
  # q = ncol(quantiles)-2-1
  
  # name in first column
  quants = ncol(quantiles)-1
  
  subClassNames_train = unique(quantiles[,1])
  
  quantTrain = data.frame(matrix(0, ncol = sampleSize*quants+1, nrow = sampleTimes*length(subClassNames_train)))

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
      # ordered = as.matrix(unordered[ do.call(order, unordered), ])
      
      # print(unordered)
      # print(rowSums(unordered))
      # print(order(rowSums(unordered)))
      
      ordered = as.matrix(unordered[order(rowSums(unordered)),])
      
      # print(ordered)
      
      quantTrain[j+(i-1)*sampleTimes,2:ncol(quantTrain)] = array_reshape(ordered, dim = c(1,sampleSize*quants), order = "C")
    }
  }
  
  
  
  return(quantTrain[shuffle(1:nrow(quantTrain)),])
}

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


getTrainAndTestOnlySurf <- function(quantiles, sampleSize, sampleTimes, fName = "NOPE", reDo = FALSE, euklid = FALSE){
  # if(!file.exists(fName) || reDo == TRUE){
  
    if(euklid == FALSE){
      quantiles = quantiles[,1:(1+(ncol(quantiles)-1)/2)]
    }
  
    subClassNames = unique(quantiles[,1])
    numObjects = length(subClassNames)
    
    numClasses = length(unique(getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")))
    
    subClassNames_test_ind = sample(c(1:numObjects), size = numObjects*0.1, replace = FALSE)
    subClassNames_test =subClassNames[subClassNames_test_ind]
    subClassNames_train =subClassNames[-subClassNames_test_ind]
    
    print("starting sampling ...")
    sampledQuantiles = sampleMultipleTimes(quantiles = quantiles,sampleSize = sampleSize, sampleTimes = sampleTimes)

    library(permute)
    sampledQuantiles = sampledQuantiles[shuffle(1:nrow(sampledQuantiles)), ]
    y_all = getClassNamesFromSubClasses(sampledQuantiles[,1],splitPattern = "_")
    
    train_indices = which(sampledQuantiles[,1] %in% subClassNames_train)
    test_indices = which(sampledQuantiles[,1] %in% subClassNames_test)
    
    y_train = y_all[train_indices]
    x_train = as.matrix(sampledQuantiles[train_indices,-1])
    
    y_test= y_all[test_indices]
    x_test = as.matrix(sampledQuantiles[test_indices,-1])
    
    classLevels = unique(y_all)
    y_train = as.numeric(as.factor(y_train))-1
    y_train <- to_categorical(y_train, numClasses)
    
    y_test = as.numeric(as.factor(y_test))-1
    y_test <- to_categorical(y_test, numClasses)
  
  return(list("x_train" = x_train, "y_train" = y_train, "x_test" = x_test, "y_test" = y_test, "numClasses" = numClasses))
}

# getTrainAndTestOnlySurf(quantiles[,1:4],sampleSize = 2,sampleTimes = 1)


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


convModel4 <- function(TrainTest, sampleSize, sampleTimes, q, epochs = 30){
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_reshape(target_shape = c(sampleSize, ncol(x_train)/sampleSize,1),input_shape = c(ncol(x_train))) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(10,1), padding = 'same', input_shape = c(sampleSize, ncol(x_train)/sampleSize,1)) %>% 
    layer_activation("relu") %>%
    layer_average_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same', input_shape = c(sampleSize, ncol(x_train)/sampleSize,1)) %>% 
    layer_activation("relu") %>%
    layer_average_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(10,3), padding = 'same', input_shape = c(sampleSize, ncol(x_train)/sampleSize,1)) %>% 
    layer_activation("relu") %>%
    layer_average_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(1,2), padding = 'same', input_shape = c(sampleSize, ncol(x_train)/sampleSize,1)) %>% 
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same', input_shape = c(sampleSize, ncol(x_train)/sampleSize,1)) %>% 
    layer_activation("relu") %>%
    layer_average_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
  

    
    layer_flatten() %>%
    layer_dense(400) %>%
    layer_dense(200) %>%
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
    epochs = epochs, batch_size = 8, 
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



#---------------------------------------------------------------------------------------------------------
#
# euclidean and geodesic
#
#---------------------------------------------------------------------------------------------------------

sampleSize = 100
sampleTimes = 1
q = 10


# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_90_m_3_q_1.csv"
#
fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDirStandard/All_ind_0_nE_1000_nD_50_n_40_m_100_q_10.csv"
quantiles = read.csv(file =fName, header = TRUE, row.names = 1)

unique(getClassNamesFromSubClasses(quantiles[,1],"_"))
# "bathtub" "chair"   "sofa"    "table"   "toilet" 

# "bathtub" "chair"
# $acc
# [1] 0.91

# sub = which(getClassNamesFromSubClasses(quantiles[,1],"_") %in% c("bathtub","chair"))
# quantiles = quantiles[sub,]

# Tr = getTrainAndTestOnlySurf(quantiles[,1:(q+3)],sampleSize = sampleSize,sampleTimes = sampleTimes)
Tr = getTrainAndTestOnlySurf(quantiles,sampleSize = sampleSize,sampleTimes = sampleTimes,euklid = FALSE)


plotQuantiles(quantiles,euklid = TRUE)

array_reshape(Tr$x_train[1,], dim = c(sampleSize,(q+2)))
# convModel4(TrainTest = Tr,sampleSize = sampleSize,sampleTimes = sampleTimes,m = m,q = q+3,epochs = 100)
convModel4(TrainTest = Tr,sampleSize = sampleSize,sampleTimes = sampleTimes,q = (q+2),epochs = 300)




