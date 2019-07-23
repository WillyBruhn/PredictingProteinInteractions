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
sampleTimes = 100


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
# example bathtub vs toilet
#

sampleSize = 20
m = 100
sampleTimes = 100


# fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_90_m_3_q_1.csv"
#
fName = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/All_ind_0_nE_1000_nD_100_n_12_m_3_q_10.csv"
quantiles = read.csv(file =fName, header = TRUE, row.names = 1)

# sub = which(getClassNamesFromSubClasses(quantiles[,1], splitPattern = "_") %in% c("bathtub", "toilet", "sofa"))
# quantiles = quantiles[sub,]

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
  layer_dense(units = 600, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 400, activation = 'relu') %>% 
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
  epochs = 50, batch_size = 30, 
  validation_split = 0.2
)

model %>% evaluate(x_test, y_test)




