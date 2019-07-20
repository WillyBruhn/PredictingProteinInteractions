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
    
    # print(paste(start,end))
    # print(out[1,])
    # print(quantiles[i,])
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


n = 100
m = 1000
q = 1

sampleSize = 20
sampleTimes = 100

path = "/home/willy/PredictingProteinInteractions/data/animals/"
fName = "quantiles"
quantiles = readQuantilesFromFile(path = path, fName = fName, n = n, m= m, q = q)






NNInput = getNNInputFromQuantiles(quantiles,m,sampleSize = sampleSize,sampleTimes = sampleTimes)

y = getClassNamesFromSubClasses(NNInput[,1],splitPattern = "-")

classLevels = unique(y)
numClasses = length(classLevels)


y_train = as.numeric(as.factor(y))-1
x_train = as.matrix(NNInput[,-1])


y_train <- to_categorical(y_train, numClasses)

# y_test <- to_categorical(y_test, 10)

#---------------------------------------------------------
model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 120, activation = 'relu', input_shape = c(ncol(NNInput)-1)) %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10, activation = 'relu') %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10, activation = 'relu') %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10, activation = 'relu') %>% 
  layer_dense(units = numClasses, activation = 'softmax')

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

history <- model %>% fit(
  x_train, y_train, 
  epochs = 100, batch_size = 10, 
  validation_split = 0.4
)


history$metrics$val_los


