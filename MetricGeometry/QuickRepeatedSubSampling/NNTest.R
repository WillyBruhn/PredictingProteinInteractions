#-----------------------------
# Willy Bruhn, 19.7.19
#
# Testing if a nn can distinguish 
# distributions in 2d.
#
#-----------------------------

library(randomForest)
library(mlbench)
library(caret)
library(doParallel)
library(doBy)
library(rfUtilities)
library(keras)

unfold2dProj <- function(mat){
  # x, y
  
  distributions = data.frame(matrix(0,ncol = 2, nrow = nrow(mat)))
  colnames(distributions) = c("q_1", "q_2")
  
  for(i in 1:nrow(mat)){
    distributions[i,1] = mat[i,1]   
    distributions[i,2] = mat[i,2]+mat[i,1]
  }
  
  return(distributions)
}

unfoldProjArbitraryDim <- function(mat){
  # x, y
  
  distributions = data.frame(matrix(0,ncol = ncol(mat), nrow = nrow(mat)))
  colnames(distributions) = c("class", seq(1:(ncol(mat)-1)))
  
  for(i in 1:nrow(mat)){
    # name
    distributions[i,1] = as.character(mat[i,1])
    
    # 1st quantile
    distributions[i,2] = abs(mat[i,2])
    for(j in 3:ncol(mat))
    distributions[i,j] = distributions[i,j-1] + mat[i,j]
  }
  
  return(distributions)
}


sampleDistributions <- function(mat1, k, classLabel, m, plot = FALSE, col = NULL, order = TRUE){
  # mat1 ... points in 2d-space
  # simulate sampling m distributions  for k times
  #---------------------------------------------------
  
  dataOut = data.frame(matrix(0, ncol = 2*m+1, nrow = k))
  colnames(dataOut)[1] = "class"
  
  dataOut[,1] = classLabel
  
  for(j in 1:k){
    mat1_sample = mat1[c(1:nrow(mat1))[sample(c(1:nrow(mat1)), size = m, replace = FALSE)],]
    if(plot == TRUE){
      points(x = mat1_sample[,1], y = mat1_sample[,2], pch = 19, cex = 1, col = col)
    }
      
    print(mat1_sample)
    distributions_ordered = mat1_sample
    if(length(mat1_sample) > 2){
      distributions = unfold2dProj(mat1_sample)
      distributions_ordered = distributions
      if(order == TRUE){
        distributions_ordered = distributions[order(distributions[,1]),]
      }
      
      
      # input-nodes
      inputNodes = matrix(0, ncol = nrow(distributions_ordered)*2)
      inputNodes
      
      for( i in 1:nrow(distributions_ordered)){
        inputNodes[(i-1)*2+1] = distributions_ordered[i,1]
        inputNodes[(i-1)*2+1+1] = distributions_ordered[i,2]
      }
      
      dataOut[j,2:ncol(dataOut)] = inputNodes
    } else {
      dataOut[j,2:ncol(dataOut)] = matrix(mat1_sample, ncol = 2)
    }
  }
  
  return(dataOut)
}



generateExample1 <- function(){
  # emd vs geo
  n1 = 1000
  n2 = 1000
  
  mean1 = c(1,1)
  sd1 = c(1,5)
  mean2 = c(1,1)
  sd2 = c(2,1)
  
  mat1 = matrix(c(rnorm(n1,mean = mean1[1], sd = sd1[1]),rnorm(n1,mean = mean1[2],sd = sd1[2])), ncol = 2, byrow = FALSE)
  mat2 = matrix(c(rnorm(n2,mean = mean2[1], sd = sd2[1]),rnorm(n2,mean = mean2[2],sd = sd2[2])), ncol = 2, byrow = FALSE)
  
  mat_both = rbind(mat1,mat2)
  plot(mat_both)
  points(x = mat1[,1], y = mat1[,2], col ="red", pch = 19)
  points(x = mat2[,1], y = mat2[,2], col ="blue", pch = 19)
  
  
  # generate training set and test set
  k = 200
  splitRatio = 0.8
  
  X = sampleDistributions(mat1 = mat1,k = k,"X",m = 10)
  Y = sampleDistributions(mat1 = mat2,k = k,"Y",m = 10)
  
  All = rbind(X,Y)
  
  return(All)
}


fitNNModel <- function(All, splitRatio = 0.8){
  a <- createDataPartition(All$class, p = splitRatio, list=FALSE)
  training <- All[a,]
  test <- All[-a,]
  
  metric <- "Accuracy"
  numFolds <- trainControl(method = 'cv', number = 10, classProbs = TRUE, verboseIter = TRUE, preProcOptions = list(thresh = 0.75, ICAcomp = 3, k = 5))
  
  cl <- makePSOCKcluster(5)
  registerDoParallel(cl)
  nn <- train(class~., data=training, method="nnet", metric=metric, trControl=numFolds, tuneGrid=expand.grid(size=c(10), decay=c(0.1)))
  stopCluster(cl)
  
  print(nn)
  
  accuracy(predict(object = nn, newdata = test),test[,1])
}


#-----------------------------------------------------------------------------------------------
# Example 1
All = generateExample1()

# fitNNModel(All)


#-----------------------------------------------------------------------------------------------
# Example2

getPoints <- function(n1,mean1,sd1){
  return(matrix(c(rnorm(n1,mean = mean1[1], sd = sd1[1]),rnorm(n1,mean = mean1[2],sd = sd1[2])), ncol = 2, byrow = FALSE))
}

X1 = getPoints(1000, c(0,0), c(1,1))
X2 = getPoints(1000, c(10,0), c(1,1))
X3 = getPoints(1000, c(0,4), c(1,3))
Y1 = getPoints(n1 = 1000, mean = c(1,1), sd = c(1,2))

mat_both = rbind(X1,X2,X3,Y1)

# pdf("/home/willy/PredictingProteinInteractions/Results/Images/NN2dExample.pdf")
par(mfrow = c(1,3))
xli = c(-6,15)
yli = c(-6,15)
plot(mat_both, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
plot(Y1, col = "blue", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
plot(mat_both, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)
legend(x = 5,y = 10, legend = c("X","Y"),col = c("red", "blue"),pch = 19)
# dev.off()


# pdf("/home/willy/PredictingProteinInteractions/Results/Images/NN2dExampleSampleX1vsY1.pdf")
par(mfrow = c(1,2))
plot(mat_both, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)
legend(x = 5,y = 10, legend = c("X","Y"),col = c("red", "blue"),pch = 19)

Y_demo = sampleDistributions(mat1 = Y1,k = 1,"Y",m = 20,TRUE, col = "green")

# par(mfrow = c(1,1))
plot(mat_both, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)
legend(x = 5,y = 10, legend = c("X","Y"),col = c("red", "blue"),pch = 19)
X1_demo = sampleDistributions(mat1 = X1,k = 1,"X",m = 20,TRUE, col ="green")
# dev.off()

# generate training set and test set
experimentWithModel <- function(k = 200, splitRatio = 0.8, sampleSize = 10, X1, X2,X3,Y1, model, order = TRUE, numClasses = 2){
  X1_s = sampleDistributions(mat1 = X1,k = k,"X",m = sampleSize,order = order)
  X2_s = sampleDistributions(mat1 = X2,k = k,"X",m = sampleSize,order = order)
  X3_s = sampleDistributions(mat1 = X3,k = k,"X",m = sampleSize,order = order)
  Y_s = sampleDistributions(mat1 = Y1,k = k*3,"Y",m = sampleSize,order = order)
  
  All = rbind(X1_s,X2_s,X3_s,Y_s)
  
  inds = sample((1:nrow(All)), size = nrow(All)*0.3, replace = FALSE)
  test = All[inds,]
  train = All[-inds,]
  
  y = train[,1]
  classLevels = unique(y)
  numClasses = length(classLevels)
  
  y_train = as.numeric(as.factor(y))-1
  x_train = as.matrix(train[,-1])
  colnames(x_train) = seq(1:ncol(x_train))
  
  y_train <- to_categorical(y_train, numClasses)
  
  
  x_test = as.matrix(test[,-1])
  y_test = as.numeric(as.factor(test[,1]))-1
  y_test <- to_categorical(y_test, numClasses)
  
  
  # model <- keras_model_sequential() 
  # model %>% 
  #   layer_dense(units = 100, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
  #   layer_dropout(rate = 0.1) %>%
  #   layer_dense(units = 10, activation = 'relu') %>% 
  #   layer_dropout(rate = 0.1) %>%
  #   layer_dense(units = 10, activation = 'relu') %>% 
  #   layer_dropout(rate = 0.1) %>%
  #   layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = 30, batch_size = 10, 
    validation_split = 0.1
  )
  
  return((model %>% evaluate(x_test, y_test))$acc)
}


sampleSizes = c(1:30)
accuracies = rep(0,length(sampleSizes))
for(i in 1:length(sampleSizes)){
  accuracies[i] = experimentWithModel(X1 = X1, X2 = X2, X3 = X3, Y1 = Y1, sampleSize = sampleSizes[i])
}

par(mfrow = c(1,1))
plot(x = sampleSizes, y = accuracies, type = "l")



# X1 = getPoints(1000, c(1,-3), c(1,1))
# X2 = getPoints(1000, c(1,1), c(1,1))
# X3 = getPoints(1000, c(1,5), c(1,1))
Y1 = getPoints(n1 = 1000, mean = c(0,0), sd = c(1,2))

a = 1.3
X1_inds = which(Y1[,2] > a)
X2_inds = which(Y1[,2] < -a )
X3_inds = c(1:nrow(Y1))[-c(X1_inds, X2_inds)]

X1 = Y1[X1_inds,]
X2 = Y1[X2_inds,]
X3 = Y1[X3_inds,]

# experimentWithModel(X1 = X1, X2 = X2, X3 = X3, Y1 = Y1, sampleSize = 40)

numClasses = 2
sampleSize = 20
model <- keras_model_sequential()
model %>%
  layer_dense(units = 100, activation = 'relu', input_shape = c(sampleSize*2)) %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10, activation = 'relu') %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10, activation = 'relu') %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = numClasses, activation = 'softmax')
experimentWithModel(X1 = X1, X2 = X2, X3 = X3, Y1 = Y1, sampleSize = sampleSize, model = model,order = TRUE)
# 0.9666666

sampleSize = 1
model <- keras_model_sequential()
model %>%
  layer_dense(units = 100, activation = 'relu', input_shape = c(sampleSize*2)) %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10, activation = 'relu') %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10, activation = 'relu') %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = numClasses, activation = 'softmax')
experimentWithModel(X1 = X1, X2 = X2, X3 = X3, Y1 = Y1, sampleSize = sampleSize, model = model)
# 0.4777778

sampleSize = 1
model <- keras_model_sequential()
model %>%
  layer_dense(units = 1, activation = 'relu', input_shape = c(2*sampleSize)) %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = numClasses, activation = 'softmax')
experimentWithModel(X1 = X1, X2 = X2, X3 = X3, Y1 = Y1, sampleSize = sampleSize, model = model)


mat_both = rbind(X1,X2,X3,Y1)

pdf("/home/willy/PredictingProteinInteractions/Results/Images/NN2dExample.pdf")
par(mfrow = c(1,4))
xli = c(-6,5)
yli = c(-6,10)
plot(Y1, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
plot(X1, col = "blue", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
plot(X3, col = "blue", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
plot(X2, col = "blue", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
# points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)
legend(x = 5,y = 10, legend = c("X","Y"),col = c("red", "blue"),pch = 19)
dev.off()


pdf("/home/willy/PredictingProteinInteractions/Results/Images/NN2dExampleSampleX1vsY1.pdf")
par(mfrow = c(1,4))
plot(Y1, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
# points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)
legend(x = 3,y = 12, legend = c("X","Y"),col = c("red", "blue"),pch = 19)

Y_demo = sampleDistributions(mat1 = Y1,k = 1,"Y",m = 20,TRUE, col = "green")

# par(mfrow = c(1,1))
plot(mat_both, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)
legend(x = 3,y = 12, legend = c("X","Y"),col = c("red", "blue"),pch = 19)
X1_demo = sampleDistributions(mat1 = X1,k = 1,"X",m = 20,TRUE, col ="green")


plot(mat_both, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)
legend(x = 3,y = 12, legend = c("X","Y"),col = c("red", "blue"),pch = 19)
X2_demo = sampleDistributions(mat1 = X2,k = 1,"X",m = 20,TRUE, col ="green")

plot(mat_both, col = "red", pch = 19, xlab = "", ylab ="", xlim = xli, ylim = yli)
points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)
legend(x = 3,y = 12, legend = c("X","Y"),col = c("red", "blue"),pch = 19)
X3_demo = sampleDistributions(mat1 = X3,k = 1,"X",m = 20,TRUE, col ="green")
dev.off()

#-----------------------------------------------------------------------
# Example 3

# install_keras(method = c("auto", "virtualenv", "conda"),
#               conda = "auto", version = "default", tensorflow = "default",
#               extra_packages = c("tensorflow-hub"))
# 


# install_keras(tensorflow = "gpu")




# 240/240 [==============================] - 0s 23us/sample - loss: 0.2200 - acc: 0.9667
# $loss
# [1] 0.2200407
# 
# $acc
# [1] 0.9666666

# pdf("/home/willy/PredictingProteinInteractions/Results/Images/NN2dExampleHistoryTraining.pdf")
plot(history)
dev.off()

# model %>% predict_classes(x_test)
#----------------------------------------------------------------------------------------------------
# 3dim
getPointsArbitraryDim <- function(n1,mean1,sd1, className){
  points = data.frame(matrix(rnorm(n1*length(mean1),mean = mean1, sd = sd1), ncol = length(mean1), byrow = TRUE))
  
  points = abs(points)
  
  points = cbind(rep(className,nrow(points)), points)
  colnames(points)[1] = "class"
  return(points)
}


sampleDistributionsArbitraryDim <- function(mat1, sampleTimes, m, sampleSize, plot = FALSE, col = NULL, pointPlotSize = 10){
  # mat1 ... points in 2d-space
  # simulate sampling m distributions  for k times
  
  dataOut = data.frame(matrix(0, ncol = sampleSize*(ncol(mat1)-1)+1, nrow = sampleTimes))
  colnames(dataOut)[1] = "class"
  
  for(j in 1:sampleTimes){
    mat1_sample = mat1[sample(nrow(mat1), size = sampleSize, replace = FALSE),]
    if(plot == TRUE && ncol(mat1) == 3+1){
      if(is.null(col)) col = "green"
      points3d(x = mat1_sample[,2], y = mat1_sample[,3], z = mat1_sample[,4], pch = 19, size = pointPlotSize, col = col, add = TRUE)
    }
    
    # print(mat1_sample[,1])
    distributions = unfoldProjArbitraryDim(mat1_sample)
    # distributions_ordered = distributions[order(distributions[,1]),]
    
    distributions_ordered = distributions[ do.call(order, distributions), ]
    
    # print(distributions_ordered)
    
    q = ncol(distributions_ordered)-1

    # input-nodes
    inputNodes = data.frame(matrix(0, ncol = q*sampleSize+1, nrow = 1))
    inputNodes[1,1] = distributions[1,1]
    colnames(inputNodes) = c("class", seq(1:(ncol(inputNodes)-1)))
    

    for( i in 1:sampleSize){
      start = (i-1)*q+1+1
      end = start+q-1

      # print(paste(start,end))
      length(distributions_ordered[i,])
      # print(inputNodes)
      inputNodes[1,start:end] = distributions_ordered[i,-1]
    }

    # print(inputNodes)
    
    dataOut[j,] = inputNodes
  }
  
  return(dataOut)
}

plotBoth <- function(all){
  classes = unique(all[,1])
  
  colMap = c("red", "blue", "yellow", "green", "pink","black")
  for(i in 1:length(classes)){
    pts = all[which(all[,1] == classes[i]),]
    
    print(pts[1,])
    
    flag = TRUE
    if(i == 1) flag == FALSE
    points3d(x = pts[,2], y = pts[,3], z = pts[,4], col = colMap[i], add = flag)
  }
}


library(rgl)



NNexampleArbitraryDim <- function(mat_both_sampled,trainTestSplit = 0.3, epochs = 300){
  
  inds = sample((1:nrow(mat_both_sampled)), size = nrow(mat_both_sampled)*trainTestSplit, replace = FALSE)
  test = mat_both_sampled[inds,]
  train = mat_both_sampled[-inds,]
  
  #----------------------------------------------------
  y = train[,1]
  
  classLevels = unique(y)
  numClasses = length(classLevels)
  
  y_train = as.numeric(as.factor(y))-1
  x_train = as.matrix(train[,-1])
  colnames(x_train) = seq(1:ncol(x_train))
  
  y_train <- to_categorical(y_train, numClasses)
  
  x_test = as.matrix(test[,-1])
  y_test = as.numeric(as.factor(test[,1]))-1
  y_test <- to_categorical(y_test, numClasses)
  
  model <- keras_model_sequential() 
  model %>% 
    layer_dense(units = 100, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
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
    epochs = epochs, batch_size = 10, 
    validation_split = 0.1
  )
  
  return(model %>% evaluate(x_test, y_test))
}

#--------------------------------------
# example1
X1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,0,0), sd1 = c(1,1,1), className = "X")
Y1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(1,0,0), sd1 = c(1,1,1), className = "Y")

sampleSize = 100
sampleTimes = 200

both_all = rbind(X1,Y1)

# plotBoth(both_all)
X1_s = sampleDistributionsArbitraryDim(mat1 = X1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
Y1_s = sampleDistributionsArbitraryDim(mat1 = Y1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)

mat_both_sampled = rbind(X1_s,Y1_s)

NNexampleArbitraryDim(mat_both_sampled = mat_both_sampled,trainTestSplit = 0.3,epochs = 300)


#--------------------------------------
# example2
X1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,0,0), sd1 = c(1,1,1), className = "X")
X2 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,1,0), sd1 = c(1,1,1), className = "X")
X3 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,0,1), sd1 = c(1,1,1), className = "X")
Y1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(1,0,0), sd1 = c(1,1,1), className = "Y")

sampleSize = 100
sampleTimes = 200

both_all = rbind(X1,Y1,X2,X3)

# plotBoth(both_all)
X1_s = sampleDistributionsArbitraryDim(mat1 = X1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
X2_s = sampleDistributionsArbitraryDim(mat1 = X2,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
X3_s = sampleDistributionsArbitraryDim(mat1 = X3,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
# X1_s = sampleDistributionsArbitraryDim(mat1 = X1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)

Y1_s = sampleDistributionsArbitraryDim(mat1 = Y1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)

mat_both_sampled = rbind(X1_s,Y1_s,X2_s,X3_s)

NNexampleArbitraryDim(mat_both_sampled = mat_both_sampled,trainTestSplit = 0.3,epochs = 30)
# $loss
# [1] 0.08430995
# 
# $acc
# [1] 0.9833333

#--------------------------------------
# example3
# this should not be possible to predict. If it separates the classes, then the method does not
# work.
X1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,0,0), sd1 = c(1,1,1), className = "X")
Y1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,0,0), sd1 = c(1,1,1), className = "Y")

sampleSize = 100
sampleTimes = 200

both_all = rbind(X1,Y1)

# plotBoth(both_all)
X1_s = sampleDistributionsArbitraryDim(mat1 = X1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
Y1_s = sampleDistributionsArbitraryDim(mat1 = Y1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)

mat_both_sampled = rbind(X1_s,Y1_s)

NNexampleArbitraryDim(mat_both_sampled = mat_both_sampled,trainTestSplit = 0.3,epochs = 30)
# $loss
# [1] 0.7091806
# 
# $acc
# [1] 0.475


#--------------------------------------
# example4
# more than two classes
X1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,0,0), sd1 = c(1,1,1), className = "X")
Y1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,0,0), sd1 = c(2,1,1), className = "Y")
Z1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(0,0,0), sd1 = c(1,2,1), className = "Z")
P1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(1,0,0), sd1 = c(1,1,1), className = "P")
Q1 = getPointsArbitraryDim(n1 = 1000, mean1 = c(1,0,0), sd1 = c(1,2,1), className = "Q")
sampleSize = 100
sampleTimes = 200

both_all = rbind(X1,Y1,Z1,P1,Q1)

plotBoth(both_all)
X1_s = sampleDistributionsArbitraryDim(mat1 = X1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
Y1_s = sampleDistributionsArbitraryDim(mat1 = Y1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
Z1_s = sampleDistributionsArbitraryDim(mat1 = Z1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
P1_s = sampleDistributionsArbitraryDim(mat1 = P1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)
Q1_s = sampleDistributionsArbitraryDim(mat1 = Q1,sampleTimes = sampleTimes,sampleSize = sampleSize, m = 10, plot = FALSE,pointPlotSize = 10)

mat_both_sampled = rbind(X1_s,Y1_s,Z1_s, P1_s, Q1_s)


NNexampleArbitraryDim(mat_both_sampled = mat_both_sampled,trainTestSplit = 0.3,epochs = 300)
# $loss
# [1] 0.4626881
# 
# $acc
# [1] 0.8866667
