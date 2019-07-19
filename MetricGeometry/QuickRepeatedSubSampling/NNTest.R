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

# install.packages("rfUtilities")
library(rfUtilities)

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

sampleDistributions <- function(mat1, k, classLabel, m){
  # mat1 ... points in 2d-space
  
  # simulate sampling m distributions  for k times
  
  dataOut = data.frame(matrix(0, ncol = m+1, nrow = k))
  colnames(dataOut)[1] = "class"
  
  dataOut[,1] = classLabel
  
  for(j in 1:k){
    mat1_sample = mat1[sample(nrow(mat1), size = m, replace = FALSE),]
    # points(x = mat1_sample[,1], y = mat1_sample[,2], col ="green", pch = 19, cex = 1)
    
    distributions = unfold2dProj(mat1_sample)
    distributions_ordered = distributions[order(distributions[,1]),]
    
    # input-nodes
    inputNodes = matrix(0, ncol = nrow(distributions_ordered)*2)
    inputNodes
    
    for( i in 1:nrow(distributions_ordered)){
      inputNodes[(i-1)*2+1] = distributions_ordered[i,1]
      inputNodes[(i-1)*2+1+1] = distributions_ordered[i,2]
    }
    
    dataOut[j,2:ncol(dataOut)] = inputNodes
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

fitNNModel(All)


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
plot(mat_both, col = "red", pch = 19)
points(x = Y1[,1], y = Y1[,2], col ="blue", pch = 19)

# generate training set and test set
k = 200
splitRatio = 0.8

X1_s = sampleDistributions(mat1 = X1,k = k,"X",m = 10)
X2_s = sampleDistributions(mat1 = X2,k = k,"X",m = 10)
X3_s = sampleDistributions(mat1 = X3,k = k,"X",m = 10)
Y_s = sampleDistributions(mat1 = Y1,k = k,"Y",m = 10)

All = rbind(X1_s,X2_s,X3_s,Y_s)

fitNNModel(All)



