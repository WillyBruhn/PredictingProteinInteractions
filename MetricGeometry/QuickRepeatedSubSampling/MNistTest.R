
s1 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R"
source(s1)

s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s2)

# library(readmnist)

# test_images <- Read.mnist("/home/willy/PredictingProteinInteractions/data/MNist/t10k-images-idx3-ubyte")
# test_labels <- Read.mnist("/home/willy/PredictingProteinInteractions/data/MNist/t10k-labels-idx1-ubyte")
# train_images <- Read.mnist("/home/willy/PredictingProteinInteractions/data/MNist/train-images-idx3-ubyte")
# train_labels <- Read.mnist("/home/willy/PredictingProteinInteractions/data/MNist/train-labels-idx1-ubyte")


# plotNumber <- function(number){
#   plot(x = 0, y = 0, xlim = c(0,28), ylim = c(0,28))
#   for(i in 1:nrow(number)){
#     for(j in 1:ncol(number)){
#       points(x = i, y = j, col = number[i,j])
#     }
#   }
# }
# 
# 
# test_labels$labels[3]
# 
# print(matrix(test_images$pic[3,], ncol = 28, byrow = FALSE))
# plotNumber(matrix(test_images$pic[3,], ncol = 28, byrow = FALSE))


# install.packages("keras")
library(keras)
mnist <- dataset_mnist()
x_train <- mnist$train$x
y_train <- mnist$train$y

# visualize the digits
par(mfcol=c(6,6))
par(mar=c(0, 0, 3, 0), xaxs='i', yaxs='i')
for (idx in 1:36) { 
  im <- x_train[idx,,]
  im <- t(apply(im, 2, rev)) 
  image(1:28, 1:28, im, col=gray((0:255)/255), 
        xaxt='n', main=paste(y_train[idx]))
}

generateF_approximations_3dModel <- function(model_points, n = 100, m = 10, q = 2, pos =TRUE){
  
  pos13_F_list = list()
  pos13_F_approx_list = list()
  
  for(i in 1:m){
    pos13_F_list[[i]] = samplePointsAndCalculateCDFofEc(all_pts = model_points, n = n,plot = FALSE)
    pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
  }
  
  return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
}

getAllModel_F_approximationsDigits <- function(model_vec, n = 100, m = 50, q = 2, maxNum = NULL){
  
  distributions_lists = list()
  
  
  if(is.null(maxNum) || maxNum > length(model_vec)) maxNum = length(model_vec)
  
  print(paste(n,m,q, length(model_vec)))
  
  for(i in 1:maxNum){
    print(i/length(maxNum))
    F_app = generateF_approximations_3dModel(model_vec[[i]]$image,n = n,q = q, m = m)
    distributions_lists[[i]] =  list("name" = model_vec[[i]]$name,"F" = F_app)
  }
  
  return(distributions_lists)
}

combine_custom <- function(list1, list2){
  return(list(list1, list2))
}

getAllModel_F_approximationsDigits_parallel <- function(model_vec, n = 100, m = 50, q = 2, maxNum = NULL){
  print("-------------------------------------------")
  print(paste("Calculating all models",sep =""))
  
  if(is.null(maxNum) || maxNum > length(model_vec)) maxNum = length(model_vec)
  
  print(maxNum)
  
  out = foreach(i = 1:length(maxNum), .combine = c) %do% {
    F_app = generateF_approximations_3dModel(model_vec[[i]]$image,n = n,q = q, m = m)
    list("name" = model_vec[[i]]$name,"F" = F_app)
  }
  
  print(length(out))
  
  return(out)
}

# ?foreach


getCoords <- function(img){
  num = length(which(img != 0))
  
  df = data.frame(matrix(0,ncol = 2, nrow = num))
  colnames(df) = c("x","y")
  
  # print(df)
  
  count = 1
  for(i in 1:nrow(img)){
    for(j in 1:ncol(img)){
      if(img[i,j] != 0){
        df[count,1] = j
        df[count,2] = nrow(img)-i
        count = count + 1
      }
    }
  }
  
  return(df)
}

n_train = length(x_train)/28/28

maxNum = 60000

imagesFile = "/home/willy/PredictingProteinInteractions/data/MNist/allImages.RData"
recreateImages = TRUE

if(!file.exists(imagesFile) || recreateImages){
  images = list()
  for(i in 1:(maxNum)){
    print(paste("coordtrans ", i/maxNum, sep =""))
    image = getCoords(x_train[i,,])
    images[[i]] = list("name" = y_train[i], "image" = image)
  }
  
  saveRDS(object = images, file = imagesFile)
}
images = readRDS(imagesFile)

numPoints = rep(0,length(images))
for(i in 1:length(images)){
  numPoints[i] = nrow(images[[i]]$image)
}

plot(images[[1001]]$image,xlim =c(0,29),ylim =c(0,29))

min(numPoints)

n = 10
m = 10
q = 8

projectionFile = paste("/home/willy/PredictingProteinInteractions/data/MNist/proj_n_", n, "_m_", m, "_q_", q,"_maxNum_",maxNum,".csv",sep ="")
maxNum = 10000
recreateProj = TRUE

if(!file.exists(projectionFile) || recreateProj){
  allDigits = getAllModel_F_approximationsDigits(model_vec = images,
                                                 n = n,
                                                 m = m,
                                                 q = q,
                                                 maxNum)
  
  df = getManhattanProjection(allDigits)
  write.csv(df,file = projectionFile, row.names = FALSE)
}

df = read.csv(projectionFile, header = TRUE)
#----------
# par(mfcol=c(1,1))
# colmap = c("steelblue1", "yellow", "pink1", "red", "green3", "blue", "aliceblue", "black", "brown", "green2")

# cols = as.vector(colmap[as.numeric(as.factor(df[,1]))])
# 
# plot(x = df[,3], y = df[,4])
# points(x = df[,3], y = df[,4], col = cols, pch = 19)
# 
# geos = getGeometricCenters(df, n = n)
# 
# for(i in 1:nrow(geos)){
#   text(x = geos[i,3], y = geos[i,4], labels = geos[i,1], col = "black", cex = 3)
# }
#----------

x = df[,3:ncol(df)]
colnames(x)

y = as.factor(df[,1])

data = cbind(x,y)
data[1:5,]

library(caret)
control <- trainControl(method="repeatedcv", number=5, repeats=10)
seed <- 7919
metric <- "Accuracy"
set.seed(seed)


mtry <- data.frame(matrix(0, ncol = 1, nrow = 3))
colnames(mtry) = c("k")
mtry[1:4,] = c(1,5,10,20)

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
knn <- train(y~., data=data, method="knn", metric=metric, tuneGrid=mtry, trControl=control)
stopCluster(cl)

saveRDS(knn,file = "/home/willy/PredictingProteinInteractions/data/MNist/knn.RData")


#--------------------------------------------------------
# translate test-set into projection

x_test = mnist$test$x
y_test = mnist$test$y

n_test = length(x_test)/28/28

maxNum_test = min(maxNum*0.4,n_test)

images_test = list()
for(i in 1:(maxNum_test)){
  print(i)
  image = getCoords(x_test[i,,])
  images_test[[i]] = list("name" = y_test[i], "image" = image)
}

y_test = y_test[c(1:maxNum_test)]


allDigits_test = getAllModel_F_approximationsDigits(images_test,
                                               n = n,
                                               m = m,
                                               q = q)

df_test = getManhattanProjection(allDigits_test)

x_test2 = df_test[,3:ncol(df_test)]
colnames(x_test2)

y_test2 = as.factor(df_test[,1])

y_pred = predict(rf_default, newdata = x_test2)
y_pred2 = matrix(y_pred, ncol = m, byrow = TRUE)

nrow(y_pred2)

y_pred_final =rep(0,length(y_test))

for(i in 1:length(y_test)){
  y_pred_final[i] = calculate_mode(y_pred2[i,])
}

calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

print(paste("accuracy:",length(which((as.numeric(y_pred_final) == as.numeric(y_test)) == TRUE))/length(y_test)))

evaluateAccuracy <- function(model, x_test, y_test){

  y_pred = predict(model, newdata = x_test)
  y_pred2 = matrix(y_pred, ncol = m, byrow = TRUE)
  
  y_pred_final =rep(0,length(y_test))
  
  for(i in 1:length(y_test)){
    y_pred_final[i] = calculate_mode(y_pred2[i,])
  }
  
  calculate_mode <- function(x) {
    uniqx <- unique(na.omit(x))
    uniqx[which.max(tabulate(match(x, uniqx)))]
  }
  
  print(paste("accuracy:",length(which((as.numeric(y_pred_final) == as.numeric(y_test)) == TRUE))/length(y_test)))
}


