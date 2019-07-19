library(rgl)

s1 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R"
source(s1)

s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s2)


s3 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
source(s3)


generateF_approximations_3dModel <- function(model_points, n = 100, m = 10, q = 2, pos =TRUE){
  
  print(nrow(model_points))
  
  pos13_F_list = list()
  pos13_F_approx_list = list()
  
  for(i in 1:m){
    pos13_F_list[[i]] = samplePointsAndCalculateCDFofEc(all_pts = model_points, n = n,plot = FALSE)
    pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
  }
  
  return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
}

getAllModel_F_approximations <- function(model_vec, num, n = 100, m = 50, q = 2){
  
  distributions_lists = list()
  
  for(i in 1:num){
    F_app = generateF_approximations_3dModel(model_vec[3+(i-1)*5][[1]],n = n,q = q, m = m)
    distributions_lists[[i]] =  list("name" = model_vec[1+(i-1)*5],"F" = F_app)
  }
  
  return(distributions_lists)
}

n = 100
m = 5
q = 8

model_vec = getMemoliModelsInPath(n_s_euclidean = 4000, n_s_dijkstra = 50)

# 5 values consecutive values form one model
length(model_vec)

models_all = getAllModel_F_approximations(model_vec, 72, n = n,m = m,q = q)

df = getManhattanProjection(models_all)

writeProjectionToFile(proj = df,n = n,m = m,q = q,path = "/home/willy/PredictingProteinInteractions/data/animals/",fName = "proj")



# grap the test names
testNames = sample(unique(df[,1]), size = 72*0.3,replace = FALSE)

testIndices = which(df[,1] %in% testNames)
df_test = df[testIndices,]
df_train = df[-testIndices,]

colnames(df_test)
colnames(df_train)

classNames = rep("",nrow(df))
for(i in 1:length(classNames)){
  classNames[i] = strsplit(df$name[i],split = "-")[[1]][1]
}

numClasses = length(unique(classNames))
unique(classNames)
# as.numeric(as.factor(classNames))

start=50
# colmap = colors()[seq(start,length(colors()),(length(colors())-start)/numClasses)]
colors()

colmap = c("steelblue1", "yellow", "pink1", "red", "green3", "blue", "aliceblue")

cols = as.vector(colmap[as.numeric(as.factor(classNames))])

plot(x = df[,3], y = df[,4], xlim = c(0,0.8), ylim = c(0,0.15))
points(x = df[,3], y = df[,4], col = cols, pch = 19)

df_relabeled = cbind(classNames,df[,2:ncol(df)])
colnames(df_relabeled) = c("name", colnames(df_relabeled)[2:ncol(df_relabeled)])

geos = getGeometricCenters(df_relabeled, n = n)
geos[,1] = unique(classNames)

for(i in 1:nrow(geos)){
  text(x = geos[i,3], y = geos[i,4], labels = geos[i,1], col = "black", cex = 3)
}

par(xpd=FALSE)

legend("topright", inset=c(-0.2,0), legend = c(unique(classNames)), col = colmap, lty = rep(19,numClasses))


library(randomForest)
library(mlbench)
library(caret)


x = df_train[,3:ncol(df_train)]
colnames(x)

y = df_train[,1]

data_train = cbind(x,y)
data_train[1:5,]

control <- trainControl(method="repeatedcv", number=5, repeats=10)
seed <- 71
metric <- "Accuracy"
set.seed(seed)

mtry <- data.frame(matrix(0, ncol = 1, nrow = 0))
colnames(mtry) = c("k")
mtry[1:3,] = c(1,10,30)

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
knn <- train(y~., data=data_train, method="knn", metric=metric, tuneGrid=mtry, trControl=control)
stopCluster(cl)

knn


# Testing on testset
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

evaluateAccuracy <- function(model, x_test, y_test, m){
  
  y_pred = predict(model, newdata = x_test)
  y_pred2 = matrix(y_pred, ncol = m, byrow = TRUE)
  
  y_pred_final =rep(0,length(y_test)/m)
  
  # print(y_pred2)
  
  y_pred2_class = matrix(rep("",length(y_pred)), ncol = m, byrow = TRUE)
  print(y_pred2_class)
  for(i in 1:nrow(y_pred2)){
    for(j in 1:ncol(y_pred2)){
      y_pred2_class[i,j] = strsplit(as.character(y_pred2[i,j]),split = "-")[[1]][1]
    }
  }
  
  for(i in 1:length(y_pred_final)){
    y_pred_final[i] = calculate_mode(y_pred2_class[i,])
  }

  
  y_test_class = rep("",length(y_test)/m)
  for(i in 1:length(y_test_class)){
    y_test_class[i] = strsplit(as.character(y_test[(i-1)*m+1]),split = "-")[[1]][1]
  }
  
  
  
  print(y_pred_final)
  print(y_test_class)
  
  dfOut = cbind(y_pred_final,y_test_class)
  
  print(paste("accuracy:",length(which((y_pred_final == y_test_class) == TRUE))/length(y_test_class)))
  
  
  list("all" = dfOut, "wrongPred" = dfOut[which(dfOut[,1] != dfOut[,2]),])
}

x_test = df_test[,3:ncol(df_test)]
colnames(x_test2)
y_test = as.factor(df_test[,1])

evaluateAccuracy(knn, x_test, y_test, m)

#-----------------------------------------------------------------------------------------
# Accuracy   Kappa    
# 0.9849722  0.9824543
#
# using the whole data-set with cv


x = df_relabeled[,3:ncol(df_relabeled)]
colnames(x)

y = df_relabeled[,1]

data = cbind(x,y)
data[1:5,]

control <- trainControl(method="repeatedcv", number=5, repeats=10)
seed <- 716
metric <- "Accuracy"
set.seed(seed)

mtry <- data.frame(matrix(0, ncol = 1, nrow = 0))
colnames(mtry) = c("k")
mtry[1:2,] = c(1,10)

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
knn_relabeled <- train(y~., data=data, method="knn", metric=metric, tuneGrid=mtry, trControl=control)
stopCluster(cl)

knn_relabeled
# k   Accuracy   Kappa    
# 1  0.9839306  0.9812382
# 10  0.9854833  0.9830515
#------------------------------------------------------------------------------------------
# nonrigid world

getAllModel_F_approximations_non_rigid_world <- function(model_vec, n = 100, m = 50, q = 2){
  
  distributions_lists = list()
  
  for(i in 1:length(model_vec)){
    F_app = generateF_approximations_3dModel(model_vec[[i]]$vert,n = n,q = q, m = m)
    distributions_lists[[i]] =  list("name" = model_vec[[i]]$name,"F" = F_app)
  }
  
  return(distributions_lists)
}

# read in files
vert_files = list.files("/home/willy/PredictingProteinInteractions/data/nonRigidWorld/nonrigid3d/",
                        pattern = ".vert",
                        full.names = TRUE)
vert_files_short = list.files("/home/willy/PredictingProteinInteractions/data/nonRigidWorld/nonrigid3d/",
                        pattern = ".vert")

vert_list = list()
for(i in 1:length(vert_files)){
  vert_list[[i]] = list("name" = vert_files_short[i], "vert" = read.table(vert_files[i]))
}


allVertModels = getAllModel_F_approximations_non_rigid_world(vert_list,n = 500,m = 500,q = 10)

dfVert = getManhattanProjection(allVertModels)
classNames = rep("",nrow(dfVert))
for(i in 1:length(classNames)){
  classNames[i] = strsplit(dfVert[i,1],split = ".vert")[[1]][1]
}

# classNames

x = dfVert[,3:ncol(dfVert)]
colnames(x)

y = classNames

data = cbind(x,y)
data[1:5,]

control <- trainControl(method="repeatedcv", number=5, repeats=10)
seed <- 79
metric <- "Accuracy"
set.seed(seed)


mtry <- data.frame(matrix(0, ncol = 1, nrow = 0))
colnames(mtry) = c("k")
mtry[1,] = c(10)

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
rf_default <- train(y~., data=data, method="knn", metric=metric, tuneGrid=mtry, trControl=control)
stopCluster(cl)

rf_default

# Accuracy   Kappa    
# 0.8474649  0.8464272



#------------------------------------------------------------------------------------------
# readMnist


