#!/usr/bin/Rscript

#---------------------------------------------------------------------
# Willy Bruhn, 5.8.19
#
# Process models of animals.
#
#---------------------------------------------------------------------


wsPath = "/home/willy/PredictingProteinInteractions/setUp/SourceLoader.R"
# wsPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/setUp/SourceLoader.R"

# mode = "onlyExperiments"
# mode = "onlyGenerateModels"

wsPath = as.character(paste(funr::get_script_path(), "/../../setUp/SourceLoader.R", sep = ""))

source(wsPath)
sourceFiles(c("helperFunctions"))
sourceFiles(c("UltraQuickRepeatedSubSampling"))
sourceFiles(c("QuickRepeatedSubSampling"))
sourceFiles(c("TriangulateIsoSurface"))
sourceFiles(c("kerasFunctions"))

animalPath = getPath("animals")

library(rgl)
library(permute)


packagesLoadedFrom("animalTest.R")

downsampleEuclideanAndGetGeodesic <- function(objPath, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE){
  
  print(objPath)
  
  model_rgl = read.obj(objPath, convert.rgl = FALSE)
  model_rgl_plot = read.obj(objPath, convert.rgl = TRUE)
  
  print("plotting")
  if(plot) shade3d(model_rgl_plot)
  points = t(model_rgl$shapes[[1]]$positions)
  edges = t(model_rgl$shapes[[1]]$indices)+1
  
  # if(checkForLargerModel(path,proteinName = name,n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra) == FALSE)
  print(paste("model has ", nrow(points),", ", nrow(edges), " points", sep ="" ))
  
  ob = preProcessMesh(points = points, edges = edges, plot = FALSE)
  
  print(paste("processed model has ", nrow(ob$points), "points", sep ="" ))
  
  graph = ob$graph
  edges = ob$edges
  
  
  library(rdist)
  print("step 1: euclidean fps ...")
  
  sampled_indices = myFarthestPointSampling(points, k = n_s_euclidean)
  if(plot) plotDownsampledPoints(model_rgl,ob$points,sampled_indices, plotModel = FALSE)
  
  # furthermore subsample with the distances on the surface
  print("step 2: surface distance of sampled points ...")
  d_surface = myShortestDistances(graph,sampled_indices)
  
  print("step 3: surface fps ...")
  # ??farthest_point_sampling
  fps_surface <- farthest_point_sampling(d_surface)
  sampled_indices2 = fps_surface[1:n_s_dijkstra]
  
  # rgl.open()
  if(plot) plotDownsampledPoints(model_rgl,points,sampled_indices[sampled_indices2],FALSE, col = "blue", size = 51)
  
  v2 = myVoronoi(points[sampled_indices[sampled_indices2],], points)
  v_n = v2/sum(v2)
  
  l = list("centers" = points[sampled_indices[sampled_indices2],], "mu" = v_n, "indices_order" = sampled_indices[sampled_indices2], "d_surface" = d_surface[sampled_indices2,sampled_indices2], "sampled_indices_geo" = sampled_indices)
  return(l)
}

# getGeoDistanceName <- function(path,ind,n_s_euclidean,n_s_dijkstra){
#   return(paste(path,"/geoDist_ind_",ind, "_nE_", n_s_euclidean, "_nD_", n_s_dijkstra, ".csv", sep =""))
# }
# 
# getGeoDistanceQuantileName <- function(path,ind,n_s_euclidean,n_s_dijkstra, n, m, q){
#   return(paste(path,"/Quantiles_ind_",ind, "_nE_", n_s_euclidean, "_nD_", n_s_dijkstra,
#                "_n_", n, "_m_", m, "_q_",q,".csv", sep =""))
# }


getAllModelsFromAnimal <- function(path = "/home/willy/PredictingProteinInteractions/data/animals/models/",
                                   name = "camel",
                                   n_s_euclidean = 100,
                                   n_s_dijkstra = 50,
                                   n = 10,
                                   m = 3,
                                   q = 1,
                                   plot = TRUE){
  
  path_final = paste(path, "/", name, "-poses/", sep ="")
  all.files = list.files(path_final, recursive = FALSE,pattern = ".obj", full.names = TRUE)
  
  all.names = list.files(path_final, recursive = FALSE,pattern = ".obj", full.names = FALSE)
  all.names_final = unlist(strsplit(all.names, split = ".obj"))
  
  distances_path = paste(path_final,"/geoDistances/", sep ="")
  if(!dir.exists(distances_path)) dir.create(distances_path)
  
  quantilesOut = data.frame(matrix(0,ncol = q+3, nrow = 0))
  colnames(quantilesOut) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")))
  
  for(i in 1:length(all.names_final)){
    quantilesName = getGeoDistanceQuantileName(distances_path,i,n_s_euclidean,n_s_dijkstra,n = n,m = m,q = q)
    if(!file.exists(quantilesName)){
      geoDistanceName = getGeoDistanceName(distances_path,i,n_s_euclidean,n_s_dijkstra)
      
      if(!file.exists(geoDistanceName)){
        mod = downsampleEuclideanAndGetGeodesic(objPath = all.files[i],
                                                n_s_euclidean = n_s_euclidean,
                                                n_s_dijkstra = n_s_dijkstra,
                                                plot = plot)
        
        write.csv(x = mod$d_surface, file = geoDistanceName,row.names = FALSE)
      }
      
      geoDist = read.csv(geoDistanceName)
      
      Fapp = generateF_approximations_3dModelWithMetric(d_surface = geoDist,n = n,m = m,q = q)
      
      quantiles = data.frame(matrix(0,ncol = q+3, nrow = m))
      colnames(quantiles) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")))
      
      quantiles[,1] = rep(all.names_final[i], m)
      for(i in 1:length(Fapp$F_app_list)){
        quantiles[i,2:ncol(quantiles)] = Fapp$F_app_list[[i]]  
      }
      
      write.csv(x = quantiles, quantilesName)
    }
    
    quantiles = read.csv(quantilesName, row.names = 1)
    
    quantilesOut = rbind(quantilesOut,quantiles)
  }
  
  return(quantilesOut)
}

createImageWithArrows <- function(quantiles, headToo = FALSE){
  df = quantiles
  #
  # df = readProjectionFromFile(n = n,m = m,q = q,path = "/home/willy/PredictingProteinInteractions/data/animals/",fName = "proj")
  classes = unique(getClassNamesFromSubClasses(df[,1]))
  colors = as.numeric(as.factor(classes))
  colorMap = c("red", "blue", "green", "yellow", "black", "pink", "orange")
  for(i in 1:length(colorMap)){
    print(paste(classes[i], colorMap[i]))
    
    if(!(classes[i] == "head") || headToo == TRUE){
      inds = which(getClassNamesFromSubClasses(df[,1]) == classes[i])
      # print(inds)
      points3d(x = df[inds,2], y = df[inds,3], z = df[inds,4],colorMap[i], add = TRUE,size = 10)
      
      geox = sum(df[inds,2])/length(inds)
      geoy = sum(df[inds,3])/length(inds)
      geoz = sum(df[inds,4])/length(inds)
      
      textCoord = c(geox-0.1,geoy-0.1,geoz)
      if(classes[i] == "horse") {
        textCoord = c(geox+0.1,geoy+0.1,geoz)
        text3d(x = textCoord[1]+0.05, y = textCoord[2]+0.05, z = textCoord[3],texts = classes[i],cex = 2)
        
        arrow3d(p0 = textCoord, c(geox+0.05,geoy+0.05,geoz), col = "black")
      }else{
        text3d(x = textCoord[1]-0.05, y = textCoord[2]-0.05, z = textCoord[3],texts = classes[i],cex = 2)
        
        arrow3d(p0 = textCoord, c(geox-0.05,geoy-0.05,geoz), col = "black")
      }
      
      
      
    }
  }
  
  # rgl.snapshot("/home/willy/PredictingProteinInteractions/Results/Images/animals3dProjectionExample.png")
}


# sampleDistancesAndCalculateCDwithD
# generateF_approximations_3dModel <- function(model_points, n = 100, m = 10, q = 2, pos =TRUE){
#   
#   # print(nrow(model_points))
#   
#   pos13_F_list = list()
#   pos13_F_approx_list = list()
#   
#   for(i in 1:m){
#     pos13_F_list[[i]] = samplePointsAndCalculateCDFofEc(all_pts = model_points, n = n,plot = FALSE)
#     pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
#   }
#   
#   return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
# }
# 
# generateF_approximations_3dModelWithMetric <- function(d, n = 100, m = 10, q = 2, pos =TRUE){
#   
#   # print(nrow(model_points))
#   
#   pos13_F_list = list()
#   pos13_F_approx_list = list()
#   
#   for(i in 1:m){
#     pos13_F_list[[i]] = sampleDistancesAndCalculateCDFofEcWith(d, n = n,plot = FALSE)
#     pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
#   }
#   
#   return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
# }

getAllModel_F_approximations <- function(model_vec, num, n = 100, m = 50, q = 2){
  
  distributions_lists = list()
  
  for(i in 1:num){
    print(i/num)
    
    F_app = generateF_approximations_3dModel(model_vec[3+(i-1)*5][[1]],n = n,q = q, m = m)
    distributions_lists[[i]] =  list("name" = model_vec[1+(i-1)*5],"F" = F_app)
  }
  
  return(distributions_lists)
}


animalModel <- function(TrainTest, batch_size = 10, epochs = 30, validation_split = 0.1){
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test

  numClasses = TrainTest$numClasses

  print("building model ...")
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 300, activation = 'relu', input_shape = c(ncol(x_train))) %>%
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
    validation_split = validation_split
  )

  model %>% evaluate(x_test, y_test)

  return(list("model" = model, "history" = history))
}


animalModelSimple <- function(TrainTest, batch_size = 10, epochs = 30, validation_split = 0.1){
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  print("building model ...")
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 9, activation = 'relu', input_shape = c(ncol(x_train))) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 3, activation = 'relu') %>%
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
    validation_split = validation_split
  )
  
  model %>% evaluate(x_test, y_test)
  
  return(list("model" = model, "history" = history))
}

#----------------------------------------------------------------------------------------------------------------------------
# stratified k-fold cv
#----------------------------------------------------------------------------------------------------------------------------

n_s_euclidean = 1000
n_s_dijkstra = 100
n = 96
m = 100
q = 1
plot = FALSE

# n_s_euclidean = 200
# n_s_dijkstra = 50
# n = 48
# m = 100
# q = 1

print("fetching quantiles ...")

camel = getAllModelsFromAnimal(name = "camel",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = m,q = q, plot = plot)
cat = getAllModelsFromAnimal(name = "cat",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = m,q = q, plot = plot)
lion = getAllModelsFromAnimal(name = "lion",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = m,q = q, plot = plot)
elephant = getAllModelsFromAnimal(name = "elephant",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = m,q = q, plot = plot)
flam = getAllModelsFromAnimal(name = "flam",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = m,q = q, plot = plot)
horse = getAllModelsFromAnimal(name = "horse",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = m,q = q, plot = plot)
head = getAllModelsFromAnimal(name = "head",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = m,q = q, plot = plot)

all_models = rbind(camel,cat,lion,elephant,flam,horse,head)

# write.csv(all_models,file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_200_nS_50_n_48_m_100_q_1.csv",
#           row.names = FALSE)

# write.csv(all_models,file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_1000_nS_100_n_99_m_100_q_1.csv",
          # row.names = FALSE)

write.csv(all_models,file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_1000_nS_100_n_96_m_100_q_1.csv",
          row.names = FALSE)

# quantiles = read.csv(file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_200_nS_50_n_48_m_100_q_1.csv")
# quantiles = read.csv(file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_1000_nS_100_n_99_m_100_q_1.csv")

quantiles = read.csv(file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_1000_nS_100_n_96_m_100_q_1.csv")

# createImageWithArrows(quantiles, FALSE)
# points3d(quantiles[1:100, 2:4], col = "green", size = 20)
# points3d(quantiles[101:200, 2:4], col = "blue", size = 20)

sampleSize = 10
sampleTimes = 100
sampleTimes_test = 100
epochs = 30
batch_size = 10
k = 10

ClassLevels = sort(unique(getClassNamesFromSubClasses(quantiles[,1])))
numClasses = length(ClassLevels)

library(caret)

# create stratified folds
# these are the actual model-names
uniqueModels = unique(quantiles[,1])

# these are the classes to which each model belongs
classNames = getClassNamesFromSubClasses(uniqueModels)

folds = createFolds(classNames, k = k,list = TRUE)

dirName = "/home/willy/PredictingProteinInteractions/data/animals/Evaluation/"
if(!dir.exists(dirName)) dir.create(dirName)


for(foldInd in 1:length(folds)){

  
    y_test_name_inds = folds[[foldInd]]
    test_inds = which(quantiles[,1] %in% uniqueModels[y_test_name_inds])
    train_inds = c(1:nrow(quantiles))[-test_inds]
    
    # specify euklid here, because otherwise only half of the quantiles is used (not nicely programmed)
    AllQuantiles_train = getSamplesSurf2(quantiles = quantiles[train_inds,],sampleSize = sampleSize,sampleTimes = sampleTimes,numPermutations = 1,reDo = TRUE,numClasses = numClasses,m = m,euklid = TRUE,splitPattern = "-", sort = TRUE)
    
    AllQuantiles_test = getSamplesSurf2(quantiles = quantiles[test_inds,],sampleSize = sampleSize,sampleTimes = sampleTimes_test,numPermutations = 1,reDo = TRUE,numClasses = numClasses,m = m,euklid = TRUE,splitPattern = "-", sort = TRUE)
    
    shuf = shuffle(1:nrow(AllQuantiles_train$X))
    shuf2 = shuffle(1:nrow(AllQuantiles_test$X))
    TrainTest = list("x_train" = AllQuantiles_train$X[shuf,], "y_train" = AllQuantiles_train$y[shuf,], "x_test" = AllQuantiles_test$X[shuf2,], "y_test" = AllQuantiles_test$y[shuf2,], numClasses = numClasses)
    
    
    doAgain = TRUE
    doAgainCount = 0
    while(doAgain && doAgainCount < 2){
      
      modelList = animalModel(TrainTest = TrainTest,epochs = epochs, batch_size = batch_size)
    
      predictions <- modelList$model %>% predict_classes(TrainTest$x_test)
      
      pred = predictions+1
      print(pred)
      gt = reverseToCategorical(TrainTest$y_test,ClassLevels)
      y_test_pred = rep("0",length(pred))
      su = 0
      for(i in 1:length(gt)){
        if(gt[i] == ClassLevels[pred[i]]) su = su + 1
        
        y_test_pred[i] = ClassLevels[pred[i]]
      }
      su/length(gt)
      y_test_pred
      
      as.numeric(as.factor(y_test_pred))
      
      confMat = table(factor(as.character(ClassLevels[as.numeric(as.factor(y_test_pred))]),
                             levels=ClassLevels),
                      factor(as.character(ClassLevels[as.numeric(as.factor(gt))]),
                             levels=ClassLevels))
      
      sum(confMat)
      confMatNormalized = confMat/colSums(confMat)[col(confMat)]
      
      print(confMat)
      print(confMatNormalized)
      
      accuracy = sum(diag(confMat)) / sum(confMat)
      
      write.table(confMat, file = paste(dirName, "/confMat_fold_", foldInd, ".txt", sep = ""))
      write.table(accuracy, file = paste(dirName, "/accuracy_fold_", foldInd, ".txt", sep = ""), row.names = FALSE)
      write.table(uniqueModels[folds[[foldInd]]], file = paste(dirName, "/names_fold_", foldInd, ".txt", sep = ""), row.names = FALSE)
      
      write.table(modelList$history, file = paste(dirName, "/history_", foldInd, ".txt", sep = ""), row.names = FALSE)
      
      if(accuracy > 0.75) {
        doAgain = FALSE
      }
      doAgainCount = doAgainCount +1
    }
}



# now read in all confusion-matrices and average over the performance
conf_all = read.table(paste(dirName, "/confMat_fold_", 1, ".txt", sep = ""))

for(i in 2:length(folds)){
  conf_i = read.table(paste(dirName, "/confMat_fold_", i, ".txt", sep = ""))
  
  conf_all = conf_all + conf_i
}

conf_all = as.matrix(conf_all)

accuracy = sum(diag(conf_all)) / sum(conf_all)
accuracy


confMatNormalized = conf_all/colSums(conf_all)[col(conf_all)]

write.table(x = signif(accuracy,2),file = paste("/home/willy/PredictingProteinInteractions/data/animals/Evaluation/Accuracy.tex", sep = ""),
            quote = FALSE, col.names = FALSE, row.names = FALSE)

print(xtable(x = conf_all,caption = paste("Confusion-matrix from the animal-test-set.",sep =""),label = "animalConfusionNN", type = "latex"),
      file = paste("/home/willy/PredictingProteinInteractions/data/animals/Evaluation/Confusion.tex", sep = ""))
#----------------------------------------------------------------------------------------------------------------------------


# getMemoliModelsInPath <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/", n_s_euclidean = 4000, n_s_dijkstra = 50){
#   obj_files = list.files(path = path, pattern = ".obj",recursive = TRUE, full.names = TRUE)
#   
#   print(obj_files)
#   
#   
#   model_vec = c()
#   
#   
#   for(i in 1:length(obj_files)){
#     obj_name = strsplit(obj_files[i],split = "/")[[1]][10]
#     f_path = strsplit(obj_files[i], split = obj_name)[[1]][1]
#     
#     obj_name = strsplit(obj_name,split = ".obj")[[1]][1]
#     
#     print(obj_name)
#     
#     while (rgl.cur() > 0) { rgl.close() }
#     model = getMemoliModel3(f_path,obj_name, n_s_euclidean, n_s_dijkstra)
#     
#     
#     
#     model_vec = c(model_vec, obj_name, model)
#   }
#   
#   return(model_vec)
# }

#-----------------------------------------------------------------------------------

# n = 100
# m = 1000
# q = 1
# 
# model_vec = getMemoliModelsInPath(n_s_euclidean = 4000, n_s_dijkstra = 50)
# 
# model = getMemoliModel()
# 
# 
# # 5 values consecutive values form one model
# length(model_vec)
# 
# points3d(model_vec[[4]])
# model_vec[[5]]
# 
# 
# models_all = getAllModel_F_approximations(model_vec, 72, n = n,m = m,q = q)
# 
# 
# # quants = readQuantilesFromFile(n = n,m=m,q=q)
# writeQuantilesToFileAnimal(model_vec = models_all,n = n,m=m,q=q)
# 
# 
# df = getManhattanProjection(models_all)
# writeProjectionToFile(proj = df,n = n,m = m,q = q,path = "/home/willy/PredictingProteinInteractions/data/animals/",fName = "proj")
# 
# 
# 
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
#     geox = sum(df[inds,2])/length(inds)+0.01
#     geoy = sum(df[inds,3])/length(inds)
#     geoz = sum(df[inds,4])/length(inds)
#     
#     text3d(x = geox, y = geoy, z = geoz,texts = classes[i],cex = 2)
#   }
# }
# 
# x_sd = sd(df[,2])
# y_sd = sd(df[,3])
# z_sd = sd(df[,4])
# 
# #
# rgl.snapshot("/home/willy/PredictingProteinInteractions/Results/Images/animals3dProjectionExample.png")
# 
# 
# 
# 
# # grab the test names
# testNames = sample(unique(df[,1]), size = 72*0.3,replace = FALSE)
# 
# testIndices = which(df[,1] %in% testNames)
# df_test = df[testIndices,]
# df_train = df[-testIndices,]
# 
# colnames(df_test)
# colnames(df_train)
# 
# classNames = rep("",nrow(df))
# for(i in 1:length(classNames)){
#   classNames[i] = strsplit(df$name[i],split = "-")[[1]][1]
# }
# 
# numClasses = length(unique(classNames))
# unique(classNames)
# # as.numeric(as.factor(classNames))
# 
# start=50
# # colmap = colors()[seq(start,length(colors()),(length(colors())-start)/numClasses)]
# colors()
# 
# colmap = c("steelblue1", "yellow", "pink1", "red", "green3", "blue", "aliceblue")
# 
# cols = as.vector(colmap[as.numeric(as.factor(classNames))])
# 
# plot(x = df[,3], y = df[,4], xlim = c(0,0.8), ylim = c(0,0.15))
# points(x = df[,3], y = df[,4], col = cols, pch = 19)
# 
# df_relabeled = cbind(classNames,df[,2:ncol(df)])
# colnames(df_relabeled) = c("name", colnames(df_relabeled)[2:ncol(df_relabeled)])
# 
# geos = getGeometricCenters(df_relabeled, n = n)
# geos[,1] = unique(classNames)
# 
# for(i in 1:nrow(geos)){
#   text(x = geos[i,3], y = geos[i,4], labels = geos[i,1], col = "black", cex = 3)
# }
# 
# par(xpd=FALSE)
# 
# legend("topright", inset=c(-0.2,0), legend = c(unique(classNames)), col = colmap, lty = rep(19,numClasses))
# 
# 
# library(randomForest)
# library(mlbench)
# library(caret)
# 
# 
# x = df_train[,3:ncol(df_train)]
# colnames(x)
# 
# y = df_train[,1]
# 
# data_train = cbind(x,y)
# data_train[1:5,]
# 
# control <- trainControl(method="repeatedcv", number=5, repeats=10)
# seed <- 71
# metric <- "Accuracy"
# set.seed(seed)
# 
# mtry <- data.frame(matrix(0, ncol = 1, nrow = 0))
# colnames(mtry) = c("k")
# mtry[1:3,] = c(1,10,30)
# 
# library(doParallel)
# cl <- makePSOCKcluster(5)
# registerDoParallel(cl)
# knn <- train(y~., data=data_train, method="knn", metric=metric, tuneGrid=mtry, trControl=control)
# stopCluster(cl)
# 
# knn
# 
# 
# # Testing on testset
# calculate_mode <- function(x) {
#   uniqx <- unique(na.omit(x))
#   uniqx[which.max(tabulate(match(x, uniqx)))]
# }
# 
# evaluateAccuracy <- function(model, x_test, y_test, m){
#   
#   y_pred = predict(model, newdata = x_test)
#   y_pred2 = matrix(y_pred, ncol = m, byrow = TRUE)
#   
#   y_pred_final =rep(0,length(y_test)/m)
#   
#   # print(y_pred2)
#   
#   y_pred2_class = matrix(rep("",length(y_pred)), ncol = m, byrow = TRUE)
#   print(y_pred2_class)
#   for(i in 1:nrow(y_pred2)){
#     for(j in 1:ncol(y_pred2)){
#       y_pred2_class[i,j] = strsplit(as.character(y_pred2[i,j]),split = "-")[[1]][1]
#     }
#   }
#   
#   for(i in 1:length(y_pred_final)){
#     y_pred_final[i] = calculate_mode(y_pred2_class[i,])
#   }
# 
#   
#   y_test_class = rep("",length(y_test)/m)
#   for(i in 1:length(y_test_class)){
#     y_test_class[i] = strsplit(as.character(y_test[(i-1)*m+1]),split = "-")[[1]][1]
#   }
#   
#   
#   
#   print(y_pred_final)
#   print(y_test_class)
#   
#   dfOut = cbind(y_pred_final,y_test_class)
#   
#   print(paste("accuracy:",length(which((y_pred_final == y_test_class) == TRUE))/length(y_test_class)))
#   
#   
#   list("all" = dfOut, "wrongPred" = dfOut[which(dfOut[,1] != dfOut[,2]),])
# }
# 
# x_test = df_test[,3:ncol(df_test)]
# colnames(x_test2)
# y_test = as.factor(df_test[,1])
# 
# evaluateAccuracy(knn, x_test, y_test, m)
# 
# #-----------------------------------------------------------------------------------------
# # Accuracy   Kappa    
# # 0.9849722  0.9824543
# #
# # using the whole data-set with cv
# 
# 
# x = df_relabeled[,3:ncol(df_relabeled)]
# colnames(x)
# 
# y = df_relabeled[,1]
# 
# data = cbind(x,y)
# data[1:5,]
# 
# control <- trainControl(method="repeatedcv", number=5, repeats=10)
# seed <- 716
# metric <- "Accuracy"
# set.seed(seed)
# 
# mtry <- data.frame(matrix(0, ncol = 1, nrow = 0))
# colnames(mtry) = c("k")
# mtry[1:2,] = c(1,10)
# 
# library(doParallel)
# cl <- makePSOCKcluster(5)
# registerDoParallel(cl)
# knn_relabeled <- train(y~., data=data, method="knn", metric=metric, tuneGrid=mtry, trControl=control)
# stopCluster(cl)
# 
# knn_relabeled
# # k   Accuracy   Kappa    
# # 1  0.9839306  0.9812382
# # 10  0.9854833  0.9830515
















#----------------------------------------------------------------------------------------------------------------------------
# # old Test
# n_s_euclidean = 200
# n_s_dijkstra = 50
# 
# n = 48
# 
# camel = getAllModelsFromAnimal(name = "camel",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = 100,q = 1)
# cat = getAllModelsFromAnimal(name = "cat",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = 100,q = 1)
# lion = getAllModelsFromAnimal(name = "lion",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = 100,q = 1)
# elephant = getAllModelsFromAnimal(name = "elephant",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = 100,q = 1)
# flam = getAllModelsFromAnimal(name = "flam",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = 100,q = 1)
# horse = getAllModelsFromAnimal(name = "horse",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = 100,q = 1)
# head = getAllModelsFromAnimal(name = "head",n_s_euclidean =n_s_euclidean,n_s_dijkstra = n_s_dijkstra, n = n,m = 100,q = 1)
# 
# all_models = rbind(camel,cat,lion,elephant,flam,horse,head)
# 
# write.csv(all_models,file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_200_nS_50_n_48_m_100_q_1.csv", row.names = FALSE)
# 
# 
# n = 40
# m = 100
# q = 1
# 
# sampleSize = 20
# sampleTimes = 200
# 
# testSplitRatio = 0.3
# 
# quantiles = read.csv(file ="/home/willy/PredictingProteinInteractions/data/animals/models/all_models_nE_200_nS_50_n_48_m_100_q_1.csv")
# 
# subClassNames = unique(quantiles[,1])
# numObjects = length(subClassNames)
# 
# numClasses = length(unique(getClassNamesFromSubClasses(quantiles[,1],splitPattern = "-")))
# 
# subClassNames_test_ind = sample(c(1:numObjects), size = numObjects*testSplitRatio, replace = FALSE)
# subClassNames_test =subClassNames[subClassNames_test_ind]
# subClassNames_train =subClassNames[-subClassNames_test_ind]
# 
# 
# NNInput_train = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_train),],
#                                         m,sampleSize = sampleSize,sampleTimes = sampleTimes)
# 
# 
# NNInput_train = NNInput_train[shuffle(1:nrow(NNInput_train)), ]
# 
# NNInput_test = getNNInputFromQuantiles(quantiles[which(quantiles[,1] %in% subClassNames_test),],
#                                        m,sampleSize = sampleSize,sampleTimes = sampleTimes)
# 
# y_train = getClassNamesFromSubClasses(NNInput_train[,1],splitPattern = "-")
# 
# classLevels = unique(y_train)
# y_train = as.numeric(as.factor(y_train))-1
# x_train = as.matrix(NNInput_train[,-1])
# y_train <- to_categorical(y_train, numClasses)
# 
# 
# # unique(NNInput_test[,1])
# 
# y_test = getClassNamesFromSubClasses(NNInput_test[,1],splitPattern = "-")
# y_test = as.numeric(as.factor(y_test))-1
# x_test = as.matrix(NNInput_test[,-1])
# y_test <- to_categorical(y_test, numClasses)
# 
# 
# # y_test <- to_categorical(y_test, 10)
# 
# #---------------------------------------------------------
# model <- keras_model_sequential()
# model %>% 
#   layer_dense(units = 300, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
#   layer_dropout(rate = 0.1) %>%
#   layer_dense(units = 10, activation = 'relu') %>% 
#   layer_dropout(rate = 0.1) %>%
#   layer_dense(units = 10, activation = 'relu') %>% 
#   layer_dropout(rate = 0.1) %>%
#   layer_dense(units = numClasses, activation = 'softmax')
# 
# model %>% compile(
#   loss = 'categorical_crossentropy',
#   optimizer = optimizer_rmsprop(),
#   metrics = c('accuracy')
# )
# 
# history <- model %>% fit(
#   x_train, y_train, 
#   epochs = 30, batch_size = 10, 
#   validation_split = 0.2
# )
# 
# model %>% evaluate(x_test, y_test)
# 

# 
# # 4200/4200 [==============================] - 0s 27us/sample - loss: 0.0453 - acc: 0.9936
# # $loss
# # [1] 0.04529072
# # 
# # $acc
# # [1] 0.9935714