#!/usr/bin/Rscript

#---------------------------------------------------------------------
# Willy Bruhn, 19.7.19
#
# Read in the model-net10 data-set and create the usual projections.
#
#---------------------------------------------------------------------

# s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
# source(s1)
# 
# s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
# source(s2)
# 
# s3 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R"
# source(s3)




wsPath = "/home/willy/PredictingProteinInteractions/setUp/SourceLoader.R"
wsPath = "../../setUp/SourceLoader.R"

wsPath = as.character(paste(funr::get_script_path(), "/../../setUp/SourceLoader.R", sep = ""))

wsPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/setUp/SourceLoader.R"

source(wsPath)
sourceFiles(c("helperFunctions"))
sourceFiles(c("UltraQuickRepeatedSubSampling"))
sourceFiles(c("TriangulateIsoSurface"))

path2Manifold = "../../Manifold/build/"
datasetPath = "../../data/ModelNet10/ModelNet10/"

# path2Manifold = "/home/willy/PredictingProteinInteractions/Manifold/build/"
# datasetPath = "/home/willy/PredictingProteinInteractions/data/ModelNet10/"

path2Manifold = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/Manifold/build/"
datasetPath = "/home/sysgen/Documents/LWB//PredictingProteinInteractions/data/ModelNet10/"


# install.packages("igraph")

# install.packages("spatstat")
# install.packages("deldir")

# installations needed
# sudo apt-get install gfortran    # deldir -> spatstat
# apt-get install mesa-common-dev  # rgl

# library(rgl)

library(lubridate)
library(rdist)

getModel10Net <- function(fName, plot = FALSE){
  
  tab = read.table(fName,skip=2, col.names = c("x","y", "z", "trash"), fill = TRUE)
  tab_xyz = tab[is.na(tab$trash),1:3]
  
  if(plot == TRUE) points3d(tab_xyz)
  return(tab_xyz)
}

getDataSet <- function(datasetPath){
  allFiles_off = list.files(datasetPath, recursive = TRUE, full.names = TRUE,pattern = ".off")
  
  dataSet = data.frame(matrix(0,nrow = length(allFiles_off), ncol = 3))
  colnames(dataSet) = c("model", "trainTest", "file")

  for(i in 1:length(allFiles_off)){
    vec = strsplit(allFiles_off[i],split = "/")[[1]]
    vecSmall = vec[(length(vec)-1):length(vec)]
    
    trainTest = vecSmall[1]
    modelName = strsplit(vecSmall[2],split = ".off")[[1]]
    
    my_print(allFiles_off[i])
    
    objNameIn = paste(strsplit(allFiles_off[i],split = ".off")[[1]][1],".objPre", sep ="")
    objNameOut = paste(strsplit(allFiles_off[i],split = ".off")[[1]][1],".obj", sep ="")
    # print(objName)
    
    if(!file.exists(objNameIn)){
      my_print(paste("Creating obj file ", objNameIn ,"...", sep =""))
      
      # from off 2 obj
      system(paste("off2obj ",allFiles_off[i] ," > ", objNameIn, sep =""))
    }
    
    if(!file.exists(objNameOut)){
      # make watertight
      # obj = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/bathtub/test/bathtub_0110.obj"
      
      manifoldCommand = "./manifold"
      args = paste(" ",objNameIn," ",objNameOut, " 2000 ", sep="")
      system(paste(path2Manifold,manifoldCommand,args, sep =""))
      
    }
    
    dataSet[i,] = c(modelName, trainTest, objNameOut)
  }  

  return(dataSet)
}

generateF_approximations_3dModel <- function(model_points, n = 100, m = 10, q = 2){
  
  pos13_F_list = list()
  pos13_F_approx_list = list()
  
  for(i in 1:m){
    pos13_F_list[[i]] = samplePointsAndCalculateCDFofEc(all_pts = model_points, n = n,plot = FALSE)
    pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
  }
  
  return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
}

getAllModel_F_approximations <- function(model_vec, n = 100, m = 50, q = 2){
  
  distributions_lists = list()
  
  for(i in 1:length(model_vec)){
    my_print(i/length(model_vec))
    F_app = generateF_approximations_3dModel(model_vec[[i]]$vert,n = n,q = q, m = m)
    distributions_lists[[i]] =  list("name" = model_vec[[i]]$name,"F" = F_app)
  }
  
  return(distributions_lists)
}

getAllModels <- function(dataSet, maxNum = NULL){
  all_models = list()
  
  if(is.null(maxNum)) maxNum = nrow(dataSet)
  
  for(i in 1:maxNum){
    my_print(paste(i/maxNum))
    m1 = getModel10Net(dataSet$file[i], FALSE)
    n1 = dataSet$model[i]
    all_models[[i]] = list("vert" = m1, "name" = n1)
  }
  
  return(all_models)
}


# model_rgl = read.obj("/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10//bathtub/train/bathtub_0034.obj", convert.rgl = FALSE)
# # model_rgl_plot = read.obj(objPath, convert.rgl = TRUE)
# 
# 
# 
# downsampleEuclideanAndGetGeodesicModel10Net("/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10//bathtub/train/bathtub_0034.obj",
#                                             n_s_euclidean = 1000,
#                                             n_s_dijkstra = 50,
#                                             plot = TRUE)

# install.packages("rdist")
library(rdist)
downsampleEuclideanAndGetGeodesicModel10Net <- function(objPath, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE, verbose = FALSE){
  model_rgl = read.obj(objPath, convert.rgl = FALSE)
  model_rgl_plot = read.obj(objPath, convert.rgl = TRUE)
  
  if(verbose) my_print("plotting")
  if(plot) shade3d(model_rgl_plot)
  
  points = t(model_rgl$shapes[[1]]$positions)
  edges = t(model_rgl$shapes[[1]]$indices)+1
  
  
  # model = read.ply(objPath)
  # points = t(model$vb)[,1:3]
  # edges = t(model$it)
  # 
  # if(plot) shade3d(model)
  
  my_print(paste("model has ", nrow(points),", points", nrow(edges), " edges", sep ="" ))
  
  ob = preProcessMesh(points = points, edges = edges, plot = FALSE)
  my_print(paste("processed model has ", nrow(ob$points), "points", sep ="" ))
  
  if(ob$numOfConComps != 1){
    print(paste(" ... and ", ob$numOfConComps, " connected components.", sep =""))
    return(NULL)
  }
  
  graph = ob$graph
  edges = ob$edges
  
  
  my_print("step 1: euclidean fps ...")
  
  sampled_indices = c(1:nrow(points))
  if(n_s_euclidean < nrow(points)){
    sampled_indices = myFarthestPointSampling(points, k = n_s_euclidean)
  }
  
  if(plot) {
    points3d(points[sampled_indices,], size = 10, col = "green")
  }
  
  sampled_indices2 = c(1:length(sampled_indices))
  # furthermore subsample with the distances on the surface
  my_print("step 2: surface distance of sampled points ...")
  d_surface = myShortestDistances(graph,sampled_indices)
  
  my_print("step 3: surface fps ...")
  # ??farthest_point_sampling
  
  sampled_indices2 = c(1:n_s_euclidean)
  if(n_s_dijkstra != n_s_euclidean){
    fps_surface <- farthest_point_sampling(d_surface)
    sampled_indices2 = fps_surface[1:n_s_dijkstra]
  }


  # rgl.open()
  if(FALSE) plotDownsampledPoints(model_rgl,points,sampled_indices[sampled_indices2],FALSE, col = "blue", size = 51)
  
  v2 = myVoronoi(points[sampled_indices[sampled_indices2],], points)
  v_n = v2/sum(v2)
  
  
  d_euclid = as.matrix(dist(x = points[sampled_indices[sampled_indices2],]))
  
  l = list("centers" = points[sampled_indices[sampled_indices2],], "mu" = v_n, "indices_order" = sampled_indices[sampled_indices2], "d_surface" = d_surface[sampled_indices2,sampled_indices2], "sampled_indices_geo" = sampled_indices, "d_euclid" = d_euclid)
  return(l)
}



getSmallDataSet <- function(dataSet,NumOfObjectsFromEachClass = 10, onlyTrain = TRUE){
  
  dataSetTrain = dataSet
  if(onlyTrain) dataSetTrain = dataSet[which(dataSet[,2] == "train"),]
  
  classNames = rep("", nrow(dataSetTrain))
  
  for(i in 1:length(dataSetTrain[,1])){
    t = strsplit(dataSetTrain[i,1], split = "_")[[1]]
    classNames[i] = paste(t[1:(length(t)-1)], collapse = '_')
  }
  
  classLevels = unique(classNames)
  numClasses = length(classLevels)
  
  print(classLevels)
  
  smallDataSet = data.frame(matrix(0,ncol = 3, nrow = numClasses*NumOfObjectsFromEachClass))
  colnames(smallDataSet) = colnames(dataSetTrain)
  
  for(i in 1:length(classLevels)){
    inds = which(classNames == classLevels[i])[1:NumOfObjectsFromEachClass]
    
    start = (i-1)*NumOfObjectsFromEachClass+1
    end = start+NumOfObjectsFromEachClass-1
    
    smallDataSet[start:end,] = dataSetTrain[inds,]
  }
  
  return(smallDataSet)
}

# install.packages("lubridate")

getSurfaceSampledModels <- function(dataSet, n_s_euclidean = 1000, n_s_dijkstra = 100, plot = TRUE){
  
  models = list()
  
  times = rep(0,nrow(dataSet))
  
  start.time <- Sys.time()
  for(i in 1:nrow(dataSet)){
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    start.time <- Sys.time()
    
    times[i] = time.taken
    eta = 0
    if(i > 10){
      avgLast10 = mean(times[(i-10):i])
      
      eta = seconds_to_period((nrow(dataSet)-i)*avgLast10)
    }
    
    my_print(paste(dataSet[i,3], i/nrow(dataSet), ", eta: ",eta ),0)
    
    
    # close all other rgl-windows
    while (rgl.cur() > 0) { rgl.close() }
    
    # get folder
    t = strsplit(dataSet[i,3],"/")[[1]]
    dir = paste(t[1:(length(t)-1)],collapse = "/")
    
    # check if distance folder is allready existent
    distancesDir = paste(dir, "/Distances/", sep ="")
    if(!dir.exists(distancesDir)) dir.create(distancesDir)
    
    distanceFile_surf = getGeoDistanceName(path = distancesDir,ind = 0,n_s_euclidean = n_s_euclidean,n_s_dijkstra = n_s_dijkstra,fname = paste(dataSet[i,1],"_surf",sep = ""))
    distanceFile_euclid = getGeoDistanceName(path = distancesDir,ind = 0,n_s_euclidean = n_s_euclidean,n_s_dijkstra = n_s_dijkstra,fname = paste(dataSet[i,1],"_euclid",sep = ""))
    
    # check if distance-file allready exists
    if(!file.exists(distanceFile_surf) || !file.exists(distanceFile_euclid) ){
      mod = downsampleEuclideanAndGetGeodesicModel10Net(objPath = dataSet[i,3], n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra, plot = plot)
      
      if(is.null(mod)){
        write.csv("",file = distanceFile_surf, row.names = FALSE)
        write.csv("",file = distanceFile_euclid, row.names = FALSE)
      } else {
        write.csv(mod$d_surface,file = distanceFile_surf, row.names = FALSE)
        write.csv(mod$d_euclid,file = distanceFile_euclid, row.names = FALSE)
        
        if(plot) {
          points3d(mod$centers, col = "blue", size = 20)
          rgl.snapshot(paste(distancesDir, "/", dataSet[i,1], ".png", sep =""))
        }
      }
    }
    
    info_surf = file.info(distanceFile_surf)
    info_euclid = file.info(distanceFile_euclid)
    
    if(info_surf$size > 1 && info_euclid$size > 1){
      distance_matrix_surf = read.csv(distanceFile_surf,header = TRUE)
      distance_matrix_euclid = read.csv(distanceFile_euclid,header = TRUE)
      
      models[[i]] = list("d_surface" = distance_matrix_surf, "d_euclid" = distance_matrix_euclid, "name" = dataSet[i,1], "path" = dataSet[i,3], "n_s_euclidean" = n_s_euclidean, "n_s_dijkstra" = n_s_dijkstra, "testTrain" = dataSet[i,2])
    }
  }
  
  
  models[sapply(models, is.null)] <- NULL
  
  return(models)
}

distributionOfDE <- function(models,
                             AllQuantilesPath = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDirStandard/",
                             n = 10,
                             m = 3,
                             q = 1,
                             mode = "Distances",
                             plot = TRUE,
                             recalculate = FALSE){
  
  quantilesOut = data.frame(matrix(0,ncol = (q+2)*2+1+1, nrow = 0))
  colnames(quantilesOut) = c("class","testTrain", as.vector(paste(c("q"),c(1:(q+2)), sep = "")), as.vector(paste(c("q_euc"),c(1:(q+2)), sep = "")))
  
  
  for(i in 1:length(models)){
    print(paste(models[[i]]$name, i/length(models)))
    
    path = strsplit(models[[i]]$path,split = models[[i]]$name)[[1]][1]
    quantilesDir = paste(path, "/Quantiles/", sep = "")
    
    if(!dir.exists(quantilesDir)) dir.create(quantilesDir)
    
    quantilesName = getGeoDistanceQuantileName(quantilesDir,mode,models[[i]]$n_s_euclidean,models[[i]]$n_s_dijkstra,n = n,m = m,q = q, fname = models[[i]]$name)
    if(!file.exists(quantilesName) || recalculate == TRUE){
      Fapp = generateF_approximations_3dModelWithMetric(d_surface = models[[i]]$d_surface, d_euclid = models[[i]]$d_euclid, n = n,m = m,q = q,mode = mode)

      # q+2 quantiles per distance + 1 name
      quantiles = data.frame(matrix(0,ncol = (q+2)*2+1+1, nrow = m))
      colnames(quantiles) = c("class","testTrain", as.vector(paste(c("q"),c(1:(q+2)), sep = "")), as.vector(paste(c("q_euc"),c(1:(q+2)), sep = "")))
      
      quantiles[,1] = rep(models[[i]]$name, m)
      quantiles[,2] = rep(models[[i]]$testTrain, m)
      for(i in 1:length(Fapp$F_app_list)){
        quantiles[i,3:(q+4)] = Fapp$F_app_list[[i]]
        quantiles[i,(q+5):ncol(quantiles)] = Fapp$F_app_euclid[[i]]
      }
      
      write.csv(x = quantiles, quantilesName)
    }
    quantiles = read.csv(quantilesName, row.names = 1)
    
    quantilesOut = rbind(quantilesOut,quantiles)
  }
  
  
  if(!dir.exists(AllQuantilesPath)) dir.create(AllQuantilesPath)
  allQuantiles = getGeoDistanceQuantileName(AllQuantilesPath,mode,models[[1]]$n_s_euclidean,models[[1]]$n_s_dijkstra,n = n,m = m,q = q, fname = "All")
  
  write.csv(quantilesOut,allQuantiles,row.names = FALSE)
  
  return(quantilesOut)
}



distributionOfDEParallel <- function(models,
                             AllQuantilesPath = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDirStandard/",
                             n = 10,
                             m = 3,
                             q = 1,
                             mode = "Distances",
                             plot = TRUE,
                             recalculate = FALSE){




  quantilesOutMat = t(apply(matrix(c(1:length(models),nrow = 1)), 1,  FUN = function(i){
                        print(paste(models[[i]]$name, i/length(models)))

                        path = strsplit(models[[i]]$path,split = models[[i]]$name)[[1]][1]
                        quantilesDir = paste(path, "/Quantiles/", sep = "")

                        if(!dir.exists(quantilesDir)) dir.create(quantilesDir)

                        quantilesName = getGeoDistanceQuantileName(quantilesDir,mode,models[[i]]$n_s_euclidean,models[[i]]$n_s_dijkstra,n = n,m = m,q = q, fname = models[[i]]$name)
                        if(!file.exists(quantilesName) || recalculate == TRUE){
                          Fapp = generateF_approximations_3dModelWithMetric(d_surface = models[[i]]$d_surface, d_euclid = models[[i]]$d_euclid, n = n,m = m,q = q,mode = mode)

                          # q+2 quantiles per distance + 1 name
                          quantiles = data.frame(matrix(0,ncol = (q+2)*2+1+1, nrow = m))
                          colnames(quantiles) = c("class","testTrain", as.vector(paste(c("q"),c(1:(q+2)), sep = "")), as.vector(paste(c("q_euc"),c(1:(q+2)), sep = "")))

                          quantiles[,1] = rep(models[[i]]$name, m)
                          quantiles[,2] = rep(models[[i]]$testTrain, m)
                          for(i in 1:length(Fapp$F_app_list)){
                            quantiles[i,3:(q+4)] = Fapp$F_app_list[[i]]
                            quantiles[i,(q+5):ncol(quantiles)] = Fapp$F_app_euclid[[i]]
                          }

                          write.csv(x = quantiles, quantilesName)
                        }
                        quantiles = read.csv(quantilesName, row.names = 1)

                        quantiles
                      }))
  
  quantilesOut = data.frame(matrix(0,ncol = (q+2)*2+1+1, nrow = m*length(models)))
  colnames(quantilesOut) = c("class","testTrain", as.vector(paste(c("q"),c(1:(q+2)), sep = "")), as.vector(paste(c("q_euc"),c(1:(q+2)), sep = "")))


  if(!dir.exists(AllQuantilesPath)) dir.create(AllQuantilesPath)
  allQuantiles = getGeoDistanceQuantileName(AllQuantilesPath,mode,models[[1]]$n_s_euclidean,models[[1]]$n_s_dijkstra,n = n,m = m,q = q, fname = "All")

  # write.csv(quantilesOut,allQuantiles,row.names = FALSE)

  return(quantilesOut)
}

distributionOfDEAllLocalities <- function( models,
                                           AllQuantilesPath = "/home/willy/PredictingProteinInteractions/data/ModelNet10/AllQuantilesDir/",
                                           n = 10,
                                           m = 3,
                                           q = 1,
                                           plot = TRUE){
  #-------------------------------------------------------------------------------
  # for each point in the model take the closest n neighbors and calculate the DE
  #-------------------------------------------------------------------------------
  
  quantilesOut = data.frame(matrix(0,ncol = q+3, nrow = 0))
  colnames(quantilesOut) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")))
  
  for(i in 1:length(models)){
    print(paste(models[[i]]$name, i/length(models)))
    
    path = strsplit(models[[i]]$path,split = models[[i]]$name)[[1]][1]
    quantilesDir = paste(path, "/QuantilesNearestNeighbors/", sep = "")
    
    quantilesName = getGeoDistanceQuantileName(quantilesDir,0,models[[i]]$n_s_euclidean,models[[i]]$n_s_dijkstra,n = n,m = m,q = q, fname = models[[i]]$name)
    
    if(!dir.exists(quantilesDir)) dir.create(quantilesDir)
    
    d_tmp = models[[i]]$d_surface
    
    if(!file.exists(quantilesName)){
      Fapp = generateF_approximations_3dModelWithMetricAllAroundSurface(d = d_tmp, n = n, q = q)
      
      quantiles = data.frame(matrix(0,ncol = q+3, nrow = nrow(d_tmp)))
      colnames(quantiles) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")))
      
      quantiles[,1] = rep(models[[i]]$name, nrow(d_tmp))
      for(i in 1:length(Fapp$F_app_list)){
        quantiles[i,2:ncol(quantiles)] = Fapp$F_app_list[[i]]  
      }
      
      write.csv(x = quantiles, quantilesName)
    }
    quantiles = read.csv(quantilesName, row.names = 1)
    
    quantilesOut = rbind(quantilesOut,quantiles)
  }
  
  if(!dir.exists(AllQuantilesPath)) dir.create(AllQuantilesPath)
  allQuantiles = getGeoDistanceQuantileName(AllQuantilesPath,0,models[[1]]$n_s_euclidean,models[[1]]$n_s_dijkstra,n = n,m = m,q = q, fname = "All")
  
  write.csv(quantilesOut,allQuantiles)
  
  return(quantilesOut)
}


plotQuantiles <- function(quantiles, euklid = FALSE){
  classes = getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")
  classLevels = unique(classes)
  numOfClasses = length(classLevels)
  
  print(classLevels)
  
  print(paste("number of classes: ", numOfClasses))
  
  colMap = c("red", "blue", "green", "black", "brown", "yellow", "pink", "lightgreen", "darkblue", "grey")
  
  
  for(i in 1:numOfClasses){
    inds = which(classes == classLevels[i])
    points3d(quantiles[inds,2:4], col = colMap[(i-1)*2+1])
    
    if(euklid) points3d(quantiles[inds,5:7], col = colMap[(i-1)*2+2])
    
  }
}


plot3Examples <- function(quantiles){
  plotQuantiles(quantiles)
  inds2 = which(quantiles[,1] == "sofa_0003")
  points3d(quantiles[inds2,2:4], col = "green", size = 20)
  points3d(quantiles[inds2,5:7], col = "green", size = 20)
  
  inds3 = which(quantiles[,1] == "sofa_0004")
  points3d(quantiles[inds3,2:4], col = "lightgreen", size = 20)
  points3d(quantiles[inds3,5:7], col = "lightgreen", size = 20)
  
  inds4 = which(quantiles[,1] == "bathtub_0001")
  points3d(quantiles[inds4,2:4], col = "pink", size = 20)
  points3d(quantiles[inds4,5:7], col = "pink", size = 20)
  
  inds5 = which(quantiles[,1] == "toilet_0019")
  points3d(quantiles[inds5,2:4], col = "brown", size = 20)
  points3d(quantiles[inds5,5:7], col = "brown", size = 20)
}


transformQuants <- function(quantiles, q = 1){
  tranf = quantiles
  
  quant = q+2
  
  for(i in 1:nrow(tranf)){
    for(j in 2:quant){
      tranf[i,1+j] = (tranf[i,1+j])/tranf[i,2]
    }
  }
  
  return(tranf)
}

#------------------------------------------------------------------------

# datasetPath = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/"
# datasetPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/ModelNet10/"


# get all the file names and information if it belongs to train or test
dataSet = getDataSet(datasetPath)
nrow(dataSet)

# get the first 20 models from each class
# smallDataSet = getSmallDataSet(dataSet,2)

smallDataSet = dataSet

smallDataSet = na.omit(smallDataSet)

nrow(smallDataSet)

GLOBAL_VERBOSITY = 1
models = getSurfaceSampledModels(smallDataSet,plot = FALSE,n_s_euclidean = 1000,n_s_dijkstra = 100)
# quit()


for(i in 1:length(models)){
  if(nrow(models[[i]]$d_surface) == 0 || nrow(models[[i]]$d_euclid) == 0) {
    print(i)
    models[[i]] <- NULL
    i = i-1
  }
}

length(models)

saveRDS(models, file = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/AllModels_ne_1000_nd_100.Rdata")


# # randomly sample and calculate DE
quantilesDist = distributionOfDE(models = models,
                             n = 40,
                             m =100,
                             mode = "Distances",
                             q = 10,
                             recalculate = FALSE)


# quantilesDist = distributionOfDEParallel(models = models[1:10],
#                                  n = 10,
#                                  m =1,
#                                  mode = "Distances",
#                                  q = 2,
#                                  recalculate = FALSE)
# quantilesDist


# quantilesEcc = distributionOfDE(models = models,
#                                  n = 1000,
#                                  m =1,
#                                  mode = "Eccentricities",
#                                  q = 1)

# quantilesDist
# 
# quantilesDist_trans = transformQuants(quantilesDist, q=1)
# plotQuantiles(quantilesDist_trans,FALSE)
# plotQuantiles(quantilesDist,FALSE)
# 
# 
# 
# rgl.open()
# plotQuantiles(quantilesDist,FALSE)
# plotQuantiles(quantilesEcc,TRUE)
# 
# cor(quantilesDist[,-1],quantilesEcc[,-1])
# 






