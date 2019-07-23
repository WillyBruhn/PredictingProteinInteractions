#---------------------------------------------------------------------
# Willy Bruhn, 19.7.19
#
# Read in the model-net10 data-set and create the usual projections.
#
#---------------------------------------------------------------------

s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
source(s1)

s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s2)

s3 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R"
source(s3)

library(rgl)

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
    
    print(allFiles_off[i])
    
    objNameIn = paste(strsplit(allFiles_off[i],split = ".off")[[1]][1],".objPre", sep ="")
    objNameOut = paste(strsplit(allFiles_off[i],split = ".off")[[1]][1],".obj", sep ="")
    # print(objName)
    
    if(!file.exists(objNameIn)){
      print(paste("Creating obj file ", objNameIn ,"...", sep =""))
      
      # from off 2 obj
      system(paste("off2obj ",allFiles_off[i] ," > ", objNameIn, sep =""))
    }
    
    if(!file.exists(objNameOut)){
      # make watertight
      # obj = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/bathtub/test/bathtub_0110.obj"
      path2Manifold = "/home/willy/Manifold/build/"
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
    print(i/length(model_vec))
    F_app = generateF_approximations_3dModel(model_vec[[i]]$vert,n = n,q = q, m = m)
    distributions_lists[[i]] =  list("name" = model_vec[[i]]$name,"F" = F_app)
  }
  
  return(distributions_lists)
}

getAllModels <- function(dataSet, maxNum = NULL){
  all_models = list()
  
  if(is.null(maxNum)) maxNum = nrow(dataSet)
  
  for(i in 1:maxNum){
    print(paste(i/maxNum))
    m1 = getModel10Net(dataSet$file[i], FALSE)
    n1 = dataSet$model[i]
    all_models[[i]] = list("vert" = m1, "name" = n1)
  }
  
  return(all_models)
}


downsampleEuclideanAndGetGeodesicModel10Net <- function(objPath, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE, doGeo = TRUE){
  model_rgl = read.obj(objPath, convert.rgl = FALSE)
  model_rgl_plot = read.obj(objPath, convert.rgl = TRUE)
  
  print("plotting")
  if(plot) shade3d(model_rgl_plot)
  
  points = t(model_rgl$shapes[[1]]$positions)
  edges = t(model_rgl$shapes[[1]]$indices)+1
  
  
  # model = read.ply(objPath)
  # points = t(model$vb)[,1:3]
  # edges = t(model$it)
  # 
  # if(plot) shade3d(model)
  
  print(paste("model has ", nrow(points),", points", nrow(edges), " endges", sep ="" ))
  
  ob = preProcessMesh(points = points, edges = edges, plot = FALSE)
  print(paste("processed model has ", nrow(ob$points), "points", sep ="" ))
  
  if(ob$numOfConComps != 1){
    print(paste(" ... and ", ob$numOfConComps, " connected components.", sep =""))
    return(NULL)
  }
  
  graph = ob$graph
  edges = ob$edges
  
  # return(ob)
  
  # points = points[unique(unlist(graph[[c(1:length(graph[1]))]])),]
  
  # return(graph)
  
  library(rdist)
  print("step 1: euclidean fps ...")
  
  sampled_indices = c(1:nrow(points))
  if(n_s_euclidean < nrow(points)){
    sampled_indices = myFarthestPointSampling(points, k = n_s_euclidean)
  }
  
  if(plot) {
    points3d(points[sampled_indices,], size = 10, col = "green")
  }
  
  sampled_indices2 = c(1:length(sampled_indices))
  if(doGeo){
    # furthermore subsample with the distances on the surface
    print("step 2: surface distance of sampled points ...")
    d_surface = myShortestDistances(graph,sampled_indices)
    
    print("step 3: surface fps ...")
    # ??farthest_point_sampling
    fps_surface <- farthest_point_sampling(d_surface)
    sampled_indices2 = fps_surface[1:n_s_dijkstra]
  } else {
    d_surface = as.matrix(dist(points[sampled_indices2,]))
  }
  # rgl.open()
  if(FALSE) plotDownsampledPoints(model_rgl,points,sampled_indices[sampled_indices2],FALSE, col = "blue", size = 51)
  
  v2 = myVoronoi(points[sampled_indices[sampled_indices2],], points)
  v_n = v2/sum(v2)
  
  
  l = list("centers" = points[sampled_indices[sampled_indices2],], "mu" = v_n, "indices_order" = sampled_indices[sampled_indices2], "d_surface" = d_surface[sampled_indices2,sampled_indices2], "sampled_indices_geo" = sampled_indices)
  return(l)
}



getSmallDataSet <- function(dataSet,NumOfObjectsFromEachClass = 10){
  
  dataSetTrain = dataSet[which(dataSet[,2] == "train"),]
  
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

getSurfaceSampledModels <- function(dataSet, n_s_euclidean = 1000, n_s_dijkstra = 100, plot = TRUE){
  
  models = list()
  
  for(i in 1:nrow(dataSet)){
    print(paste(dataSet[i,3], i/nrow(dataSet)))
    
    # close all other rgl-windows
    while (rgl.cur() > 0) { rgl.close() }
    
    # get folder
    t = strsplit(dataSet[i,3],"/")[[1]]
    dir = paste(t[1:(length(t)-1)],collapse = "/")
    
    # check if distance folder is allready existent
    distancesDir = paste(dir, "/Distances/", sep ="")
    if(!dir.exists(distancesDir)) dir.create(distancesDir)
    
    distanceFile = getGeoDistanceName(path = distancesDir,ind = 0,n_s_euclidean = n_s_euclidean,n_s_dijkstra = n_s_dijkstra,fname = dataSet[i,1])
    
    # check if distance-file allready exists
    if(!file.exists(distanceFile)){
      mod = downsampleEuclideanAndGetGeodesicModel10Net(objPath = dataSet[i,3], n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra, plot = plot)
      
      if(is.null(mod)){
        # write.csv("NumberOfConCompsTooLarge",file = distanceFile, row.names = FALSE)
      } else {
        write.csv(mod$d_surface,file = distanceFile, row.names = FALSE)
        
        points3d(mod$centers, col = "blue", size = 20)
        rgl.snapshot(paste(distancesDir, "/", dataSet[i,1], ".png", sep =""))
      }
    }
    
    info = file.info(distanceFile)
    
    if(info$size > 1){
      distance_matrix = read.csv(distanceFile,header = TRUE)
      
      models[[i]] = list("d_surface" = distance_matrix, "name" = dataSet[i,1], "path" = dataSet[i,3], "n_s_euclidean" = n_s_euclidean, "n_s_dijkstra" = n_s_dijkstra)
    }
  }
  
  models[sapply(models, is.null)] <- NULL
  
  return(models)
}

distributionOfDE <- function(models,
                             n = 10,
                             m = 3,
                             q = 1,
                             plot = TRUE){
  
  quantilesOut = data.frame(matrix(0,ncol = q+3, nrow = 0))
  colnames(quantilesOut) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")))

  for(i in 1:length(models)){
    print(paste(models[[i]]$name, i/length(models)))
    
    path = strsplit(models[[i]]$path,split = models[[i]]$name)[[1]][1]
    quantilesDir = paste(path, "/Quantiles/", sep = "")
    
    if(!dir.exists(quantilesDir)) dir.create(quantilesDir)
    
    quantilesName = getGeoDistanceQuantileName(quantilesDir,0,models[[i]]$n_s_euclidean,models[[i]]$n_s_dijkstra,n = n,m = m,q = q, fname = models[[i]]$name)
    if(!file.exists(quantilesName)){
      Fapp = generateF_approximations_3dModelWithMetric(d = models[[i]]$d_surface, n = n,m = m,q = q)

      quantiles = data.frame(matrix(0,ncol = q+3, nrow = m))
      colnames(quantiles) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")))
      
      quantiles[,1] = rep(models[[i]]$name, m)
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
#------------------------------------------------------------------------
n = 30
m = 100
q = 10

pathToProjection = "/home/willy/PredictingProteinInteractions/data/ModelNet10/projections/"

datasetPath = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/"


# get all the file names and information if it belongs to train or test
# dataSet = getDataSet(datasetPath)

# get the first 20 models from each class
smallDataSet = getSmallDataSet(dataSet,20)

# apply farthest point sampling and store the geodesic distances
models = getSurfaceSampledModels(smallDataSet)

# randomly sample and calculate DE
quantiles = distributionOfDE(models = models,
                             n = 10,
                             m =1000,
                             q = 1)



classes = getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")
classLevels = unique(classes)
numOfClasses = length(classLevels)

colMap = c("red", "blue", "green", "black", "brown", "yellow", "pink", "lightgreen", "darkblue", "grey")


for(i in 1:10){
  inds = which(classes == classLevels[i])
  points3d(quantiles[inds,-1], col = colMap[i])  
}














