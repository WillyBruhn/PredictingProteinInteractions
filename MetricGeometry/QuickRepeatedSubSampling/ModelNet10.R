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


# model_rgl = read.obj("/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10//bathtub/train/bathtub_0034.obj", convert.rgl = FALSE)
# # model_rgl_plot = read.obj(objPath, convert.rgl = TRUE)
# 
# 
# 
# downsampleEuclideanAndGetGeodesicModel10Net("/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10//bathtub/train/bathtub_0034.obj",
#                                             n_s_euclidean = 1000,
#                                             n_s_dijkstra = 50,
#                                             plot = TRUE)


downsampleEuclideanAndGetGeodesicModel10Net <- function(objPath, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE){
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
  
  print(paste("model has ", nrow(points),", points", nrow(edges), " edges", sep ="" ))
  
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
  # furthermore subsample with the distances on the surface
  print("step 2: surface distance of sampled points ...")
  d_surface = myShortestDistances(graph,sampled_indices)
  
  print("step 3: surface fps ...")
  # ??farthest_point_sampling
  fps_surface <- farthest_point_sampling(d_surface)
  sampled_indices2 = fps_surface[1:n_s_dijkstra]

  # rgl.open()
  if(FALSE) plotDownsampledPoints(model_rgl,points,sampled_indices[sampled_indices2],FALSE, col = "blue", size = 51)
  
  v2 = myVoronoi(points[sampled_indices[sampled_indices2],], points)
  v_n = v2/sum(v2)
  
  
  d_euclid = as.matrix(dist(x = points[sampled_indices[sampled_indices2],]))
  
  l = list("centers" = points[sampled_indices[sampled_indices2],], "mu" = v_n, "indices_order" = sampled_indices[sampled_indices2], "d_surface" = d_surface[sampled_indices2,sampled_indices2], "sampled_indices_geo" = sampled_indices, "d_euclid" = d_euclid)
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
    
    if(info_surf$size > 1 && info_euclid$size){
      distance_matrix_surf = read.csv(distanceFile_surf,header = TRUE)
      distance_matrix_euclid = read.csv(distanceFile_euclid,header = TRUE)
      
      models[[i]] = list("d_surface" = distance_matrix_surf, "d_euclid" = distance_matrix_euclid, "name" = dataSet[i,1], "path" = dataSet[i,3], "n_s_euclidean" = n_s_euclidean, "n_s_dijkstra" = n_s_dijkstra)
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
      Fapp = generateF_approximations_3dModelWithMetric(d_surface = models[[i]]$d_surface, d_euclid = models[[i]]$d_euclid, n = n,m = m,q = q)

      # q+2 quantiles per distance + 1 name
      quantiles = data.frame(matrix(0,ncol = (q+2)*2+1, nrow = m))
      colnames(quantiles) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")), as.vector(paste(c("q_euc"),c(1:(q+2)), sep = "")))
      
      quantiles[,1] = rep(models[[i]]$name, m)
      for(i in 1:length(Fapp$F_app_list)){
        quantiles[i,2:(q+3)] = Fapp$F_app_list[[i]]
        quantiles[i,(q+4):ncol(quantiles)] = Fapp$F_app_euclid[[i]]
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


plotQuantiles <- function(quantiles){
  classes = getClassNamesFromSubClasses(quantiles[,1],splitPattern = "_")
  classLevels = unique(classes)
  numOfClasses = length(classLevels)
  
  print(classLevels)
  
  print(paste("number of classes: ", numOfClasses))
  
  colMap = c("red", "blue", "green", "black", "brown", "yellow", "pink", "lightgreen", "darkblue", "grey")
  
  
  for(i in 1:numOfClasses){
    inds = which(classes == classLevels[i])
    points3d(quantiles[inds,2:4], col = colMap[(i-1)*2+1])
    points3d(quantiles[inds,5:7], col = colMap[(i-1)*2+2])
  }
}

#------------------------------------------------------------------------
n = 30
m = 100
q = 10

pathToProjection = "/home/willy/PredictingProteinInteractions/data/ModelNet10/projections/"

datasetPath = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/"


# get all the file names and information if it belongs to train or test
dataSet = getDataSet(datasetPath)

# get the first 20 models from each class
smallDataSet = getSmallDataSet(dataSet,5)


smallDataSet = na.omit(smallDataSet)

nrow(smallDataSet)

sub = which(getClassNamesFromSubClasses(smallDataSet[,1], splitPattern = "_") %in% c("bathtub", "toilet", "sofa", "table"))
# sub = which(getClassNamesFromSubClasses(smallDataSet[,1], splitPattern = "_") %in% c("bathtub", "toilet"))
smallDataSet = smallDataSet[sub,]


# smallDataSet[1,]

# apply farthest point sampling and store the geodesic distances
models = getSurfaceSampledModels(smallDataSet,plot = FALSE,n_s_euclidean = 1000,n_s_dijkstra = 50)


for(i in 1:length(models)){
  if(nrow(models[[i]]$d_surface) == 0 || nrow(models[[i]]$d_euclid) == 0) {
    print(i)
    models[[i]] <- NULL
    i = i-1
  }
}

# length(models)


# # randomly sample and calculate DE
quantiles = distributionOfDE(models = models,
                             n = 48,
                             m =2,
                             q = 1)





plotQuantiles(quantiles)

# i = 10
# points3d(x = quantiles[i,2], y = quantiles[i,3], z = quantiles[i,4], col = "black", size = 10)
# points3d(x = quantiles[i,5], y = quantiles[i,6], z = quantiles[i,7], col = "black", size = 10)
# 
# i = 11
# points3d(x = quantiles[i,2], y = quantiles[i,3], z = quantiles[i,4], col = "red", size = 10)
# points3d(x = quantiles[i,5], y = quantiles[i,6], z = quantiles[i,7], col = "red", size = 10)

plot3Examples(quantiles)


quantilesLocal = distributionOfDEAllLocalities(models = models,
                                             n = 10,
                                             q = 1)

plotQuantiles(quantilesLocal)

quantilesLocal2 = distributionOfDEAllLocalities(models = models,
                                               n = 20,
                                               q = 1)
plotQuantiles(quantilesLocal2)


quantilesLocal3 = distributionOfDEAllLocalities(models = models,
                                                n = 10,
                                                q = 1)

plotQuantiles(quantilesLocal3)
plot3Examples(quantilesLocal3)

plot3Examples(quantilesLocal2)
plot3Examples(quantiles)

plot3Examples(quantilesLocal)


quantiles[which.max(quantiles[,2]),1]

models[12487]
12487/200

# points3d(models[[63]]$d_surface)
models[[63]]$path

mod = downsampleEuclideanAndGetGeodesicModel10Net(objPath = models[[63]]$path, n_s_euclidean = 1000, n_s_dijkstra = 100, plot = TRUE)

points3d(mod$centers, col ="blue", size = 20)


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




sub = which(getClassNamesFromSubClasses(quantilesLocal[,1], splitPattern = "_") %in% c("bathtub", "toilet", "sofa"))
subQuantiles = quantilesLocal[sub,]

plotQuantiles(subQuantiles)
rgl.open()
plotQuantiles(quantilesLocal)


quantOut = quantiles

for(i in 1:nrow(quantiles)){
  print(i/nrow(quantiles))
  quantOut[i,2] = quantiles[i,2]
  quantOut[i,3] = quantiles[i,3]-quantiles[i,2]
  quantOut[i,4] = quantiles[i,4]-quantiles[i,3]
}

quantiles = quantOut


# sub = which(getClassNamesFromSubClasses(quantiles[,1], splitPattern = "_") %in% c("sofa"))
# plotQuantiles(quantiles[sub,])
plotQuantiles(quantiles)

pca <- prcomp(quantiles[,2:4], center = FALSE,scale. = FALSE)

pca$sdev
pca$rotation

pca_coords = pca$x

pca_coords[,1] = pca_coords[,1] - min(pca_coords[,1])+1
pca_coords[,2] = pca_coords[,2] - min(pca_coords[,2])+1
pca_coords[,3] = pca_coords[,3] - min(pca_coords[,3])+1

logSubQuantiles = quantiles
logSubQuantiles[,2] = log(pca_coords[,1])
logSubQuantiles[,3] = log(pca_coords[,2])
logSubQuantiles[,4] = log(pca_coords[,3])

plotQuantiles(logSubQuantiles)
# plotQuantiles(quantiles)

# which(quantiles[,1] %in% c("sofa_0001"))
# 
# points3d(logSubQuantiles[which(logSubQuantiles[,1] %in% c("sofa_0001")),2:4], col = "blue", size = 20)
# points3d(logSubQuantiles[which(logSubQuantiles[,1] %in% c("table_0001")),2:4], col = "red", size = 20)
# 

# points3d(x = logSubQuantiles[,2], y = logSubQuantiles[,3], z = logSubQuantiles[,4], col = "red", alpha = 0.7, aspect = c(1, 1, 0.5))

x_dim = c(min(logSubQuantiles[,2]),max(logSubQuantiles[,2]))
y_dim = c(min(logSubQuantiles[,3]),max(logSubQuantiles[,3]))
z_dim = c(min(logSubQuantiles[,4]),max(logSubQuantiles[,4]))

x_width = x_dim[2]-x_dim[1]
y_width = y_dim[2]-y_dim[1]
z_width = z_dim[2]-z_dim[1]

# nZ = 100
# nY = nZ*d_y/d_z
# nX = nZ*d_x/d_z



gridSections = 5
dims = ceil(pca$sdev*gridSections)
# dims = ceil(c(100,100,100)*gridSections)
ceil(dims[1])*ceil(dims[2])*ceil(dims[3])

d_x = x_width/dims[1]
d_y = y_width/dims[2]
d_z = z_width/dims[3]

x_grid_vals = x_dim[1]+(seq(1:(dims[1]+1))-1)*d_x
y_grid_vals = y_dim[1]+(seq(1:(dims[2]+1))-1)*d_y
z_grid_vals = z_dim[1]+(seq(1:(dims[3]+1))-1)*d_z

x_grid_vals
y_grid_vals
z_grid_vals

# g = grid3d(c("x", "y+", "z"), n = dims)


# given one set of points that come from one object we want 
# to get the voxels in which these points are

numberOfDistributionsPerObject = 100
outGrid = data.frame(matrix(0, ncol = 3, nrow = numberOfDistributionsPerObject))
colnames(outGrid) = c("x","y","z")

whereInGrid <- function(x_grid_vals,val){
  if(val > x_grid_vals[length(x_grid_vals)]){
    return(length(x_grid_vals))
  }
  if(val < x_grid_vals[1]){
    return(1)
  }
  
  smaller = which(x_grid_vals <= val)
  greater = which(x_grid_vals >= val)-1
  inter = intersect(smaller,greater)
  return(inter)
}

# x_grid_vals
# length(x_grid_vals)
# whereInGrid(x_grid_vals,-700)

# currentObj[1,2]
# whereInGrid(x_grid_vals,currentObj[2,2])



getGridRepResentation <- function(quantiles,name){
  
  inds = which(quantiles[,1] == name)
  currentObj = quantiles[inds,]
  # dims specifies the number of borders in each dimension
  # for each distribution find the location in the outGrid
  for(i in 1:nrow(currentObj)){
    outGrid[i,1] = whereInGrid(x_grid_vals = x_grid_vals, val = currentObj[i,2])
    outGrid[i,2] = whereInGrid(x_grid_vals = y_grid_vals, val = currentObj[i,3])
    outGrid[i,3] = whereInGrid(x_grid_vals = z_grid_vals, val = currentObj[i,4])
  }
  return(outGrid)
}

p1 = getGridRepResentation(logSubQuantiles,"bathtub_0013")
p2 = getGridRepResentation(logSubQuantiles,"chair_0013")
# p3 = getGridRepResentation(logSubQuantiles,"chair_0014")
p4 = getGridRepResentation(logSubQuantiles,"chair_0015")

p5 = getGridRepResentation(logSubQuantiles,"chair_0018")

points3d(p1)
points3d(p2, col = "red")
# points3d(p3, col = "blue")
points3d(p4, col = "blue")
points3d(p5, col = "green")




