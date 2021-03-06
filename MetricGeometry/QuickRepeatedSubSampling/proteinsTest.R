library(rgl)

s1 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R"
source(s1)

s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s2)


s3 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
source(s3)

library(keras)


downsampleEuclideanAndGetGeodesicProtein <- function(objPath, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE){
  model_rgl = read.obj(objPath, convert.rgl = FALSE)
  model_rgl_plot = read.obj(objPath, convert.rgl = TRUE)
  
  print("plotting")
  if(plot) shade3d(model_rgl_plot)
  
  # pos at 2, neg at 3
  points = t(model_rgl$shapes[[2]]$positions)
  edges = t(model_rgl$shapes[[2]]$indices)+1
  
  
  # if(checkForLargerModel(path,proteinName = name,n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra) == FALSE)
  print(paste("model has ", nrow(points),", ", nrow(edges), " points", sep ="" ))
  
  ob = preProcessMesh(points = points, edges = edges, plot = FALSE)
  print(paste("processed model has ", nrow(ob$points), "points", sep ="" ))
  
  graph = ob$graph
  edges = ob$edges
  
  
  library(rdist)
  print("step 1: euclidean fps ...")
  
  sampled_indices = myFarthestPointSampling(points, k = n_s_euclidean)
  
  
  if(plot) {
    points3d(points[sampled_indices,], size = 10, col = "green")
  }
  
  return(graph)
  
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

getAllModelsProtein <- function(path = "/home/willy/Schreibtisch/106Test/Output/",
                                   name_prot = "001",
                                   n_s_euclidean = 100,
                                   n_s_dijkstra = 50,
                                   n = 10,
                                   m = 3,
                                   q = 1,
                                   positive = TRUE,
                                   plot = TRUE){
  
  path_final = paste(path, "/", name_prot, "/", sep ="")
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
      geoDistanceName_pos = getGeoDistanceName(distances_path,i,n_s_euclidean,n_s_dijkstra,fname = "pos")
      geoDistanceName_neg = getGeoDistanceName(distances_path,i,n_s_euclidean,n_s_dijkstra,fname = "neg")
      
      if(!file.exists(geoDistanceName_pos) || !file.exists(geoDistanceName_neg)){
        # mod = downsampleEuclideanAndGetGeodesicProtein(objPath = all.files[i],
        #                                         n_s_euclidean = n_s_euclidean,
        #                                         n_s_dijkstra = n_s_dijkstra,
        #                                         plot = plot)
        
        p = getProtein(path = path, name_prot, n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra,
                       plot = plot, reCalculate = TRUE)
        
        write.csv(x = p$prot_pos_sampled$geoDistances, file = geoDistanceName_pos,row.names = FALSE)
        write.csv(x = p$prot_neg_sampled$geoDistances, file = geoDistanceName_neg,row.names = FALSE)
      }
      
      geoDist_pos = read.csv(geoDistanceName_pos)
      geoDist_neg = read.csv(geoDistanceName_neg)
      
      Fapp_pos = generateF_approximations_3dModelWithMetric(d = geoDist_pos,n = n,m = m,q = q)
      Fapp_neg = generateF_approximations_3dModelWithMetric(d = geoDist_neg,n = n,m = m,q = q)
      
      quantiles = data.frame(matrix(0,ncol = q+3, nrow = 2*m),stringsAsFactors = FALSE)
      colnames(quantiles) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")))
      
      quantiles[1:m,1] = rep(paste(name_prot,"_pos",sep=""), m)
      print(name_prot)
      print(quantiles[,1])
      for(i in 1:length(Fapp_pos$F_app_list)){
        quantiles[i,2:ncol(quantiles)] = Fapp_pos$F_app_list[[i]]  
      }
      
      quantiles[(m+1):(2*m),1] = rep(paste(name_prot,"_neg",sep=""), m)
      for(i in 1:length(Fapp_neg$F_app_list)){
        quantiles[i+length(Fapp_pos$F_app_list),2:ncol(quantiles)] = Fapp_neg$F_app_list[[i]]  
      }
      
      write.csv(x = quantiles, quantilesName)
    }
    
    quantiles = read.csv(quantilesName, row.names = 1)
    
    quantilesOut = rbind(quantilesOut,quantiles)
  }
  
  return(quantilesOut)
}

getAllModelsProteinWrapper <- function(path = "/home/willy/Schreibtisch/106Test/Output/", prot_names = c("000_Trx"),
                                       n_s_euclidean = n_s_euclidean,
                                       n_s_dijkstra = n_s_dijkstra,
                                       n = n,
                                       m = m,
                                       q = q,
                                       positive = TRUE,
                                       plot = TRUE){
  
  quantilesOut = data.frame(matrix(0,ncol = q+3, nrow = 0), stringsAsFactors = FALSE)
  colnames(quantilesOut) = c("class", as.vector(paste(c("q"),c(1:(q+2)), sep = "")))
  
  QuantPath = paste(strsplit(path, "/Output/")[[1]][1],"/QuantileDistances/",sep = "")
  if(!dir.exists(QuantPath)) dir.create(QuantPath)
  
  quantFile = getGeoDistanceQuantileName(path = QuantPath,ind = 0,n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra,n = n,m = m,q = q,fname = "quant")
  
  for(i in 1:length(prot_names)){
    print(prot_names[i])
    
    quant  = getAllModelsProtein(path = path,
                                      name = prot_names[i],
                                      n_s_euclidean = n_s_euclidean,
                                      n_s_dijkstra = n_s_dijkstra,
                                      n = n,
                                      m = m,
                                      q = q,
                                      plot = plot)
    
    quantilesOut = rbind(quantilesOut,quant)
  }
  
  write.csv(quantilesOut,quantFile)
  
  return(quantilesOut)
}


n_s_euclidean = 500
n_s_dijkstra = 50
n = 48
m = 100
q = 1

path = "/home/willy/Schreibtisch/106Test/Output/"
name = "000_Trx"
# quantiles = getAllModelsProtein(path = path, name = "000_Trx",n_s_euclidean = 900,n_s_dijkstra = 500,n = 450,m = 500, q = 1,plot = TRUE)
# points3d(trx)

proteins = list.dirs(path = path, recursive = FALSE, full.names = FALSE)[1:15]

# proteins = c("046", "016","013", "027", proteins)

quantiles = getAllModelsProteinWrapper(path = path,
                           prot_names =  proteins,
                           n_s_euclidean = n_s_euclidean,
                           n_s_dijkstra = n_s_dijkstra,
                           n = n,
                           m = m,
                           q = q,
                           FALSE)


while (rgl.cur() > 0) { rgl.close() }
labels = readLabels("/home/willy/PredictingProteinInteractions/data/labels.txt")
functionals = c(getFunctionalProteins(), "000_Trx")

inds = which(labels$name %in% functionals)

# unique(quantiles[which(quantiles[,1] %in% functionals),1])


classes = unique(getClassNamesFromSubClasses(quantiles[,1]))
colors = as.numeric(as.factor(classes))

t = strsplit(classes[10], split = "_")[[1]]
paste(t[1:(length(t)-1)], collapse = '_')

for(i in 1:length(classes)){
  print(paste(classes[i]))
  
    inds = which(getClassNamesFromSubClasses(quantiles[,1]) == classes[i])
    
    t = strsplit(classes[i], split = "_")[[1]]
    name = paste(t[1:(length(t)-1)], collapse = '_')
    
    if(name %in% functionals) {
      points3d(x = quantiles[inds,2], y = quantiles[inds,3], z = quantiles[inds,4],col = "red", add = TRUE,size = 5)
    } else {
      points3d(x = quantiles[inds,2], y = quantiles[inds,3], z = quantiles[inds,4],col = "blue", add = TRUE,size = 5)
    }

    geox = sum(quantiles[inds,2])/length(inds)+0.01
    geoy = sum(quantiles[inds,3])/length(inds)
    geoz = sum(quantiles[inds,4])/length(inds)
    
    text3d(x = geox, y = geoy, z = geoz,texts = classes[i],cex = 2)
}

#---------------------------------------------------------------


#------------------------------------------------------------------
# proteins
#------------------------------------------------------------------
getClassNamesFromSubClassesProteins <- function(subClasses, splitPattern = "_"){
  # gets the className from a subclass
  # e.g.
  # Lion-01 --> Lion
  #-------------------------------------------------------
  
  classNames = rep("", length(subClasses))
  for(i in 1:length(subClasses)){
    t = strsplit(as.character(subClasses[i]),split = splitPattern)[[1]]
    classNames[i] = paste(t[1:(length(t)-1)], collapse = '_')
  }
  
  return(classNames)
}


transformToNodesProtein <- function(quantiles, sSize){
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

getNNInputFromQuantilesProtein <- function(quantiles,m, sampleSize=m, sampleTimes=1){
  # m number of distributions
  # sampleSize how many of the m distributions to select
  # sampleTimes how many times to select
  
  quantilesAsNNInput_pos = data.frame(matrix(0,ncol = (ncol(quantiles)-1)*sampleSize+1, nrow = sampleTimes*nrow(quantiles)/m))
  for(i in 1:(nrow(quantiles)/m)){
    print(i*sampleTimes/(nrow(quantilesAsNNInput)))
    
    start = (i-1)*m+1
    end = start+m-1
    
    for(j in 1:sampleTimes){
      inds = sample(c(start:end),size = sampleSize,replace = FALSE)
      quantilesAsNNInput_pos[(i-1)*(sampleTimes)+j,] = transformToNodesProtein(quantiles[inds,],sampleSize)
    }
  }
  
  quantilesAsNNInput_neg = data.frame(matrix(0,ncol = (ncol(quantiles)-1)*sampleSize+1, nrow = sampleTimes*nrow(quantiles)/m))
  for(i in 1:(nrow(quantiles)/m)){
    print(i*sampleTimes/(nrow(quantilesAsNNInput)))
    
    start = (i-1)*m+1
    end = start+m-1
    
    for(j in 1:sampleTimes){
      inds = sample(c(start:end),size = sampleSize,replace = FALSE)
      quantilesAsNNInput_neg[(i-1)*(sampleTimes)+j,] = transformToNodesProtein(quantiles[inds,],sampleSize)
    }
  }
  
  return(quantilesAsNNInput)
}


subClassNames = unique(quantiles[,1])
numObjects = length(subClassNames)
ClassNames = unique(getClassNamesFromSubClassesProteins(quantiles[,1]))

numClasses = length(unique(getClassNamesFromSubClassesProteins(quantiles[,1])))

ClassNames_test_ind = sample(c(1:numObjects), size = numObjects*0.3, replace = FALSE)
ClassNames_test =ClassNames[ClassNames_test_ind]
ClassNames_train =ClassNames[-ClassNames_test_ind]

quantiles[,1] %in% subClassNames_train

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


#----------------------------------------
fName = "/home/willy/Schreibtisch/106Test/Output/000_Trx/000_TrxModified.obj"

rglMod = read.obj(fName)

rglMod2 = read.obj(fName,convert.rgl = TRUE)

# shade3d(rglMod2)

points = t(rglMod$shapes[[1]]$positions)
edges = t(rglMod$shapes[[1]]$indices)+1

points3d(points)

rglModOrig = read.obj("/home/willy/Schreibtisch/106Test/Output/000_Trx/000_Trx.obj")

pointsOrig = t(rglModOrig$shapes[[2]]$positions)
points3d(pointsOrig, col = "red")

pointsOrigNeg = t(rglModOrig$shapes[[3]]$positions)
points3d(pointsOrigNeg, col = "blue")



ob = preProcessMesh(points = points, edges = edges, plot = FALSE)
print(paste("processed model has ", nrow(ob$points), "points", sep ="" ))

# nrow(points)
graph = ob$graph
edges = ob$edges

library(rdist)

sampled_indices = myFarthestPointSampling(points, k = 1000)

points3d(points[sampled_indices,], col = "green", size = 20)

d_surface = myShortestDistances(graph,sampled_indices)

fps_surface <- farthest_point_sampling(d_surface)
sampled_indices2 = fps_surface[1:500]
sampled_indices2


centers = points[sampled_indices[sampled_indices2],]
points3d(centers, col = "blue", size = 25)

# find out for each point if it belongs to positive or negative
# for each point get the distance to the closest point in positive and the closest distance to negative
# whichever is smaller determines if the point belongs to that potential

pos = knnx.dist(data = pointsOrig, query = t(centers[1,]),k = 1)
neg = knnx.dist(data = pointsOrigNeg, query = t(centers[1,]),k = 1)
points3d(x = centers[1,1], y = centers[1,2], z = centers[1,3], col = "green", size = 50)


pos = knnx.dist(data = pointsOrig, query = centers,k = 1)
neg = knnx.dist(data = pointsOrigNeg, query = centers,k = 1)


i = 100
pos[i]
neg[i]
points3d(x = centers[i,1], y = centers[i,2], z = centers[i,3], col = "green", size = 50)

posNeg = cbind(pos,neg)
posNegVector = posNeg[,1] < posNeg[,2]

# specifies for each point if it belongs to positive
posNegVector

length(which(posNegVector == TRUE))/length(posNegVector)



prot = getProteinModelStichedSurface(plot = TRUE)

prot2 = getProteinModelStichedSurface(protName = "003", plot = TRUE)

getAllProteinModels()

li = list(prot2)

models = li

li[[1]]$d_surface

mod = li[[1]]

sampleEccentricitiesAndGetQuantiles(mod, n = 500, q = 1)

sampleEccentricitiesAndGetQuantiles <- function(model,
                                                n = 10,
                                                m = 2,
                                                q = 1)
{
  F_app = getFApproximationsSumMethod(mod = model, n = n, m = m, q = q)
  
  quantiles = data.frame(matrix(0,ncol = (q+2)*3+1, nrow = m))
  colnames(quantiles) = c("name",   as.vector(paste(c("q_pos"),c(1:(q+2)), sep = "")), 
                                    as.vector(paste(c("q_neg"),c(1:(q+2)), sep = "")),
                                    as.vector(paste(c("q_pos_neg"),c(1:(q+2)), sep = "")))
  
  quantiles[,1] = rep(models[[1]]$name, m)
  for(i in 1:m){
    quantiles[i,(1:(q+2))+1] = F_app[[i]]$F_pos_approx
    quantiles[i,(1:(q+2))+1+(q+2)] = F_app[[i]]$F_neg_approx
    quantiles[i,(1:(q+2))+1+(q+2)*2] = F_app[[i]]$F_pos_neg_approx
  }
  
  return(quantiles)
}

getFApproximationsSumMethod <- function(mod, n = 10,m = 100,q=1){
  pos_indices = which(mod$posNegVector == TRUE)
  
  d_pos = mod$d_surface[pos_indices,pos_indices]
  d_neg = mod$d_surface[-pos_indices,-pos_indices]
  
  F_approximations = list()
  for(i in 1:m){
    indices = sample(c(1:nrow(mod$centers)), size = n, replace = FALSE)
    
    pos_indices_samp = which(indices %in% pos_indices)
    neg_indices_samp = pos_indices[-pos_indices_samp]
    
    d_pos = mod$d_surface[pos_indices_samp,pos_indices_samp]
    F_pos = DistributionOfEccentricities(d_pos)
    F_pos_approx = approximateCDF(F_pos, q)
    
    d_neg = mod$d_surface[neg_indices_samp,neg_indices_samp]
    F_neg = DistributionOfEccentricities(d_neg)
    F_neg_approx = approximateCDF(F_neg, q)
    
    F_pos_neg = DistributionOfEccentricities(mod$d_surface[indices,indices])
    F_pos_neg_approx = approximateCDF(F_pos_neg, q)
    
    F_approximations[[i]] = list("F_pos_approx" = F_pos_approx, "F_neg_approx" = F_neg_approx, "F_pos_neg_approx" = F_pos_neg_approx)
  }
  return(F_approximations)
}




generateF_approximations_Protein <- function(model_points, posNegVector, n = 100, m = 10, q = 2){
  
  pos13_F_list = list()
  pos13_F_approx_list = list()
  
  for(i in 1:m){
    pos13_F_list[[i]] = samplePointsAndCalculateCDFofEc(all_pts = model_points, n = n,plot = FALSE)
    pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
  }
  
  return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
}

getAllModel_F_approximations <- function(model_vec, n = 100, m = 50, q = 2, c_vec = c(1,1,1)){
  
  distributions_lists = list()
  
  for(i in 1:length(model_vec)){
    print(i/length(model_vec))
    F_app = generateF_approximations_Protein(model_vec[[i]]$centers, model_vec[[i]]$posNegVector, n = n,q = q, m = m)
    distributions_lists[[i]] =  list("name" = model_vec[[i]]$name,"F" = F_app)
  }
  
  return(distributions_lists)
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







getAllProteinModels = function(path = "/home/willy/Schreibtisch/106Test/Output/",
                               n_s_euclidean = 1000,
                               n_s_dijkstra = 500){
  
  dirs = list.dirs(path = path, recursive = FALSE, full.names = FALSE)
  
  proteinModels = list()
  for(i in 1:length(dirs)){
    proteinModels[[i]] = getProteinModelStichedSurface(path = path, protName = dirs[i], n_s_euclidean = n_s_euclidean,n_s_dijkstra = n_s_dijkstra, plot = FALSE)
  }
  
  return(proteinModels)
}


getProteinModelStichedSurface <- function(path = "/home/willy/Schreibtisch/106Test/Output/",
                                          protName = "000_Trx",
                                          n_s_euclidean = 1000,
                                          n_s_dijkstra = 500,
                                          plot = FALSE){
  # first stich the model together
  fNameOrigingal = paste(path, "/", protName, "/", protName, ".obj", sep ="")
  fNameStitched = paste(path, "/", protName, "/", protName, "_stitched.obj", sep = "")
  
  fNameModelDownsampled = paste(path, "/", protName, "/", protName, "_model_downsampled.rData", sep ="")
  
  if(!file.exists(fNameModelDownsampled)){
    if(!file.exists(fNameStitched)){
      # make watertight
      # obj = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/bathtub/test/bathtub_0110.obj"
      path2Manifold = "/home/willy/Manifold/build/"
      manifoldCommand = "./manifold"
      args = paste(" ",fNameOrigingal," ",fNameStitched, " 2000 ", sep="")
      system(paste(path2Manifold,manifoldCommand,args, sep =""))
      
    }
    
    rglModStitched = read.obj(fNameStitched)
    points = t(rglModStitched$shapes[[1]]$positions)
    edges = t(rglModStitched$shapes[[1]]$indices)+1
    
    rglModOrig = read.obj(fNameOrigingal)
    
    pointsOrig_pos = t(rglModOrig$shapes[[2]]$positions)
    pointsOrig_neg = t(rglModOrig$shapes[[3]]$positions)
    
    ob = preProcessMesh(points = points, edges = edges, plot = FALSE)
    print(paste("processed model has ", nrow(ob$points), "points", sep ="" ))
    
    # nrow(points)
    graph = ob$graph
    edges = ob$edges
    
    library(rdist)
    
    sampled_indices = myFarthestPointSampling(points, k = n_s_euclidean)
    
    # points3d(points[sampled_indices,], col = "green", size = 20)
    
    d_surface = myShortestDistances(graph,sampled_indices)
    
    fps_surface <- farthest_point_sampling(d_surface)
    sampled_indices2 = fps_surface[1:n_s_dijkstra]
    sampled_indices2
    
    d_surface = d_surface[sampled_indices2, sampled_indices2]
    
    
    centers = points[sampled_indices[sampled_indices2],]
    
    # find out for each point if it belongs to positive or negative
    # for each point get the distance to the closest point in positive and the closest distance to negative
    # whichever is smaller determines if the point belongs to that potential
    
    pos = knnx.dist(data = pointsOrig_pos, query = centers,k = 1)
    neg = knnx.dist(data = pointsOrig_neg, query = centers,k = 1)
    
    posNeg = cbind(pos,neg)
    posNegVector = posNeg[,1] < posNeg[,2]
    
    # specifies for each point if it belongs to positive
    posNegVector
    
    length(which(posNegVector == TRUE))/length(posNegVector)
    
    if(plot){
      # points3d(points)
      points3d(centers[which(posNegVector == FALSE),], col = "red")
      points3d(centers[which(posNegVector == TRUE),], col = "blue")
    }
    
    lis = list("centers" = centers, "posNegVector" = posNegVector, "d_surface" = d_surface, "name" = protName)
    
    saveRDS(lis,file = fNameModelDownsampled)
  }
  
  lis = readRDS(fNameModelDownsampled)
  
  return(lis)
}
