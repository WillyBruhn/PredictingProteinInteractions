

if(WS == FALSE){
  s1 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R"
  s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
  s3 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
}

source(s1)
source(s2)
source(s3)


library(keras)
library(readobj)
library(FNN)

install.packages(file.choose(), repos=NULL)

install.packages("/home/sysgen/Documents/LWB/PredictingProteinInteractions/packages/rgl/")
setwd("/home/sysgen/Documents/LWB/");
packages<-dir();
install.packages("rgl", repos=NULL)

library(rgl)


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

sampleEccentricitiesAndGetQuantiles <- function(model,
                                                n = 10,
                                                m = 2,
                                                q = 1)
{
  F_app = getFApproximationsSumMethod(mod = model, n = n, m = m, q = q)
  
  quantiles = data.frame(matrix(0,ncol = (q+2)*6+1, nrow = m))
  colnames(quantiles) = c("name",
                          as.vector(paste(c("q_pos"),c(1:(q+2)), sep = "")), 
                          as.vector(paste(c("q_neg"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_neg"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_euclid"),c(1:(q+2)), sep = "")), 
                          as.vector(paste(c("q_neg_euclid"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_neg_euclid"),c(1:(q+2)), sep = "")))
  
  quantiles[,1] = rep(model$name, m)
  for(i in 1:m){
    quantiles[i,(1:(q+2))+1] = F_app[[i]]$F_pos_approx
    quantiles[i,(1:(q+2))+1+(q+2)] = F_app[[i]]$F_neg_approx
    quantiles[i,(1:(q+2))+1+(q+2)*2] = F_app[[i]]$F_pos_neg_approx
    
    quantiles[i,(1:(q+2))+1+(q+2)*3] = F_app[[i]]$F_pos_approx_euclid
    quantiles[i,(1:(q+2))+1+(q+2)*4] = F_app[[i]]$F_neg_approx_euclid
    quantiles[i,(1:(q+2))+1+(q+2)*5] = F_app[[i]]$F_pos_neg_approx_euclid
  }
  
  return(quantiles)
}

getFApproximationsSumMethod <- function(mod, n = 10,m = 100,q=1){
  pos_indices = which(mod$posNegVector == TRUE)
  neg_indices = which(mod$posNegVector == FALSE)
  
  F_approximations = list()
  for(i in 1:m){
    pos_indices_samp = sample(pos_indices, size = n*length(pos_indices), replace = FALSE)
    neg_indices_samp = sample(neg_indices, size = n*length(neg_indices), replace = FALSE)
    
    #----------------------------------------------------------------------------
    # surface
    d_pos = mod$d_surface[pos_indices_samp,pos_indices_samp]
    F_pos = DistributionOfEccentricities(d_pos, mod$measure[pos_indices_samp])
    F_pos_approx = approximateCDF(F_pos, q)
    
    d_neg = mod$d_surface[neg_indices_samp,neg_indices_samp]
    F_neg = DistributionOfEccentricities(d_neg, mod$measure[neg_indices_samp])
    F_neg_approx = approximateCDF(F_neg, q)
    
    merged = c(pos_indices_samp,neg_indices_samp)
    F_pos_neg = DistributionOfEccentricities(mod$d_surface[merged,merged], mod$measure[merged])
    F_pos_neg_approx = approximateCDF(F_pos_neg, q)
    
    #----------------------------------------------------------------------------
    # euclid
    d_pos_euclid = as.matrix(dist(mod$centers[pos_indices_samp,]))
    F_pos_euclid  = DistributionOfEccentricities(d_pos_euclid, mod$measure[pos_indices_samp])
    F_pos_approx_euclid = approximateCDF(F_pos_euclid, q)
    
    d_neg_euclid = as.matrix(dist(mod$centers[neg_indices_samp,]))
    F_neg_euclid = DistributionOfEccentricities(d_neg_euclid, mod$measure[neg_indices_samp])
    F_neg_approx_euclid = approximateCDF(F_neg_euclid, q)
    
    merged = c(pos_indices_samp,neg_indices_samp)
    F_pos_neg_euclid = DistributionOfEccentricities(as.matrix(dist(mod$centers[merged,])), mod$measure[merged])
    F_pos_neg_approx_euclid = approximateCDF(F_pos_neg_euclid, q)
    
    F_approximations[[i]] = list("F_pos_approx" = F_pos_approx,
                                 "F_neg_approx" = F_neg_approx,
                                 "F_pos_neg_approx" = F_pos_neg_approx,
                                 "F_pos_approx_euclid" = F_pos_approx_euclid,
                                 "F_neg_approx_euclid" = F_neg_approx_euclid,
                                 "F_pos_neg_approx_euclid" = F_pos_neg_approx_euclid)
  }
  return(F_approximations)
}

# 
# hist(models5[[1]]$measure)
# hist(models10[[1]]$measure)
# models5[[1]]$measure - models20[[1]]$measure
# getFApproximationsSumMethod(models5[[1]], n = 1, m = 1)



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
                               n_s_euclidean = 5000,
                               n_s_dijkstra = 5000,
                               stitchNum = 10000,
                               measureNearestNeighbors = 20,
                               onlyTheseIndices = NULL,
                               recalculate = FALSE){
  
  dirs = list.dirs(path = path, recursive = FALSE, full.names = FALSE)
  if(!is.null(onlyTheseIndices)) dirs = dirs[onlyTheseIndices]
  
  proteinModels = list()
  for(i in 1:length(dirs)){
    proteinModels[[i]] = getProteinModelStichedSurface(path = path, protName = dirs[i], n_s_euclidean = n_s_euclidean,
                                                       n_s_dijkstra = n_s_dijkstra, plot = FALSE,
                                                       recalculate = recalculate,
                                                       stitchNum =stitchNum,
                                                       measureNearestNeighbors = measureNearestNeighbors)
  }
  
  return(proteinModels)
}

plotProteinModel <- function(fNameOrigingal = NULL,lis, openNew = TRUE, sz = 1){
  # points3d(points)
  
  if(openNew) rgl.open()
  points3d(lis$centers[which(lis$posNegVector == FALSE),], col = "red", size = sz)
  points3d(lis$centers[which(lis$posNegVector == TRUE),], col = "blue", size = sz)
  
  sc = sz*(1/min(lis$measure) +1)
  
  for(i in 1:nrow(lis$centers)){
    points3d(x = lis$centers[i,1], y = lis$centers[i,2], z = lis$centers[i,3], col = "green", size = lis$measure[i]*sc)
  }
  
  if(!is.null(fNameOrigingal)){
    rglModOrigPlot = read.obj(fNameOrigingal, convert.rgl = TRUE)
    shade3d(rglModOrigPlot)
  }

  # ?rgl.bg
  rgl.bg( sphere = FALSE, color = "white")
}

getProteinModelStichedSurface <- function(path = "/home/willy/Schreibtisch/106Test/Output/",
                                          protName = "000_Trx",
                                          n_s_euclidean = 5000,
                                          n_s_dijkstra = 5000,
                                          stitchNum = 5000,
                                          measureNearestNeighbors = 20,
                                          plot = FALSE,
                                          recalculate = FALSE){
  # first stich the model together
  fNameOrigingal = paste(path, "/", protName, "/", protName, ".obj", sep ="")
  fNameStitched = paste(path, "/", protName, "/", protName, "_",stitchNum,"_stitched.obj", sep = "")
  
  fNameModelDownsampled = paste(path, "/", protName, "/", protName, "_sitch_",stitchNum, "_nE_",n_s_euclidean, "_nD_", n_s_dijkstra,
                                "_model_downsampled.rData", sep ="")
  
  if(!file.exists(fNameModelDownsampled) || recalculate == TRUE){
    if(!file.exists(fNameStitched) || recalculate == TRUE){
      # make watertight
      # obj = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/bathtub/test/bathtub_0110.obj"
      path2Manifold = "/home/willy/Manifold/build/"
      manifoldCommand = "./manifold"
      args = paste(" ",fNameOrigingal," ",fNameStitched, " ",stitchNum," ", sep="")
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
    
    lis = list("centers" = centers, "posNegVector" = posNegVector, "d_surface" = d_surface, "name" = protName)
    
    mu = getProteinMeasure(lis, neighbors = measureNearestNeighbors)
    lis = list("centers" = centers, "posNegVector" = posNegVector, "d_surface" = d_surface, "name" = protName, "measure" = mu, "measureNearestNeighbors" = measureNearestNeighbors)

    saveRDS(lis,file = fNameModelDownsampled)
  }
  
  lis = readRDS(fNameModelDownsampled)
  
  my_print(fNameModelDownsampled, 1)
  
  if(is.null(lis$measureNearestNeighbors) || lis$measureNearestNeighbors != measureNearestNeighbors) {
    my_print(paste("recalculating measure for ", fNameModelDownsampled, sep = ""), 1)
    
    mu = getProteinMeasure(lis, neighbors = measureNearestNeighbors)
    lis = list("centers" = lis$centers, "posNegVector" = lis$posNegVector, "d_surface" = lis$d_surface, "name" = lis$name, "measure" = mu, "measureNearestNeighbors" = measureNearestNeighbors)
    saveRDS(lis,file = fNameModelDownsampled)
  }
  
    
  
  if(plot){
    plotProteinModel(fNameOrigingal, lis)
  }
  
  return(lis)
}


get_quantile_name_protein <- function(name, n,m,q, measureNN){
  paste(name,"_n_",n,"_m_",m,"_q_",q,"_muNN_",measureNN, ".csv",sep ="")
}

get_quantiles_protein <- function(path = "/home/willy/Schreibtisch/106Test/Output/", model, n, m, q, measureNN, recalculate = FALSE){
  quantFolder = paste(path = path, "/", model$name, "/Quantiles/", sep ="")
  if(!dir.exists(quantFolder)) dir.create(quantFolder)
  
  quantFileName = paste(quantFolder, get_quantile_name_protein("quantiles", n, m, q, measureNN), sep = "")
  
  if(!file.exists(quantFileName) || recalculate == TRUE){
    quantiles = sampleEccentricitiesAndGetQuantiles(model = model, n = n, m = m, q = q)
    
    write.csv(quantiles,file = quantFileName, row.names = FALSE)
  }
  
  quantiles = read.csv(quantFileName, header = TRUE, colClasses=c("character",rep("numeric",(q+2)*3)))
  return(quantiles)
}


get_quantiles_all_proteins <- function(model_vec, path, n, m ,q, recalculate = FALSE){
  # quantiles_all = data.frame(matrix(0,ncol = (q+2)*3+1, nrow = m*length(model_vec)))
  # colnames(quantiles_all) = c("name",   as.vector(paste(c("q_pos"),c(1:(q+2)), sep = "")), 
  #                             as.vector(paste(c("q_neg"),c(1:(q+2)), sep = "")),
  #                             as.vector(paste(c("q_pos_neg"),c(1:(q+2)), sep = "")))
  
  
  quantiles_all = data.frame(matrix(0,ncol = (q+2)*6+1, nrow = m))
  colnames(quantiles_all) = c("name",
                          as.vector(paste(c("q_pos"),c(1:(q+2)), sep = "")), 
                          as.vector(paste(c("q_neg"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_neg"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_euclid"),c(1:(q+2)), sep = "")), 
                          as.vector(paste(c("q_neg_euclid"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_neg_euclid"),c(1:(q+2)), sep = "")))
  
  vec = strsplit(path, "/")
  QuantFolder = paste(paste(vec[[1]][1:(length(vec[[1]])-1)],collapse = "/"),"/Quantiles/", sep ="")
  
  measureNN = model_vec[[1]]$measureNearestNeighbors
  QuantName = paste(QuantFolder,get_quantile_name_protein(name = "All",n = n, m = m, q = q, measureNN = measureNN), sep ="")
  
  if(!file.exists(QuantName) || recalculate == TRUE){
    if(!dir.exists(QuantFolder)) dir.create(QuantFolder)
    
    for(i in 1:length(model_vec)){
      print(paste(model_vec[[i]]$name, i/length(model_vec)))
      
      quant = get_quantiles_protein(path = path, model = model_vec[[i]],n = n, m = m, q = q, recalculate = recalculate, measureNN = measureNN)
  
      start_ind = (i-1)*m+1
      end_ind = start_ind + m-1
  
      # print(quant)    
      quantiles_all[start_ind:end_ind,] = quant
      quantiles_all[start_ind:end_ind,1] = as.character(quant[,1])
    }
    
    write.csv(x = quantiles_all,file = QuantName, row.names = FALSE)
  } else {
    quantiles_all = read.csv(file = QuantName, header = TRUE)
  }
  
  return(quantiles_all)
}



plot_one_prot_quant <- function(quantiles,q, col = "yellow", size = 15){
  
  quants = (1:(q+2))
  pos = quants+1
  neg = quants+1 + length(quants)
  pos_neg = quants+1 + length(quants)*2
  
  pos_euclid = quants+1 + length(quants)*3
  neg_euclid = quants+1 + length(quants)*4
  pos_neg_euclid = quants+1 + length(quants)*5
  
  
  points3d(quantiles[,pos], col = col, size = size)
  points3d(quantiles[,neg], col = col, size = size)
  points3d(quantiles[,pos_neg], col = col, size = size)
  
  points3d(quantiles[,pos_euclid], col = col, size = size)
  points3d(quantiles[,neg_euclid], col = col, size = size)
  points3d(quantiles[,pos_neg_euclid], col = col, size = size)
  
  geo = c(sum(quantiles[,2]), sum(quantiles[,3]),sum(quantiles[,4]))/nrow(quantiles)
  print(geo)
  
  print(paste(min(quantiles[,2]), max(quantiles[,2])))
  
  text3d(quantiles[1,1], x = geo[1], y = geo[2], z = geo[3], cex = 5, col = "yellow")
}

plot_prot_quants <- function(quantiles, q, functionals, plotMode = "Booth", withEuclid = FALSE){

  functionalInds = which(quantiles[,1] %in% functionals)
  
  
  nonFunctionalInds = c(1:nrow(quantiles))[-functionalInds]
  
  quants = (1:(q+2))
  pos = quants+1
  neg = quants+1 + length(quants)
  pos_neg = quants+1 + length(quants)*2
  
  pos_euclid = quants+1 + length(quants)*3
  neg_euclid = quants+1 + length(quants)*4
  pos_neg_euclid = quants+1 + length(quants)*5
  
  print(pos)
  print(neg)
  print(pos_neg)
  
  if(plotMode == "Booth" || plotMode == "onlyNonFunctional"){
    points3d(quantiles[nonFunctionalInds,pos], col = "red")
    points3d(quantiles[nonFunctionalInds,neg], col = "blue")
    points3d(quantiles[nonFunctionalInds,pos_neg], col = "green")

    if(withEuclid){
      points3d(quantiles[nonFunctionalInds,pos_euclid], col = "brown")
      points3d(quantiles[nonFunctionalInds,neg_euclid], col = "pink")
      points3d(quantiles[nonFunctionalInds,pos_neg_euclid], col = "lightblue")
    }
  }
  
  if(plotMode == "Booth" || plotMode == "onlyFunctional") {

    points3d(quantiles[functionalInds,pos], col = "red", size = 10)
    points3d(quantiles[functionalInds,neg], col = "blue", size = 10)
    points3d(quantiles[functionalInds,pos_neg], col = "green", size = 10)
    
    if(withEuclid){
      points3d(quantiles[functionalInds,pos_euclid], col = "brown",size = 10)
      points3d(quantiles[functionalInds,neg_euclid], col = "pink", size = 10)
      points3d(quantiles[functionalInds,pos_neg_euclid], col = "lightblue", size = 10)
    }
  }
}

# plot_prot_quants(quantiles, q=1, functionals, plotMode = "Booth", withEuclid = TRUE)

getProteinMeasure <- function(model,neighbors = 20){
  
  measure = rep(1,nrow(model$centers))
  for(i in 1:nrow(model$centers)){
    nearestPoints = which.minn(model$d_surface[i,], n = neighbors)
    if(model$posNegVector[i] == FALSE && sum(model$posNegVector[nearestPoints]) > 0){
      measure[i] = sum(model$posNegVector[nearestPoints])
    }
    
    if(model$posNegVector[i] == TRUE && sum(model$posNegVector[nearestPoints]) < neighbors){
      measure[i] = (neighbors - sum(model$posNegVector[nearestPoints]))
    }
  }
  
  measure = normalize(measure)
  
  return(measure)
}

#-------------------------------------------------------------------------------------------------------------

n_s_euclidean = 1000
n_s_dijkstra = 500
n = 10
m = 100
q = 1
# path = "/home/willy/Schreibtisch/106Test/Output/"
path = "/home/willy/PredictingProteinInteractions/data/120Experiment/Output/"

GLOBAL_VERBOSITY = 2
models1 = getAllProteinModels(path = path, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 1, recalculate = FALSE)
quantiles1 = get_quantiles_all_proteins(model_vec = models1, path = path, n = 1, m = 1, q = 1, recalculate = FALSE)

models5 = getAllProteinModels(path = path, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 5, recalculate = FALSE)
quantiles5 = get_quantiles_all_proteins(model_vec = models5, path = path, n = 1, m = 1, q = 1, recalculate = FALSE)


models10 = getAllProteinModels(path = path, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 10, recalculate = FALSE)
quantiles10 = get_quantiles_all_proteins(model_vec = models10, path = path, n = 1, m = 1, q = 1, recalculate = FALSE)

models20 = getAllProteinModels(path = path, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 20, recalculate = FALSE)
quantiles20 = get_quantiles_all_proteins(model_vec = models20, path = path, n = 1, m = 1, q = 1, recalculate = FALSE)

models500 = getAllProteinModels(path = path, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 500, recalculate = FALSE)
quantiles500 = get_quantiles_all_proteins(model_vec = models500, path = path, n = 1, m = 1, q = 1, recalculate = FALSE)

models10[[1]]$measureNearestNeighbors
models5[[1]]$measureNearestNeighbors
models20[[1]]$measureNearestNeighbors

quantiles10[1:1,1:5]
quantiles5[1:1,1:5]
quantiles20[1:1,1:5]

quantiles10 + quantiles5

t = read.table("/home/willy/PredictingProteinInteractions/data/labels120.txt", header = TRUE)
functionals = tolower(as.character(t[which(t[,2] == "glutathionineBinding"),1]))

inds = which(quantiles[,1] %in% functionals)
notInds = c(1:nrow(quantiles))[-inds]


plot_prot_quants(quantiles5, q=1, functionals, plotMode = "Booth", withEuclid = TRUE)
plot_prot_quants(quantiles10, q=1, functionals, plotMode = "Booth", withEuclid = TRUE)
plot_prot_quants(quantiles20, q=1, functionals, plotMode = "Booth", withEuclid = TRUE)
plot_prot_quants(quantiles500, q=1, functionals, plotMode = "Booth", withEuclid = TRUE)
plot_prot_quants(quantiles1, q=1, functionals, plotMode = "Booth", withEuclid = TRUE)
# plot_prot_quants



colnames(quantiles)[11:13]
df = cbind(as.character(quantiles[,1]), rep(0.1,nrow(quantiles)),quantiles[,c(11:13)])


getGeometricCenters(df, 0.1)
getGeometricCenters


m = as.matrix(dist(quantiles[,8:10], method = "manhattan"))
colnames(m) = quantiles[,1]
rownames(m) = quantiles[,1]


colnames(m)[which.minn(m[7,],n = 106)]

colnames(quantiles)

distances = list()
for(i in 1:106){
  distances[[i]] = dist(matrix(quantiles[i,-1], nrow = 3), method = "manhattan")
  
  if(colnames(m)[i] %in% functionals) points3d(distances[[i]][1], distances[[i]][2], distances[[i]][3], col ="red", size = 10)
  else points3d(distances[[i]][1], distances[[i]][2], distances[[i]][3], col ="blue")
}


print(distances[[1]])



functionals = c(getFunctionalProteins(), "000_Trx")





prot_indices = which(quantiles[,1] == "000_Trx")

plot_prot_quants(quantiles[prot_indices,],functionals = functionals,q)

q = 1
plot_one_prot_quant(quantiles[which(quantiles[,1] == "000_Trx"),],q = q,col = "black", size = 20)
plot_one_prot_quant(quantiles[which(quantiles[,1] == "027"),],q = q,col = "black", size = 20)
plot_one_prot_quant(quantiles[which(quantiles[,1] == "013"),],q = q,col = "black", size = 20)
plot_one_prot_quant(quantiles[which(quantiles[,1] == "016"),],q = q,col = "black", size = 20)
plot_one_prot_quant(quantiles[which(quantiles[,1] == "053"),],q = q,col = "black", size = 20)

plot_one_prot_quant(quantiles[which(quantiles[,1] == "002"),],q = q,col = "black", size = 20)

functionals

