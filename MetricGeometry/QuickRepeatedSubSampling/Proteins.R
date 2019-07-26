library(rgl)

s1 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R"
source(s1)

s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s2)


s3 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
source(s3)

library(keras)
library(readobj)
library(FNN)


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
  
  quantiles = data.frame(matrix(0,ncol = (q+2)*3+1, nrow = m))
  colnames(quantiles) = c("name",   as.vector(paste(c("q_pos"),c(1:(q+2)), sep = "")), 
                          as.vector(paste(c("q_neg"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_neg"),c(1:(q+2)), sep = "")))
  
  quantiles[,1] = rep(model$name, m)
  for(i in 1:m){
    quantiles[i,(1:(q+2))+1] = F_app[[i]]$F_pos_approx
    quantiles[i,(1:(q+2))+1+(q+2)] = F_app[[i]]$F_neg_approx
    quantiles[i,(1:(q+2))+1+(q+2)*2] = F_app[[i]]$F_pos_neg_approx
  }
  
  return(quantiles)
}

getFApproximationsSumMethod <- function(mod, n = 10,m = 100,q=1){
  pos_indices = which(mod$posNegVector == TRUE)
  neg_indices = which(mod$posNegVector == FALSE)
  
  F_approximations = list()
  for(i in 1:m){
    pos_indices_samp = sample(pos_indices, size = n/2, replace = FALSE)
    neg_indices_samp = sample(neg_indices, size = n/2, replace = FALSE)
    
    d_pos = mod$d_surface[pos_indices_samp,pos_indices_samp]
    
    F_pos = DistributionOfEccentricities(d_pos)
    F_pos_approx = approximateCDF(F_pos, q)
    
    d_neg = mod$d_surface[neg_indices_samp,neg_indices_samp]
    F_neg = DistributionOfEccentricities(d_neg)
    F_neg_approx = approximateCDF(F_neg, q)
    
    F_pos_neg = DistributionOfEccentricities(mod$d_surface[c(pos_indices_samp,neg_indices_samp),c(pos_indices_samp,neg_indices_samp)])
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
                               n_s_dijkstra = 500,
                               onlyTheseIndices = NULL,
                               recalculate = FALSE){
  
  dirs = list.dirs(path = path, recursive = FALSE, full.names = FALSE)
  if(!is.null(onlyTheseIndices)) dirs = dirs[onlyTheseIndices]
  
  proteinModels = list()
  for(i in 1:length(dirs)){
    proteinModels[[i]] = getProteinModelStichedSurface(path = path, protName = dirs[i], n_s_euclidean = n_s_euclidean,n_s_dijkstra = n_s_dijkstra, plot = FALSE, recalculate = recalculate)
  }
  
  return(proteinModels)
}


getProteinModelStichedSurface <- function(path = "/home/willy/Schreibtisch/106Test/Output/",
                                          protName = "000_Trx",
                                          n_s_euclidean = 1000,
                                          n_s_dijkstra = 500,
                                          plot = FALSE,
                                          recalculate = FALSE){
  # first stich the model together
  fNameOrigingal = paste(path, "/", protName, "/", protName, ".obj", sep ="")
  fNameStitched = paste(path, "/", protName, "/", protName, "_stitched.obj", sep = "")
  
  fNameModelDownsampled = paste(path, "/", protName, "/", protName, "_model_downsampled.rData", sep ="")
  
  if(!file.exists(fNameModelDownsampled) || recalculate == TRUE){
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


get_quantile_name_protein <- function(name, n,m,q){
  paste(name,"_n_",n,"_m_",m,"_q_",q,".csv",sep ="")
}

get_quantiles_protein <- function(path = "/home/willy/Schreibtisch/106Test/Output/", model, n, m, q, recalculate = FALSE){
  quantFolder = paste(path = path, "/", model$name, "/Quantiles/", sep ="")
  if(!dir.exists(quantFolder)) dir.create(quantFolder)
  
  quantFileName = paste(quantFolder, get_quantile_name_protein("quantiles", n, m, q), sep = "")
  
  if(!file.exists(quantFileName) || recalculate == TRUE){
    quantiles = sampleEccentricitiesAndGetQuantiles(model = model, n = n, m = m, q = q)
    
    write.csv(quantiles,file = quantFileName, row.names = FALSE)
  }
  quantiles = read.csv(quantFileName, header = TRUE)
  return(quantiles)
}


get_quantiles_all_proteins <- function(model_vec, path, n, m ,q, recalculate = FALSE){
  quantiles_all = data.frame(matrix(0,ncol = (q+2)*3+1, nrow = m*length(model_vec)))
  colnames(quantiles_all) = c("name",   as.vector(paste(c("q_pos"),c(1:(q+2)), sep = "")), 
                              as.vector(paste(c("q_neg"),c(1:(q+2)), sep = "")),
                              as.vector(paste(c("q_pos_neg"),c(1:(q+2)), sep = "")))
  
  for(i in 1:length(model_vec)){
    print(paste(model_vec[[i]]$name, i/length(model_vec)))
    
    quant = get_quantiles_protein(path = path, model = model_vec[[i]],n = n, m = m, q = q, recalculate = recalculate)

    start_ind = (i-1)*m+1
    end_ind = start_ind + m-1

    # print(quant)    
    quantiles_all[start_ind:end_ind,] = quant
    quantiles_all[start_ind:end_ind,1] = as.character(quant[,1])
  }
  
  return(quantiles_all)
}

library(rgl)

# plot_prot_quants(quantiles)

plot_prot_quants <- function(quantiles){
  
  points3d(quantiles[,((1:(q+2))+1)], col = "red")
  points3d(quantiles[,((1:(q+2))+1 +(q+2))], col = "blue")
  points3d(quantiles[,((1:(q+2))+1 +(q+2)*2)], col = "green")
}


#-------------------------------------------------------------------------------------------------------------

n_s_euclidean = 1000
n_s_dijkstra = 500
n = 10
m = 100
q = 1
path = "/home/willy/Schreibtisch/106Test/Output/"

prot2 = getProteinModelStichedSurface(protName = "003", plot = TRUE)

models = getAllProteinModels()

quantiles = get_quantiles_all_proteins(model_vec = models, path = path, n = 100, m = 100, q = 1, recalculate = TRUE)


plot_prot_quants(quantiles)


