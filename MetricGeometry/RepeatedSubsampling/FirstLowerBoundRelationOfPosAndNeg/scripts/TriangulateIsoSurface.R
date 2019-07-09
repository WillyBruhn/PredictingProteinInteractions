#!/usr/bin/Rscript

source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/isoFaces.R")
source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/extrinsicDistances.R")

library(misc3d)
library(rgl)

library(geometry)
adjacencyMatrixFromMesh <- function(mesh, points, point_indices = c(1:nrow(points))){
  #-----------------------------------------------
  # mesh ...   triangulated mesh
  # points ... point positions
  # point_indices ... indices of points to keep
  # returns
  # The weighted adjacencymatrix
  # of the reduced model with only the points in point_indices
  #-----------------------------------------------
  
  # mat = matrix(0,ncol = max(mesh), nrow = max(mesh))
  mat = matrix(0,ncol = length(point_indices), nrow = length(point_indices))
  d = as.matrix(dist(points[point_indices,], method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
  
  for(i in 1:nrow(mesh)){
    # in each row all points are connected
    for(j in 1:ncol(mesh)){
      for(k in 1:ncol(mesh)){
        ind1 = mesh[i,j]
        ind2 = mesh[i,k]
        
        # we have to check if the points are within the selected ones
        if(ind1 %in% point_indices && ind2 %in% point_indices){
          # we have to get the index of the points in the sampled version
          # that means we need to kow how many indices have been thrown out
          # before ind1
          
          # example:
          # point_indices = (1,2,3,4,10)
          # ind1 := 3
          # which(point_indices <= ind1) # = 1 2 3
          #
          # ind_sampled = 3
          #
          # point_indices = (3,4,10)
          # ind1 := 3
          # which(point_indices <= ind1) # = 3
          #
          # ind_sampled = 1
          
          ind1_sampled = which(point_indices == ind1)
          ind2_sampled = which(point_indices == ind2)
          
          # print(paste(ind1,ind2))
          dist = d[ind1_sampled,ind2_sampled]
          # print(dist)
          
          
          mat[ind1_sampled,ind2_sampled] = dist
          mat[ind2_sampled,ind1_sampled] = dist
        }
      }
    }
  }
  return(mat)
}

# install.packages("igraph")
library(igraph)

plotMesh <- function(points,edges){
  points3d(points, col = "red")
  triangles3d(points[t(edges),])
}

edgesOfPoints <- function(pointIndices, edges){
  maxP = max(pointIndices)
  
  edges_keep = rep(FALSE, nrow(edges))
  
  for(i in 1:nrow(edges)){
    if(max(edges[i,]) <= maxP){
      edges_keep[i] = TRUE
    }
  }
  
  edges_subset = edges[which(edges_keep == TRUE),]
  
  return(edges_subset)
}

plotDownsampledPoints <- function(model,points, indices, plotModel = TRUE){
  if(plotModel == TRUE) shade3d(model)
  points3d(points[indices,], col = "red")
}

plotMemoliModel <- function(model,centers, measure, scale = 100, plotModel = TRUE){
  if(plotModel == TRUE) shade3d(model)
  
  for(i in 1:nrow(centers)){
    # points3d(centers[i,], col = "red", size = measure[i]*scale)
    plot3d(x = centers[i,1], y = centers[i,2], z = centers[i,3], add = TRUE, col = "green", size = measure[i]*scale, xlab = "", ylab = "", zlab = "")
    
  }

}

plotMemoliModel2 <- function(model,centers1, measure1, centers2, measure2, scale = 100, plotModel = TRUE){
  if(plotModel == TRUE) shade3d(model)
  
  for(i in 1:nrow(centers1)){
    # points3d(centers[i,], col = "red", size = measure[i]*scale)
    plot3d(x = centers1[i,1], y = centers1[i,2], z = centers1[i,3], add = TRUE, col = "green", size = measure1[i]*scale, xlab = "", ylab = "", zlab = "")
    
  }
  
  for(i in 1:nrow(centers2)){
    # points3d(centers[i,], col = "red", size = measure[i]*scale)
    plot3d(x = centers2[i,1], y = centers2[i,2], z = centers2[i,3], add = TRUE, col = "green", size = measure2[i]*scale, xlab = "", ylab = "", zlab = "")
    
  }
  
}

myVoronoi <- function(centers,points){
  measureCounts = rep(0,nrow(centers))
  
  d = rep(-1,nrow(centers))
  for(i in 1:nrow(points)){
    for(j in 1:nrow(centers)){
      d[j] = dist(rbind(points[i,], centers[j,]))
    }
    nearestCenter = which(d == min(d))
    
    measureCounts[nearestCenter] = measureCounts[nearestCenter] + 1 
  }
  
  return(measureCounts)
}

downsampleMemoliMethod <- function(model_rgl,points, adjacencyMatrix, n_s_euclidean = 4000, n_s_dijkstra = 50){
  
  # print("step 1: Setting up graph ...")
  # m = adjacencyMatrixFromMesh(edges, points, indices)
  
  library(igraph)
  # create graph from adjacency matrix
  g <- graph.adjacency(adjacencyMatrix, weighted=TRUE)
  
  # Get all path distances
  s.paths <- shortest.paths(g, algorithm = "dijkstra")
  
  # install.packages("rdist")
  library(rdist)
  
  print("step 2 Calculating distance-matrix ...")
  d_total = dist(points)
  
  # sample 2000 points in the first step 
  # n_s1 = 2000
  
  print("step 3: euclidean fps ...")
  fps <- farthest_point_sampling(d_total)
  sampled_indices = fps[1:n_s_euclidean]
  
  plotDownsampledPoints(model_rgl,points,sampled_indices)
  
  # furthermore subsample with the distances on the surface
  # n_s2 = 50
  
  print("step 4: surface fps ...")
  d_surface = s.paths[sampled_indices,sampled_indices]
  fps_surface <- farthest_point_sampling(d_surface)
  sampled_indices2 = fps_surface[1:n_s_dijkstra]
  
  plotDownsampledPoints(model_rgl,points,sampled_indices[sampled_indices2],TRUE)
  
  v2 = myVoronoi(points[sampled_indices[sampled_indices2],], points)
  
  v_n = v2/sum(v2)
  
  l = list("centers" = points[sampled_indices[sampled_indices2],], "mu" = v_n)
  return(l)
}

getMemoliModel <- function(model_rgl,points, indices, edges, memoliPtsPath, memoliMuPath, n_s_euclidean = 4000, n_s_dijkstra = 50){
  if(!file.exists(memoliPtsPath) || !file.exists(memoliMuPath)){
    print("step 1: Setting up graph ...")
    m = adjacencyMatrixFromMesh(edges, points, indices)
    
    mem = downsampleMemoliMethod(model_rgl, points = points[indices,], adjacencyMatrix = m, n_s_euclidean, n_s_dijkstra)
    
    # mem = downsampleMemoliMethod(model_rgl,points, indices, edges, n_s_euclidean, n_s_dijkstra)
    
    write.table(mem$centers,file=memoliPtsPath, row.names = FALSE, col.names = c("x","y","z"))
    write.table(mem$mu,file=memoliMuPath, row.names = FALSE, col.names = c("mu"))
  }
  
  pts = read.table(file = memoliPtsPath, header = TRUE, colClasses = c("numeric", "numeric", "numeric"), stringsAsFactors=FALSE)
  mu2 = unlist(read.table(file = memoliMuPath, header = TRUE, colClasses = c("numeric"), stringsAsFactors=FALSE))

  mu = as.numeric(mu2)
    
  l = list("centers" = pts, "mu" = mu)
  return(l)
}

getMemoliModel2 <- function(memoliObjPath, memoliPtsPath, memoliMuPath, n = 1){
  horse1 = read.obj(memoliObjPath, convert.rgl = FALSE)
  horse1_rgl = read.obj(memoliObjPath, convert.rgl = TRUE)
  shade3d(horse1_rgl)
  horse1_points = t(horse1$shapes[[1]]$positions)
  horse1_edges = t(horse1$shapes[[1]]$indices)+1
  
  horse1_points = myDownsampleIndices(horse1_points, n = n)
  
  horse_sampled = getMemoliModel(horse1_rgl,points = horse1_points, edges = horse1_edges, memoliPtsPath, memoliMuPath)
  
  plotMemoliModel(horse1_rgl,horse_sampled$centers, measure = horse_sampled$mu, scale = 500, plotModel = FALSE)
  
  l = list("rgl_model" = horse1_rgl, "points" = horse1_points, "centers" = horse_sampled$centers, "mu" = horse_sampled$mu)
  return(l)
}

getMemoliModel3 <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/horse-poses/", name){
  memoliModelPath = path
  model_name = name
  memoliPtsPath = paste(memoliModelPath, model_name, ".pts",sep ="")
  memoliMuPath = paste(memoliModelPath, model_name, ".mu",sep = "")
  memoliObjPath = paste(path, model_name, ".obj", sep ="")
  
  horse_sampled = getMemoliModel2(memoliObjPath = memoliObjPath, memoliPtsPath = memoliPtsPath, memoliMuPath = memoliMuPath)
  
  return(horse_sampled)
}
#--------------------------------------------------------------------------
# FLB method

get_s <- function(model_X){
  n_X = nrow(model_X$centers)
  d_X = as.matrix(dist(model_X$centers))
  s_X = rep(0,n_X)
  for(i in 1:n_X){
    s = 0
    for(k in 1:n_X){
      s = s + model_X$mu[k] * d_X[i,k]
    }
    
    s_X[i] = s
  }
  
  
  
  return(s_X)
}

FLB <- function(model_X,model_Y,s_X,s_Y,E,L = 10){
  # print(seq(0,max(E)+0.0001,1/(L-1)))
  
  # print(min(E))
  # print(max(E))
  
  h = hist(E, breaks = seq(0,max(E)+1,1/(L-1)), plot = FALSE)
  u_vals = h$breaks
  
  
  S_X = rep(0,L)
  S_Y = rep(0,L)
  for(u in 1:L){
    B_X_u = which(s_X <= u_vals[u])
    S_X[u] = sum(model_X$mu[B_X_u])
    
    B_Y_u = which(s_Y <= u_vals[u])
    S_Y[u] = sum(model_Y$mu[B_Y_u])
  }  
  
  flb = 0
  for(i in 1:(L-1)){
    u_1 = u_vals[i]
    u_2 = u_vals[i+1]
    
    flb = flb + abs(u_2 - u_1) * abs(S_X[i] - S_Y[i]) 
  }
  flb = flb*0.5
  
  return(flb)
}

FLB_method <- function(model_X,model_Y, L = 10){
  s_X = get_s(model_X)
  s_Y = get_s(model_Y)
  
  E_X = unique(s_X)
  E_Y = unique(s_Y)
  E = sort(unique(c(E_X,E_Y)))
  
  return(FLB(model_X,model_Y,s_X,s_Y,E,L))
}


FLB_protein <- function(model_X,model_Y,s_X_pos,s_Y_pos,E_pos, s_X_neg,s_Y_neg,E_neg, L = 10){
  
  h_pos = hist(E_pos, breaks = seq(0,max(E_pos)+1,1/(L-1)), plot = FALSE)
  u_vals_pos = h_pos$breaks
  
  
  S_X_pos = rep(0,L)
  S_Y_pos = rep(0,L)
  for(u in 1:L){
    B_X_u = which(s_X <= u_vals[u])
    S_X[u] = sum(model_X$mu[B_X_u])
    
    B_Y_u = which(s_Y <= u_vals[u])
    S_Y[u] = sum(model_Y$mu[B_Y_u])
  }  
  
  flb = 0
  for(i in 1:(L-1)){
    u_1 = u_vals[i]
    u_2 = u_vals[i+1]
    
    flb = flb + abs(u_2 - u_1) * abs(S_X[i] - S_Y[i]) 
  }
  flb = flb*0.5
  
  return(flb)
}

FLB_method_protein <- function(model_X,model_Y, L = 10){
  s_X_pos = get_s(model_X$prot_pos_sampled)
  s_Y_pos = get_s(model_Y$prot_pos_sampled)
  
  E_X_pos = unique(s_X_pos)
  E_Y_pos = unique(s_Y_pos)
  E_pos = sort(unique(c(E_X_pos,E_Y_pos)))
  
  s_X_neg = get_s(model_X$prot_neg_sampled)
  s_Y_neg = get_s(model_Y$prot_neg_sampled)
  
  E_X_neg = unique(s_X_neg)
  E_Y_neg = unique(s_Y_neg)
  E_neg = sort(unique(c(E_X_neg,E_Y_neg)))
  
  
  
  return(FLB(model_X,model_Y,s_X,s_Y,E,L))
}


flb_distances <- function(m1,m2,m3,m4,m5,m6,m7,m8, names){
  d = matrix(0,nrow = 8, ncol = 8)
  
  d[1,1]=FLB_method(m1,m1)
  d[1,2]=FLB_method(m1,m2)
  d[1,3]=FLB_method(m1,m3)
  d[1,4]=FLB_method(m1,m4)
  d[1,5]=FLB_method(m1,m5)
  d[1,6]=FLB_method(m1,m6)
  d[1,7]=FLB_method(m1,m7)
  d[1,8]=FLB_method(m1,m8)
  
  d[2,1]=FLB_method(m2,m1)
  d[2,2]=FLB_method(m2,m2)
  d[2,3]=FLB_method(m2,m3)
  d[2,4]=FLB_method(m2,m4)
  d[2,5]=FLB_method(m2,m5)
  d[2,6]=FLB_method(m2,m6)
  d[2,7]=FLB_method(m2,m7)
  d[2,8]=FLB_method(m2,m8)
  
  d[3,1]=FLB_method(m3,m1)
  d[3,2]=FLB_method(m3,m2)
  d[3,3]=FLB_method(m3,m3)
  d[3,4]=FLB_method(m3,m4)
  d[3,5]=FLB_method(m3,m5)
  d[3,6]=FLB_method(m3,m6)
  d[3,7]=FLB_method(m3,m7)
  d[3,8]=FLB_method(m3,m8)
  
  d[4,1]=FLB_method(m4,m1)
  d[4,2]=FLB_method(m4,m2)
  d[4,3]=FLB_method(m4,m3)
  d[4,4]=FLB_method(m4,m4)
  d[4,5]=FLB_method(m4,m5)
  d[4,6]=FLB_method(m4,m6)
  d[4,7]=FLB_method(m4,m7)
  d[4,8]=FLB_method(m4,m8)
  
  d[5,1]=FLB_method(m5,m1)
  d[5,2]=FLB_method(m5,m2)
  d[5,3]=FLB_method(m5,m3)
  d[5,4]=FLB_method(m5,m4)
  d[5,5]=FLB_method(m5,m5)
  d[5,6]=FLB_method(m5,m6)
  d[5,7]=FLB_method(m5,m7)
  d[5,8]=FLB_method(m5,m8)
  
  d[6,1]=FLB_method(m6,m1)
  d[6,2]=FLB_method(m6,m2)
  d[6,3]=FLB_method(m6,m3)
  d[6,4]=FLB_method(m6,m4)
  d[6,5]=FLB_method(m6,m5)
  d[6,6]=FLB_method(m6,m6)
  d[6,7]=FLB_method(m6,m7)
  d[6,8]=FLB_method(m6,m8)
  
  d[7,1]=FLB_method(m7,m1)
  d[7,2]=FLB_method(m7,m2)
  d[7,3]=FLB_method(m7,m3)
  d[7,4]=FLB_method(m7,m4)
  d[7,5]=FLB_method(m7,m5)
  d[7,6]=FLB_method(m7,m6)
  d[7,7]=FLB_method(m7,m7)
  d[7,8]=FLB_method(m7,m8)
  
  d[8,1]=FLB_method(m8,m1)
  d[8,2]=FLB_method(m8,m2)
  d[8,3]=FLB_method(m8,m3)
  d[8,4]=FLB_method(m8,m4)
  d[8,5]=FLB_method(m8,m5)
  d[8,6]=FLB_method(m8,m6)
  d[8,7]=FLB_method(m8,m7)
  d[8,8]=FLB_method(m8,m8)
 
  r = data.frame(d)
  names(r) = names
  
  rownames(r) = names
  
  return(r)
}
#--------------------------------------------------------------
# install.packages("readobj")
library(readobj)
# 
# # memoliModelPath = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/"
# # horse-01.obj
# memoliModelPath = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/horse-poses/"
# model_name = "horse-01"
# 
# horse01 = getMemoliModel3(memoliModelPath,"horse-01")
# horse02 = getMemoliModel3(memoliModelPath,"horse-02")
# horse03 = getMemoliModel3(memoliModelPath,"horse-03")
# horse04 = getMemoliModel3(memoliModelPath,"horse-04")
# 
# head01 = getMemoliModel3(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/head-poses/","head-01-anger")
# head02 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/head-poses/","head-02-cry")
# 
# cat01 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/cat-poses/","cat-01")
# cat02 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/cat-poses/","cat-02")
# cat03 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/cat-poses/","cat-03")
# cat04 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/cat-poses/","cat-04")
# 
# 
# 
# FLB_method(model_X=horse01, model_Y=horse02)
# FLB_method(model_X=horse04, model_Y=horse02)
# 
# FLB_method(model_X=head02, model_Y=horse02)
# FLB_method(model_X=cat01, model_Y=horse01)
# FLB_method(model_X=cat01, model_Y=horse01)
# FLB_method(model_X=cat01, model_Y=horse01)
# 
# d = flb_distances(horse01,horse02,horse03,horse04,cat01,cat02,cat03,head02, c("horse01","horse02","horse03","horse04","cat01","cat02","cat03","head02"))
# d
# 
# heatmap(as.matrix(d))
# 
# 
# getCentroid(cat01$points)
# 
# d_m = dist(cat01$points)
# 
# 
# hc = hclust(d_m)
# n = 500
# cut = cutree(hc, k = n)
# 
# centroids = matrix(0,ncol = 3, nrow = n)
# for(i in 1:n){
#   centroids[i,] = getCentroid(cat01$points[which(cut == i),])
# }
# 
# centroids
# 
# points3d(cat01$points)
# points3d(cat01$centers, col = "blue")
# 
# # install.packages("FNN")
# library(FNN)
# 
# centroids_nn = cat01$points[get.knnx(data = cat01$points, query = centroids, k = 1, algo = "kd_tree")$nn.index,]
# 
# points3d(centroids_nn, col = "red")
# 
# 



#----------------------------------------------------------------------------------------
library(readobj)
# myToyPath = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/000_Trx/000_Trx.obj"
# toy = read.obj(myToyPath, convert.rgl = FALSE)
# 
# toy_rgl = read.obj(myToyPath, convert.rgl = TRUE)
# shade3d(toy_rgl)
# 
# # toy$shapes[[2]]
# 
# toy_points_pos = t(toy$shapes[[2]]$positions)
# points3d(toy_points_pos, col = "blue")
# 
# toy_points_neg = t(toy$shapes[[3]]$positions)
# points3d(toy_points_neg, col = "red")
#----------------------------------------------------------------------------------------

myDownsampleIndices <- function(points, n = 1){
  library(FNN)
  points_sample = seq(1, nrow(points), n)
  points_downsampled2 = points[points_sample,]
  # toy_points_pos_downsampled = points[get.knnx(data = points, query = points_downsampled2, k = 1, algo = "kd_tree")$nn.index,]
  
  sampled_points_knn = get.knnx(data = points, query = points_downsampled2, k = 1, algo = "kd_tree")$nn.index
  
  return(sampled_points_knn)
}

getProtein <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", name, n = 10, n_s_euclidean = 4000, n_s_dijkstra = 50){
 
  if(!dir.exists(paste(path, "/", name,"/PtsAndMu/", sep = ""))) dir.create(paste(path, "/", name,"/PtsAndMu/", sep = ""))
  
  ptsPath_pos = paste(path, "/", name,"/PtsAndMu/",name, "_",n,"_pos.pts",sep ="")
  muPath_pos = paste(path, "/", name,"/PtsAndMu/",name, "_",n,"_pos.mu",sep ="")
  
  ptsPath_neg = paste(path, "/", name,"/PtsAndMu/",name, "_",n,"_neg.pts",sep ="")
  muPath_neg = paste(path, "/", name,"/PtsAndMu/",name, "_",n,"_neg.mu",sep ="")
  objPath = paste(path, "/", name,"/",name, ".obj",sep ="")
  
  # horse_sampled = getMemoliModel2(memoliObjPath = objPath, memoliPtsPath = ptsPath, memoliMuPath = muPath)
  
  # getMemoliModel2 <- function(memoliObjPath, memoliPtsPath, memoliMuPath, n = 1){
  prot = read.obj(objPath, convert.rgl = FALSE)
  prot_rgl = read.obj(objPath, convert.rgl = TRUE)
  shade3d(prot_rgl)
  prot_points_pos = t(prot$shapes[[2]]$positions)
  prot_edges_pos = t(prot$shapes[[2]]$indices)+1
  
  prot_points_neg = t(prot$shapes[[3]]$positions)
  prot_edges_neg = t(prot$shapes[[3]]$indices)+1
  
  
  print(paste("model has ", nrow(prot_points_pos),", ", nrow(prot_points_neg), " points", sep ="" ))
  
  prot_points_pos_sampled_ind = myDownsampleIndices(prot_points_pos, n = n)
  prot_points_neg_sampled_ind = myDownsampleIndices(prot_points_neg, n = n)
  
  print(paste("downsampled model has ", length(prot_points_pos_sampled_ind),", ", length(prot_points_neg_sampled_ind), " points", sep ="" ))
  
  # print("step 1: Setting up graph ...")
  # m = adjacencyMatrixFromMesh(prot_edges_pos, prot_points_pos, prot_points_pos_sampled_ind)
  # 
  # downsampleMemoliMethod(prot_rgl, points = prot_points_pos[prot_points_pos_sampled_ind,], adjacencyMatrix = m)
  
  prot_pos_sampled = getMemoliModel(prot_rgl, points = prot_points_pos, indices = prot_points_pos_sampled_ind, edges = prot_edges_pos, ptsPath_pos, muPath_pos, n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra)
  prot_neg_sampled = getMemoliModel(prot_rgl, points = prot_points_neg, indices = prot_points_neg_sampled_ind, edges = prot_edges_neg, ptsPath_neg, muPath_neg, n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra)
  
  
  plotMemoliModel2(prot_rgl, centers1 = prot_pos_sampled$centers, measure1 = prot_pos_sampled$mu,
                   centers2 = prot_neg_sampled$centers, measure2 = prot_neg_sampled$mu, scale = 500, plotModel = FALSE)
  
  
  l = list("rgl_model" = prot_rgl, "prot_pos_sampled" = prot_pos_sampled, "prot_neg_sampled" = prot_neg_sampled)
  return(l)
}

getProteinsInPath <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", n = 100, n_s_euclidean = 100, n_s_dijkstra = 15){

  proteinNames = list.dirs(path, full.names = FALSE, recursive = FALSE)

  print("------------------------------------------------------------------")
  print(paste("Found the following proteins in ", path, ": "))
  print(proteinNames)
  print("------------------------------------------------------------------")

  p_vec = c()
  for(pName in proteinNames){
    p = getProtein(path = path, pName, n = n, n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra)
    p_vec = c(p_vec,pName,p)
  }

  return(p_vec)
}

unfoldProteinVector <- function(p_vec, ind){
  name = p_vec[1 + (ind-1)*4]
  rgl_model = p_vec[2 + (ind-1)*4]

  model_positive = p_vec[3 + (ind-1)*4]
  model_negative = p_vec[4 + (ind-1)*4]
  
  return(list("name" = name, "rgl_model" = rgl_model, "model_positive" = model_positive, "model_negative" = model_negative))
}

calculateFlbDistancesOfProteins(p_vec)

calculateFlbDistancesOfProteins <- function(p_vec){
  numberOfProteins = length(p_vec)/4
  
  d_pos = matrix(0,ncol = numberOfProteins, nrow = numberOfProteins)
  d_neg = matrix(0,ncol = numberOfProteins, nrow = numberOfProteins)
  
  pNames = rep(" ", numberOfProteins)
  for(i in 1:numberOfProteins){
    for(j in 1:numberOfProteins){
      
      unfolded =  unfoldProteinVector(p_vec,i)
      unfolded2 =  unfoldProteinVector(p_vec,j)

      d_pos[i,j] = FLB_method(unfolded$model_positive$prot_pos_sampled, unfolded2$model_positive$prot_pos_sampled)
      d_neg[i,j] = FLB_method(unfolded$model_negative$prot_neg_sampled, unfolded2$model_negative$prot_neg_sampled)
      
      pNames[i] = unfolded$name
    }
  }
  
  row.names(d_pos) = pNames
  colnames(d_pos) = pNames
  
  row.names(d_neg) = pNames
  colnames(d_neg) = pNames
  
  return(list("pos" = d_pos, "neg" = d_neg))
}

calcFlbAndPlot <- function(p_vec){
  d = calculateFlbDistancesOfProteins(p_vec, path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/")
  
  pdf(paste(path,"positive.pdf", sep = ""))
  heatmap(as.matrix(d$pos), main = "pos")
  dev.off()
  
  pdf(paste(path,"positive.pdf", sep = ""))
  heatmap(as.matrix(d$neg), main = "neg")
  dev.off()
  
}

# p_vec is a list, 4 consecutive elements are one protein
p_vec = getProteinsInPath()
calcFlbAndPlot(p_vec)



#------------------------------------------------------------------------
# more RAM?
# install.packages("openssl")
# install.packages("httr")
# install.packages("devtools")
# devtools::install_github("krlmlr/ulimit")
# 
# ulimit::memory_limit(140000)
#------------------------------------------------------------------------
