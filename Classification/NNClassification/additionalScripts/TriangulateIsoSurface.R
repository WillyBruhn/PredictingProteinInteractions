#!/usr/bin/Rscript

# source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/isoFaces.R")
# source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/extrinsicDistances.R")


is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}

if(!is.installed("misc3d")){install.packages("misc3d")}
library("misc3d")


# tell horst to install GL library
# if(!is.installed("rgl")){install.packages("rgl")}
# library("rgl")


library(misc3d)
# library(rgl)

if(!is.installed("doBy")){install.packages("doBy")}
library("doBy")

if(!is.installed("geometry")){install.packages("geometry")}
library("geometry")

if(!is.installed("FNN")){install.packages("FNN")}
library("FNN")

if(!is.installed("gplots")){install.packages("gplots")}
library("gplots")

if(!is.installed("readobj")){install.packages("readobj")}
library("readobj")
if(!is.installed("xtable")){install.packages("xtable")}
library("xtable")

# if(!is.installed("spatstat")){install.packages("spatstat")}
library("spatstat")

# if(!is.installed("igraph")){install.packages("igraph")}
library("igraph")




# citation(package = "readobj")
# citation(package = "igraph")


# higher value means more is printed
MY_PRINT_SWITCH = 0

printMyMessage <- function(message, val, debug_lvl = 5){
  # debug_lvl ... this message is printed if debug_lvl <= MY_PRINT_SWITCH
  
  if(debug_lvl <= MY_PRINT_SWITCH) print(paste(message, ":",val,sep = " "))
}

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

plotDownsampledPoints <- function(model,points, indices, plotModel = TRUE, size = 20, col = "green"){
  if(plotModel == TRUE){
    shade3d(model)
  } else {
    points3d(points, col = "red")
  }
  points3d(points[indices,], col = col, size  = size)
}

plotMemoliModel <- function(model,points, centers, measure, scale = 100, plotModel = TRUE, plotMu = TRUE){
  # if(plotModel == TRUE) shade3d(model)
  
  if(plotModel == TRUE){
    shade3d(model)
  } else {
    # points3d(points, col = "black", add = FALSE)
  }
  
  if(plotMu == TRUE){
    for(i in 1:nrow(centers)){
      # points3d(centers[i,], col = "red", size = measure[i]*scale)
      plot3d(x = centers[i,1], y = centers[i,2], z = centers[i,3], add = TRUE, col = "green", size = measure[i]*scale, xlab = "", ylab = "", zlab = "")
      
    }
    
  }

}

plotMemoliObj <- function(obj, plotModel = TRUE, scale = 100){
  if(plotModel == TRUE){
    shade3d(obj$rgl_model)
  } else {
    points3d(obj$points, col = "black", add = FALSE)
  }
  
  for(i in 1:nrow(obj$centers)){
    print(obj$centers[i,1])
    plot3d(x = obj$centers[i,1], y = obj$centers[i,2], z = obj$centers[i,3], add = TRUE, col = "green", size = obj$measure[i]*scale, xlab = "", ylab = "", zlab = "")
  }
}


plotMemoliModel2 <- function(model,centers1, measure1, centers2, measure2, scale = 5000, plotModel = TRUE, plotMu = TRUE){
  rgl.open()
  # rgl.bg(color = "white")
  
  par3d(windowRect = c(0, 0, 1000, 1000))
  
  if(plotModel == TRUE) {
    shade3d(model$vmd_mol0_rep2)
    shade3d(model$vmd_mol0_rep3)
  }
  
  if(plotMu == TRUE){
  for(i in 1:nrow(centers1)){
    # points3d(centers[i,], col = "red", size = measure[i]*scale)
    plot3d(x = centers1[i,1], y = centers1[i,2], z = centers1[i,3], add = TRUE, col = "green", size = measure1[i]*scale, xlab = "", ylab = "", zlab = "")
    
  }
  
  for(i in 1:nrow(centers2)){
    # points3d(centers[i,], col = "red", size = measure[i]*scale)
    plot3d(x = centers2[i,1], y = centers2[i,2], z = centers2[i,3], add = TRUE, col = "yellow", size = measure2[i]*scale, xlab = "", ylab = "", zlab = "")
    
  }
  }
}

myVoronoi <- function(centers,points){
  measureCounts = rep(0,nrow(centers))
  
  # d = rep(-1,nrow(centers))
  # for(i in 1:nrow(points)){
  #   for(j in 1:nrow(centers)){
  #     d[j] = dist(rbind(points[i,], centers[j,]))
  #   }
  #   nearestCenter = which(d == min(d))
  #   
  #   measureCounts[nearestCenter] = measureCounts[nearestCenter] + 1 
  # }
  
  indices = knnx.index(data = centers, query = points, k = 1)
  
  for(c in 1:length(measureCounts)){
    measureCounts[c] = length(which(indices == c))
  }
  
  return(measureCounts)
}

# knnx.index(data = centers, query = points, k = 1)
# 
# ?aggregate()


downsampleMemoliMethod <- function(model_rgl,points, edges, graph, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE){
  
  library(rdist)
  
  print("step 1: euclidean fps ...")
  
  sampled_indices = myFarthestPointSampling(points, k = n_s_euclidean)
  
  if(plot) plotDownsampledPoints(model_rgl,points,sampled_indices, plotModel = FALSE)
  
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
  
  l = list("centers" = points[sampled_indices[sampled_indices2],], "mu" = v_n, "indices_order" = sampled_indices[sampled_indices2], "d_surface" = d_surface[sampled_indices2,sampled_indices2], "sampled_indices_large" = sampled_indices)
  return(l)
}

# rgl.open()
# downsampleMemoliMethod(model_rgl = prot_rgl$vmd_mol0_rep2, points = ob$points, edges = ob$edges, graph = ob$graph,plot = TRUE, n_s_euclidean = 20, n_s_dijkstra = 6)



downsampleMemoliMethodWithoutSurface <- function(model_rgl, points, n_s_euclidean_1 = 4000, n_s_euclidean_2 = 50, plot = TRUE){
  print("step 1: euclidean fps ...")
  sampled_indices_1 = myFarthestPointSampling(points, k = n_s_euclidean_1)
  
  print(paste("downsampled model has ", length(sampled_indices_1), " points", sep = ""))
  
  # print(paste("further downsampled model has ", length(sampled_indices_2), " points", sep = ""))
  
  if(plot) plotDownsampledPoints(model_rgl,points,sampled_indices_1,TRUE)
  
  print("calculating voronoi ...")
  v2 = myVoronoi(points[sampled_indices_1,], points)
  
  print("... done")
  
  v_n = v2/sum(v2)
  
  # v_n = v2
  
  l = list("centers" = points[sampled_indices_1,], "mu" = v_n)
  return(l)
}



# 
# ob = p_vec = getProteinsInPath(n_s_euclidean = 100, n_s_dijkstra = 50)
# 
# vertices = c(1,2,3000)
# myShortestDistances(ob$graph,vertices)
# distances(ob$graph,v = vertices, to = vertices)



myShortestDistances <- function(g,vertices){
  d = matrix(0,ncol = length(vertices), nrow = length(vertices))
  # print(d)
  
  for(i in 1:length(vertices)){
    print(paste(100 * i / length(vertices), " %", sep = ""))
    
    d_i = distances(graph = g,v = vertices[i], to = vertices, algorithm = "dijkstra")

    d[i,]  = d_i
    
    # d_i = distances(graph = g,v = vertices[i], to = vertices[-c(1:i)], algorithm = "dijkstra")
    # 
    # d[i,]  = d_i
    # d[,i] = d_i
  }

  return(d)
}

# 
# # example graph
# set.seed(1)
# g <- erdos.renyi.game(20, 1/20)
# V(g)$name <- letters[1:20]
# par(mar=rep(0,4))
# plot(g)
# 
# cl = clusters(g)
# 
# cl$membership 
# wm = which.max(cl$csize)
# 
# which(cl$membership ==wm)
# 
# g[wm,]
# 
# cl$no


getBiggestConnectedComponent <- function(g){
  cl = clusters(g)
  # print(paste("no: ", cl$no))
  
  wm = which.max(cl$csize)
  return(which(cl$membership ==wm))
}


# getBiggestConnectedComponent(g)

# 
# getProtein(pathname = "105", n_s_euclidean = 100)
# p = getProtein(name = "000_Trx", n_s_euclidean = 10, n_s_dijkstra = 3, plot = TRUE)
# 
# p$prot_pos_sampled$centers
# 
# myVoronoi(p$prot_pos_sampled$centers, prot_points_pos)
# 
# points3d(prot_points_pos, col = "blue")
# points3d(p$prot_pos_sampled$centers, col = "red", size = 10)
# 
# 
# 
# 

# prot = read.obj("/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/000_Trx/000_Trx.obj", convert.rgl = FALSE)
# prot_rgl = read.obj("/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/000_Trx/000_Trx.obj", convert.rgl = TRUE)
# 
# # prot_rgl$vmd_mol0_rep2
# # 
# prot_points_pos = t(prot$shapes[[2]]$positions)
# prot_edges_pos = t(prot$shapes[[2]]$indices)+1
# 
# prot_points_neg = t(prot$shapes[[3]]$positions)
# prot_edges_neg = t(prot$shapes[[3]]$indices)+1
# # 
# # prot_edges_pos[which.max(prot_edges_pos),]
# # 
# # prot_edges_pos[which.min(prot_edges_pos),]
# # 
# # prot_edges_neg[1:5,]
# # 
# # 
# # shade3d(prot_rgl$vmd_mol0_rep2)
# # shade3d(prot_rgl$vmd_mol0_rep3)
# # 
# # points3d(prot_points_neg, col = "blue")
# # points3d(prot_points_pos, col = "red")
# # 
# ob = preProcessMesh(prot_points_pos, prot_edges_pos, plot = TRUE)


# 
# 
# downsampleMemoliMethod(model_rgl = prot_rgl$vmd_mol0_rep2, points = ob$points, edges = ob$edges, graph = ob$graph,plot = TRUE, n_s_euclidean = 100, n_s_dijkstra = 20)
# 
# 



# ob$graph
# myDrawShortestPath(ob,1,294)
# 1   2  23  61 133 294 303 156 155  72  35  10

# ob$graph[133,294]
# rgl.snapshot(filename = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/shortestPath", fmt = "png", top = TRUE )


edgesOfPoint <- function(p, edges){
  return(rbind(edges[which(edges[,1] == p),],
               edges[which(edges[,2] == p),],
               edges[which(edges[,3] == p),]))
}

edgesOfPoint2 <- function(p, edges){
  
  ret = c()
  ind = which(edges == p)
  for(i in ind){
    if(i %% 2 == 0){
      ret = c(ret, i-1)
    }
    else {
      ret = c(ret, i+1)
    }
  }
  
  # print(edges[min(ind):max(ind)])
  
  return(edges[ret])
}


# edgesOfPoint(133, prot_edges_pos)
# 
# edgesOfPoint(294, prot_edges_pos)


myDrawShortestPath <- function(ob, start, end, rgl = c(""), col = "blue"){
  p = shortest_paths(ob$graph, start, end)$vpath[[1]]
  
  print(p)
  
  if(rgl !=c("")){
    shade3d(rgl)
  } else {
    points3d(ob$points[], col = "red", size = 10)
  }

  
  lines3d(ob$points[p,], col = col, lwd = 10)
  points3d(ob$points[p,], col = col, size = 20)
}

testGeodesic <- function(ob, start,end, eps = 0.00001){

  is = distances(ob$graph,start,end)[1]
  should = myCalcLengthOfPath(ob$points,shortest_paths(ob$graph, start, end)$vpath[[1]])[[1]]
  
  if(abs(is - should) > eps){
    print(paste("is: ", is, ", should: ", should), sep ="")
    return(FALSE)
  } else {
    return(TRUE)
  }

}

tetGeodesicMultiple <- function(ob, n = 100){
  for(i in 1:n){
    print(i/n*100)
    start = sample(c(1:nrow(ob$points)),replace = TRUE,size = 1)
    end = sample(c(1:nrow(ob$points)),replace = TRUE,size = 1)
    if(FALSE == testGeodesic(ob,start,end)) return(FALSE)
  }
  return(TRUE)
}

# tetGeodesicMultiple(ob,n=1000)
# 
# # 
# # testGeodesic(ob,1,3)
# # testGeodesic(ob,1,2277)
# # testGeodesic(ob,1885,294)
# 
# # 1    3   28   59  115  260  259  460  672  659  658  923 1196 1493 1886 1885 2277
# testGeodesic(ob,1885,2277)
# 
# 
# 
# myDrawShortestPath(ob,2000,10000, col = "green")
# rgl.snapshot(filename = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/TrxShortestPath", fmt = "png", top = TRUE )


# l = myMakeUniquePointsAndEdges(prot_points_pos,prot_edges_pos)
# ob = makeMyEdges(l$points,l$edges)
# 
# edgesOfPoint2(1,ob$edges)
# 
# # 1   2  23  61 133 294
# edgesOfPoint(1,l$edges)
# edgesOfPoint(2,l$edges)
# edgesOfPoint(23,l$edges)
# edgesOfPoint(61,l$edges)
# edgesOfPoint(133,l$edges)
# edgesOfPoint(294,l$edges)
# 
# edgesOfPoint2(1,ob$edges)
# edgesOfPoint2(2,ob$edges)
# edgesOfPoint2(23,ob$edges)
# edgesOfPoint2(61,ob$edges)
# edgesOfPoint2(133,ob$edges)
# edgesOfPoint2(294,ob$edges)


myCalcLengthOfPath <- function(points, path){
  s = 0
  for(i in 1:(length(path)-1)){
    # print(path[i])
    s = s + euklid_dist(points[path[i],], points[path[i+1],])
  }
  return(s)
}


# points3d(ob$points[c(1,9000),], col = "blue", size = 100)

# 
# preProcessMesh(prot_points_neg, prot_edges_neg, plot = TRUE)


makeMyEdges <- function(points,edges){
  edges_out = rep(0,nrow(edges)*2*3)
  
  weights = rep(0,nrow(edges))
  
  # edges2 = matrix(0,nrow())
  
  for(i in 1:nrow(edges)){
    edges_out[(i-1)*6+1] = edges[i,1]
    edges_out[(i-1)*6+2] = edges[i,2]
    
    edges_out[(i-1)*6+3] = edges[i,1]
    edges_out[(i-1)*6+4] = edges[i,3]
    
    edges_out[(i-1)*6+5] = edges[i,2]
    edges_out[(i-1)*6+6] = edges[i,3]
    
    
    
    weights[(i-1)*3+1] = euklid_dist(points[edges[i,1],],points[edges[i,2],])

    weights[(i-1)*3+2] = euklid_dist(points[edges[i,1],],points[edges[i,3],])

    weights[(i-1)*3+3] = euklid_dist(points[edges[i,2],],points[edges[i,3],])
  }
  
  return(list("edges" = edges_out, "weights" = weights))
}


preProcessMesh <- function(points, edges, plot = FALSE){
  print("downsampling to unique points ...")
  l = myMakeUniquePointsAndEdges(points,edges)
  
  print("creating edges and edge-weights ...")
  ob = makeMyEdges(l$points,l$edges)
  
  print("creating graph ... ")
  g <- make_graph(ob$edges, directed = FALSE)
  
  g = set_edge_attr(g, "weight", index=E(g), ob$weights)
  
  weights2 = get.edge.attribute(g,"weight",index = E(g))

  multiples = count_multiple(g)
  for(i in 1:length(multiples)){
    # if(multiples[i] == TRUE){
    #   weights2[i] = weights2[i] / 2
    # }
    
    weights2[i] = weights2[i]/multiples[i]
  }
  
  g = set_edge_attr(g, "weight", index=E(g), weights2)
  
  g = simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  numOfConComps = clusters(g)$no
  bg = getBiggestConnectedComponent(g)
  
  if(plot){
    points3d(l$points, col = "red")
    points3d(l$points[bg,], col = "green", size = 10)
  }
  
  points = l$points[bg,]
  # points = l$points
  
  g = induced_subgraph(g,bg)
  # ob = makeMyEdges(points,l$edges)

  return(list("points" = points, "edges" = ob$edges, "graph" = g, "numOfConComps" = numOfConComps))
}

# 
# n = 5
# dim = 2
# m = matrix(rnorm(dim*n), ncol = dim)
# 
# ed = c(1,3, 3,4, 4,5)
# d = as.matrix(dist(m, upper = TRUE))
# 
# weights = c(d[1,3], d[3,4], d[5,4])
# 
# 
# g = make_graph(ed, directed = FALSE)
# 
# g = set_edge_attr(g, "weight", index = E(g), weights)
# 
# 
# distances(g, 1, 5)
# myCalcLengthOfPath(m,shortest_paths(g, 1, 5)$vpath[[1]])
# 
# plot(g)
# plot(m)
# text(x = m[,1]+0.05, y = m[,2], c(1:5))


myMakeUniquePointsAndEdges <- function(points,edges){
  # Points are appearing multiple times with the same coordinates
  # we have to reduce it to have one connected mesh instead of
  # many small ones
  
  points_unique = unique(points)
  indices = knnx.index(points_unique,query = points, k = 1)
  
  edges_correct_indices = edges

  for(i in 1:nrow(edges)){
    edges_correct_indices[i,1] = indices[edges[i,1]]
    edges_correct_indices[i,2] = indices[edges[i,2]]
    edges_correct_indices[i,3] = indices[edges[i,3]]
  }

  return(list("points" = points_unique, "edges" = edges_correct_indices))
}



checkForLargerModel <- function(path, proteinName, n_s_euclidean, n_s_dijkstra, positive = TRUE){
  # first try to reuse the points sampled with larger parameters, that means n_s_euclidean has to be same
  # but n_s_dijkstra has to be larger in another model. In that case we can use only the first few points
  # of that model.
  # Uses the largest model that can be found.
  
  dir = paste(path, "/", proteinName, "/PtsAndMu/", sep ="")
  s = list.files(path = dir,pattern = "*.order")
  
  dijkstraVals = rep(-1,length(s))
  
  printMyMessage("s",s)
  
  posNegExt = "pos"
  if(positive == FALSE) posNegExt = "neg"
  
  files = strsplit(s, split = ".order")
  for(i in 1:length(files)){
    if(length(files) == 0) return(FALSE)
    
    # print(files[[i]])
    # s3 = strsplit(files[[i]], split = proteinName) # bad for names that are also numbers
    
    s4 = substr(files[[i]], nchar(proteinName)+1,nchar(files[[i]]))
    printMyMessage("s4",s4)
    
    s5 = strsplit(s4, split = "_")
    
    printMyMessage("s5[[1]][1]",s5[[1]][1])
    printMyMessage("s5[[1]][1]",s5[[1]][2])
    printMyMessage("s5[[1]][1]",s5[[1]][3])

    if(as.numeric(s5[[1]][1]) == n_s_euclidean && as.numeric(s5[[1]][2]) >= n_s_dijkstra && s5[[1]][3] == posNegExt){
      # found a file
      # return(paste(dir,files[i][[1]],".order", sep =""))
      
      dijkstraVals[i] = as.numeric(s5[[1]][2])
    }
  }
  
  if(max(dijkstraVals != -1)) return(paste(dir,files[which.max(dijkstraVals)][[1]],".order", sep =""))
  
  return(FALSE)
}

# p1 = getProtein(name = "060",n_s_euclidean = 4000, n_s_dijkstra = 4000)
# p2 = getProtein(name = "060",n_s_euclidean = 4000, n_s_dijkstra = 60)
# 
# 
# all.equal(p1$prot_pos_sampled$centers[1:60,]
# ,p2$prot_pos_sampled$centers[1:60,])
# 
# all.equal(p1$prot_neg_sampled$centers[1:60,]
#           ,p2$prot_neg_sampled$centers[1:60,])
# checkForLargerModel(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", proteinName = "060",n_s_euclidean = 4000, n_s_dijkstra = 60, positive = FALSE)

# MY_PRINT_SWITCH = 5
# checkForLargerModel(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", proteinName = "099",n_s_euclidean = 4000, n_s_dijkstra = 100)



#

# s = checkForLargerModel("/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", "000_Trx", 4000,500)
# 
# q = read.table(s, header = TRUE)
# q$pointInd


getModelFromLargerModel <- function(orderFile, points, n_s_dijkstra){
  
  indices = read.table(orderFile, header = TRUE)$pointInd[1:n_s_dijkstra]
  
  v2 = myVoronoi(points[indices,], points)
  
  v_n = v2/sum(v2)
  
  l = list("centers" = points[indices,], "mu" = v_n, "indices_order" = indices)
}

getMemoliModel <- function(model_rgl,points, path = "", name="", edges, graph, memoliPtsPath, memoliMuPath, memoliOrderPath, n_s_euclidean = 4000, n_s_dijkstra = 50, reCalculate = FALSE){
  if(!file.exists(memoliPtsPath) || !file.exists(memoliMuPath) || reCalculate == TRUE){
    
    # first try to reuse the points sampled with larger parameters, that means n_s_euclidean has to be same
    # but n_s_dijkstra has to be larger in another model. In that case we can use only the first few points
    # of that model
    
    # largerModelFile = c()
    # if(length(strfind(memoliOrderPath, "pos")) > 0){
    #   largerModelFile = checkForLargerModel(path, name,n_s_euclidean, n_s_dijkstra, positive = TRUE)
    # } else {
    #   largerModelFile = checkForLargerModel(path, name,n_s_euclidean, n_s_dijkstra, positive = FALSE)
    # }

    mem = c()
    # if(largerModelFile == FALSE){
      mem = downsampleMemoliMethod(model_rgl, points = points, edges = edges, graph = graph,  n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra)
    # } else {
    #   mem = getModelFromLargerModel(largerModelFile,points,n_s_dijkstra)
    # }
    
    
    print("writing model to file ...")
    write.table(mem$centers,file=memoliPtsPath, row.names = FALSE, col.names = c("x","y","z"))
    write.table(mem$mu,file=memoliMuPath, row.names = FALSE, col.names = c("mu"))
    
    write.table(mem$indices_order, file = memoliOrderPath, row.names = FALSE, col.names = c("pointInd"))
    
    # write.table(mem$d_surface, file = memoliDSurfacePath)
  }
  
  pts = read.table(file = memoliPtsPath, header = TRUE, colClasses = c("numeric", "numeric", "numeric"), stringsAsFactors=FALSE)
  mu2 = unlist(read.table(file = memoliMuPath, header = TRUE, colClasses = c("numeric"), stringsAsFactors=FALSE))


  mu = as.numeric(mu2)
    
  l = list("centers" = pts, "mu" = mu, "geoDistances" = mem$d_surface)
  return(l)
}

getMemoliModel2 <- function(memoliObjPath, memoliPtsPath, memoliMuPath, memoliOrderPath, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE){
  horse1 = read.obj(memoliObjPath, convert.rgl = FALSE)
  horse1_rgl = read.obj(memoliObjPath, convert.rgl = TRUE)
  shade3d(horse1_rgl)
  horse1_points = t(horse1$shapes[[1]]$positions)
  horse1_edges = t(horse1$shapes[[1]]$indices)+1
  
  
  ob = preProcessMesh(points = horse1_points, edges = horse1_edges, plot = FALSE)
  print(paste("original model has: ", nrow(horse1_points), ",processed model has ", nrow(ob$points)))

  horse_sampled = getMemoliModel(model_rgl = prot_rgl, points = ob$points, edges = ob$edges, graph = ob$graph, memoliPtsPath = memoliPtsPath, memoliMuPath = memoliMuPath, memoliOrderPath = memoliOrderPath,n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra)
  
  # getMemoliModel()
  
  # horse_sampled = getMemoliModel(horse1_rgl,points = horse1_points, edges = horse1_edges, memoliPtsPath, memoliMuPath)
  
  if(plot) plotMemoliModel(horse1_rgl,horse_sampled$centers, measure = horse_sampled$mu, scale = 500, plotModel = FALSE)
  
  l = list("rgl_model" = horse1_rgl, "points" = horse1_points, "centers" = horse_sampled$centers, "mu" = horse_sampled$mu)
  return(l)
}

getMemoliModel3 <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/horse-poses/", name, n_s_euclidean = 4000, n_s_dijkstra = 50){
  memoliModelPath = path
  model_name = name
  
  ext = paste("_", n_s_euclidean, "_", n_s_dijkstra, "_")
  
  if(!dir.exists(paste(path,"/PtsAndMu",sep=""))) dir.create(paste(path,"/PtsAndMu",sep=""))
  memoliPtsPath = paste(memoliModelPath, "/PtsAndMu/", model_name, ext,  ".pts",sep ="")
  memoliMuPath = paste(memoliModelPath, "/PtsAndMu/", model_name, ext, ".mu",sep = "")
  
  
  memoliOrderPath = paste(memoliModelPath, "/PtsAndMu/", model_name, ext, ".order",sep = "")
  memoliObjPath = paste(path, model_name, ".obj", sep ="")
  
  
  horse_sampled = getMemoliModel2(memoliObjPath = memoliObjPath, memoliPtsPath = memoliPtsPath, memoliMuPath = memoliMuPath, memoliOrderPath = memoliOrderPath, n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra)
  
  return(horse_sampled)
}

# cat01 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/cat-poses/","cat-01")

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

FLB_deprecated <- function(model_X,model_Y,s_X,s_Y,E,L = 10){
  # TODO: do like FLB_2d with distribution
  
  # print(seq(0,max(E)+0.0001,1/(L-1)))
  
  # print(min(E))
  # print(max(E))
  
  b = seq(min(E),max(E),(max(E)-min(E))/(L))
  h = hist(E, breaks = b, plot = FALSE)
  u_vals = h$breaks
  
  # h = hist(E, breaks = seq(0,max(E)+1,1/(L-1)), plot = FALSE)
  # u_vals = h$breaks
  
  
  S_X = rep(0,length(u_vals)-1)
  S_Y = rep(0,length(u_vals)-1)
  for(u in 1:(length(u_vals)-1)){
    B_X_u = which(s_X <= u_vals[u])
    
    s = 0
    if(length(B_X_u) > 0) s = sum(model_X$mu[B_X_u[1]])
    S_X[u] = s
    
    s = 0
    B_Y_u = which(s_Y <= u_vals[u])
    
    if(length(B_Y_u) > 0) s = sum(model_Y$mu[B_Y_u[1]])
    S_Y[u] = s
  } 
  
  # print("s_X")
  # print(s_X)
  # 
  # print("S_X")
  # print(S_X)
  
  flb = 0
  for(i in 1:(length(u_vals)-1)){
    u_1 = u_vals[i]
    u_2 = u_vals[i+1]
    
    flb = flb + abs(u_2 - u_1) * abs(S_X[i] - S_Y[i]) 
  }
  flb = flb*0.5
  
  return(flb)
}


FLB <- function(model_X,model_Y,s_X,s_Y,E,L){
  # model_X ... combined model of pos and neg
  # model_Y ... combined model of pos and neg
  # s_X ... is $s_X_pos or $s_X_neg of all points
  # s_Y ... is $s_X_pos or $s_X_neg of all points
  # E ... occuring vals in s_X$s_X_pos or s_Y$s_X_pos
  
  b = seq(min(E),max(E),(max(E)-min(E))/length(unique(E)))
  h = hist(E, breaks = b, plot = FALSE)
  u_vals = h$breaks
  
  # print(model_X$mu)
  
  library(spatstat)
  f = ewcdf(x = s_X, weights = model_X$mu, normalise = FALSE)
  f2 = ewcdf(x = s_Y, weights = model_Y$mu, normalise = FALSE)
  
  f1_vals = f(u_vals)
  f2_vals = f2(u_vals)
  
  S_diff = matrix(0,nrow = length(u_vals)-1, ncol = length(u_vals)-1)
  flb =0
  
  flb_vals = rep(0,nrow(S_diff)*ncol(S_diff))
  for(i in 1:(length(u_vals)-1)){
    u_1 = u_vals[i]
    u_2 = u_vals[i+1]
    
    # it is sorted so we don't need abs
    # u_p_d = abs(u_2_pos - u_1_pos)
    u_p_d = u_2 - u_1
    
    flb = flb + u_p_d * abs(f1_vals[i]-f2_vals[i])
  }
  
  flb = flb*0.5
  return(flb)
}

#---------------------------
# p = getProteinsInPath(subset_names = c("000_Trx", "016", "034"),n_s_euclidean = 4000,n_s_dijkstra = 5)
# calculateFlbDistancesOfProteins(p_vec = p,sumMethod = TRUE)



#---------------------------

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


FLB_protein_sum_method_deprecated <- function(model_X_pos, model_X_neg,model_Y_pos, model_Y_neg, L = 10){
  
  # positive
  s_X_pos = get_s(model_X_pos$prot_pos_sampled)
  s_Y_pos = get_s(model_Y_pos$prot_pos_sampled)
  
  E_X_pos = unique(s_X_pos)
  E_Y_pos = unique(s_Y_pos)
  E_pos = sort(unique(c(E_X_pos,E_Y_pos)))
  
  # negative
  s_X_neg = get_s(model_X_neg$prot_neg_sampled)
  s_Y_neg = get_s(model_Y_neg$prot_neg_sampled)
  
  E_X_neg = unique(s_X_neg)
  E_Y_neg = unique(s_Y_neg)
  E_neg = sort(unique(c(E_X_neg,E_Y_neg)))
  
  # pos and neg
  model_X_pos_and_neg = list("centers", "mu")
  model_X_pos_and_neg$centers = rbind(model_X_pos$prot_pos_sampled$centers, model_X_neg$prot_neg_sampled$centers)
  model_X_pos_and_neg$mu = rbind(model_X_pos$prot_pos_sampled$mu, model_X_neg$prot_neg_sampled$mu)
  
  model_Y_pos_and_neg = list("centers", "mu")
  model_Y_pos_and_neg$centers = rbind(model_Y_pos$prot_pos_sampled$centers, model_Y_neg$prot_neg_sampled$centers)
  model_Y_pos_and_neg$mu = rbind(model_Y_pos$prot_pos_sampled$mu, model_Y_neg$prot_neg_sampled$mu)
  
  s_X_pos_and_neg = get_s(model_X_pos_and_neg)
  s_Y_pos_and_neg = get_s(model_Y_pos_and_neg)
  
  E_X_pos_and_neg = unique(s_X_pos_and_neg)
  E_Y_pos_and_neg = unique(s_Y_pos_and_neg)
  E_pos_and_neg = sort(unique(c(E_X_pos_and_neg,E_Y_pos_and_neg)))
  
  
  f =   FLB(model_X_pos$prot_pos_sampled,model_Y_pos$prot_pos_sampled,s_X_pos,s_Y_pos,E_pos,L)
  + FLB(model_X_neg$prot_neg_sampled,model_Y_neg$prot_neg_sampled,s_X_neg,s_Y_neg,E_neg,L)
  + FLB(model_X_pos_and_neg,model_Y_pos_and_neg,s_X_pos_and_neg,s_Y_pos_and_neg,E_pos_and_neg,L)
  
  return(f)
}

get_s_2d <-function(model_X){
  
  n_X = nrow(model_X$centers)
  
  middle = (n_X/2)
  
  d_X = as.matrix(dist(model_X$centers))
  
  s_X_pos = rep(0,n_X)
  s_X_neg = rep(0,n_X)
  
  for(i in 1:n_X){
    s_pos = 0
    for(k in 1:middle){
      s_pos = s_pos + model_X$mu[k] * d_X[i,k]
    }
    s_X_pos[i] = s_pos
    
    s_neg = 0
    for(k in (middle+1):n_X){
      s_neg = s_neg + model_X$mu[k] * d_X[i,k]
    }
    s_X_neg[i] = s_neg
  }
  
  return(list("s_X_pos" = s_X_pos, "s_X_neg" = s_X_neg))
}


# E = c(1,3,4,5,6,6,9,10)
# L = 4
# b = seq(min(E),max(E),(max(E)-min(E))/(L))
# v = hist(E,breaks= b)
# v$breaks

findNext <- function(vals,bound, measure, ind_old = 0,sum_old = 0){
  for(i in (ind_old+1):length(vals)){
    if(bound < vals[i]){
      return(list("sum" = sum_old, "ind" = i-1))
    }
    sum_old = sum_old + measure[i]
  }
  
  return(list("sum" = sum_old, "ind" = length(vals)))
}

# v=c(1:10)
# measure = rep(0.1,10)
# q = findNext(v, 4,measure)
# findNext(v,5,measure,q$ind,q$sum)
# 
# findNext(v,10,measure,q$ind,q$sum)


get_S_2d <- function(model, s,u_vals_pos, u_vals_neg){
  S = matrix(0,nrow = length(u_vals_pos)-1, ncol = length(u_vals_neg)-1)
  
  for(u_pos in 1:(length(u_vals_pos)-1)){
    # s_pos = 0

    for(u_neg in 1:(length(u_vals_pos)-1)){
      B_u_pos = which(s$s_X_pos <= u_vals_pos[u_pos])
      B_u_neg = which(s$s_X_neg <= u_vals_neg[u_neg])
      
      s1 = 0
      if(length(B_u_pos) != 0) s1 = sum(model$mu[B_u_pos[1]])
      s2 = 0
      if(length(B_u_neg) != 0) s2 = sum(model$mu[B_u_neg[1]])
      S[u_pos,u_neg] = s1 + s2
    }
  }
  return(S)
}

# get_S_2d_fast <- function(model, s,model2, s2,u_vals_pos, u_vals_neg){
#   S = matrix(0,nrow = length(u_vals_pos)-1, ncol = length(u_vals_neg)-1)
#   S2 = matrix(0,nrow = length(u_vals_pos)-1, ncol = length(u_vals_neg)-1)
#   
#   library(spatstat)
#   f_pos = ewcdf(x = s$s_X_pos, weights = model$mu, normalise = FALSE)
#   f_neg = ewcdf(x = s$s_X_neg, weights = model$mu, normalise = FALSE)
#   
#   f_pos2 = ewcdf(x = s2$s_X_pos, weights = model2$mu, normalise = FALSE)
#   f_neg2 = ewcdf(x = s2$s_X_neg, weights = model2$mu, normalise = FALSE)
#   
#   
#   for(u_pos in 1:(length(u_vals_pos)-1)){
#     for(u_neg in 1:(length(u_vals_neg)-1)){
#       
#       S[u_pos,u_neg] = f_pos(u_vals_pos[u_pos]) + f_neg(u_vals_neg[u_neg])
#       S2[u_pos,u_neg] = f_pos2(u_vals_pos[u_pos]) + f_neg2(u_vals_neg[u_neg])
#     }
#   }
#   
#   
#   # S1 = apply(u1,c(1,2), 
#   #            function(x) {
#   #              f(x)
#   #            }
#   # )
#   # 
#   # S2 = S1 + apply(u2,c(1,2), 
#   #                 function(x) {
#   #                   f2(x)
#   #                 }
#   # )
#   
#   return(list("S_X" = S, "S_Y" = S2))
# }

get_S_2d_fast <- function(model, s,model2, s2,u_vals_pos, u_vals_neg){
  S = matrix(0,nrow = length(u_vals_pos)-1, ncol = length(u_vals_neg)-1)
  S2 = matrix(0,nrow = length(u_vals_pos)-1, ncol = length(u_vals_neg)-1)
  
  library(spatstat)
  f_pos = ewcdf(x = s$s_X_pos, weights = model$mu, normalise = FALSE)
  f_neg = ewcdf(x = s$s_X_neg, weights = model$mu, normalise = FALSE)
  
  f_pos2 = ewcdf(x = s2$s_X_pos, weights = model2$mu, normalise = FALSE)
  f_neg2 = ewcdf(x = s2$s_X_neg, weights = model2$mu, normalise = FALSE)
  
  S_diff = matrix(0,nrow = length(u_vals_pos)-1, ncol = length(u_vals_neg)-1)
  
  for(u_pos in 1:(length(u_vals_pos)-1)){
    for(u_neg in 1:(length(u_vals_neg)-1)){
      
      S[u_pos,u_neg] = f_pos(u_vals_pos[u_pos]) + f_neg(u_vals_neg[u_neg])
      S2[u_pos,u_neg] = f_pos2(u_vals_pos[u_pos]) + f_neg2(u_vals_neg[u_neg])
      
      S_diff[u_pos,u_neg]  = abs(S[u_pos,u_neg]-S2[u_pos,u_neg] )
    }
  }
  
  
  # S1 = apply(u1,c(1,2), 
  #            function(x) {
  #              f(x)
  #            }
  # )
  # 
  # S2 = S1 + apply(u2,c(1,2), 
  #                 function(x) {
  #                   f2(x)
  #                 }
  # )
  
  return(list("S_X" = S, "S_Y" = S2, "S_diff" = S_diff))
}


# s = c(1.155726, 1.168095, 1.143569, 1.145343, 1.352365, 1.348140, 1.171360, 1.301305)
# sUn = sort(unique(s))
# 
# mu = c(0.3023864, 0.2329294 ,0.2384817, 0.2262026, 0.2989116, 0.2681929, 0.2278942, 0.2050013)
# f = ecdf(s)
# f(sUn)
# plot(f)

# install.packages("spatstat")

# f = ewcdf(x = s, weights = mu, normalise = FALSE)
# f2 = ewcdf(x = s, weights = mu, normalise = FALSE)
# plot(f)
# 
# S = matrix(0,2,2)
# f(c(1.25,20.0))
# 
# l = c(1.25,20.0)
# u1 = matrix(c(l,l),byrow = TRUE,ncol = 2,nrow = 2)
# u2 = matrix(c(l,l),byrow = TRUE,ncol = 2,nrow = 2)
# 
# u[1,]
# f(u[1,])
# 
# S[1,] = t(as.matrix(f(u[1,])))
# 
# ?apply
# 
# ind = matrix()

# S1 = apply(u1,c(1,2), 
#   function(x) {
#     f(x)
#   }
# )
# 
# S2 = S1 + apply(u2,c(1,2), 
#            function(x) {
#              f2(x)
#            }
# )




FLB_2d <- function(model_X,model_Y,s_X,s_Y,E_pos, E_neg,L){
  # model_X ... combined model of pos and neg
  # model_Y ... combined model of pos and neg
  # s_X ... has $s_X_pos and $s_X_neg of all points
  # s_Y ... has $s_X_pos and $s_X_neg of all points
  # E_pos ... occuring vals in s_X$s_X_pos and s_Y$s_X_pos
  # E_neg ... occuring vals in s_X$s_X_neg and s_Y$s_X_neg
  
  b_pos = seq(min(E_pos),max(E_pos),(max(E_pos)-min(E_pos))/length(unique(E_pos)))
  h_pos = hist(E_pos, breaks = b_pos, plot = FALSE)
  u_vals_pos = h_pos$breaks
  
  b_neg = seq(min(E_neg),max(E_neg),(max(E_neg)-min(E_neg))/length(unique(E_neg)))
  h_neg = hist(E_neg, b_neg, plot = FALSE)
  u_vals_neg = h_neg$breaks
  
  library(spatstat)
  f_pos = ewcdf(x = s_X$s_X_pos, weights = model_X$mu, normalise = FALSE)
  f_neg = ewcdf(x = s_X$s_X_neg, weights = model_X$mu, normalise = FALSE)
  
  f_pos2 = ewcdf(x = s_Y$s_X_pos, weights = model_Y$mu, normalise = FALSE)
  f_neg2 = ewcdf(x = s_Y$s_X_neg, weights = model_Y$mu, normalise = FALSE)
  
  
  f1_vals_pos = f_pos(u_vals_pos)
  f1_vals_neg = f_neg(u_vals_neg)
  f2_vals_pos = f_pos2(u_vals_pos)
  f2_vals_neg = f_neg2(u_vals_neg)
  
  S_diff = matrix(0,nrow = length(u_vals_pos)-1, ncol = length(u_vals_neg)-1)
  flb =0
  
  flb_vals = rep(0,nrow(S_diff)*ncol(S_diff))
  for(i in 1:(length(u_vals_pos)-1)){
    u_1_pos = u_vals_pos[i]
    u_2_pos = u_vals_pos[i+1]
    
    # it is sorted so we don't need abs
    # u_p_d = abs(u_2_pos - u_1_pos)
    u_p_d = u_2_pos - u_1_pos
    for(j in 1:(length(u_vals_neg)-1)){
      u_1_neg = u_vals_neg[j]
      u_2_neg = u_vals_neg[j+1]

      flb = flb + u_p_d * (u_2_neg - u_1_neg) * abs((f1_vals_pos[i] + f1_vals_neg[j])-(f2_vals_pos[i]+f2_vals_neg[j]))
      }
  }
  
  flb = flb*0.5
  return(flb)
}


# p_vec = getProteinsInPath(subset_names = c("000_Trx","013","017", "081", "082"),n_s_euclidean = 1000, n_s_dijkstra = 500)
# p_vec = getProteinsInPath(n_s_euclidean = 4000, n_s_dijkstra = 500)

writeDistanceMatrixFLB <- function(d,path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/",euclidean,dijkstra, extra = ""){
  write.csv(d, file = paste(path,"/d_",euclidean,"_",dijkstra,"_",extra,".csv", sep =""))
}

readDistanceMatrixFLB <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/",euclidean,dijkstra, extra = ""){
  d = as.matrix(read.csv(file = paste(path,"/d_",euclidean,"_",dijkstra,"_",extra,".csv", sep =""), header = TRUE, check.names = FALSE, row.names = 1))
  return(d)
}

#------------------------------------------------------
# d = calcFlbAndPlot(p_vec,plotToFile = FALSE)
# write.csv(d, file = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/d_1000_500_small.csv")
# d2 = as.matrix(read.csv(file = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/d_1000_500_small.csv", header = TRUE, check.names = FALSE, row.names = 1))
# 
# writeDistanceMatrixFLB(d,euclidean = 1000, dijkstra = 500, extra = "small")
# d2 = readDistanceMatrixFLB(euclidean = 1000, dijkstra = 500, extra = "small")
# 
# memoliNNclassificationErrorEstimate(d,classlabels = getProteinLabels(d), normalized = TRUE)
# memoliNNclassificationErrorEstimate(d2,classlabels = getProteinLabels(d2), normalized = TRUE)

# plotModelsSideBySide(p_vec, "000_Trx", FALSE)

# trx = unfoldProteinVector(p_vec,1)
# p2 = unfoldProteinVector(p_vec,1)
# FLB_protein_combined_method(model_X_pos = trx$model_positive$prot_pos_sampled, model_X_neg = trx$model_negative$prot_neg_sampled,
#                             model_Y_pos = p2$model_positive$prot_pos_sampled, model_Y_neg = p2$model_negative$prot_neg_sampled)
#------------------------------------------------------

FLB_protein_combined_method <- function(model_X_pos, model_X_neg,model_Y_pos, model_Y_neg, L = 10){
  # pos and neg
  model_X_pos_and_neg = list("centers", "mu")
  model_X_pos_and_neg$centers = rbind(model_X_pos$centers, model_X_neg$centers)
  model_X_pos_and_neg$mu = c(model_X_pos$mu, model_X_neg$mu)
  
  
  model_Y_pos_and_neg = list("centers", "mu")
  model_Y_pos_and_neg$centers = rbind(model_Y_pos$centers, model_Y_neg$centers)
  model_Y_pos_and_neg$mu = c(model_Y_pos$mu, model_Y_neg$mu)
  
  s_X_pos_and_neg = get_s_2d(model_X_pos_and_neg)
  s_Y_pos_and_neg = get_s_2d(model_Y_pos_and_neg)
  
  E_pos = sort((c(s_X_pos_and_neg$s_X_pos, s_Y_pos_and_neg$s_X_pos)))
  E_neg = sort((c(s_X_pos_and_neg$s_X_neg, s_Y_pos_and_neg$s_X_neg)))
  
  
  f =   FLB_2d(model_X = model_X_pos_and_neg, model_Y = model_Y_pos_and_neg, 
               s_X = s_X_pos_and_neg, s_Y = s_Y_pos_and_neg,
               E_pos = E_pos, E_neg = E_neg, L = L)
  
  return(f)
}

# p_vec = getProteinsInPath(subset_names = c("000_Trx", "013"), n = 100, n_s_euclidean = 500, n_s_dijkstra = 2)
# trx = subsetOfProteins(p_vec, c("000_Trx"))
# p2 = subsetOfProteins(p_vec, c("013"))
# get_s_2d(trx$prot_pos_sampled)
# 
# FLB_protein_combined_method(model_X_pos = trx$prot_pos_sampled,
#                             model_X_neg = trx$prot_neg_sampled,
#                             model_Y_pos = p2$prot_pos_sampled,
#                             model_Y_neg = p2$prot_neg_sampled, L = 2)
# 
# FLB_protein_combined_method(model_Y_pos = trx$prot_pos_sampled,
#                             model_Y_neg = trx$prot_neg_sampled,
#                             model_X_pos = p2$prot_pos_sampled,
#                             model_X_neg = p2$prot_neg_sampled, L = 2)




FLB_protein_sum_method <- function(model_X_pos, model_X_neg,model_Y_pos, model_Y_neg, L = 10, c_vec = c(1,1,1)){
  
  # positive
  s_X_pos = get_s(model_X_pos$prot_pos_sampled)
  s_Y_pos = get_s(model_Y_pos$prot_pos_sampled)
  
  E_X_pos = unique(s_X_pos)
  E_Y_pos = unique(s_Y_pos)
  E_pos = sort(unique(c(E_X_pos,E_Y_pos)))
  
  # negative
  s_X_neg = get_s(model_X_neg$prot_neg_sampled)
  s_Y_neg = get_s(model_Y_neg$prot_neg_sampled)
  
  E_X_neg = unique(s_X_neg)
  E_Y_neg = unique(s_Y_neg)
  E_neg = sort(unique(c(E_X_neg,E_Y_neg)))
  
  # pos and neg
  model_X_pos_and_neg = list("centers", "mu")
  model_X_pos_and_neg$centers = rbind(model_X_pos$prot_pos_sampled$centers, model_X_neg$prot_neg_sampled$centers)
  model_X_pos_and_neg$mu = c(model_X_pos$prot_pos_sampled$mu, model_X_neg$prot_neg_sampled$mu)
  
  model_Y_pos_and_neg = list("centers", "mu")
  model_Y_pos_and_neg$centers = rbind(model_Y_pos$prot_pos_sampled$centers, model_Y_neg$prot_neg_sampled$centers)
  model_Y_pos_and_neg$mu = c(model_Y_pos$prot_pos_sampled$mu, model_Y_neg$prot_neg_sampled$mu)
  
  s_X_pos_and_neg = get_s(model_X_pos_and_neg)
  s_Y_pos_and_neg = get_s(model_Y_pos_and_neg)
  
  E_X_pos_and_neg = unique(s_X_pos_and_neg)
  E_Y_pos_and_neg = unique(s_Y_pos_and_neg)
  E_pos_and_neg = sort(unique(c(E_X_pos_and_neg,E_Y_pos_and_neg)))
  
  
  f =   c_vec[1] * FLB(model_X_pos$prot_pos_sampled,model_Y_pos$prot_pos_sampled,s_X_pos,s_Y_pos,E_pos,L)
  + c_vec[2] * FLB(model_X_neg$prot_neg_sampled,model_Y_neg$prot_neg_sampled,s_X_neg,s_Y_neg,E_neg,L)
  + c_vec[3] * FLB(model_X_pos_and_neg,model_Y_pos_and_neg,s_X_pos_and_neg,s_Y_pos_and_neg,E_pos_and_neg,L)
  
  return(f)
}

#---------------------------
# p = getProteinsInPath(subset_names = c("000_Trx", "016", "034"),n_s_euclidean = 4000,n_s_dijkstra = 5)
# calculateFlbDistancesOfProteins(p_vec = p,sumMethod = TRUE)



#---------------------------

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
# library(readobj)
# 
# # memoliModelPath = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/"
# # horse-01.obj
# memoliModelPath = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/horse-poses/"
# # model_name = "horse-01"
# # 
# horse01 = getMemoliModel3(memoliModelPath,"horse-01")
# horse02 = getMemoliModel3(memoliModelPath,"horse-02")
# horse03 = getMemoliModel3(memoliModelPath,"horse-03")
# horse04 = getMemoliModel3(memoliModelPath,"horse-04")
# # 
# head01 = getMemoliModel3(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/head-poses/","head-01-anger")
# head02 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/head-poses/","head-02-cry")
# # 
# cat01 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/cat-poses/","cat-01", n_s_euclidean = 500, n_s_dijkstra = 50)
# 
# 
# plotMemoliModel(model = cat01$rgl_model, points = cat01$points, centers = cat01$centers, measure = cat01$mu, plotModel = FALSE, scale = 1000)
# 
# 
# 
# cat02 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/cat-poses/","cat-02")
# cat03 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/cat-poses/","cat-03")
# cat04 = getMemoliModel3("/home/willy/RedoxChallenges/MasterThesis/memoliModels/cat-poses/","cat-04")
# 
#

#-----------------------------------------------------------------------------
# # shortest path
# horse1 = read.obj("/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/cat-poses/cat-02.obj", convert.rgl = FALSE)
# horse1_rgl = read.obj("/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/cat-poses/cat-02.obj", convert.rgl = TRUE)
# shade3d(horse1_rgl)
# horse1_points = t(horse1$shapes[[1]]$positions)
# horse1_edges = t(horse1$shapes[[1]]$indices)+1
# 
# ob = preProcessMesh(points = horse1_points, edges = horse1_edges, plot = FALSE)
# 
# myDrawShortestPath(ob,200,5000, horse1_rgl)
# rgl.snapshot(filename = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/cat_shortestPath", fmt = "png", top = TRUE )
#-----------------------------------------------------------------------------

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

# heatmap.2(as.matrix(d), trace="none", main = "posAndNeg",Rowv = FALSE,dendrogram ="column", key = FALSE)


getMemoliModelsInPath <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/", n_s_euclidean = 4000, n_s_dijkstra = 50){
  obj_files = list.files(path = path, pattern = ".obj",recursive = TRUE, full.names = TRUE)
  
  print(obj_files)
  
  
  model_vec = c()
  
  
  for(i in 1:length(obj_files)){
    obj_name = strsplit(obj_files[i],split = "/")[[1]][10]
    f_path = strsplit(obj_files[i], split = obj_name)[[1]][1]
    
    obj_name = strsplit(obj_name,split = ".obj")[[1]][1]
    
    print(obj_name)
    
    while (rgl.cur() > 0) { rgl.close() }
    model = getMemoliModel3(f_path,obj_name, n_s_euclidean, n_s_dijkstra)

    model_vec = c(model_vec, obj_name, model)
  }

  return(model_vec)
}

# model_vec = getMemoliModelsInPath(n_s_euclidean = 500, n_s_dijkstra = 50)
# 
# plotModelAtInd(model_vec, 1, FALSE)
# 
# 
# o = getModelAtInd(model_vec,1)





getModelAtInd <- function(model_vec, ind){

  m = (ind-1)*5
  
  obj = list("name" =  model_vec[m+1][[1]][1], "rgl_model" = model_vec[m+2]$rgl_model, "points" = model_vec[m+3]$points, "centers" = model_vec[m+4]$centers, "mu" = model_vec[m+5]$mu)
  return(obj)
}


plotModelAtInd <- function(model_vec, ind, plotModel = TRUE, plotMu = FALSE){
  mod = getModelAtInd(model_vec,ind)
  
  plotMemoliModel(model = mod$rgl_model, points = mod$points, centers = mod$centers, measure = mod$mu, plotModel = plotModel, scale = 1000, plotMu = plotMu)
}



calc_flb_distances <- function(model_vec){
  n = length(model_vec)/5
  d = matrix(0,nrow = n, ncol = n)
  
  # print(n)
  
  modelnames = rep("",n)
  
  for(i in 1:n){
    for(j in 1:n){
      model_X = getModelAtInd(model_vec = model_vec, ind = i)
      model_Y = getModelAtInd(model_vec = model_vec, ind = j)
      
      # d[i][j] = 10
      d[i,j] = FLB_method(model_X=model_X, model_Y=model_Y)
      
      modelnames[i] = model_X$name
      modelnames[j] = model_Y$name
    }
  }
  
  colnames(d) = modelnames
  rownames(d) = modelnames

  return(d)
}

myMakeSnapshot <- function(modelVec, ind, fileName, fileExt = "png", theta = -70, phi = 25, plotModel = TRUE, plotMu = TRUE){
  # close all other rgl-windows
  while (rgl.cur() > 0) { rgl.close() }
  
  par3d(windowRect = c(0, 0, 1000, 1000))
  
  plotModelAtInd(modelVec,ind,plotMu = plotMu, plotModel = plotModel)
  
  rgl.viewpoint(theta = theta, phi = phi)

  
  fName = paste(fileName,".",fileExt, sep ="")
  
  # for example svg
  if(fileExt != "png") rgl.postscript(filename = fName, fmt = fileExt)
  else rgl.snapshot(filename = fName, fmt = fileExt, top = TRUE )
  
  while (rgl.cur() > 0) { rgl.close() }
}

myMakeSnapshotProtein <- function(modelVec, names, fileName, fileExt = "png", theta = -70, phi = 25, plotModel = TRUE, plotMu = TRUE){
  # close all other rgl-windows
  while (rgl.cur() > 0) { rgl.close() }
  
  plotModelsSideBySide(modelVec,names,plotMu = plotMu, plotModel = plotModel)

  rgl.viewpoint(theta = theta, phi = phi)
  
  fName = paste(fileName,".",fileExt, sep ="")
  
  # for example svg
  if(fileExt != "png") rgl.postscript(filename = fName, fmt = fileExt)
  else rgl.snapshot(filename = fName, fmt = fileExt, top = TRUE )
  
  while (rgl.cur() > 0) { rgl.close() }
}

# ---------------------------------------------------------------------------------------------------------
# # memoli models
# model_vec = getMemoliModelsInPath(n_s_euclidean = 4000, n_s_dijkstra = 50)
# d = calc_flb_distances(model_vec)
# 
# # write.csv(d, file = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/dist_matrix_4000_50.csv")
# 
# # pdf("/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/heatmap.pdf")
# # heatmap.2(as.matrix(d), trace="none", main = "models",Rowv = FALSE,dendrogram ="column", key = FALSE, breaks = seq(0,0.7,0.01))
# heatmap.2(as.matrix(d), trace="none", main = "models",Rowv = FALSE,dendrogram ="column", key = FALSE)
# dev.off()
# 
# 
# plotModelAtInd(model_vec,4,plotMu = FALSE)
# rgl.snapshot(filename = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/camel4", fmt = "png", top = TRUE )
# 
# plotModelAtInd(model_vec,17)
# rgl.snapshot(filename = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/cat", fmt = "png", top = TRUE )
# 
# 
# myMakeSnapshot(modelVec = model_vec, ind = 57, fileName = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/horse1", fileExt = "png")
# myMakeSnapshot(modelVec = model_vec, ind = 60, fileName = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/horse2", fileExt = "png")
# myMakeSnapshot(modelVec = model_vec, ind = 62, fileName = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/horse3", fileExt = "png")
# myMakeSnapshot(modelVec = model_vec, ind = 58, fileName = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/horse4", fileExt = "png")
# myMakeSnapshot(modelVec = model_vec, ind = 61, fileName = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/horse5", fileExt = "png")
# 
# 
# # animalLabels = c(rep("camel", 11), rep("cat", 10), rep("elephant", 11), rep("flamingo", 11), rep("head", 10), rep("horse", 10), rep("lion", 10))
# # animalLabels
# # 
# # 
# # conf = memoliNNclassificationErrorEstimate(d,animalLabels,10000)
# # write.csv(conf, file = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/Models/conf.csv")

# ---------------------------------------------------------------------------------------------------------

# memoliNNclassification(d,animalLabels, c(3,13,19), returnClassName = FALSE)

memoliNNclassification_old <- function(d, classlabels, sample_to_predict, n=1, returnClassName = TRUE){
  
  colnames(d) = classlabels
  
  # print(colnames(d))
  
  predictions = rep("", length(sample_to_predict))
  for(i in 1:length(sample_to_predict)){
    
    d_s = d[sample_to_predict[i], -sample_to_predict]
    
    # if(returnClassName){
    # predictions[i] = names(d_s[which.min(d_s)])
    
    names = names(d_s[which.minn(d_s, n=n)])
    if(n > 1){
      if(length(which(names == "functional")) > 0){
        predictions[i] = "functional"
      } else {
        predictions[i] = "not_functional"
      }
    } else {
      predictions[i] = names
    }
    # predictions[i] = names(d_s[which.minn(d_s, n=n)])
    # } else {
    #   predictions[i] = names(d_s[which.min(d_s)])
    # }
    
    # ?which.minn() <-------------
    
  }
  
  # predictions = strsplit(predictions,"\\.")[[1]][1]
  
  # ?strsplit()
  
  return(predictions)
}

assignWeights <- function(d){
  
  names = names(d)
  uniqueLabels = unique(names)
  
  # print(d)
  
  labelWeights = rep(0,length(uniqueLabels))
  
  weights = rep(0,length(names))
  for(i in 1:length(uniqueLabels)){
    currentSet = which(names == uniqueLabels[i])
    labelWeights[i] = length(currentSet)/length(names)
    
    weights[currentSet] = 1/labelWeights[i]
  }
  
  # print(weights)
  
  return(weights)
  # return(list("labelWeights" = labelWeights, "weights" = weights))
}

memoliNNclassification <- function(d, classlabels, sample_to_predict, n=1, returnClassName = TRUE){
  # with weights
  
  colnames(d) = classlabels
  
  predictions = rep("", length(sample_to_predict))
  for(i in 1:length(sample_to_predict)){
    
    d_s = d[sample_to_predict[i], -sample_to_predict]

    w = assignWeights(d_s)

      
      # names = names(d_s[which.minn(d_s, n=n)])
      # if(n > 1){
      #   if(length(which(names == "functional")) > 0){
      #     predictions[i] = "functional"
      #   } else {
      #     predictions[i] = "not_functional"
      #   }
      # } else {
      #   predictions[i] = names
      # }
    
    NN = which.minn(d_s,n=n)
    candidateLabels = names(d_s)[NN]
    
    # print(candidateLabels)
    # print(w[NN])
    
    df = data.frame(candidateLabels,w[NN])
    names(df) = c("label", "weight")
    
    # df_avg_dist = aggregate(weight ~ label, df, sum)

    df = aggregate(weight ~ label, df, sum)
    
    
    
    
    # return(df)
    
    # print(df$label[2])
    
    # print(which.max(df$weight))
    # print(df$label[which.max(df$weight)])
    
    predictions[i] = as.character(df$label[which.max(df$weight)])
  }
  
  return(predictions)
}

# memoliNNclassification(d_felix,classlabels = getProteinLabels(d_felix), 12, n=105)
# 
# 
# aggregate(weight ~ label, df, sum)
# 
# 
# d_labeled = d_felix
# colnames(d_labeled) = getProteinLabels(d_felix)
# 
# assignWeights(d_labeled[1,])
# 
# 
# d = readDistanceMatrixFLB(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output_distances/", euclidean = 4000, dijkstra = 5)
# # memoliNNclassification(d,classlabels = getProteinLabels(d), 12, n=1)
# memoliNNclassificationErrorEstimate(d, getProteinLabels(d), n = 1000, kNN = 10,normalized = TRUE)



memoliNNclassificationErrorEstimate <- function(d,classlabels,n = 10000, kNN = 1, normalized = FALSE){
  labels_un = sort(unique(classlabels))
  random_sample = matrix(0,nrow = n, ncol = length(labels_un))
  colnames(random_sample) = labels_un
  for(i in 1:length(labels_un)){
    random_camel = which(classlabels == labels_un[i])
    # print(random_camel)
    random_sample[,i] = sample(random_camel,n, replace = TRUE)
  }
  
  conf = matrix(0,ncol = length(labels_un), nrow = length(labels_un))
  colnames(conf) = paste(labels_un,"_pred",sep="")
  rownames(conf) = labels_un
  
  for(i in 1:n){
    pred = memoliNNclassification(d,classlabels,random_sample[i,],n=kNN)
    
    for(j in 1:nrow(conf)){
      for(k in 1:nrow(conf)){
        # print(names(conf[,1])[k])
        if(names(conf[,1])[k] == pred[j]) conf[j,k] = conf[j,k] + 1
      }
    }
  }
  
  # test if rowsum = n
  for(i in 1:nrow(conf)){
    if(sum(conf[i,]) != n) {
      print(conf)
      
      print(paste("sum(conf[i,]) = ", sum(conf[i,]), " != ", n))
      return(FALSE)
    }
  }
  
  if(normalized == TRUE){
    conf = conf/n
  }
  
  return(conf)
}


optimizeNNclassification <- function(d,classlabels,n=10000,k_values = c(1:25), normalized = FALSE){
  
  k_stats = data.frame(matrix(0,nrow = length(k_values), ncol = 5))
  names(k_stats) = c("k", "TP", "TN", "FN", "FP")
  for(i in 1:length(k_values)){
    
    printMyMessage(message = "optimizeNN: ",val = paste(i/length(k_values)*100,"%"), debug_lvl = 1)
    t = memoliNNclassificationErrorEstimate(d = d,classlabels = classlabels, n = n, normalized = normalized,kNN = k_values[i])
    
    k_stats[i,] = c(k_values[i],t[1,1], t[2,2], t[1,2],t[2,1])
  }
  
  return(k_stats)
}

plotKNNOpt <- function(fName = "", k_stats, plotToFile = FALSE, method = ""){
  if(plotToFile == TRUE) pdf(file = fName)
  
  colors = c("red", "blue", "green")
  
  accuracy = (k_stats$TP+k_stats$TN)/(k_stats$TP+k_stats$TN+k_stats$FP+k_stats$FN)
  k_max = which.max(accuracy)
  acc_max = max(accuracy)
  
  # Accuracy = (TP+TN)/(TP+TN+FP+FN) 
  plot(y = k_stats$TP, x = k_stats$k, type = "l", col = colors[1], ylim = c(0,1), xlab = "k", ylab = "rate", main = paste(method,": maximal Accuracy at k = ", k_max, "(", acc_max, ")", sep = ""))
  points(y = k_stats$TN, x = k_stats$k, type = "l", col = colors[2])
  

  points(y = accuracy, x = k_stats$k, type = "l", col = colors[3])
  
  abline(v = k_max)
  
  # ?abline
  
  
  # legend(10,0.2, legend = c("TP", "TN", "Acc"), col = colors,pch = rep("l",3))
  legend(10, 0.2, legend=c("TP", "TN", "Acc"),
         col=colors, lty=rep(1,4), cex=0.8)
  
  
  if(plotToFile == TRUE) dev.off()
  
  k_stats = cbind(k_stats,accuracy)
  return(k_stats[k_max,])
}

# memoliNNclassificationErrorEstimate(d_felix, getProteinLabels(d_felix), n = 100, kNN = 12,normalized = TRUE)
# conf


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

calcFarthestPoints <- function(data, k = 10){
  # data ... points, columns are the x,y,z-position, each row one point
  # k ... number of farthest points to calculate
  # farthest ... matrix of size nrow(data)Xk with the distances to the farthest k points
  # farthest_indices ... matrix of the farthest k points 
  
  farthest = matrix(0,nrow = nrow(data), ncol = k)
  farthest_indices = matrix(0,nrow = nrow(data), ncol = k)
  
  for(i in 1:nrow(data)){
    d_1 <- apply(data, 1, function(x){sqrt(sum((x-data[i,])^2))})
    
    
    farthest_indices1 = which.maxn(d_1, n = k)
    farthest[i,] = d_1[farthest_indices1]
    farthest_indices[i,] = farthest_indices1
  }
  
  return(list("farthest" = farthest, "farthest_indices" = farthest_indices))
}

myDirectionalHausdorfDistance <- function(open_points,selected_points){
  # print(open_points)
  # print(selected_points)
  
  knn_index = knnx.index(data = selected_points,query = open_points, k = 1)
  knn_dist = knnx.dist(data = selected_points,query = open_points, k = 1)
  
  return(list("knn_index" = knn_index, "knn_dist" = knn_dist))
}

myFarthestPointSampling <- function(data, start_ind = 1, k = 10){
  # data ... points, columns are the x,y,z-position, each row one point
  # k ... number of farthest points to calculate
  # farthest ... matrix of size nrow(data)Xk with the distances to the farthest k points
  # farthest_indices ... matrix of the farthest k points 
  
  library(pracma)
  
  available_indices = c(1:nrow(data))
  
  farthest_index = start_ind
  
  sample_indices = rep(0,k)
  
  sample_indices[1] = farthest_index

  
  available_indices = available_indices[-which(available_indices == farthest_index)]
  
  d_1 <- apply(data[available_indices,], 1, function(x){sqrt(sum((x-data[farthest_index,])^2))})

  farthest_indices1 = which.maxn(d_1, n = 1)
  
  sample_indices[2] = available_indices[farthest_indices1]
  
  available_indices = available_indices[-farthest_indices1]
  
  sampleCount = 2
  
  if(k >2){
    # smaple k points
    for(sampleCount in 3:k){
      
      print(paste(sampleCount/k*100, " %", sep = ""))
      
      # gets the minimal distances in the still available points
      f = myDirectionalHausdorfDistance(data[available_indices,],data[-available_indices,])
      
      # off all the minimal distances take the biggest one
      # point_to_sample = f$knn_index[which(f$knn_dist == max(f$knn_dist), arr.ind = TRUE)]
      
      ps = which(f$knn_dist == max(f$knn_dist))[1]
      point_to_sample = available_indices[ps]
      
      # print(f)
      
      # print("ps")
      # print(ps)
      
      # sample[sampleCount,] = data[point_to_sample,]
      # data = data[-point_to_sample,]
      
      # print(point_to_sample)
      
      available_indices = available_indices[-which(available_indices == point_to_sample)]
      
      sample_indices[sampleCount] = point_to_sample
    }
  }
  
  return(sample_indices)
}

# #-------------------------------------------
# # ## example farthest point sampling
# fpsPlot(n = 100,k = 10)
# n = 50
# dim = 2
# m = matrix(rnorm(n*dim),nrow=n, ncol = dim)
# ind = c()
# for(i in 1:10){
#   ind = fpsPlot2d(m, k = i, plotToFile = "/home/willy/RedoxChallenges/MasterThesis/WillyMT/images/FarthestPointSampling/")
# }
# 
# ind
# 
# x = m[,1]
# y = m[,2]
# points <- data.frame(x[ind], y[ind], distance = sqrt((x-100)^2 + (y-100)^2))
# 
# library(ggplot2)
# # install.packages("ggvoronoi")
# library(ggvoronoi)
# 
# pdf(file = "/home/willy/RedoxChallenges/MasterThesis/WillyMT/images/FarthestPointSampling/voronoi10.pdf")
# p <- ggplot(points[ind,],aes(x[ind],y[ind]), col = "green") + stat_voronoi(geom="path") + geom_point(color="green") + geom_point(x = c(1,2),y = c(1,1))
# plot(p)
# 
# 
# 
# dev.off()
# 
# # fpsPlot2d(m, k = 10, plotToFile = "/home/willy/RedoxChallenges/MasterThesis/WillyMT/images/FarthestPointSampling/")
# 
# 
# # ?aes
# 
# #-------------------------------------------

fpsPlot <- function(n,k){
  # n = 10
  dim = 3
  m = matrix(rnorm(n*dim),nrow=n, ncol = dim)
  
  points3d(m, col = "red", size = 10)
  ind = myFarthestPointSampling(m, 1, k)
  print(ind)
  points3d(m[ind,], col = "green", size = 15)
  texts3d(m[ind,], texts = as.character(c(1:length(ind))), cex = 4)
}

fpsPlot2d <- function(m,k, plotToFile="empty"){

  startPoint = 1
  

  # ?pdf
  if(plotToFile != "empty"){
    pdf(file = paste(plotToFile,"FarthestPointSampling_",k, ".pdf", sep = ""))
  }
  
  l = c(min(min(m[,1]), min(m[,2])), max(max(m[,1]), max(m[,2])))
  
  plot(m, col = "red", xlim = l, ylim = l,  ylab = "", xlab = "")
    
  ind = myFarthestPointSampling(m, startPoint, k)
  
  if(k == 1) {
    print(m)
    points(m[startPoint,1], m[startPoint,2], col = "green")
    text(m[startPoint,1] + c(0.05), m[startPoint,2], as.character(startPoint), cex = 1)
  }
  # print(ind)
  else {
    points(m[ind,], col = "green")
    text(m[ind,] + c(0.05), as.character(c(1:length(ind))), cex = 1)
  }

  if(plotToFile != "empty"){
    dev.off()
  }
  
    
  # text(m + c(0.05), as.character(c(1:nrow(m))), cex = 1, col = "red")
  
  # print(dist(m), upper = TRUE)
  
  return(ind)
}


plot2dFarthest <- function(data, ind, n){
  # obj = calcFarthestPoints(m, k = n)
  
  obj = myFarthestPointSampling(data, 1, k = n)
  
  plot(x = data[,1], y = data[,2])
  
  points(x = obj[,1], y = obj[,2],col="red")
  
  points(x = data[ind,1], y = data[ind,2], col ="green")
}

myDownsampleIndices <- function(points, n = 1){
  # n ... number of points to keep
  
  points_sample = seq(1, nrow(points), nrow(points)/n)
  points_downsampled2 = points[points_sample,]
  # toy_points_pos_downsampled = points[get.knnx(data = points, query = points_downsampled2, k = 1, algo = "kd_tree")$nn.index,]
  
  sampled_points_knn = get.knnx(data = points, query = points_downsampled2, k = 1, algo = "kd_tree")$nn.index
  
  return(sampled_points_knn)
}

# myDownsampleIndices(points = m, n = 3)

# myFarthestPointSampling

getOrderPath <- function(path, name, n_s_euclidean, n_s_dijkstra, positive = TRUE){
  ext = paste(n_s_euclidean, "_", n_s_dijkstra, sep = "")
  
  if(positive == TRUE){
    return(paste(path, "/", name,"/PtsAndMu/",name, ext,"_pos.order",sep =""))
  } else {
    return(paste(path, "/", name,"/PtsAndMu/",name, ext,"_neg.order",sep =""))
  }

}

getProtein <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", name, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE, reCalculate = FALSE){
 
  if(!dir.exists(paste(path, "/", name,"/PtsAndMu/", sep = ""))) dir.create(paste(path, "/", name,"/PtsAndMu/", sep = ""))
  
  ext = paste(n_s_euclidean, "_", n_s_dijkstra, sep = "")
  
  ptsPath_pos = paste(path, "/", name,"/PtsAndMu/",name, ext,"_pos.pts",sep ="")
  muPath_pos = paste(path, "/", name,"/PtsAndMu/",name, ext, "_pos.mu",sep ="")
  
  ptsPath_neg = paste(path, "/", name,"/PtsAndMu/",name, ext,"_neg.pts",sep ="")
  muPath_neg = paste(path, "/", name,"/PtsAndMu/",name, ext,"_neg.mu",sep ="")
  objPath = paste(path, "/", name,"/",name, ".obj",sep ="")
  
  orderPathNeg = getOrderPath(path, name, n_s_euclidean, n_s_dijkstra, positive = FALSE)
  orderPathPos = getOrderPath(path, name, n_s_euclidean, n_s_dijkstra, positive = TRUE)
  
  ob_pos = list("points" = c(), "edges" = c(), "graph" = c())
  ob_neg = list("points" = c(), "edges" = c(), "graph" = c())
  
  prot_rgl = read.obj(objPath, convert.rgl = TRUE)
  # preprocessing mesh
  if(!file.exists(ptsPath_pos) || !file.exists(ptsPath_neg) || !file.exists(muPath_pos) || !file.exists(muPath_neg) || reCalculate == TRUE){
    prot = read.obj(objPath, convert.rgl = FALSE)
    if(plot) shade3d(prot_rgl)
    prot_points_pos = t(prot$shapes[[2]]$positions)
    prot_edges_pos = t(prot$shapes[[2]]$indices)+1
    
    prot_points_neg = t(prot$shapes[[3]]$positions)
    prot_edges_neg = t(prot$shapes[[3]]$indices)+1
    
    
    # if(checkForLargerModel(path,proteinName = name,n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra) == FALSE)
    print(paste("model has ", nrow(prot_points_pos),", ", nrow(prot_points_neg), " points", sep ="" ))
    
    ob_pos = preProcessMesh(points = prot_points_pos, edges = prot_edges_pos, plot = FALSE)
    ob_neg = preProcessMesh(points = prot_points_neg, edges = prot_edges_neg, plot = FALSE)
    
    print(paste("processed model has ", nrow(ob_pos$points),", ", nrow(ob_neg$points), " points", sep ="" ))
  }
  
  prot_pos_sampled = getMemoliModel(prot_rgl, path = path, name = name,  points = ob_pos$points, edges = ob_pos$edges, graph = ob_pos$graph, memoliPtsPath = ptsPath_pos, memoliMuPath = muPath_pos, memoliOrderPath = orderPathPos,n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra, reCalculate = reCalculate)
  prot_neg_sampled = getMemoliModel(prot_rgl, path = path, name = name, points = ob_neg$points, edges = ob_neg$edges, graph = ob_neg$graph, memoliPtsPath = ptsPath_neg, memoliMuPath = muPath_neg, memoliOrderPath = orderPathNeg,n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra, reCalculate = reCalculate)
  
  
if(plot) plotMemoliModel2(prot_rgl, centers1 = prot_pos_sampled$centers, measure1 = prot_pos_sampled$mu,
                          centers2 = prot_neg_sampled$centers, measure2 = prot_neg_sampled$mu, scale = 500, plotModel = FALSE)
  
  
  l = list("rgl_model" = prot_rgl, "prot_pos_sampled" = prot_pos_sampled, "prot_neg_sampled" = prot_neg_sampled)
  return(l)
}

# p_vec = getProteinsInPath(n_s_euclidean = 1000, n_s_dijkstra = 50)

getProteinsInPath <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", n_s_euclidean = 100, n_s_dijkstra = 15, subset_names = FALSE, reCalculate = FALSE){

  proteinNames = list.dirs(path, full.names = FALSE, recursive = FALSE)

  print("------------------------------------------------------------------")
  print(paste("Found the following proteins in ", path, ": "))
  print(proteinNames)
  print("------------------------------------------------------------------")

  p_vec = c()
  for(pName in proteinNames){
    
    if(subset_names == FALSE || pName %in% subset_names){
      print(paste("reading ", pName, " ...", sep =""))
      
      p = getProtein(path = path, pName, n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra, reCalculate = reCalculate)
      p_vec = c(p_vec,pName,p)
    }

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

calculateFlbDistancesOfProteins <- function(p_vec, sumMethod = FALSE, L = 10, c_vec = c(1,1,1)){
  numberOfProteins = length(p_vec)/4
  
  d = matrix(0,ncol = numberOfProteins, nrow = numberOfProteins)
  
  pNames = rep(" ", numberOfProteins)
  for(i in 1:(numberOfProteins-1)){
    for(j in (i+1):numberOfProteins){
      
      print(paste(i,j,sep = " "))
      print(paste(((i-1)*numberOfProteins+j) / ((numberOfProteins*numberOfProteins)/2), " %", sep = ""))

      
      unfolded =  unfoldProteinVector(p_vec,i)
      unfolded2 =  unfoldProteinVector(p_vec,j)

      # d_pos[i,j] = FLB_method(unfolded$model_positive$prot_pos_sampled, unfolded2$model_positive$prot_pos_sampled)
      # d_neg[i,j] = FLB_method(unfolded$model_negative$prot_neg_sampled, unfolded2$model_negative$prot_neg_sampled)
      
      if(sumMethod == TRUE){
        d[i,j] = FLB_protein_sum_method(unfolded$model_positive, unfolded$model_negative,unfolded2$model_positive, unfolded2$model_negative, L = L, c_vec = c_vec)
      } else {
        d[i,j] = FLB_protein_combined_method(unfolded$model_positive$prot_pos_sampled, unfolded$model_negative$prot_neg_sampled
                                             ,unfolded2$model_positive$prot_pos_sampled, unfolded2$model_negative$prot_neg_sampled,L = L)
      }
      
      d[j,i] = d[i,j]
      
      pNames[i] = unfolded$name
      pNames[j] = unfolded2$name
    }
  }
  
  row.names(d) = pNames
  colnames(d) = pNames

  
  return(d)
}

plotHeatMap <- function(d, fName, plotToFile = TRUE){
  if(plotToFile) pdf(fName)
  
  f = getFunctionalProteins()
  library(gplots)
  colors = colnames(d) %in% f
  colors = replace(colors, colors==TRUE, "red")
  colors = replace(colors, colors==FALSE, "black")
  colors[which(colnames(d) == "000_Trx")] = "green"
  
  heatmap.2(d, trace="none", colCol = colors, main = "Heatmap",Rowv = FALSE,dendrogram ="column", key = FALSE)
  if(plotToFile) dev.off()
}

calcFlbAndPlot <- function(p_vec,path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", plotToFile = TRUE, L = 10){
  d = calculateFlbDistancesOfProteins(p_vec, L = L)
  
  if(plotToFile) pdf(paste(path,"posAndNeg.pdf", sep = ""))
  
  f = getFunctionalProteins()
  library(gplots)
  colors = colnames(d) %in% f
  colors = replace(colors, colors==TRUE, "red")
  colors = replace(colors, colors==FALSE, "black")
  colors[which(colnames(d) == "000_Trx")] = "green"
  
  # ?heatmap.2()
  
  heatmap.2(d, trace="none", colCol = colors, main = "posAndNeg",Rowv = FALSE,dendrogram ="column", key = FALSE)
  if(plotToFile) dev.off()
  
  return(d)
}

getKnnOfProtein <- function(d, k = 20, name = "000_Trx"){
  query_ind = as.numeric(which(colnames(d) == name)[1])
  
  # indices = get.knn(data = d[query_ind,], k = k, algo = "kd_tree")$nn.index[query_ind,]
  
  d = d[,-query_ind]
  
  indices = sort(d[query_ind,])
  # colnames(d[get.knnx(data = d, query = ind, k = k, algo = "kd_tree")$nn.index,])
  
  # print(indices)
  
  return(names(indices[1:k]))
}

# getKnnOfProtein(d[1:7,1:7],k=2)



evalKnn <- function(d, k = 20){
  q = getKnnOfProtein(d, k = k)
  
  f = getFunctionalProteins()
  
  # print(length(intersect(q,f)))
  
  return(list("n" = length(intersect(q,f)), "proteins" = intersect(q,f)))
}



plotKnnOfTrx <- function(d){
  # makes a plot
  # x-axis ... the number of nearest neighbours
  # y-axis ... the number of functional proteins
  
  f = getFunctionalProteins()
  
  trx_row = which(rownames(d) == "000_Trx")
  
  knns = rep(0,ncol(d)-1)
  knns_proteins = c()
  
  new_proteins_order = c()
  new_proteins_ind_x = c()
  
  for(k in 1:length(knns)){
    q = evalKnn(d,k=k)
    knns[k] = q$n
    
    new_protein = setdiff(q$proteins,knns_proteins)
    
    # print(new_protein)
    
    knns_proteins = q$proteins
    
    if(length(new_protein) != 0){
      new_proteins_order = c(new_proteins_order, new_protein)
      new_proteins_ind_x = c(new_proteins_ind_x,k)
    }
  }
  
  plot(c(1:length(knns)),knns, type = "l", xlab = "kNN", ylab = "functional proteins", ylim = c(0,length(f)))
  
  # print(length(new_proteins_ind_x))
  
  for(i in 1:length(new_proteins_ind_x)){
    text(labels = new_proteins_order[i], x = new_proteins_ind_x[i]-3, y = i, col = "red")
  }
  
}

# plotKnnOfTrx(d)
# d
# getKnnOfProtein(d[1:7,1:7],k = 2)
# 
# d[1:7,1:7]

subsetOfProteins <- function(p_vec, names){
  p_vec_out = c()
  
  # indices = which(p_vec == names)
  
  for(i in seq(1,length(p_vec), 4)){
    for(j in 1:length(names)){
      if(p_vec[i] == names[j]) p_vec_out = c(p_vec_out, p_vec[i], p_vec[i+1], p_vec[i+2], p_vec[i+3])
    }
  }
  
  return(p_vec_out)
}

plotModelsSideBySide <- function(p_vec, names, plotModel = TRUE, plotMu = TRUE){
  
  
  for(i in names){
    p_vec_small2 = subsetOfProteins(p_vec, i)
    
    # I don't understand why, but if we don't call shade3d
    # once before our model, the model will not be plotted
    # with a solid surface
    #---------
    shade3d(p_vec_small2$rgl_model$vmd_mol0_rep3)
    while (rgl.cur() > 0) { rgl.close() }  
    #---------
    
    plotMemoliModel2( model = p_vec_small2$rgl_model,
                      centers1 = p_vec_small2$prot_pos_sampled$centers, 
                      centers2 = p_vec_small2$prot_neg_sampled$centers,
                      measure1 = p_vec_small2$prot_pos_sampled$mu,
                      measure2 = p_vec_small2$prot_neg_sampled$mu,
                      scale = 100000/50,
                      plotModel = plotModel, plotMu = plotMu)
  }
  rgl.bg(color = "white")
}


# plotModelsSideBySide(trx,c("000_Trx"), plotModel = TRUE, plotMu = TRUE)
# 
# p_vec_small2 = subsetOfProteins(trx, c("000_Trx"))
# 
# shade3d(p_vec_small2$rgl_model$vmd_mol0_rep3)
# 
# shade3d(p_vec_small2$rgl_model)

# p_vec is a list, 4 consecutive elements are one protein
# p_vec = getProteinsInPath()
# p_vec = getProteinsInPath(n_s_euclidean = 1000, n_s_dijkstra = 500)

# p_vec = getProteinsInPath(n_s_euclidean = 4000, n_s_dijkstra = 1000)


# p_vec = getProteinsInPath(n = 20, n_s_euclidean = 1000, n_s_dijkstra = 500)

# proteinNames = list.dirs("/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output/", full.names = FALSE, recursive = FALSE)

# p_vec = getProteinsInPath(subset_names = proteinNames[1:60], n_s_euclidean = 4000, n_s_dijkstra = 25)
# p_vec = getProteinsInPath(subset_names = c("000_Trx", "013", "026", "010", "016", "017"), n_s_euclidean = 4000, n_s_dijkstra = 25)
# p_vec = getProteinsInPath(subset_names = c("000_Trx"), n_s_euclidean = 1000, n_s_dijkstra = 250)
# p_vec = getProteinsInPath(subset_names = c("000_Trx"),n_s_euclidean = 4000, n_s_dijkstra = 250)

getProteinLabels <- function(d){
  f = c(getFunctionalProteins(), "000_Trx")
  proteinLabels = colnames(d)
  
  for(i in f){
    ind = which(proteinLabels == i)
    # print(ind)
    if(length(ind) > 0) proteinLabels[ind[1]] = "functional"
  }
  
  proteinLabels = replace(proteinLabels,which(proteinLabels != "functional"), "not_functional")
  
  return(proteinLabels)
}

createConfMatrixAndHeatmap <- function(d,outPath, conf_name, pseudoRoc_Name, heatmap_name, dropNames = c("000_Ars", "000_DsbA", "000_DsbC", "000_Grx1", "000_Grx4_1", "000_NrdH")){
  namesToKeep = setdiff(colnames(d),dropNames)
  
  d = d[namesToKeep,namesToKeep]
  
  conf = memoliNNclassificationErrorEstimate(d,getProteinLabels(d),10000, kNN = 1,normalized = TRUE)
  write.csv(conf, file = conf_name)
  
  pdf(pseudoRoc_Name)
  plotKnnOfTrx(d)
  dev.off()
  
  plotHeatMap(d, heatmap_name,plotToFile = TRUE)
}

calcDistances <- function(outPath = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output_distances/", n_s_euclidean_vec = c(4000),n_s_dijkstra_vec=c(50), recalcConf = FALSE, reCalculate = FALSE, reCalculateSampling = FALSE, sumMethod = FALSE, c_vec = c(1,1,1)){
  
  # if(sumMethod == TRUE) outPath = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output_distances_sumMethod/"
  
  if(!dir.exists(outPath)) dir.create(outPath)
  
  for(i in 1:length(n_s_euclidean_vec)){
    for(j in 1:length(n_s_dijkstra_vec)){
      d_name = paste(outPath,"d_",n_s_euclidean_vec[i],"_", n_s_dijkstra_vec[j],"_.csv", sep = "")
      
      conf_name = paste(outPath,"conf_",n_s_euclidean_vec[i],"_", n_s_dijkstra_vec[j],".csv",  sep = "")
      pseudoRoc_Name = paste(outPath,"pseudoRoc_",n_s_euclidean_vec[i],"_", n_s_dijkstra_vec[j],".pdf",  sep = "")
      
      heatmap_name = paste(outPath,"heatmap_",n_s_euclidean_vec[i],"_", n_s_dijkstra_vec[j],".pdf",  sep = "")
      
      
      if(!file.exists(d_name) || !file.exists(conf_name) || !file.exists(pseudoRoc_Name) || !file.exists(heatmap_name) || recalcConf == TRUE || reCalculate == TRUE || reCalculateSampling == TRUE){
        
        p_vec = getProteinsInPath(n_s_euclidean = n_s_euclidean_vec[i], n_s_dijkstra = n_s_dijkstra_vec[j], reCalculate = reCalculateSampling)
        
        if(!file.exists(d_name) || reCalculate == TRUE){
          d = calculateFlbDistancesOfProteins(p_vec, sumMethod = sumMethod, c_vec = c_vec)
          # write.csv(d, file = d_name)
          writeDistanceMatrixFLB(d=d, path = outPath, euclidean = n_s_euclidean_vec[i], dijkstra = n_s_dijkstra_vec[j])
        }

        d = readDistanceMatrixFLB(path = outPath, euclidean = n_s_euclidean_vec[i], dijkstra = n_s_dijkstra_vec[j])
        
        
        # namesToKeep = setdiff(colnames(d),c("000_Ars", "000_DsbA", "000_DsbC", "000_Grx1", "000_Grx4_1", "000_NrdH"))
        # 
        # d = d[namesToKeep,namesToKeep]
        # 
        # conf = memoliNNclassificationErrorEstimate(d,getProteinLabels(d),10000, kNN = 1,normalized = TRUE)
        # write.csv(conf, file = conf_name)
        # 
        # pdf(pseudoRoc_Name)
        # plotKnnOfTrx(d)
        # dev.off()
        # 
        # plotHeatMap(d, heatmap_name,plotToFile = TRUE)
        createConfMatrixAndHeatmap(d = d,outPath = outPath, conf_name = conf_name, pseudoRoc_Name = pseudoRoc_Name, heatmap_name = heatmap_name)
      }
    }
  }
}

mergeConfusionMatrices <- function(path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output_distances/"){
  
  files = list.files(path = path, pattern = "conf_")

  if(length(files) == 0) return(FALSE)
  
  
  conf_matrices = data.frame(matrix(0,ncol = 6, nrow = length(files)))
  colnames(conf_matrices) = c("euclid", "surface", "TP", "TN", "FP", "FN")
  
  
  for(i in 1:length(files)){
    f2 = strsplit(files[i],split = ".csv")[[1]] 
    print(f2)
    
    f3 = strsplit(f2,split = "conf_")[[1]][2]
    print(f3)
    
    euclidean_n = as.numeric(strsplit(f3, split = "_")[[1]][1])
    print(euclidean_n)
    
    dijkstra_n = as.numeric(strsplit(f3, split = "_")[[1]][2])
    print(dijkstra_n)
    
    t = read.csv(file = paste(path,files[i], sep=""), header = TRUE, row.names = 1)
    
    print(t)
    conf_matrices[i,] = c(euclidean_n,dijkstra_n, t[1,1], t[2,2], t[1,2],t[2,1])
  }

  # conf_matrices = conf_matrices[order(conf_matrices$surface),]
  
  return(conf_matrices)
}

plotConfusionMatrices <- function(conf_matrices, path = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/myVmdToys/Output_distances/", plotToFile = TRUE, main =""){
  conf_matrices = conf_matrices[order(conf_matrices$surface),]
  
  colors = c("red", "green")
  maxSurf = max(conf_matrices$surface)
  
  # print(maxSurf)
  
  if(plotToFile == TRUE) pdf(paste(path,"PlotConfusionMatrices.pdf", sep =""))
  plot(conf_matrices$surface, conf_matrices$TP, ylim = c(0,1), xlim = c(0,maxSurf), col = colors[1], type = "l", main = main, ylab = "rate", xlab = "n_surface")
  points(conf_matrices$surface, conf_matrices$TN, ylim = c(0,1), xlim = c(0,maxSurf), col = colors[2], type = "l")
  # points(conf_matrices$surface, conf_matrices$FN, ylim = c(0,1), xlim = c(0,maxSurf), col = colors[3], type = "l")
  # points(conf_matrices$surface, conf_matrices$FP, ylim = c(0,1), xlim = c(0,maxSurf), col = colors[4], type = "l")
  
  legend(0, 0.6, legend=c("TP (=1-FP)", "TN (=1-FN)"),
         col=colors, lty=rep(1,4), cex=0.8)
  if(plotToFile == TRUE) dev.off()
  
  
  bestConf_ind = which.max(conf_matrices$TP+conf_matrices$TN)
  bestConf = matrix(c(conf_matrices$TP[bestConf_ind], conf_matrices$FP[bestConf_ind], 
                    conf_matrices$FN[bestConf_ind], conf_matrices$TN[bestConf_ind]), byrow = TRUE, nrow = 2, ncol = 2)
  
  colnames(bestConf) = c("functional_pred", "not_functional_pred")
  rownames(bestConf) = c("functional", "not_functional")
  
  print(xtable(bestConf, type = "latex", caption = paste("confusion matrix with $s_{\\text{dijkstra}} = $", conf_matrices$surface[bestConf_ind], " of mehtod ", main, sep = "")
               , label = "bestConfMatrix"), file = paste(path,"bestConfusion.tex", sep = ""))
  
  
  bestConf_ind_tpr = which.max(conf_matrices$TP)
  bestConf_tpr = matrix(c(conf_matrices$TP[bestConf_ind_tpr], conf_matrices$FP[bestConf_ind_tpr], 
                      conf_matrices$FN[bestConf_ind_tpr], conf_matrices$TN[bestConf_ind_tpr]), byrow = TRUE, nrow = 2, ncol = 2)
  
  colnames(bestConf_tpr) = c("functional_pred", "not_functional_pred")
  rownames(bestConf_tpr) = c("functional", "not_functional")
  
  print(xtable(bestConf_tpr, type = "latex", caption = paste("confusion matrix with the highest True-positive-rate ($s_{\\text{dijkstra}} = $",conf_matrices$surface[bestConf_ind], ") of method ", main, sep = "")
               , label = "bestConfMatrixTruePositiveRate"), file = paste(path,"bestConfusionTruePositiveRate.tex", sep = ""))


  return(conf_matrices[bestConf_ind,])  
}

# calcDistances(n_s_euclidean_vec = c(4000), n_s_dijkstra_vec = c(5,10,20,50,60,70,80,90), recalcConf = TRUE, reCalculate = TRUE)
# calcDistances(n_s_euclidean_vec = c(4000), n_s_dijkstra_vec = c(5,10,20,30,40,50,90,85,75,81,82,83,84,79,78,77,76,74,200), recalcConf = FALSE, reCalculate = FALSE)



readInConfmatrix <- function(fName, method_id){
  
  t = read.csv(file = fName, header = TRUE, row.names = 1)
  
  frame = data.frame(method_id,t[1,1],t[2,2],t[1,2],t[2,1])
  
  colnames(frame) = c("methods", "TP", "TN", "FP", "FN")
  
  return(frame)
}



#------------------------------------------------------------------------
# remove all pts-files with this command in the shell
# find -type d -name PtsAndMu -exec rm -rf {} \;
#------------------------------------------------------------------------


