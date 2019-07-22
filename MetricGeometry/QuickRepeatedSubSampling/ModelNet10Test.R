library(rgl)

s1 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R"
source(s1)

s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s2)


s3 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
source(s3)




downsampleEuclideanAndGetGeodesic <- function(objPath, n_s_euclidean = 4000, n_s_dijkstra = 50, plot = FALSE){
  model_rgl = read.obj(objPath, convert.rgl = FALSE)
  model_rgl_plot = read.obj(objPath, convert.rgl = TRUE)
  
  print("plotting")
  if(plot) shade3d(model_rgl_plot)
  points = t(model_rgl$shapes[[1]]$positions)
  edges = t(model_rgl$shapes[[1]]$indices)+1

  print("extracted points and edges")
  
  # if(checkForLargerModel(path,proteinName = name,n_s_euclidean = n_s_euclidean, n_s_dijkstra = n_s_dijkstra) == FALSE)
  print(paste("model has ", nrow(points),", ", nrow(edges), " points", sep ="" ))
  
  ob = preProcessMesh(points = points, edges = edges, plot = FALSE)
  
  print(paste("processed model has ", nrow(ob$points), "points", sep ="" ))
  
  graph = ob$graph
  edges = ob$edges
  
  return(graph)
  
  
  library(rdist)
  print("step 1: euclidean fps ...")
  
  sampled_indices = myFarthestPointSampling(points, k = n_s_euclidean)
  # if(plot) plotDownsampledPoints(model_rgl,ob$points,sampled_indices, plotModel = plot)
  
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



obj = read.obj("/home/willy/Schreibtisch/bathtub_0001.off")

#------------------------------------------------------------------
# off2obj bathtub_0001.off > bath.obj
# to convert to object wave format that we can read in with "read.obj"
obj = read.obj("/home/willy/Schreibtisch/bath.obj", convert.rgl = TRUE)

path = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/bathtub/train/"
name = "bathtub_0001.off"

name_final = paste(path,"/",name, sep ="")
name_final = "/home/willy/Schreibtisch/bath.obj"

g = downsampleEuclideanAndGetGeodesic(objPath = name_final,
                                  n_s_euclidean = 10,
                                  n_s_dijkstra = 5,
                                  plot = TRUE)


g[[3]]





