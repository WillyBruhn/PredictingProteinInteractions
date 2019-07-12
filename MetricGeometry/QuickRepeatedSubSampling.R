#---------------------------------------------------------------------
# Willy Bruhn, 10.7.2019
# Precalculating the Distributions F,G
# Then not the points but the distributions are sampled
#---------------------------------------------------------------------
print("loading libraries")

# install.packages("parallel")
library(parallel)

# install.packages("rbenchmark")
library(rbenchmark)

# install.packages("doBy")
library(doBy)

# install.packages("smacpod")
# library(smacpod)

# install.packages("raster")
library(raster)

# install.packages("plotrix")
# library(plotrix)

# install.packages("foreach")
library(foreach)

# install.packages("doMC")
library(doMC)

registerDoMC(cores=Sys.getenv('LSB_MAX_NUM_PROCESSORS'))
registerDoMC(cores= 6)

# install.packages("rlist")
# library(rlist)

# install.packages("ade4")
# library(ade4)

# install.packages("emdist")
library(emdist)

print("done lodaing libraries")

read_pts_file <- function(OutputPath,protName, pos = TRUE, ext_pos = "_pot_positive.pts", ext_neg = "_pot_negative.pts"){
  pts_file_name = ""
  if(pos == TRUE) {
    pts_file_name = paste(OutputPath,"/",protName,"/",protName,ext_pos, sep = "")
  } else{
    pts_file_name = paste(OutputPath,"/",protName,"/",protName,ext_neg, sep = "")
  } 
  
  return(read.csv2(pts_file_name,header = TRUE,colClasses = rep("numeric",3),dec = "."))
}

euclideanDistance <- function(a,b){
  d = a-b
  
  # print(d)
  return(sqrt(sum(d*d)))
}
# euclideanDistance(c(1,0),c(0,1))
# c(1,2)*c(1,2)
# euclideanDistance(pts_pos[1,],pts_pos[1,])


euclidNorm <- function(a){
  return(a/euclideanDistance(a,rep(0,length(a))))
}

calcEuclideanDistances <- function(all_points, ind_a, ind_b){
  # all_points ... all points; a whole pts_file
  # ind_a ... start indices
  # ind_b ... end_indices 
  # all distances from ind_a to ind_b are calculated
  # could use this to calculate DE of all points
  #----------------------------------------------------
  
  distances = matrix(0,nrow = length(ind_a), ncol = length(ind_b))
  # print(distances)
  for(i in 1:length(ind_a)){
    for(j in 1:length(ind_b)){
      print(paste(i,j))
      distances[i,j] = euclideanDistance(all_points[ind_a[i],], all_points[ind_b[j],])
    }
  }

  return(distances)
    
}
# calcEuclideanDistances(pts_pos, c(1), c(1:nrow(pts_pos)))

eccentricities <- function(d, mu = rep(1/nrow(d),nrow(d))){
  return(colSums(x = d*mu))
}

DistributionOfEccentricities <- function(d, mu = rep(1/nrow(d),nrow(d))){
  F_ = ecdf(eccentricities(d))
  return(F_)
}

samplePointsAndCalculateCDFofEc <- function(all_pts,n, plot = FALSE){
  sample_ind = sample(1:nrow(all_pts), size = n, replace = FALSE)

  d = as.matrix(dist(all_pts[sample_ind,],upper = TRUE, diag = TRUE))
  
  F_ = DistributionOfEccentricities(d)
  
  if(plot){
    plot(F_)
  }
  
  return(F_)
}

DifferenceOfIntegral <- function(F_,G_){
  # F_ ... ecdf
  # G_ ... ecdf
  
  # print(knots(F_))
  # print(knots(G_))
  
  vals = sort(unique(c(knots(F_),knots(G_))))
  
  # print(vals)
  
  val = 0
  for(i in 1:(length(vals)-1)){
    val = val +abs(vals[i+1]-vals[i])*abs(F_(vals[i]) - G_(vals[i]))
  }
  
  return(val)
}

# willy
# DifferenceOfIntegral(F_ = all_protein_models[[ind2]]$distributions[[1]], G_ = all_protein_models[[ind2]]$distributions[[3]])
# DifferenceOfIntegral(F_ = all_protein_models[[ind2]]$distributions[[3]], G_ = all_protein_models[[ind2]]$distributions[[1]])


DistributionOfDifferenceOfIntegral <- function(points_a, points_b, n, m = 100){
  # points ... a content of pts-file
  # n ... number of points to sample
  #-------------------------------------
  
  vals = rep(0,m)
  for(i in 1:m){
    F_ = samplePointsAndCalculateCDFofEc(points_a,n)
    G_ = samplePointsAndCalculateCDFofEc(points_b,n)
    
    vals[i] = DifferenceOfIntegral(F_,G_)
  }

  return(vals)
}


checkTriangleInEq <- function(dxy,dyz,dxz,eps = 0.01){
  v = c(dxy,dyz,dxz)
  m = max(v)
  ind = which.maxn(v,n=1)
  
  # print(m)
  # print(ind)
  # 
  # sum(m[-ind])
    
  if(m <= sum(v[-ind]) + eps) return(TRUE)
  
  return(FALSE)
  # if(!(dxy <= dyz + dxz + eps)) return(FALSE)
  # if(!(dyz <= dxy + dxz + eps)) return(FALSE)
  # if(!(dxz <= dxy + dyz + eps)) return(FALSE)
  # 
  # return(TRUE)
}

# checkTriangleInEq(1,0,1)

checkIfTriangleHolds <- function(d_exact){
  
  for(i in 1:nrow(d_exact)){
    for(j in 1:nrow(d_exact)){
      for(k in 1:nrow(d_exact)){
        if(!checkTriangleInEq(d_exact[i,j], d_exact[i,k], d_exact[k,j])){
          print(paste(d_exact[i,j], d_exact[i,k], d_exact[k,j]))
          return(FALSE)
        }
      } 
    }
  }
  return(TRUE)
}



# install.packages("smacpod")


# t = circles.intersect(coords = rbind(xy_coords[1:3,]), r = c(d_1,d_2,d_3))

# install.packages("raster")

cirlceIntersect <- function(A.x, A.y, A.r, B.x, B.y, B.r, plot = FALSE, epsA = 0.01, epsB = 0.01, plotConst = 1){
  
  A.r = A.r + epsA
  B.r = B.r + epsB
  
  if(plot == TRUE){
    plot(A.x,A.y, xlim = c(min(A.x,B.x)-A.r-plotConst,max(A.x,B.x)+A.r+plotConst), ylim = c(min(A.y,B.y)-A.r-plotConst,max(A.y,B.y)+A.r+plotConst))
    # points(B.x,B.y)
    
    draw.circle(A.x,A.y,A.r,nv=100,border=NULL,col=NA,lty=1,density=NULL, angle=45,lwd=1)
    draw.circle(B.x,B.y,B.r,nv=100,border=NULL,col=NA,lty=1,density=NULL, angle=45,lwd=1)
    
    # ?draw.circle()
  }
  
  d = euclideanDistance(c(A.x,A.y),c(B.x,B.y))
  
  if (d <= A.r + B.r && d >= abs(B.r - A.r)) {
    
    ex = (B.x - A.x) / d
    ey = (B.y - A.y) / d
    
    # print(ex)
    # print(ey)
    # 
    x = (A.r * A.r - B.r * B.r + d * d) / (2 * d)
    y = sqrt(A.r * A.r - x * x)
    
    P1.x = A.x + x * ex - y * ey
    P1.y = A.y + x * ey + y * ex
    
    P2.x = A.x + x * ex + y * ey
    P2.y = A.y + x * ey - y * ex
    
    if(plot == TRUE){
      points(x =c(P1.x,P2.x), y = c(P1.y,P2.y), col = "green", cex = 1)
    }
    
    return(list("P1" = c(P1.x,P1.y), "P2" = c(P2.x,P2.y)))
    
  } else {
    # No Intersection, far outside or one circle within the other
    return(NULL)
  }
}
# c= 1
# p = cirlceIntersect(c,0,2, c+2,0,2)
# loc = cirlceIntersect(0.0192962892826289, 0.0207007353886564, 0.855259416940968, 0.0787789890704529, -0.144024990918186, 0.629172807132114, TRUE)
# loc

# no Solution example, one in other
# cirlceIntersect(A.x = 0, A.y = 0, A.r = 64.2580032348633, B.x = 1.28894641328613, B.y = 0, B.r = 62.3499984741211,plot = TRUE)


locatePoint <- function(d_1,d_2, d_3, point1, point2, point3, eps = 0.1, circleEps = 0.01){
  loc = cirlceIntersect(A.x = point1[1], A.y = point1[2], A.r = d_1, B.x = point2[1], B.y = point2[2], B.r = d_2,plot = FALSE, plotConst = 2, epsA = 0, epsB = 0)
  circleEpsA = 0
  circleEpsB = 0
  while(is.null(loc)){
    print(circleEpsA)
    print(circleEpsB)
    if(d_1 < d_2){
      circleEpsA = circleEpsA+ 0.01
    } else {
      circleEpsB = circleEpsB+ 0.01
    }
    
    loc = cirlceIntersect(A.x = point1[1], A.y = point1[2], A.r = d_1, B.x = point2[1], B.y = point2[2], B.r = d_2,plot = FALSE, plotConst = 2, epsA = circleEpsA, epsB = circleEpsB)
  }
  
  # points(x = loc$P1[1], y = loc$P1[2], col = "red")
  # points(x = loc$P2[1], y = loc$P2[2], col = "green")
  
  # print("points")
  # print(point3)
  # print(loc$P1)
  # print(loc$P2)
  # print(d_3)
  # print(euclideanDistance(loc$P1,point3))
  # print(euclideanDistance(loc$P2,point3))
  # 
  # 
  # print("loc")
  # print(min(abs(euclideanDistance(loc$P1,point3)-d_3), abs(euclideanDistance(loc$P2,point3)-d_3)))
  if(abs(euclideanDistance(loc$P1,point3)-d_3) < abs(euclideanDistance(loc$P2,point3)-d_3)){
    return(loc$P1)
  } else {
    return(loc$P2)
  } 
}


# d_1 = calcAllDistributionPairs(all_protein_models[[1]]$distributions,all_protein_models[[1]]$distributions)
# projection1 = myProjection(d_1, plot = FALSE)

# checkIfTriangleHolds(d_1)
# 
# m = matrix(rnorm(20),ncol = 2)
# d = as.matrix(dist(m,upper = TRUE, diag = TRUE))
# pro = myProjection(d)
# 
# plot(rbind(pro,m))
# points(m, col = "blue")

# project_geometry(all_protein_models = all_protein_models,ind1 = 1,ind2 = 20)


projectExample <- function(){
  matrix(rnorm(n))
  
  n = 20
  d = as.matrix(dist(matrix(rnorm(n*2), nrow = n),diag = TRUE, upper = TRUE))

  print(d)  
  xy_coords = cmdscale(d)
  plot(xy_coords)
  
  # new point
  p = matrix(rnorm(2),nrow=1)
  
  points(p, col = "red")
  
  d_1 = euclideanDistance(p,xy_coords[1,])
  d_2 = euclideanDistance(p,xy_coords[2,])
  d_3 = euclideanDistance(p,xy_coords[3,])
  
  p_predicted = locatePoint(d_1 = d_1, d_2 = d_2, d_3 = d_3, xy_coords[1,], xy_coords[2,], xy_coords[3,])
  points(x = p_predicted[1], y = p_predicted[2] , col = "green")
}

# projectExample()

calculateAllDistancesExact <- function(all_protein_models_with_distances,m){
  # Calculate all distances
  #-------------------------
  
  DE_Big_exact = matrix(0, nrow = length(all_protein_models_with_distances)*m, ncol = length(all_protein_models_with_distances)*m)
  
  rownames(DE_Big_exact) = rep("empty", length(all_protein_models_with_distances)*m)
  colnames(DE_Big_exact) = rep("empty", length(all_protein_models_with_distances)*m)
  
  for(i in 1:length(all_protein_models_with_distances)){
    for(j in 1:length(all_protein_models_with_distances)){
      print(paste(all_protein_models_with_distances[[i]]$model$name,all_protein_models_with_distances[[j]]$model$name))
      
      DE_exact = calcAllDistributionPairs(all_protein_models_with_distances[[i]]$model$distributions,all_protein_models_with_distances[[j]]$model$distributions)
      
      startI = ((i-1)*m+1)
      endI = startI+m-1
      startJ = ((j-1)*m+1)
      endJ = startJ+m-1
      DE_Big_exact[startI:endI,startJ:endJ] = DE_exact
      
      rownames(DE_Big_exact)[startI:endI] = rep(all_protein_models_with_distances[[i]]$model$name, m)
      colnames(DE_Big_exact)[startJ:endJ] = rep(all_protein_models_with_distances[[j]]$model$name, m)
    }
  }
  
  return(DE_Big_exact)
}




#-------------------------------------------------------------------------------
# n = 100
# m = 22
# times = 1
# OutputPath = "/home/willy/Schreibtisch/106Test/Output/"
# all_protein_models_with_distances = not_vectorized_get_allModels(OutputPath = OutputPath, n = n,m = m, times = times)
# DE_parallel = computeAllDistancesParallel(all_protein_models_with_distances)
# Emd_distance_matrix = emd_parallel(DE_parallel = DE_parallel,m = m)

#-------------------------------------------------------------------------------



computeAllDistancesParallel <- function(all_protein_models_with_distances,small = NULL){
  print("---------------------------------------------------")
  print("computing all pairwise DE of all distributions ...")
  
  allDistributions = list()
  
  names = list()
  
  for(i in 1:length(all_protein_models_with_distances)){
    allDistributions = list.append(allDistributions,all_protein_models_with_distances[[i]]$model$distributions)
    names = list.append(names,rep(all_protein_models_with_distances[[i]]$model$name, length(all_protein_models_with_distances[[i]]$model$distributions)))
  }
  allDistributions = do.call(c, unlist(allDistributions, recursive=FALSE))
  names = unlist(names, recursive=TRUE)

  num = length(names)
  
  if(!is.null(small)) num = small
  
  mat = foreach(i = 1:num, .combine = rbind) %:% foreach(j = 1:num) %do% {
    if(i>=j) return(as.numeric(0))
    DifferenceOfIntegral(allDistributions[[i]], allDistributions[[j]])
  }

  m <- mapply(as.matrix(mat), FUN=as.numeric)
  mat <- matrix(data=m, ncol=num, nrow=num)
  
  mat = mat + t(mat)
  
  # print(as.vector(names)[[1]])
  colnames(mat) = as.vector(names[1:num])
  rownames(mat) = as.vector(names[1:num])
  return(mat)
}

#-------------------------
# names = unique(colnames(DE_parallel))
# emds = matrix(0, nrow = length(names), ncol = length(names))
# rownames(emds) = names
# colnames(emds) = names
# emds
# 
# m = 3
# for(i in 1:length(names)){
#   for(j in 1:length(names)){
#     emds[i,j] = histDistParallel(DE_parallel,i,j,m)
#     # print(paste(i,j, emds[i,j]))
#   }
# }
# m = 3
#-------------------------

emd_parallel = function(DE_parallel, m){
  print("-------------------------------------------")
  print(paste("Calculating emd-matrix",sep =""))
  
  names = unique(colnames(DE_parallel))
  
  out = foreach(i = 1:length(names), .combine = rbind) %:% foreach(j = 1:length(names)) %do% {
    histDistParallel(DE_parallel,i,j,m)
  }
  
  mat <- mapply(as.matrix(out), FUN=as.numeric)
  mat <- matrix(data=mat, ncol=length(names), nrow=length(names), byrow = FALSE)
  
  colnames(mat) = as.vector(names)
  rownames(mat) = as.vector(names)
  
  return(mat)
}


histDistParallel <- function(DE_parallel, ind1, ind2,m){
  startI = (ind1-1)*m+1
  endI = startI-1+m
  
  startJ = (ind2-1)*m+1
  endJ = startJ-1+m
  indices1 = c(startI:endI)
  indices2 = c(startJ:endJ)
  
  vals1 = hist(DE_parallel[indices1,indices1], breaks = 100, plot = F)$counts
  vals3 = hist(DE_parallel[indices1,indices2], breaks = 100, plot = F)$counts
  
  P <- t(as.matrix(hist(vals1 , breaks = seq(from = 0,to = max(vals1,vals3)+0.05, by = 0.05), plot = F)$counts))
  Q <- t(as.matrix(hist(vals3 , breaks = seq(from = 0,to = max(vals1,vals3)+0.05, by = 0.05), plot = F)$counts))
  ed1 = emd2d(Q,P) 
  return(ed1)
}



# benchmark("parallel" = computeAllDistancesParallel(all_protein_models_with_distances),
#           "notParallel" = calculateAllDistancesExact(all_protein_models_with_distances,m=3),replications = 1,
#           columns = c("test", "replications", "elapsed",
#                       "relative", "user.self", "sys.self"))

# test replications elapsed relative user.self sys.self
# 2 notParallel            1  70.556    1.000    70.704    0.213
# 1    parallel            1  70.701    1.002    70.703    0.007




#---------------
# install.packages("RANN")
# library("RANN")
# x1 <- runif(100, 0, 2*pi)
# x2 <- runif(100, 0,3)
# DATA <- data.frame(x1, x2)
# m = matrix(c(1,1), nrow = 1)
# nearest <- nn2(data = DATA, query = m,k = 1)
# nearest$nn.idx
# 
# plot(DATA)
# points(m, col = "red")
# points(DATA[nearest$nn.idx,], col = "green")

# 
# 
# ?nn2
#---------------

# list = approximateNumberOfDistributions(points = pts_pos, n = 100, maxIt = 50, eps = 0.1, withKdTree = TRUE,plotOnTheFly = TRUE)

# ## nn 6
# length(list$distributions)
# list$distances[19,]
# 
# library(doBy)
# which.minn(list$distances[47,], n = 2)
# min(list$distances[47,-47])
# list$distances[47,18]
# which(sort(list$distances[47,]) == list$distances[47,14])




# max(list$distances)

# [1] "1.25550494618571 1.63984714880655 1.9232737790152"
# [,1]       [,2]
# [1,] -0.7304874 0.04022986
# [2,] -1.0554606 0.22199953
# [3,] -1.2497279 0.32619170

# m = matrix(c(-0.7304874, 0.04022986, -1.0554606, 0.22199953,-1.2497279, 0.32619170), nrow = 3)
# 
# locatePoint(1.25550494618571, 1.63984714880655, 1.9232737790152, m[1,], m[2,], m[3,])
# 
# 
# loc = cirlceIntersect(0.0192962892826289, 0.0207007353886564, 0.855259416940968, 0.0787789890704529, -0.144024990918186, 0.629172807132114)
# 

#-------------------------------------------------------

isThatMetric <- function(points, n = 1000, Trials = 1000, eps = 0.1){
  for(i in 1:Trials){
    print(i)
    X = samplePointsAndCalculateCDFofEc(pts_pos,n)
    Y = samplePointsAndCalculateCDFofEc(pts_pos,n)
    Z = samplePointsAndCalculateCDFofEc(pts_pos,n)
    
    dxy = DifferenceOfIntegral(X,Y)
    dyz = DifferenceOfIntegral(Y,Z)
    dxz = DifferenceOfIntegral(X,Z)
    
    if(!checkTriangleInEq(dxy,dyz,dxz, eps)) {
      print(paste(dxy, dyz, dxz))
      return(FALSE)
    }
  }
  return(TRUE)
}

# isThatMetric(pts_pos3)


#-------------------------------------------------------

precalCulateDistributionsModel <- function(name, pts, n, m){
  # n ... number of points
  # m ... number of samples
  #----------------------------
  
  distributions = list()
  for(i in 1:m){
    # print(i)
    model = samplePointsAndCalculateCDFofEc(pts,n)
    distributions[[i]] = model
  }

  return(list("name" = name, "distributions" = distributions))
}
# OutputPath = "/home/willy/Schreibtisch/106Test/Output/"
# pts = read_pts_file(OutputPath = OutputPath, protName = "013")
# precalCulateDistributionsModel("013",pts,100,500)


calcAllDistributionPairs <- function(distributions1, distributions2, asMatrix = TRUE){
  
  vals = matrix(0, nrow = length(distributions1), ncol = length(distributions2))
  
  for(i in 1:(length(distributions1))){
    for(j in (i):length(distributions2)){
      # print(paste(i,j))
      vals[i,j] =  DifferenceOfIntegral(distributions1[[i]],distributions2[[j]])
      vals[j,i] = vals[i,j]
    }
  }
  if(asMatrix) return(vals)
  
  return(as.vector(vals))
}

calcDistSampled <- function(distributions1, distributions2, m){
  
  inds_dist1 = sample(length(distributions1), size = m, replace = TRUE)
  inds_dist2 = sample(length(distributions2), size = m, replace = TRUE)
  
  vals = rep(0,m)
  for(i in 1:m){
    vals[i] = DifferenceOfIntegral(distributions1[[inds_dist1[i]]],distributions2[[inds_dist2[i]]])
  }
  
  return(as.vector(vals))
}


calcDistSampledSimple <- function(distributions1, distributions2){
  # just the pairs in the same order
  #----------------------------------

  vals = rep(0,length(distributions1))
  for(i in 1:length(distributions1)){
    vals[i] = DifferenceOfIntegral(distributions1[[i]],distributions2[[i]])
  }
  
  return(as.vector(vals))
}

getAllDistributions <- function(OutputPath, n = 100, m_distributions = 50, positive = TRUE){
  
  print(paste("Loading from ", OutputPath, sep =""))
  protein_names = list.dirs(OutputPath, recursive = FALSE, full.names = FALSE)
  
  distributions_lists = list()
  
  for(i in 1:length(protein_names)){
    # print(protein_names[i])
    pts = read_pts_file(OutputPath = OutputPath,protName = protein_names[i], pos = positive)
    model = precalCulateDistributionsModel(protein_names[i],pts, n, m_distributions)

    distributions_lists[[i]] = model
  }
  
  return(distributions_lists)
}


# all_protein_models = getAllDistributions(OutputPath = OutputPath,n = 100,m_distributions = 30)
# all_protein_models[[1]]$name

mirrorAlongLine <- function(a,b,c, x_1, y_1){
  t = -2*(a*x_1+b*y_1+c) / (a*a+b*b)
  x = t*a +x_1
  y = t*b +y_1
  return(c(x,y))
}

lineFromPoints <- function(P,Q){
  a = Q[2]-P[2]
  b = P[1]-Q[1]
  c = a*P[1]+b*P[2]
  
  return(list("a" = a, "b" = b, "c" = c))
}

# install.packages("ade4")
# library(ade4)
# 
# as.matrix(cailliez(as.dist(d_1)))

getDistancesAndPoints <- function(distributions_obj, docailiez = TRUE, doProjection = FALSE){
  name = distributions_obj$name
  DE = calcAllDistributionPairs(distributions_obj$distributions,distributions_obj$distributions)
  distances_1 = DE
  if(docailiez) {
    distances_1 = as.matrix(cailliez(as.dist(distances_1)))
  }
  # print(is.euclid(as.dist(distances_1)))
  projection1 = cmdscale(distances_1, add = TRUE, k = 2)$points
  rownames(distances_1) = rep(name,nrow(distances_1))
  colnames(distances_1) = rep(name,nrow(distances_1))
  rownames(projection1) = rownames(distances_1)
  
  return(list("DE" = DE, "distances" = distances_1, "projection" = projection1))
}

# getDistancesAndPoints(all_protein_models[[ind1]])

project_geometry <- function(all_protein_models, ind1, ind2, plot = FALSE){

  # R1 = getDistancesAndPoints(all_protein_models[[ind1]])
  # projection1 = R1$projection
  # distances_1 = R1$distances
  
  
  model1 = all_protein_models[[ind1]]$model
  projection1 = all_protein_models[[ind1]]$projection
  distances_1 = all_protein_models[[ind1]]$euclidDist
  
  # R2 = getDistancesAndPoints(all_protein_models[[ind2]])
  # projection2 = R2$projection
  # distances_2 = R2$distances
  
  model2 = all_protein_models[[ind2]]$model
  projection2 = all_protein_models[[ind2]]$projection
  distances_2 = all_protein_models[[ind2]]$euclidDist
  
  print(paste(model1$name, model2$name))
  
    
  p_1 = projection1[1,]
  p_2 = projection1[2,]
  p_3 = projection1[3,]
  
  target_origin_index = 1
  
  d_1 = DifferenceOfIntegral(F_ = model1$distributions[[1]], G_ = model2$distributions[[target_origin_index]])
  d_2 = DifferenceOfIntegral(F_ = model1$distributions[[2]], G_ = model2$distributions[[target_origin_index]])
  d_3 = DifferenceOfIntegral(F_ = model1$distributions[[3]], G_ = model2$distributions[[target_origin_index]])
  
  target_origin = locatePoint(d_1 = d_1, d_2 = d_2, d_3 = d_3, point1 = p_1, point2 = p_2, point3 = p_3)
  
  # print(target_origin)
  
  if(plot){
    plotProjections(list(projection1,projection2))
    points(x = target_origin[1], y = target_origin[2], col = "green", cex = 2)
    # points(x = projection2[target_origin_index,1], y = projection2[target_origin_index,2], col = "yellow", cex = 2)
    
  }
  
  print("hi")

  trans_vec = target_origin-projection2[1,]
  
  # print(target_origin)
  # print(trans_vec)
  
  
  projection2_translated = cbind(projection2[,1] + trans_vec[1], projection2[,2] + trans_vec[2])
  
  if(plot){
    plotProjections(list(projection1,projection2_translated))
  }

  rot_index = which.maxn(distances_2[target_origin_index,],n = 1)
  
  # print(rot_index)
  
  
  # now rotate accordingly
  p_4 = projection1[4,]
  p_5 = projection1[5,]
  p_6 = projection1[6,]
  
  d_4 = DifferenceOfIntegral(F_ = model1$distributions[[4]], G_ = model2$distributions[[rot_index]])
  d_5 = DifferenceOfIntegral(F_ = model1$distributions[[5]], G_ = model2$distributions[[rot_index]])
  d_6 = DifferenceOfIntegral(F_ = model1$distributions[[6]], G_ = model2$distributions[[rot_index]])
  
  target_rot = locatePoint(d_1 = d_4, d_2 = d_5, d_3 = d_6, point1 = p_4, point2 = p_5, point3 = p_6)
  
  # target_rot = target_rot + trans_vec
  # flipped 
  # target_rot = c(-target_rot[1], target_rot[2])
  
  line = lineFromPoints(P = projection1[1,],Q = projection1[2,])
  p_mirrored = mirrorAlongLine(a = line$a, b = line$b, c = line$c,x_1 = target_rot[1], y_1 = target_rot[2])
  print(p_mirrored)
  
  if(plot){
    plotProjections(list(projection1,projection2_translated))
    points(x = target_rot[1], y = target_rot[2], col ="green")
    points(x = projection2_translated[rot_index,1], y = projection2_translated[rot_index,2], col ="black", cex = 2)
    
    print(projection2_translated[rot_index,])
    print(target_rot)
    
    points(x = p_mirrored[1], y = p_mirrored[2], col ="green")
    points(x = target_origin[1], y = target_origin[2], col ="red", cex = 2)
    
    points(x = projection2[,1], y = projection2[,2], col ="blue")
  }

  
  d_origin_rot = euclideanDistance(target_rot, projection2_translated[target_origin_index,])
  d_origin_rot_Distributions = DifferenceOfIntegral(F_ = model2$distributions[[target_origin_index]], G_ = model2$distributions[[rot_index]])
  
  #----------------------------------
  # print("down")
  # print(d_origin_rot)
  # print(d_origin_rot_Distributions)
  # print(trans_vec)
  # print(distances_2[target_origin_index,rot_index])
  #----------------------------------
  
  eps = 0.7
  if(abs(euclideanDistance(target_rot, projection2_translated[rot_index,])) > eps 
     && abs(euclideanDistance(p_mirrored, projection2_translated[rot_index,])) > eps){
    # need to rotate
    print("need to rotate")
    # print(abs(euclideanDistance(target_rot, projection2_translated[rot_index,])))
    minDist = 10000
    
    print(target_rot)
    print(projection2_translated[rot_index,])
    
    p1 = target_rot - projection2_translated[target_origin_index,]
    p2 = projection2_translated[rot_index,] - projection2_translated[target_origin_index,]
    
    ang1 = acos((p1 %*% p2) / (norm_vec(p1) * norm_vec(p2)))
    print(ang1)
    
    p1 = p_mirrored- projection2_translated[target_origin_index,]
    ang2 = acos((p1 %*% p2) / (norm_vec(p1) * norm_vec(p2)))
    print(ang2)
    
    p_rotated = rotateByAngleAtPoint(projection2_translated, projection2_translated[target_origin_index,], -ang1)
    p_Mirror_rotated = rotateByAngleAtPoint(projection2_translated, projection2_translated[target_origin_index,], -ang2)
    
    v = c(abs(euclideanDistance(target_rot, p_rotated[rot_index,])), abs(euclideanDistance(p_mirrored, p_Mirror_rotated[rot_index,])))

    if(v[1] < v[2]){
      projection2_translated = p_rotated
    } else {
      projection2_translated = p_Mirror_rotated
    }
    # 
    # 
    # return()
    # for(i in seq(-3.2,3.2,0.1)){
    #   # print(i)
    #   p2 = rotateByAngleAtPoint(projection2_translated, projection2_translated[target_origin_index,], i)
    #   # points(p2, col = "green")
    #   # print(abs(euclideanDistance(target_rot, p2[rot_index,])))
    #   
    #   mind = min(abs(euclideanDistance(target_rot, p2[rot_index,])), abs(euclideanDistance(p_mirrored, p2[rot_index,])))
    #   
    #   if(minDist > mind) minDist = mind
    #   
    #   if(plot){
    #     plotProjections(list(projection1,p2))
    #     points(x = target_rot[1], y = target_rot[2], col ="green")
    #     points(x = p2[rot_index,1], y = p2[rot_index,2], col ="black", cex = 2)
    #     
    #     # print(p2[rot_index,])
    #     # print(target_rot)
    #     
    #     points(x = p_mirrored[1], y = p_mirrored[2], col ="green")
    #     points(x = target_origin[1], y = target_origin[2], col ="red", cex = 2)
    #     
    #     # points(x = projection2[,1], y = projection2[,2], col ="blue")
    #   }
    #   
    #   if(mind < eps){
    #     print("found good rotation")
    #     print(minDist)
    #     projection2_translated = p2
    #     
    #     print(paste(ang1,ang2, i))
    #     break;
    #   }
    # }
  }

  
  ret = rbind(projection1,projection2_translated)
  
  return(ret)
}

# install.packages("Directional")
# library("Directional")
# rot.matrix(ksi, theta, rads = FALSE)


# willey
# p = project_geometry(all_protein_models = all_protein_models,ind1 = 1,ind2 = 2)
# p2 = project_geometry(all_protein_models = all_protein_models,ind1 = 1,ind2 = 4)

preCalculateAllDistances <- function(all_protein_models,docailiez = TRUE, doProjection = TRUE){
  all_proteins_out = list()
  
  projection1 = c()
  distances_1 = c()
  Err = c()
  
  for(i in 1:length(all_protein_models)){
    R1 = getDistancesAndPoints(all_protein_models[[i]],docailiez = docailiez, doProjection)

    if(doProjection) {
      projection1 = R1$projection
      distances_1 = R1$distances
      Err = R1$DE - distances_1
    }
    
    all_proteins_out[[i]] = list("model" = all_protein_models[[i]], "euclidDist" = distances_1, "projection" = projection1, "DE" = R1$DE, "Error" = Err)
  }

  return(all_proteins_out)
}

# all_protein_models_with_distances = preCalculateAllDistances(all_protein_models, docailiez = TRUE)

# all_protein_models_with_distances = preCalculateAllDistances(all_protein_models, docailiez = FALSE)


# all_protein_models_with_distances = preCalculateAllDistances(all_protein_models, docailiez = FALSE)

mergeProjections <- function(projections, mergeInd, number){
  # ret = projections[[mergeInd]]
  ret = c()
  
  for(i in 1:length(projections)){
    # if(mergeInd != i) {
      ret = rbind(projections[[i]], ret)
    # }
  }
  return(ret)
}



not_vectorized_get_allModels <- function(OutputPath, n,m,times, positive = TRUE){
  all_protein_models_with_distances = c()
  for(i in 1:times){
    all_protein_models1 = getAllDistributions(OutputPath = OutputPath,n = n,m_distributions = m, positive = TRUE)
    all_protein_models_with_distances1 = preCalculateAllDistances(all_protein_models1, docailiez = TRUE)
  
    all_protein_models_with_distances = c(all_protein_models_with_distances, all_protein_models_with_distances1)
  
    # all_protein_models1 = getAllDistributions(OutputPath = OutputPath,n = 100,m_distributions = m, positive = FALSE)
    # all_protein_models_with_distances1 = preCalculateAllDistances(all_protein_models1, docailiez = TRUE)
    #
    # all_protein_models_with_distances = c(all_protein_models_with_distances, all_protein_models_with_distances1)
  }
  
  return(all_protein_models_with_distances)
}

vectorized_get_allModels <- function(n,m,times, positive = TRUE, MC = TRUE){
  if(MC){
    return(mclapply(c(1:times),FUN = function(x){
      all_protein_models1 = getAllDistributions(OutputPath = OutputPath,n = n,m_distributions = m, positive = positive)
      all_protein_models_with_distances1 = preCalculateAllDistances(all_protein_models1, docailiez = TRUE)
      return(all_protein_models_with_distances1)
    }))
  }else{
    return(sapply(c(1:times),FUN = function(x){
      all_protein_models1 = getAllDistributions(OutputPath = OutputPath,n = n,m_distributions = m, positive = positive)
      all_protein_models_with_distances1 = preCalculateAllDistances(all_protein_models1, docailiez = TRUE)
      return(all_protein_models_with_distances1)
    }))
  }

}

#----------------------------------
# 
# times = 1
# # all_protein_models_with_distances2 = vectorized_get_allModels(n = 100, m= 10, times = times, MC = TRUE)
# 
# all_protein_models_with_distances = rep(0,m*times)
# all_protein_models_with_distances = list()
# for(i in 1:times){
#   all_protein_models_with_distances= c(all_protein_models_with_distances, all_protein_models_with_distances2[[i]])
# }
# 
# all_protein_models_with_distances = not_vectorized_get_allModels(n = 100,m = 3, times = times)
# 
# DE_big = calculateAllDistancesExact(all_protein_models_with_distances,m=3)
# 
# (DE_parallel - DE_big)[1:6,1:6]
# cor(c(DE_parallel),c(DE_big))
# 
# is.numeric(DE_big)
# is.numeric(q)
# q
# 
# 
# names = unique(colnames(DE_big))
# emds = matrix(0, nrow = length(names), ncol = length(names))
# rownames(emds) = names
# colnames(emds) = names
# emds
# 
# for(i in 1:length(names)){
#   for(j in 1:length(names)){
#     print(paste(i,j))
#     emds[i,j] = histDist2(DE_big,names[i], names[j])
#   }
# }
# 
# felix = readDistanceMatrix2()
# 
# cor(c(emds[-106,-106]),c(felix))
# 
# 
# # devtools::install_github("eddelbuettel/rbenchmark")
# library(rbenchmark)
# 
# benchmark("mclapply" = vectorized_get_allModels(n = 10, m= 10, 5, MC = TRUE),
#           "sapply" = vectorized_get_allModels(n = 10, m= 10, 5, MC = FALSE),replications = 1,
#           columns = c("test", "replications", "elapsed",
#                       "relative", "user.self", "sys.self"))
# 
# # test replications elapsed relative user.self sys.self
# # 1 mclapply            1  17.301    1.000     0.060    0.121
# # 2   sapply            1  28.303    1.636    28.187    0.228
# 
# 
# rotateByAngleAtPoint <- function(points, center, angle){
#   p2 = cbind(points[,1] - center[1], points[,2] - center[2])
#   p2 = rotateAroundZ(cbind(p2,0), angle)[,1:2]
#   
#   p2 = cbind(p2[,1] + center[1], p2[,2] + center[2])
#   return(p2)
# }
# 
# mergeTheseIndices <- function(mergeInd, otherIndices, all_protein_models_with_distances, m, xli = NULL, yli = NULL, plotPairs = FALSE){
#   ps2 = list()
#   for(i in 1:length(otherIndices)){
#     ps = project_geometry(all_protein_models = all_protein_models_with_distances,ind1 = mergeInd,ind2 = otherIndices[i])
#     ps2[[i]] = ps
#     if(plotPairs) plotPrettyProjection(ps,getFunctionalProteins(), onlyGeomCenters = FALSE, xli = xli, yli = yli)
#   }
#   
#   p_mergeds = mergeProjections(ps2, mergeInd = mergeInd, number = m)
#   plotPrettyProjection(p_mergeds,getFunctionalProteins(), onlyGeomCenters = FALSE, xli = xli, yli = yli)
# }
# 
# mergeTheseIndices(3,c(1,3,4,8,20,30),all_protein_models_with_distances, m)
# mergeTheseIndices(1,c(1,3,4,8,20,30),all_protein_models_with_distances, m)
# mergeTheseIndices(20,c(1,3,4,8,20,30),all_protein_models_with_distances, m)
# 
# 
# ps = project_geometry(all_protein_models = all_protein_models_with_distances,ind1 = 1,ind2 = 40)
# ps2 = project_geometry(all_protein_models = all_protein_models_with_distances,ind1 = 1,ind2 = 50)
# ps3 = project_geometry(all_protein_models = all_protein_models_with_distances,ind1 = 1,ind2 = 97)
# p_mergeds = mergeProjections(list(ps,ps2,ps3), mergeInd = 1, number = m)
# plotPrettyProjection(ps,getFunctionalProteins(), onlyGeomCenters = FALSE)
# plotPrettyProjection(ps2,getFunctionalProteins(), onlyGeomCenters = FALSE)
# 
# plotPrettyProjection(p_mergeds,getFunctionalProteins(), onlyGeomCenters = FALSE)
# 
# all_protein_models_with_distances[[1]]
# 
# mergeInd = 7
# projections_list = list()
# for(i in 1:length(all_protein_models_with_distances)){
#   # for(j in 1:length(all_protein_models_with_distances[[i]])){
#     print(i)
#     # if(i != mergeInd){
#     projections_list[[i]] = project_geometry(all_protein_models = all_protein_models_with_distances, ind1 = mergeInd, ind2 = i)[(m+1):(2*m),]
#     # }
#   }
# # }
# 
# p_merged = mergeProjections(projections_list, mergeInd = mergeInd, number = m)
# 
# geoms = plotPrettyProjection(p_merged,getFunctionalProteins(), onlyGeomCenters = FALSE)
# 
# geoDists = as.matrix(dist(geoms))
# neighbors = 20
# closest = names(geoDists[which.minn(geoDists[which(rownames(geoDists) == "000_Trx"),],n = neighbors),1])
# length(which(closest %in% getFunctionalProteins() == TRUE))
# 
# 
# 
# all_protein_models_with_distances[[mergeInd]]$model$distributions
# 
# calcErrorsToInd <- function(all_protein_models_with_distances, mergeInd){
#   meanError = rep(0,length(all_protein_models_with_distances))
#   meanVarError = rep(0,length(all_protein_models_with_distances))
#   
#   meanDistances = rep(0,length(all_protein_models_with_distances))
#   meanVar = rep(0,length(all_protein_models_with_distances))
#   
#   for(i in 1:length(all_protein_models_with_distances)){
#     print(i)
#     DE_exact = calcAllDistributionPairs(all_protein_models_with_distances[[i]]$model$distributions,all_protein_models_with_distances[[mergeInd]]$model$distributions)
#     
#     dists = as.matrix(dist(projections_list[[i]]))
#     # DE_exact
#     
#     Error = DE_exact - dists
#     meanError[i] = mean(Error)
#     meanVarError[i] = var(as.vector(DE_exact - dists))
#     
#     meanDistances[i] = mean(DE_exact)
#     meanVar[i] = var(DE_exact)
#   }
#   
#   mE = mean(meanError)
#   mVE = mean(meanVarError)
#   
#   mD = mean(meanDistances)
#   mV = mean(meanVar)
#   
#   hist(meanError)
#   hist(meanDistances)
#   
#   return(list("mean(meanError)" = mE, "mean(meanVarError)" =  mVE, "mean(meanDistances)" = mD, "mean(meanVar)" = mV))
# }
# 
# 
# calcErrorsToInd(all_protein_models_with_distances, 7)
# 
# deletePointsWithHighError <- function(all_protein_models_with_distances, th = 10){
#   m2 = length(all_protein_models_with_distances[[2]])*2
#   print(m2)
#   
#   Errors = matrix(0,ncol = m2, nrow = length(all_protein_models_with_distances)*m2)
#   for(i in 1:length(all_protein_models_with_distances)){
#     print(i)
#     DE_exact = calcAllDistributionPairs(all_protein_models_with_distances[[i]]$model$distributions,all_protein_models_with_distances[[mergeInd]]$model$distributions)
#     
#     dists = as.matrix(dist(projections_list[[i]]))
#     # DE_exact
#     
#     print(dim(dists))
#     print(dim(DE_exact))
#     
#     Error = DE_exact - dists
#     
#     start = ((i-1)*m2+1)
#     end = start+m2-1
#     
#     print(paste(start,end))
#     Errors[start:end, ] = Error
#   }
#   return(Errors)
# }
# 
# Errors = deletePointsWithHighError(all_protein_models_with_distances)
# 
# nrow(Errors)
# 
# th = 1
# mean(Errors[which(rowMeans(abs(Errors)) <= th), ])
# length(which(rowMeans(abs(Errors)) <= th)) 
# p_merged_cleaned = p_merged[-which(rowMeans(abs(Errors)) > th), ]
# 
# p_merged_cleaned
# 
# calcErrorsToInd(all_protein_models_with_distances, 66)
# 
# 
# dist_projected = as.matrix(dist(p_merged))
# dist_projected
# 
# unique(colnames(dist_projected))
# 
# 
# 
# 
# 
# 
# Error_big = DE_Big_exact- as.matrix(dist(p_merged))
# 
# mean(Error_big)
# var(as.vector(Error_big))
# 
# var(as.vector(DE_Big_exact))
# 
# 
# # dist(p_merged[1:10,])
# 
# 
# 
# 
# geoms = plotPrettyProjection(p_merged_cleaned,getFunctionalProteins(), onlyGeomCenters = FALSE)
# 
# geoDists = as.matrix(dist(geoms))
# neighbors = 15
# closest = names(geoDists[which.minn(geoDists[which(rownames(geoDists) == "000_Trx"),],n = neighbors),1])
# length(which(closest %in% getFunctionalProteins() == TRUE))
# 
# 
# alldist_vs_all = as.matrix(dist(p_merged))
# 
# names = unique(colnames(alldist_vs_all))
# emds = matrix(0, nrow = length(names), ncol = length(names))
# rownames(emds) = names
# colnames(emds) = names
# emds
# 
# for(i in 1:length(names)){
#   for(j in 1:length(names)){
#     print(paste(i,j))
#     emds[i,j] = histDist2(alldist_vs_all,names[i], names[j])
#   }
# }
# 
# 
# length(which(colnames(emds)[which.minn(emds[which(rownames(emds) == "000_Trx"),], n = 20)] %in% c(getFunctionalProteins(), "000_Trx")) == TRUE)
# 
# 
# # emds_proj = cmdscale(as.matrix(cailliez(as.dist(emds))))
# # plot(emds_proj)
# # inds = which(rownames(emds_proj) %in% c(getFunctionalProteins(), "000_Trx") == TRUE)
# # points(x = emds_proj[inds,1], y = emds_proj[inds,2], col = "green")
# 
# 
# length(intersect(rownames(demd2), rownames(emds)))
# demd2 = readDistanceMatrix2()
# 
# demd2 = demd2[order(rownames(demd2)),order(rownames(demd2))]
# min(demd2)
# 
# emds2 = emds[order(rownames(emds)),order(rownames(emds))]
# emds2 = emds2[,-106]
# emds2 = emds2[-106,]
# abs(emds2 -demd2)
# 
# heatmap(cor(emds2,demd2))
# cor(c(emds2), c(demd2))
# 
# 
# # d1 = readDistanceMatrix3("/home/willy/Schreibtisch/106Test/RepSubOutput/EMD_10_5000_1.0000_0.0000_0.0000_0.1000_id_opt_NNact_0.csv")
# # d2 = readDistanceMatrix3("/home/willy/Schreibtisch/106Test/RepSubOutput/EMD_50_5000_1.0000_0.0000_1.0000_1.0000_id_opt_NNact_0.csv")
# # 
# # cor(c(d1[-106,-106]), c(demd2))
# 
# 
# library(emdist)
# histDist2 <- function(alldist_vs_all, name1,name2){
#   ind1 = which(rownames(alldist_vs_all) == name1)
#   vals1 = hist(alldist_vs_all[ind1,ind1], breaks = 100, plot = F)$counts
#   
#   ind2 = which(colnames(alldist_vs_all) == name2)
#   # vals2 = hist(alldist_vs_all[ind2,ind2], breaks = 100)$counts
# 
#   vals3 = hist(alldist_vs_all[ind1,ind2], breaks = 100, plot = F)$counts
#   
#   P <- t(as.matrix(hist(vals1 , breaks = seq(from = 0,to = max(vals1,vals3)+0.05, by = 0.05), plot = F)$counts))
#   Q <- t(as.matrix(hist(vals3 , breaks = seq(from = 0,to = max(vals1,vals3)+0.05, by = 0.05), plot = F)$counts))
#   ed1 = emd2d(Q,P) 
#   
#   # P <- t(as.matrix(hist(vals2 , breaks = seq(from = 0,to = max(vals2,vals3)+0.05, by = 0.05), plot = F)$counts))
#   # Q <- t(as.matrix(hist(vals3 , breaks = seq(from = 0,to = max(vals2,vals3)+0.05, by = 0.05), plot = F)$counts))
#   # ed2 = emd2d(Q,P) 
#   
#   return(ed1)
# }
# 
# 
# 
# 
# 
# plotProjections <- function(projections, const = 6,xli = NULL, yli = NULL){
#   
#   x_ = matrix(0,nrow = length(projections), ncol = nrow(projections[[1]]))
#   y_ = matrix(0,nrow = length(projections), ncol = nrow(projections[[1]]))
#   for(i in 1:length(projections)){
#    x_[i,] =  projections[[i]][,1]
#    y_[i,] =  projections[[i]][,2]
#   }
#   
#   xli = c(min(x_)-const,max(x_)+const)
#   yli = c(min(y_)-const,max(y_)+const)
#   
#   plot(projections[[1]], col = "red", xlim = xli, ylim = yli)
#   if(length(projections) > 1){
#     for(i in 2:length(projections)){
#       points(projections[[i]], col = "blue", xlim = xli, ylim = yli)
#     }
#   }
# 
# 
# }
# 
# 
# 
# 
# 
# 
# 
# all_protein_models = getAllDistributions(OutputPath = OutputPath,n = 100,m_distributions = 10)
# ind1 = 1
# ind2 = 7
# d_1 = calcAllDistributionPairs(all_protein_models[[ind1]]$distributions,all_protein_models[[ind1]]$distributions)
# projection1 = myProjection(d_1)
# rownames(projection1) = rep(all_protein_models[[ind1]]$name,nrow(projection1))
# 
# 
# d_2 = calcAllDistributionPairs(all_protein_models[[ind2]]$distributions,all_protein_models[[ind2]]$distributions)
# projection2 = myProjection(d_2)
# rownames(projection2) = rep(all_protein_models[[ind2]]$name,nrow(projection2))
# 
# p_1 = projection1[1,]
# p_2 = projection1[2,]
# p_3 = projection1[3,]
# 
# target_origin_index = 1
# 
# d_1_1 = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[1]], G_ = all_protein_models[[ind2]]$distributions[[target_origin_index]])
# d_1_2 = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[1]], G_ = all_protein_models[[ind2]]$distributions[[2]])
# d_1_3 = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[1]], G_ = all_protein_models[[ind2]]$distributions[[3]])
# 
# plotProjections(list(projection1,projection2))
# 
# 
# 
# 
# 
# pythagoras <-function(distance, target_row, target_column, help_column){
#   # distance between target_row,target_column
#   # we know the distance between (target_row,help_column)  and (target_column,help_column)
#   
#   dist = sqrt(distance[target_column,help_column]^2 + distance[target_row,help_column]^2)
# 
#   return(dist)
# }
# 
# 
# #-----------------------------------------------------------------------------------
# 
# getSharedDistances <- function(d_1_1, d_1_2, d_1_3,projection1,projection2, distances1 = NULL, distances2 = NULL){
#   # d_1_1 ... distance of 1 to 1
#   #---------------------------------
#   n1 = nrow(projection1)
#   n2 = nrow(projection2)
#   
#   dim = n1 + n2
#   distance_shared = matrix(0, nrow = dim, ncol = dim)
#   
#   if(is.null(distances1)){
#     distances1 = as.matrix(dist(projection1, upper = TRUE, diag = TRUE))
#   }
#   
#   if(is.null(distances2)){
#     distances2 = as.matrix(dist(projection2, upper = TRUE, diag = TRUE))
#   }
#   
#   print(checkIfTriangleHolds(distances1))
#   print(checkIfTriangleHolds(distances2))
#   print(isSymmetric(distances1))
#   print(isSymmetric(distances2))
#   
#   distance_shared[1:n1,1:n1] = distances1
#   distance_shared[(n1+1):dim,(n1+1):dim] = distances2
#   
#   
#   distance_shared[1,1+n1] = d_1_1
#   distance_shared[1+n1,1] = d_1_1
#   
#   distance_shared[1,2+n1] = d_1_2
#   distance_shared[2+n1,1] = d_1_2
#   
#   distance_shared[1,3+n1] = d_1_3
#   distance_shared[3+n1,1] = d_1_3
#   
#   for(j in 4:(n1)){
#     start = projection1[1,]
#     
#     p1 = projection2[1,]
#     p2 = projection2[2,]
#     p3 = projection2[3,]
#     
#     d1 = distance_shared[1+n1,j+n1]
#     d2 = distance_shared[2+n1,j+n1]
#     d3 = distance_shared[3+n1,j+n1]
#       
#     loc = locatePoint(d_1 = d1, d_2 = d2, d_3 = d3, point1 = p1, point2 = p2, point3 = p3)
#     
#     distance_shared[1,j] = euclideanDistance(loc,start)
#   }
#   
#   print(distance_shared)
#   # for(j in 2:n2){
#   #   for(i in 1:n1){
#   #     # pythagoras
#   #     dist = pythagoras(distance = distance_shared, target_row = i, target_column = j+n1, help_column = j+n1-1)
#   #     distance_shared[i,j+n1] = dist
#   #     distance_shared[j+n1,i] = dist
#   #   }
#   # }
#   
#   colnames(distance_shared) = c(rownames(projection1),rownames(projection2))
#   rownames(distance_shared) = c(rownames(projection1),rownames(projection2))
#   
#   # print(distance_shared)
#   
#   # return(distance_shared)
#   
#   print(isSymmetric(distance_shared))
#   print(checkIfTriangleHolds(distance_shared))
#   
#   projection_shared = myProjection(distance_shared)
#   rownames(projection_shared) = c(rownames(projection1),rownames(projection2))
#   
#   projection1_s = projection_shared[1:n1,]
#   projection2_s = projection_shared[(n1+1):dim,]
#   
#   return(list("distance_shared" = distance_shared, "projection1_s" = projection1_s, "projection2_s" = projection2_s))
# }
# 
# shared = getSharedDistances(d_1_1, d_1_2, d_1_3, projection1 = projection1, projection2)
# 
# 
# 
# DifferenceOfIntegral(all_protein_models[[1]]$distributions[[2]], all_protein_models[[7]]$distributions[[1]])
# 
# shared[1,5]+shared[5,6]
# 
# sqrt((shared[1,5])^2 + (shared[5,6])^2)
# 
# shared
# 
# distances = shared
# 
# for(k in 1:nrow(distances)){
#   for(i in 1:nrow(distances)){
#     for(j in 1:ncol(distances)){
#       
#       d1 = distances[i,k]
#       d2 = distances[i,j] 
#       d3 = distances[j,k]
#       
#       vec= c(d1,d2,d3)
#       
#       v = which.maxn(vec, n = 1)
# 
#       
#       if(!(vec[v] == sqrt((vec[-v][1])^2 - (vec[-v][1])^2))){
#         print("no")
#         print(paste(i,j,k))
#         
#         print(distances[i,k])
#         print(sqrt((distances[i,j])^2 + (distances[j,k])^2))
#         print(paste("d[",i,",",k,"] != sqrt(d[",i,",",j,"]^2 + d[",j,",",k,"])^2)",sep = ""  ))
#         print(paste(d[i,k], d[i,j], d[j,k]))
#       }
#     }
#   }
# }
# 
# distances
# 
# 
# # 0.519230772459262 +7.6700944341509
# #8.38045636368756
# 
# 
# shared$distance_shared
# 
# plotProjections(list(shared$projection1_s,shared$projection2_s))
# 
# plotProjections(list(projection1,projection2))
# 
# 
# # check if distances correct
# dim = length(all_protein_models[[ind1]]$distributions) + length(all_protein_models[[ind2]]$distributions)
# d_exact = matrix(0,nrow = dim,ncol = dim)
# 
# rownames(d_exact) = c(rownames(projection1), rownames(projection2))
# colnames(d_exact) = c(rownames(projection1), rownames(projection2))
# 
# n1 = length(all_protein_models[[ind1]]$distributions)
# n2 = length(all_protein_models[[ind2]]$distributions)
# dim = length(all_protein_models[[ind1]]$distributions) + length(all_protein_models[[ind2]]$distributions)
# for(i in 1:dim){
#   for(j in 1:dim){
#     if(i <= length(all_protein_models[[ind1]]$distributions) && j <= length(all_protein_models[[ind1]]$distributions)){
#       d_exact[i,j] = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[i]], G_ = all_protein_models[[ind1]]$distributions[[j]])
#     }  else if(i > n1 && j > n1){
#       print(paste(i-n1,j-n1))
#       d_exact[i,j] = DifferenceOfIntegral(F_ = all_protein_models[[ind2]]$distributions[[i-n1]], G_ = all_protein_models[[ind2]]$distributions[[j-n1]])
#     } else if(i > n1 && j <= n1){
#       # print(paste(i-n1,j-n1))
#       d_exact[i,j] = DifferenceOfIntegral(F_ = all_protein_models[[ind2]]$distributions[[i-n1]], G_ = all_protein_models[[ind1]]$distributions[[j]])
#     } else if(i <= n1 && j > n1){
#       # print(paste(i-n1,j-n1))
#       d_exact[i,j] = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[i]], G_ = all_protein_models[[ind2]]$distributions[[j-n1]])
#     }
#     
#     # print(paste(all_protein_models[[ind1]]$name, all_protein_models[[ind2]]$name))
#     # print(colnames(shared$distance_shared[,c(i,length(all_protein_models[[ind1]]$distributions)+j)]))
#     # 
#     # print(d_1 - shared$distance_shared[i,length(all_protein_models[[ind1]]$distributions)+j])
#   }
# }
# 
# d_exact
# 
# projection = myProjection(d_exact)
# 
# myProjection <- function(distances, plot = FALSE){
# 
#   points = matrix(0,nrow = nrow(distances), ncol = 2)
#   points[1,] = c(0,0)
#   # second_ind = which.maxn(distances[1,],n = 1)
#   # second_ind = 2
#     
#   points[2,] = c(distances[1,2],0)
#   
#   p = cirlceIntersect(A.x = points[1,1], A.y = points[1,2],A.r = distances[1,3],
#                   B.x = points[2,1], B.y = points[2,2],B.r = distances[2,3])
#   
#   points[3,] = p$P1
#   
#   if(plot){
#     plot(points[c(1,2,3),], xlim = c(-distances[1,2], distances[1,2]), ylim = c(-distances[1,2], distances[1,2]))
#   }
# 
#   for(i in 4:nrow(points)){
#     # print(i)
#       print(paste(distances[1,i], distances[2,i], points[1,1], points[1,2], points[2,1], points[2,2], points[3,1], points[3,2]))
#       q = locatePoint(distances[i,1], distances[i,2], distances[i,3], points[1,], points[2,], points[3,])
#       print(q)
#       points[i,] = q
#   }
#   
#   if(plot){
#     points(x = points[,1], y = points[,2])
#   }
#   
#   rownames(points) = rownames(distances)
#   
#   return(points)
# }
# 
# d_1 = calcAllDistributionPairs(all_protein_models[[1]]$distributions,all_protein_models[[1]]$distributions)
# projection1 = myProjection(d_1, plot = TRUE)
# 
# 
# # n = 400
# # ps = matrix(rnorm(n*2), ncol = 2)
# # d_example = as.matrix(dist(ps,upper = TRUE,diag = TRUE))
# # # checkIfTriangleHolds(d_example)
# # 
# # myProjection(d_example)
# 
# 
# plotPrettyProjection <- function(projection, functionals = NULL, onlyGeomCenters = FALSE, constPlot = 1, fontSize = 1, cex = 1, xli = NULL, yli = NULL){
#   
#   if(is.null(xli)) xli = c(min(projection[,1])-constPlot, max(projection[,1])+constPlot)
#   if(is.null(yli)) yli = c(min(projection[,2])-constPlot, max(projection[,2])+constPlot)
#   
#   if(onlyGeomCenters == FALSE){
#     plot(projection, col = "green", xlim = xli, ylim = yli, pch = 19, cex = cex, ylab = "DE", xlab = "DE")
#   } 
#   else {
#     plot(NULL, col = "green", xlim = xli, ylim = yli, pch = 19, cex = cex, ylab = "DE", xlab = "DE")
#   }
# 
#   
#   names = unique(rownames(projection))
#   geoms = matrix(0,nrow = length(names), ncol = 2)
#   for(i in 1:length(names)){
#     inds = which(rownames(projection) == names[i])
#     
#     geom = c(sum(projection[inds,1]), sum(projection[inds,2]))/length(inds)
#     
#     functional = FALSE
#     if(names[i] %in% functionals){
#       functional = TRUE
#     }
#     if(functional) {
#       text(x = geom[1], y = geom[2], labels = c(names[i]), col = "blue", cex = fontSize)
#       if(onlyGeomCenters == FALSE) points(projection[inds,], add = TRUE, col = "blue")
#     }
#     else{
#       text(x = geom[1], y = geom[2], labels = c(names[i]), col = "red", cex = fontSize)
#     } 
#     
#     geoms[i,] = geom
#   }
#   
#   rownames(geoms) = names
#   
#   return(geoms)
# }
# 
# 
# plotPrettyProjection(projection = projection, constPlot = 3, fontSize = 1, cex = 0.5)
# 
#          
# histDist <- function(model, model2){
#   print(paste(model$name,model2$name))
#   
#   vals_AA = calcDistSampled(model$distributions,model$distributions, 50)
#   vals_AB = calcDistSampled(model$distributions,model2$distributions, 50)
#   # vals_BB = calcDistSampled(all_protein_models[[j]]$distributions,all_protein_models[[j]]$distributions, 500)
#   
#   P <- t(as.matrix(hist(vals_AA , breaks = seq(from = 0,to = max(vals_AA,vals_AB)+0.05, by = 0.05), plot = F)$counts))
#   Q <- t(as.matrix(hist(vals_AB , breaks = seq(from = 0,to = max(vals_AA,vals_AB)+0.05, by = 0.05), plot = F)$counts))
#   return(emd2d(Q,P))
# }
# 
# calculateAllPairWiseEmds <- function(all_protein_models,m = 500){
#   
#   d = matrix(0,nrow = length(all_protein_models), ncol= length(all_protein_models))
#   
#   for(i in 1:length(all_protein_models)){
#     vals_AA = calcDistSampled(all_protein_models[[i]]$distributions,all_protein_models[[i]]$distributions, m)
#     for(j in 1:length(all_protein_models)){
#       print(paste(i,j))
#       vals_AB = calcDistSampled(all_protein_models[[i]]$distributions,all_protein_models[[j]]$distributions, m)
#       # vals_BB = calcDistSampled(all_protein_models[[j]]$distributions,all_protein_models[[j]]$distributions, 500)
#       
#       P <- t(as.matrix(hist(vals_AA , breaks = seq(from = 0,to = max(vals_AA,vals_AB)+0.05, by = 0.05), plot = F)$counts))
#       Q <- t(as.matrix(hist(vals_AB , breaks = seq(from = 0,to = max(vals_AA,vals_AB)+0.05, by = 0.05), plot = F)$counts))
#       d[i,j] = emd2d(Q,P) 
#     }
#   }
# 
#   df = as.data.frame(d)
#   return(df)
# }
# 
# df = calculateAllPairWiseEmds(all_protein_models,m = 50)
# 
# getFunctionalProteins()
# 
# hist(as.numeric(df[7,]),breaks = seq(from = 0, to = 160, by =0.05))
# df[7,20]
# df[7,23]
# df[7,94]
# 
# all_protein_models[[99]]$name
# which.minn(as.numeric(df[7,]), n = 20)
# 
# 
# 
# histDist(all_protein_models[[20]], all_protein_models[[7]])
# 
# 
# 
# 
# mylist <- list(a=1,b=2,c=3)
# myfxn <- function(var1,var2){
#   var1*var2
# }
# var2 <- 2
# 
# sapply(mylist,myfxn,var2=var2)
# 
# mapply(rep, times = 1:4, MoreArgs = list(x = 42))
# 
# mapply(histDist,all_protein_models,all_protein_models)
# 
# all_protein_models[1:2][1]
# 
# 
# distributions1 = precalCulateDistributionsModel(pts_pos, 100, 500)
# distributions2 = precalCulateDistributionsModel(pts_pos, 100, 50)
# 
# vals2 = calcAllDistributionPairs(all_protein_models[[7]]$distributions,all_protein_models[[1]]$distributions)
# vals = calcAllDistributionPairs(all_protein_models[[7]]$distributions,all_protein_models[[7]]$distributions)
# vals3 = calcAllDistributionPairs(all_protein_models[[20]]$distributions,all_protein_models[[7]]$distributions)
# 
# vals2 = calcDistSampled(all_protein_models[[7]]$distributions,all_protein_models[[10]]$distributions, 2500)
# 
# 
# hist(vals, xlim = c(0,max(vals2,vals,vals3)), col = "red", breaks = 100)
# integral2 = hist(vals2, add = TRUE, col = "blue", breaks = 100)
# hist(vals3, add = TRUE, col = "yellow", breaks = 100)
# 
# integral = hist(vals, breaks = 100)
# integral$breaks
# length(integral2$breaks)
# 
# library(emdist)
# 
# 
# P <- t(as.matrix(hist(vals , breaks = seq(from = 0,to = max(vals,vals2)+0.05, by = 0.05), plot = F)$counts))
# Q <- t(as.matrix(hist(vals2 , breaks = seq(from = 0,to = max(vals,vals2)+0.05, by = 0.05), plot = F)$counts))
# emd2d(Q,P) 
# 
# 
# A <- matrix(1:6 / sum(1:6), 1)
# B <- matrix(c(0, 0, 0, 0, 0, 1), 1)
# emd2d(A, B)
# # if we bring the rows closer, the distance will be reduced
# # since mass from the first row has to move down
# # emd2d(A, B,, 0.1)
# 
# 
# 
# vals = calcDistSampledSimple(distributions1,distributions2)
# vals2 = calcDistSampledSimple(distributions1,distributions1)
# 
# hist(vals, xlim = c(0,max(vals2,vals)), col = "red", breaks = 100)
# hist(vals2, add = TRUE, col = "blue", breaks = 100)


