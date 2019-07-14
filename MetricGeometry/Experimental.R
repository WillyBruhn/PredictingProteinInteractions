s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling.R"
source(s1)

approximateNumberOfDistributions <- function(points, n, eps = 1, convergenceParam = 100, maxIt = 1000, withKdTree = TRUE, plotOnTheFly = FALSE){
  # points ... a content of pts-file
  # n ... number of points to sample
  # convergenceParam ... number of times to generate new distributions when it has allready
  #                         been seen
  #
  # Randomly generate a distribution
  # store that distribution in a list
  # calculate new distributions untill
  # a new distribution is allready in the list
  #-------------------------------------
  xy_coords = c()
  
  distributions = list()
  
  # all pairwise emds of all obtained distributions
  d = matrix(0,nrow = 0, ncol = 0)
  
  
  # number of times consecutively the distribution was allready seen
  currentStreak = 0
  
  locatingPoints = NULL
  
  # minimal emd to allready seen distances
  minimalDistance = 100000000
  for(i in 1:maxIt){
    F_ = samplePointsAndCalculateCDFofEc(points,n)
    found = FALSE
    print(paste(i, currentStreak, minimalDistance))
    if(length(distributions) > 0){
      minimalDistance = 100000000
      
      if(length(distributions) < 5 || withKdTree == FALSE){
        for(j in 1:length(distributions)){
          emd = DifferenceOfIntegral(F_,distributions[[j]])
          
          if(minimalDistance > emd) minimalDistance = emd
          
          if(emd < eps) {
            # allready been seen
            found = TRUE
            break;
          }
        }
      } else if(withKdTree == TRUE){
        # use 3 points to locate
        
        whileCounter = 0
        p_predicted = NULL
        while(is.null(p_predicted)){
          
          whileCounter = whileCounter + 1
          print(whileCounter)
          if(is.null(locatingPoints) || whileCounter > 1){
            locatingPoints = sample(c(1:length(distributions)), size = 3,replace = FALSE)
          }
          
          d_1 = DifferenceOfIntegral(F_,distributions[[locatingPoints[1]]])
          d_2 = DifferenceOfIntegral(F_,distributions[[locatingPoints[2]]])
          d_3 = DifferenceOfIntegral(F_,distributions[[locatingPoints[3]]])
          
          
          p_predicted = locatePoint(d_1 = d_1, d_2 = d_2, d_3 = d_3, xy_coords[locatingPoints[1],], xy_coords[locatingPoints[2],], xy_coords[locatingPoints[3],])
        }
        
        # p_predicted = c(0,0)
        nearest <- nn2(data = xy_coords, query = matrix(p_predicted, nrow = 1),k = 1, radius = eps)
        emd = nearest$nn.dists
        
        if(plotOnTheFly){
          plot(xy_coords)
          points(xy_coords[locatingPoints,], col = "red")
          points(x = p_predicted[1], y = p_predicted[2] , col = "green")
          points(xy_coords[nearest$nn.idx,], col = "blue")
          
          print(paste("nn",nearest$nn.idx, length(distributions)+1))
        }
        
        # ?nn2
        
        minimalDistance = emd
        if(emd < eps) {
          # allready been seen
          found = TRUE
        }
      }
    }
    
    # dist has not been seen yet, insert into list
    if(found == FALSE){
      distributions[[length(distributions)+1]] = F_
      
      # recalc distance matrix 
      if(length(distributions) > 1){
        d = matrix(0, nrow = length(distributions), ncol = length(distributions))
        for(i in 1:length(distributions)){
          for(j in 1:length(distributions)){
            d[i,j] = DifferenceOfIntegral(distributions[[i]],distributions[[j]])
          }
        }
      }
      
      # print(d)
      if(nrow(d) >= 3){
        xy_coords = cmdscale(as.matrix(d))
        
        if(plotOnTheFly){
          # plot(xy_coords)
        }
        
      }
      
      
      currentStreak = 0
    } else {
      currentStreak = currentStreak + 1
    }
  }
  
  return(list("distributions" = distributions, "distances"  = d))
}


times = 1
# all_protein_models_with_distances2 = vectorized_get_allModels(n = 100, m= 10, times = times, MC = TRUE)

all_protein_models_with_distances = rep(0,m*times)
all_protein_models_with_distances = list()
for(i in 1:times){
  all_protein_models_with_distances= c(all_protein_models_with_distances, all_protein_models_with_distances2[[i]])
}

all_protein_models_with_distances = not_vectorized_get_allModels(n = 100,m = 3, times = times)

DE_big = calculateAllDistancesExact(all_protein_models_with_distances,m=3)

(DE_parallel - DE_big)[1:6,1:6]
cor(c(DE_parallel),c(DE_big))

is.numeric(DE_big)
is.numeric(q)
q


names = unique(colnames(DE_big))
emds = matrix(0, nrow = length(names), ncol = length(names))
rownames(emds) = names
colnames(emds) = names
emds

for(i in 1:length(names)){
  for(j in 1:length(names)){
    print(paste(i,j))
    emds[i,j] = histDist2(DE_big,names[i], names[j])
  }
}

felix = readDistanceMatrix2()

cor(c(emds[-106,-106]),c(felix))


# devtools::install_github("eddelbuettel/rbenchmark")
library(rbenchmark)

benchmark("mclapply" = vectorized_get_allModels(n = 10, m= 10, 5, MC = TRUE),
          "sapply" = vectorized_get_allModels(n = 10, m= 10, 5, MC = FALSE),replications = 1,
          columns = c("test", "replications", "elapsed",
                      "relative", "user.self", "sys.self"))

# test replications elapsed relative user.self sys.self
# 1 mclapply            1  17.301    1.000     0.060    0.121
# 2   sapply            1  28.303    1.636    28.187    0.228


rotateByAngleAtPoint <- function(points, center, angle){
  p2 = cbind(points[,1] - center[1], points[,2] - center[2])
  p2 = rotateAroundZ(cbind(p2,0), angle)[,1:2]
  
  p2 = cbind(p2[,1] + center[1], p2[,2] + center[2])
  return(p2)
}

mergeTheseIndices <- function(mergeInd, otherIndices, all_protein_models_with_distances, m, xli = NULL, yli = NULL, plotPairs = FALSE){
  ps2 = list()
  for(i in 1:length(otherIndices)){
    ps = project_geometry(all_protein_models = all_protein_models_with_distances,ind1 = mergeInd,ind2 = otherIndices[i])
    ps2[[i]] = ps
    if(plotPairs) plotPrettyProjection(ps,getFunctionalProteins(), onlyGeomCenters = FALSE, xli = xli, yli = yli)
  }
  
  p_mergeds = mergeProjections(ps2, mergeInd = mergeInd, number = m)
  plotPrettyProjection(p_mergeds,getFunctionalProteins(), onlyGeomCenters = FALSE, xli = xli, yli = yli)
}

mergeTheseIndices(3,c(1,3,4,8,20,30),all_protein_models_with_distances, m)
mergeTheseIndices(1,c(1,3,4,8,20,30),all_protein_models_with_distances, m)
mergeTheseIndices(20,c(1,3,4,8,20,30),all_protein_models_with_distances, m)


ps = project_geometry(all_protein_models = all_protein_models_with_distances,ind1 = 1,ind2 = 40)
ps2 = project_geometry(all_protein_models = all_protein_models_with_distances,ind1 = 1,ind2 = 50)
ps3 = project_geometry(all_protein_models = all_protein_models_with_distances,ind1 = 1,ind2 = 97)
p_mergeds = mergeProjections(list(ps,ps2,ps3), mergeInd = 1, number = m)
plotPrettyProjection(ps,getFunctionalProteins(), onlyGeomCenters = FALSE)
plotPrettyProjection(ps2,getFunctionalProteins(), onlyGeomCenters = FALSE)

plotPrettyProjection(p_mergeds,getFunctionalProteins(), onlyGeomCenters = FALSE)

all_protein_models_with_distances[[1]]

mergeInd = 7
projections_list = list()
for(i in 1:length(all_protein_models_with_distances)){
  # for(j in 1:length(all_protein_models_with_distances[[i]])){
  print(i)
  # if(i != mergeInd){
  projections_list[[i]] = project_geometry(all_protein_models = all_protein_models_with_distances, ind1 = mergeInd, ind2 = i)[(m+1):(2*m),]
  # }
}
# }

p_merged = mergeProjections(projections_list, mergeInd = mergeInd, number = m)

geoms = plotPrettyProjection(p_merged,getFunctionalProteins(), onlyGeomCenters = FALSE)

geoDists = as.matrix(dist(geoms))
neighbors = 20
closest = names(geoDists[which.minn(geoDists[which(rownames(geoDists) == "000_Trx"),],n = neighbors),1])
length(which(closest %in% getFunctionalProteins() == TRUE))


# wally


all_protein_models_with_distances[[mergeInd]]$model$distributions

calcErrorsToInd <- function(all_protein_models_with_distances, mergeInd){
  meanError = rep(0,length(all_protein_models_with_distances))
  meanVarError = rep(0,length(all_protein_models_with_distances))
  
  meanDistances = rep(0,length(all_protein_models_with_distances))
  meanVar = rep(0,length(all_protein_models_with_distances))
  
  for(i in 1:length(all_protein_models_with_distances)){
    print(i)
    DE_exact = calcAllDistributionPairs(all_protein_models_with_distances[[i]]$model$distributions,all_protein_models_with_distances[[mergeInd]]$model$distributions)
    
    dists = as.matrix(dist(projections_list[[i]]))
    # DE_exact
    
    Error = DE_exact - dists
    meanError[i] = mean(Error)
    meanVarError[i] = var(as.vector(DE_exact - dists))
    
    meanDistances[i] = mean(DE_exact)
    meanVar[i] = var(DE_exact)
  }
  
  mE = mean(meanError)
  mVE = mean(meanVarError)
  
  mD = mean(meanDistances)
  mV = mean(meanVar)
  
  hist(meanError)
  hist(meanDistances)
  
  return(list("mean(meanError)" = mE, "mean(meanVarError)" =  mVE, "mean(meanDistances)" = mD, "mean(meanVar)" = mV))
}


calcErrorsToInd(all_protein_models_with_distances, 7)

deletePointsWithHighError <- function(all_protein_models_with_distances, th = 10){
  m2 = length(all_protein_models_with_distances[[2]])*2
  print(m2)
  
  Errors = matrix(0,ncol = m2, nrow = length(all_protein_models_with_distances)*m2)
  for(i in 1:length(all_protein_models_with_distances)){
    print(i)
    DE_exact = calcAllDistributionPairs(all_protein_models_with_distances[[i]]$model$distributions,all_protein_models_with_distances[[mergeInd]]$model$distributions)
    
    dists = as.matrix(dist(projections_list[[i]]))
    # DE_exact
    
    print(dim(dists))
    print(dim(DE_exact))
    
    Error = DE_exact - dists
    
    start = ((i-1)*m2+1)
    end = start+m2-1
    
    print(paste(start,end))
    Errors[start:end, ] = Error
  }
  return(Errors)
}

Errors = deletePointsWithHighError(all_protein_models_with_distances)

nrow(Errors)

th = 1
mean(Errors[which(rowMeans(abs(Errors)) <= th), ])
length(which(rowMeans(abs(Errors)) <= th)) 
p_merged_cleaned = p_merged[-which(rowMeans(abs(Errors)) > th), ]

p_merged_cleaned

calcErrorsToInd(all_protein_models_with_distances, 66)


dist_projected = as.matrix(dist(p_merged))
dist_projected

unique(colnames(dist_projected))






Error_big = DE_Big_exact- as.matrix(dist(p_merged))

mean(Error_big)
var(as.vector(Error_big))

var(as.vector(DE_Big_exact))


# dist(p_merged[1:10,])




geoms = plotPrettyProjection(p_merged_cleaned,getFunctionalProteins(), onlyGeomCenters = FALSE)

geoDists = as.matrix(dist(geoms))
neighbors = 15
closest = names(geoDists[which.minn(geoDists[which(rownames(geoDists) == "000_Trx"),],n = neighbors),1])
length(which(closest %in% getFunctionalProteins() == TRUE))


alldist_vs_all = as.matrix(dist(p_merged))

names = unique(colnames(alldist_vs_all))
emds = matrix(0, nrow = length(names), ncol = length(names))
rownames(emds) = names
colnames(emds) = names
emds

for(i in 1:length(names)){
  for(j in 1:length(names)){
    print(paste(i,j))
    emds[i,j] = histDist2(alldist_vs_all,names[i], names[j])
  }
}


length(which(colnames(emds)[which.minn(emds[which(rownames(emds) == "000_Trx"),], n = 20)] %in% c(getFunctionalProteins(), "000_Trx")) == TRUE)


# emds_proj = cmdscale(as.matrix(cailliez(as.dist(emds))))
# plot(emds_proj)
# inds = which(rownames(emds_proj) %in% c(getFunctionalProteins(), "000_Trx") == TRUE)
# points(x = emds_proj[inds,1], y = emds_proj[inds,2], col = "green")


length(intersect(rownames(demd2), rownames(emds)))
demd2 = readDistanceMatrix2()

demd2 = demd2[order(rownames(demd2)),order(rownames(demd2))]
min(demd2)

emds2 = emds[order(rownames(emds)),order(rownames(emds))]
emds2 = emds2[,-106]
emds2 = emds2[-106,]
abs(emds2 -demd2)

heatmap(cor(emds2,demd2))
cor(c(emds2), c(demd2))


# d1 = readDistanceMatrix3("/home/willy/Schreibtisch/106Test/RepSubOutput/EMD_10_5000_1.0000_0.0000_0.0000_0.1000_id_opt_NNact_0.csv")
# d2 = readDistanceMatrix3("/home/willy/Schreibtisch/106Test/RepSubOutput/EMD_50_5000_1.0000_0.0000_1.0000_1.0000_id_opt_NNact_0.csv")
# 
# cor(c(d1[-106,-106]), c(demd2))


library(emdist)
histDist2 <- function(alldist_vs_all, name1,name2){
  ind1 = which(rownames(alldist_vs_all) == name1)
  vals1 = hist(alldist_vs_all[ind1,ind1], breaks = 100, plot = F)$counts
  
  ind2 = which(colnames(alldist_vs_all) == name2)
  # vals2 = hist(alldist_vs_all[ind2,ind2], breaks = 100)$counts
  
  vals3 = hist(alldist_vs_all[ind1,ind2], breaks = 100, plot = F)$counts
  
  P <- t(as.matrix(hist(vals1 , breaks = seq(from = 0,to = max(vals1,vals3)+0.05, by = 0.05), plot = F)$counts))
  Q <- t(as.matrix(hist(vals3 , breaks = seq(from = 0,to = max(vals1,vals3)+0.05, by = 0.05), plot = F)$counts))
  ed1 = emd2d(Q,P) 
  
  # P <- t(as.matrix(hist(vals2 , breaks = seq(from = 0,to = max(vals2,vals3)+0.05, by = 0.05), plot = F)$counts))
  # Q <- t(as.matrix(hist(vals3 , breaks = seq(from = 0,to = max(vals2,vals3)+0.05, by = 0.05), plot = F)$counts))
  # ed2 = emd2d(Q,P) 
  
  return(ed1)
}





plotProjections <- function(projections, const = 6,xli = NULL, yli = NULL){
  
  x_ = matrix(0,nrow = length(projections), ncol = nrow(projections[[1]]))
  y_ = matrix(0,nrow = length(projections), ncol = nrow(projections[[1]]))
  for(i in 1:length(projections)){
    x_[i,] =  projections[[i]][,1]
    y_[i,] =  projections[[i]][,2]
  }
  
  xli = c(min(x_)-const,max(x_)+const)
  yli = c(min(y_)-const,max(y_)+const)
  
  plot(projections[[1]], col = "red", xlim = xli, ylim = yli)
  if(length(projections) > 1){
    for(i in 2:length(projections)){
      points(projections[[i]], col = "blue", xlim = xli, ylim = yli)
    }
  }
  
  
}







all_protein_models = getAllDistributions(OutputPath = OutputPath,n = 100,m_distributions = 10)
ind1 = 1
ind2 = 7
d_1 = calcAllDistributionPairs(all_protein_models[[ind1]]$distributions,all_protein_models[[ind1]]$distributions)
projection1 = myProjection(d_1)
rownames(projection1) = rep(all_protein_models[[ind1]]$name,nrow(projection1))


d_2 = calcAllDistributionPairs(all_protein_models[[ind2]]$distributions,all_protein_models[[ind2]]$distributions)
projection2 = myProjection(d_2)
rownames(projection2) = rep(all_protein_models[[ind2]]$name,nrow(projection2))

p_1 = projection1[1,]
p_2 = projection1[2,]
p_3 = projection1[3,]

target_origin_index = 1

d_1_1 = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[1]], G_ = all_protein_models[[ind2]]$distributions[[target_origin_index]])
d_1_2 = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[1]], G_ = all_protein_models[[ind2]]$distributions[[2]])
d_1_3 = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[1]], G_ = all_protein_models[[ind2]]$distributions[[3]])

plotProjections(list(projection1,projection2))





pythagoras <-function(distance, target_row, target_column, help_column){
  # distance between target_row,target_column
  # we know the distance between (target_row,help_column)  and (target_column,help_column)
  
  dist = sqrt(distance[target_column,help_column]^2 + distance[target_row,help_column]^2)
  
  return(dist)
}


#-----------------------------------------------------------------------------------

getSharedDistances <- function(d_1_1, d_1_2, d_1_3,projection1,projection2, distances1 = NULL, distances2 = NULL){
  # d_1_1 ... distance of 1 to 1
  #---------------------------------
  n1 = nrow(projection1)
  n2 = nrow(projection2)
  
  dim = n1 + n2
  distance_shared = matrix(0, nrow = dim, ncol = dim)
  
  if(is.null(distances1)){
    distances1 = as.matrix(dist(projection1, upper = TRUE, diag = TRUE))
  }
  
  if(is.null(distances2)){
    distances2 = as.matrix(dist(projection2, upper = TRUE, diag = TRUE))
  }
  
  print(checkIfTriangleHolds(distances1))
  print(checkIfTriangleHolds(distances2))
  print(isSymmetric(distances1))
  print(isSymmetric(distances2))
  
  distance_shared[1:n1,1:n1] = distances1
  distance_shared[(n1+1):dim,(n1+1):dim] = distances2
  
  
  distance_shared[1,1+n1] = d_1_1
  distance_shared[1+n1,1] = d_1_1
  
  distance_shared[1,2+n1] = d_1_2
  distance_shared[2+n1,1] = d_1_2
  
  distance_shared[1,3+n1] = d_1_3
  distance_shared[3+n1,1] = d_1_3
  
  for(j in 4:(n1)){
    start = projection1[1,]
    
    p1 = projection2[1,]
    p2 = projection2[2,]
    p3 = projection2[3,]
    
    d1 = distance_shared[1+n1,j+n1]
    d2 = distance_shared[2+n1,j+n1]
    d3 = distance_shared[3+n1,j+n1]
    
    loc = locatePoint(d_1 = d1, d_2 = d2, d_3 = d3, point1 = p1, point2 = p2, point3 = p3)
    
    distance_shared[1,j] = euclideanDistance(loc,start)
  }
  
  print(distance_shared)
  # for(j in 2:n2){
  #   for(i in 1:n1){
  #     # pythagoras
  #     dist = pythagoras(distance = distance_shared, target_row = i, target_column = j+n1, help_column = j+n1-1)
  #     distance_shared[i,j+n1] = dist
  #     distance_shared[j+n1,i] = dist
  #   }
  # }
  
  colnames(distance_shared) = c(rownames(projection1),rownames(projection2))
  rownames(distance_shared) = c(rownames(projection1),rownames(projection2))
  
  # print(distance_shared)
  
  # return(distance_shared)
  
  print(isSymmetric(distance_shared))
  print(checkIfTriangleHolds(distance_shared))
  
  projection_shared = myProjection(distance_shared)
  rownames(projection_shared) = c(rownames(projection1),rownames(projection2))
  
  projection1_s = projection_shared[1:n1,]
  projection2_s = projection_shared[(n1+1):dim,]
  
  return(list("distance_shared" = distance_shared, "projection1_s" = projection1_s, "projection2_s" = projection2_s))
}

shared = getSharedDistances(d_1_1, d_1_2, d_1_3, projection1 = projection1, projection2)



DifferenceOfIntegral(all_protein_models[[1]]$distributions[[2]], all_protein_models[[7]]$distributions[[1]])

shared[1,5]+shared[5,6]

sqrt((shared[1,5])^2 + (shared[5,6])^2)

shared

distances = shared

for(k in 1:nrow(distances)){
  for(i in 1:nrow(distances)){
    for(j in 1:ncol(distances)){
      
      d1 = distances[i,k]
      d2 = distances[i,j] 
      d3 = distances[j,k]
      
      vec= c(d1,d2,d3)
      
      v = which.maxn(vec, n = 1)
      
      
      if(!(vec[v] == sqrt((vec[-v][1])^2 - (vec[-v][1])^2))){
        print("no")
        print(paste(i,j,k))
        
        print(distances[i,k])
        print(sqrt((distances[i,j])^2 + (distances[j,k])^2))
        print(paste("d[",i,",",k,"] != sqrt(d[",i,",",j,"]^2 + d[",j,",",k,"])^2)",sep = ""  ))
        print(paste(d[i,k], d[i,j], d[j,k]))
      }
    }
  }
}

distances


# 0.519230772459262 +7.6700944341509
#8.38045636368756


shared$distance_shared

plotProjections(list(shared$projection1_s,shared$projection2_s))

plotProjections(list(projection1,projection2))


# check if distances correct
dim = length(all_protein_models[[ind1]]$distributions) + length(all_protein_models[[ind2]]$distributions)
d_exact = matrix(0,nrow = dim,ncol = dim)

rownames(d_exact) = c(rownames(projection1), rownames(projection2))
colnames(d_exact) = c(rownames(projection1), rownames(projection2))

n1 = length(all_protein_models[[ind1]]$distributions)
n2 = length(all_protein_models[[ind2]]$distributions)
dim = length(all_protein_models[[ind1]]$distributions) + length(all_protein_models[[ind2]]$distributions)
for(i in 1:dim){
  for(j in 1:dim){
    if(i <= length(all_protein_models[[ind1]]$distributions) && j <= length(all_protein_models[[ind1]]$distributions)){
      d_exact[i,j] = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[i]], G_ = all_protein_models[[ind1]]$distributions[[j]])
    }  else if(i > n1 && j > n1){
      print(paste(i-n1,j-n1))
      d_exact[i,j] = DifferenceOfIntegral(F_ = all_protein_models[[ind2]]$distributions[[i-n1]], G_ = all_protein_models[[ind2]]$distributions[[j-n1]])
    } else if(i > n1 && j <= n1){
      # print(paste(i-n1,j-n1))
      d_exact[i,j] = DifferenceOfIntegral(F_ = all_protein_models[[ind2]]$distributions[[i-n1]], G_ = all_protein_models[[ind1]]$distributions[[j]])
    } else if(i <= n1 && j > n1){
      # print(paste(i-n1,j-n1))
      d_exact[i,j] = DifferenceOfIntegral(F_ = all_protein_models[[ind1]]$distributions[[i]], G_ = all_protein_models[[ind2]]$distributions[[j-n1]])
    }
    
    # print(paste(all_protein_models[[ind1]]$name, all_protein_models[[ind2]]$name))
    # print(colnames(shared$distance_shared[,c(i,length(all_protein_models[[ind1]]$distributions)+j)]))
    # 
    # print(d_1 - shared$distance_shared[i,length(all_protein_models[[ind1]]$distributions)+j])
  }
}

d_exact

projection = myProjection(d_exact)

myProjection <- function(distances, plot = FALSE){
  
  points = matrix(0,nrow = nrow(distances), ncol = 2)
  points[1,] = c(0,0)
  # second_ind = which.maxn(distances[1,],n = 1)
  # second_ind = 2
  
  points[2,] = c(distances[1,2],0)
  
  p = cirlceIntersect(A.x = points[1,1], A.y = points[1,2],A.r = distances[1,3],
                      B.x = points[2,1], B.y = points[2,2],B.r = distances[2,3])
  
  points[3,] = p$P1
  
  if(plot){
    plot(points[c(1,2,3),], xlim = c(-distances[1,2], distances[1,2]), ylim = c(-distances[1,2], distances[1,2]))
  }
  
  for(i in 4:nrow(points)){
    # print(i)
    print(paste(distances[1,i], distances[2,i], points[1,1], points[1,2], points[2,1], points[2,2], points[3,1], points[3,2]))
    q = locatePoint(distances[i,1], distances[i,2], distances[i,3], points[1,], points[2,], points[3,])
    print(q)
    points[i,] = q
  }
  
  if(plot){
    points(x = points[,1], y = points[,2])
  }
  
  rownames(points) = rownames(distances)
  
  return(points)
}

d_1 = calcAllDistributionPairs(all_protein_models[[1]]$distributions,all_protein_models[[1]]$distributions)
projection1 = myProjection(d_1, plot = TRUE)


# n = 400
# ps = matrix(rnorm(n*2), ncol = 2)
# d_example = as.matrix(dist(ps,upper = TRUE,diag = TRUE))
# # checkIfTriangleHolds(d_example)
# 
# myProjection(d_example)


plotPrettyProjection <- function(projection, functionals = NULL, onlyGeomCenters = FALSE, constPlot = 1, fontSize = 1, cex = 1, xli = NULL, yli = NULL){
  
  if(is.null(xli)) xli = c(min(projection[,1])-constPlot, max(projection[,1])+constPlot)
  if(is.null(yli)) yli = c(min(projection[,2])-constPlot, max(projection[,2])+constPlot)
  
  if(onlyGeomCenters == FALSE){
    plot(projection, col = "green", xlim = xli, ylim = yli, pch = 19, cex = cex, ylab = "DE", xlab = "DE")
  } 
  else {
    plot(NULL, col = "green", xlim = xli, ylim = yli, pch = 19, cex = cex, ylab = "DE", xlab = "DE")
  }
  
  
  names = unique(rownames(projection))
  geoms = matrix(0,nrow = length(names), ncol = 2)
  for(i in 1:length(names)){
    inds = which(rownames(projection) == names[i])
    
    geom = c(sum(projection[inds,1]), sum(projection[inds,2]))/length(inds)
    
    functional = FALSE
    if(names[i] %in% functionals){
      functional = TRUE
    }
    if(functional) {
      text(x = geom[1], y = geom[2], labels = c(names[i]), col = "blue", cex = fontSize)
      if(onlyGeomCenters == FALSE) points(projection[inds,], add = TRUE, col = "blue")
    }
    else{
      text(x = geom[1], y = geom[2], labels = c(names[i]), col = "red", cex = fontSize)
    } 
    
    geoms[i,] = geom
  }
  
  rownames(geoms) = names
  
  return(geoms)
}


plotPrettyProjection(projection = projection, constPlot = 3, fontSize = 1, cex = 0.5)


histDist <- function(model, model2){
  print(paste(model$name,model2$name))
  
  vals_AA = calcDistSampled(model$distributions,model$distributions, 50)
  vals_AB = calcDistSampled(model$distributions,model2$distributions, 50)
  # vals_BB = calcDistSampled(all_protein_models[[j]]$distributions,all_protein_models[[j]]$distributions, 500)
  
  P <- t(as.matrix(hist(vals_AA , breaks = seq(from = 0,to = max(vals_AA,vals_AB)+0.05, by = 0.05), plot = F)$counts))
  Q <- t(as.matrix(hist(vals_AB , breaks = seq(from = 0,to = max(vals_AA,vals_AB)+0.05, by = 0.05), plot = F)$counts))
  return(emd2d(Q,P))
}

calculateAllPairWiseEmds <- function(all_protein_models,m = 500){
  
  d = matrix(0,nrow = length(all_protein_models), ncol= length(all_protein_models))
  
  for(i in 1:length(all_protein_models)){
    vals_AA = calcDistSampled(all_protein_models[[i]]$distributions,all_protein_models[[i]]$distributions, m)
    for(j in 1:length(all_protein_models)){
      print(paste(i,j))
      vals_AB = calcDistSampled(all_protein_models[[i]]$distributions,all_protein_models[[j]]$distributions, m)
      # vals_BB = calcDistSampled(all_protein_models[[j]]$distributions,all_protein_models[[j]]$distributions, 500)
      
      P <- t(as.matrix(hist(vals_AA , breaks = seq(from = 0,to = max(vals_AA,vals_AB)+0.05, by = 0.05), plot = F)$counts))
      Q <- t(as.matrix(hist(vals_AB , breaks = seq(from = 0,to = max(vals_AA,vals_AB)+0.05, by = 0.05), plot = F)$counts))
      d[i,j] = emd2d(Q,P) 
    }
  }
  
  df = as.data.frame(d)
  return(df)
}

df = calculateAllPairWiseEmds(all_protein_models,m = 50)

getFunctionalProteins()

hist(as.numeric(df[7,]),breaks = seq(from = 0, to = 160, by =0.05))
df[7,20]
df[7,23]
df[7,94]

all_protein_models[[99]]$name
which.minn(as.numeric(df[7,]), n = 20)



histDist(all_protein_models[[20]], all_protein_models[[7]])




mylist <- list(a=1,b=2,c=3)
myfxn <- function(var1,var2){
  var1*var2
}
var2 <- 2

sapply(mylist,myfxn,var2=var2)

mapply(rep, times = 1:4, MoreArgs = list(x = 42))

mapply(histDist,all_protein_models,all_protein_models)

all_protein_models[1:2][1]


distributions1 = precalCulateDistributionsModel(pts_pos, 100, 500)
distributions2 = precalCulateDistributionsModel(pts_pos, 100, 50)

vals2 = calcAllDistributionPairs(all_protein_models[[7]]$distributions,all_protein_models[[1]]$distributions)
vals = calcAllDistributionPairs(all_protein_models[[7]]$distributions,all_protein_models[[7]]$distributions)
vals3 = calcAllDistributionPairs(all_protein_models[[20]]$distributions,all_protein_models[[7]]$distributions)

vals2 = calcDistSampled(all_protein_models[[7]]$distributions,all_protein_models[[10]]$distributions, 2500)


hist(vals, xlim = c(0,max(vals2,vals,vals3)), col = "red", breaks = 100)
integral2 = hist(vals2, add = TRUE, col = "blue", breaks = 100)
hist(vals3, add = TRUE, col = "yellow", breaks = 100)

integral = hist(vals, breaks = 100)
integral$breaks
length(integral2$breaks)

library(emdist)


P <- t(as.matrix(hist(vals , breaks = seq(from = 0,to = max(vals,vals2)+0.05, by = 0.05), plot = F)$counts))
Q <- t(as.matrix(hist(vals2 , breaks = seq(from = 0,to = max(vals,vals2)+0.05, by = 0.05), plot = F)$counts))
emd2d(Q,P) 


A <- matrix(1:6 / sum(1:6), 1)
B <- matrix(c(0, 0, 0, 0, 0, 1), 1)
emd2d(A, B)
# if we bring the rows closer, the distance will be reduced
# since mass from the first row has to move down
# emd2d(A, B,, 0.1)



vals = calcDistSampledSimple(distributions1,distributions2)
vals2 = calcDistSampledSimple(distributions1,distributions1)

hist(vals, xlim = c(0,max(vals2,vals)), col = "red", breaks = 100)
hist(vals2, add = TRUE, col = "blue", breaks = 100)


#----------------------------------------------------------------------------------------------------------------
n = 100
m = 10
times = 10
OutputPath = "/home/willy/Schreibtisch/106Test/Output/"
all_protein_models_with_distances = not_vectorized_get_allModels(OutputPath = OutputPath, n = n,m = m, times = times)
DE_parallel = computeAllDistancesParallel(all_protein_models_with_distances)
Emd_distance_matrix = emd_parallel(DE_parallel = DE_parallel,m = m)


#--------------------------------------------------------------------------------
# projectionMethod
#--------------------------------------------------------------------------------


mergeInd = 1
projections_list = list()
for(i in 1:length(all_protein_models_with_distances)){
  projections_list[[i]] = project_geometry(all_protein_models = all_protein_models_with_distances, ind1 = mergeInd, ind2 = i)[(m+1):(2*m),]
}

p_merged = mergeProjections(projections_list, mergeInd = mergeInd, number = m)

geoms = plotPrettyProjection(p_merged,getFunctionalProteins(), onlyGeomCenters = FALSE)

geoDists = as.matrix(dist(geoms))
neighbors = 20
closest = names(geoDists[which.minn(geoDists[which(rownames(geoDists) == "000_Trx"),],n = neighbors),1])
length(which(closest %in% getFunctionalProteins() == TRUE))

DE_projected = as.matrix(dist(p_merged))
Emd_distance_matrix_approx2 = emd_parallel(DE_parallel = DE_projected,m = m)

cor(c(Emd_distance_matrix_approx2),c(Emd_distance_matrix_approx))

Err = deletePointsWithHighError(all_protein_models_with_distances = all_protein_models_with_distances)
Errors = Err
nrow(Errors)

th = 1
mean(Errors[which(rowMeans(abs(Errors)) <= th), ])
length(which(rowMeans(abs(Errors)) <= th)) 
p_merged_cleaned = p_merged[-which(rowMeans(abs(Errors)) > th), ]


geoms = plotPrettyProjection(p_merged_cleaned,getFunctionalProteins(), onlyGeomCenters = FALSE)
geoDists = as.matrix(dist(geoms))
neighbors = 20
closest = names(geoDists[which.minn(geoDists[which(rownames(geoDists) == "000_Trx"),],n = neighbors),1])
length(which(closest %in% getFunctionalProteins() == TRUE))




Emd_distance_matrix2 = readDistanceMatrix1("/home/willy/Schreibtisch/106Test/Output/quickEmd_n_100_m_3.csv")


#--------------------------------------------------------------------------------
# projectionMethod
#--------------------------------------------------------------------------------

# pts = read_pts_file(OutputPath = OutputPath,protName = "000_Trx")
# CDF = samplePointsAndCalculateCDFofEc(all_pts = pts,n = 100,plot = FALSE)
# 
# integralOfDistribution(CDF)
# x = integrate(CDF, lower = 0, upper = 100)
# x$value

integrateStepFun <- function(F_, a){
  
  sup = 0
  ind = 0
  for(i in 1:length(knots(F_))){
    if(knots(F_)[i] <= a){
      sup = knots(F_)[i]
      ind = i
    } else{
      break
    }
  }
  
  su = 0
  
  if(ind > 0){
    for(i in 1:(ind-1)){
      su = su + (knots(F_)[i+1]-knots(F_)[i])*F_(knots(F_)[i])
    }
    
    su = su + (a-knots(F_)[ind])*F_(knots(F_)[ind])
  }

  return(su)
}

integrateStepFunInverse <- function(F_,b){
  val_old = 0
  for(i in seq(knots(F_)[1],knots(F_)[length(knots(F_))],0.1)){
    val = integrateStepFun(F_,i)
    if(val > b){
      return(val_old)
    }
    val_old = val
  }
  return(NULL)
}

approximateCDF <- function(F_, q){
  # F_ ... cumulitive distribution function
  # q ... number of equal heights to approximate F_
  #-------------------------------------------------
  knots = knots(F_)
  maxF = max(knots)
  n_F = length(knots)

  # print(knots)
    
  F_app = rep(0,q)
  for(i in 1:q){
    
    # ind = n_F*i/(q+1)
    # 
    # if(ind - floor(ind) > 0){
    #   F_app[i] = (knots[n_F*i/(q+1)] + knots[(n_F*i+1)/(q+1)])/2
    # } else{
    #   F_app[i] = knots[n_F*i/(q+1)]
    # }
    
    print((i/q)*(maxF))
    F_app[i]  = integrateStepFun(F_,(i/q)*(maxF))
    
    # F_app[i] = integrate(CDF, lower = 0, upper = (i/q)*(maxF-0.1))$value
    
  }

  return(F_app)  
}

pts = read_pts_file(OutputPath = OutputPath,protName = "000_Trx")
CDF = samplePointsAndCalculateCDFofEc(all_pts = pts,n = 100,plot = FALSE)

CDF_approx = approximateCDF(CDF,q = 2)

integrateStepFunInverse(CDF,3.5)
integrateStepFun(CDF,25)


plotDistWithApproximation(CDF,CDF_approx)

plotDistWithApproximation <- function(F_,F_app,plotToFile = NULL){

  if(!is.null(plotToFile)){
    pdf(file = plotToFile)
  }
    
    n = length(knots(F_))
    q = length(F_app)
    sfun = stepfun(F_app , y = c(0,seq(1/q,1,1/q)))
    
    print(knots(F_)[n/2])
    
    
    plot(F_, main = paste("n = ",n , ", q= ", q, sep = ""), xlim = c(0,knots(F_)[n/2]*2))
    plot(sfun, add = TRUE)
  if(!is.null(plotToFile)){
    dev.off()
  }  
}

plotDistWithApproximationDifferentQ <- function(OutputPath, protName, n = 100, qlist = c(1:9), pos =TRUE, plotToFile = NULL){
  
  pts = read_pts_file(OutputPath = OutputPath,protName = protName,pos = pos)
  pos13_F1 = samplePointsAndCalculateCDFofEc(all_pts = pts,n = n,plot = FALSE)
  
  if(!is.null(plotToFile)){
    pdf(file = plotToFile)
  }
  par(mfrow=c(sqrt(length(qlist)),sqrt(length(qlist))))
  
  for(i in 1:length(qlist)){
    F_app1 = approximateCDF(pos13_F1,qlist[i])
    plotDistWithApproximation(pos13_F1,F_app1)
  }

  if(!is.null(plotToFile)){
    dev.off()
  }  
  
}


DifferenceOfIntegral_F_and_Approx <- function(F_, F_app){
  q = length(F_app)
  sfun = stepfun(F_app , y = c(0,seq(1/q,1,1/q)))
  
  DifferenceOfIntegral(F_,sfun)
}

manhattan_dist <- function(rating1, rating2){
  distance <- abs(rating1-rating2)
  distance <- sum(distance)
  return(distance)
}

manhattan_dist_distribution_approx <- function(rating1, rating2){
  return(manhattan_dist(rating1,rating2)/(length(rating1)))
}

plot2d_dist_approx <-function(pos13_F_approx_list, name=NULL, col = "red", xli = NULL, yli = NULL, add = FALSE, onlyGeo = FALSE){
  
  if(add == FALSE){
    x_val = pos13_F_approx_list[[1]][1]
    y_val = pos13_F_approx_list[[1]][2]
    
    if(!is.null(xli) && !is.null(yli)) {
      plot(x = x_val, y = y_val-x_val, xlim = xli, ylim = yli,  col = col)
    } else {
      plot(x = x_val, y = y_val-x_val, col = col)
    }
  }


  x_vals = rep(0,length(pos13_F_approx_list))
  y_vals = rep(0,length(pos13_F_approx_list))
  
  for(i in 1:length(pos13_F_approx_list)){
    x_val = pos13_F_approx_list[[i]][1]
    y_val = pos13_F_approx_list[[i]][2]
    
    if(onlyGeo == FALSE) points(x = x_val, y = y_val-x_val, col = col)
    
    
    x_vals[i] = x_val
    y_vals[i] = y_val-x_val
  }
  
  # plot the name
  if(!is.null(name)){
    geo = c(sum(x_vals)/length(x_vals), sum(y_vals)/length(y_vals))
    
    print(geo)
    if(onlyGeo == TRUE) {
      text(x = geo[1], y = geo[2], labels = c(name), cex = 2, col = col)
    } else {
      text(x = geo[1], y = geo[2], labels = c(name), cex = 2, col = "black")
    }

  }
}

generateF_approximations <- function(OutputPath, protName, n = 100, m = 10, q = 2, pos =TRUE){
  
  pos13 = read_pts_file(OutputPath = OutputPath,protName = protName, pos = pos)
  
  pos13_F_list = list()
  pos13_F_approx_list = list()
  
  for(i in 1:m){
    pos13_F_list[[i]] = samplePointsAndCalculateCDFofEc(all_pts = pos13,n = n,plot = FALSE)
    pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
  }
  
  return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
}

getDistanceMatrixFromProjection <- function(df){
  
  q= ncol(df)-1
  
  distances = as.matrix(dist(x = df[,2:q], diag = TRUE, upper = TRUE),method ="manhattan")
  colnames(distances) = df[,1]
  rownames(distances) = df[,1]
  
  return(distances)
}

getGeometricCenters <- function(df){

  names = unique(df[,1])
  
  geos = data.frame(matrix(0, nrow = length(names), ncol = 3))
  colnames(geos) = c("name", "x", "y")
  
  for(i in 1:length(names)){
    inds = which(df[,1] == names[i])
    
    x = sum(df[inds,2])/length(df[inds,2])
    y = sum(df[inds,3])/length(df[inds,3])
    geos[i,1] = names[i]
    geos[i,2] = x
    geos[i,3] = y
  }
  
  return(geos)
}
#-------------------------------------------------------------------
## Example

calcError <- function(all_pts, n = 100, q = 2, plot = TRUE, verbose = TRUE, plotToFile = NULL){
  # q = 2
  # n = 100
  pos13_F1 = samplePointsAndCalculateCDFofEc(all_pts = all_pts,n = n,plot = FALSE)
  F_app1 = approximateCDF(pos13_F1,q)
  plotDistWithApproximation(pos13_F1,F_app1,plotToFile)
  DifferenceOfIntegral_F_and_Approx(pos13_F1,F_app1)
  
  pos13_F2 = samplePointsAndCalculateCDFofEc(all_pts = all_pts,n = n,plot = FALSE)
  F_app2 = approximateCDF(pos13_F2,q)
  plotDistWithApproximation(pos13_F2,F_app2,plotToFile)
  DifferenceOfIntegral_F_and_Approx(pos13_F2,F_app2)
  
  DE_exact = DifferenceOfIntegral(pos13_F1,pos13_F2)
  DE_approx = manhattan_dist_distribution_approx(F_app1,F_app2)
  
  if(verbose){
    print(paste("DE_exact = ", DE_exact, ", DE_approx = ", DE_approx, sep =""))
    print(paste("Error =", abs(DE_exact- DE_approx)))
  }

  return(abs(DE_exact- DE_approx))
}


calcErrors <- function(all_pts, times = 100, n = 100, q = 2, plot = FALSE){
  errors = rep(0,times)
  for(i in 1:times){
    errors[i] = calcError(all_pts, n, q, plot = plot, verbose = FALSE)
  }
  return(errors)
}

caclErrorCurve <- function(all_pts, times = 100, n = 100, maxQ = 99, plot = FALSE){
  meanErrors = rep(0,maxQ)
  for(i in 1:maxQ){
    print(i)
    meanErrors[i] = mean(calcErrors(pts, times = times, n = n, q = i))
  }
  
  return(meanErrors)
}

pts = read_pts_file(OutputPath = OutputPath, protName = "000_Trx")

calcError(pts, n = 100, q = 1, plot = TRUE, verbose = TRUE, 
          plotToFile = "/home/willy/PredictingProteinInteractions/Results/QuickRepSampling/Examples/ecdf_example1.pdf")


plotDistWithApproximationDifferentQ(OutputPath = OutputPath, protName = "000_Trx",n = 100, 
                                    qlist = c(1,2,3,10),
                                    plotToFile = "/home/willy/PredictingProteinInteractions/Results/QuickRepSampling/Examples/ecdf_exampleMultipleq1.pdf")


errors = calcErrors(pts, times = 100, n = 100, q = 1)
boxplot(errors, ylim = c(0,1))

errors = calcErrors(pts, times = 100, n = 100, q = 10)
boxplot(errors, ylim = c(0,1))

errors = calcErrors(pts, times = 100, n = 100, q = 99)
boxplot(errors, ylim = c(0,1))

meanErrorsQ = caclErrorCurve(pts, times = 100, n = 3, maxQ = 2)

plot(x = c(1:length(meanErrorsQ)), y = meanErrorsQ, type ="l", xlab = "q", ylab = "mean(Error)", main = "abs(DE_exact - DE_approx)")



calcError(pts,q = 10)

#-------------------------------------------------------------------



getAll_protein_F_approximations <- function(OutputPath,  n = 100, m = 50, q = 2, pos = TRUE){
  print(paste("Loading from ", OutputPath, sep =""))
  protein_names = list.dirs(OutputPath, recursive = FALSE, full.names = FALSE)
  
  distributions_lists = list()
  
  for(i in 1:length(protein_names)){
    F_app = generateF_approximations(OutputPath = OutputPath, protName = protein_names[i], n = n, m = m, q = q, pos =pos)
    distributions_lists[[i]] =  list("name" = protein_names[i],"F" = F_app)
  }
  
  return(distributions_lists)
}

get_x_and_y_vals <- function(all_protein_F_approximations,m){
  x_vals = rep(0,length(all_protein_F_approximations)*m)
  y_vals = rep(0,length(all_protein_F_approximations)*m)
  
  for(i in 1:length(all_protein_F_approximations)){
    for(j in 1:m){
      ind = (i-1)*m
      x_vals[ind+j] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
      y_vals[ind+j] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][2]-all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
    }
  }
  
  return(list("x_vals" = x_vals, "y_vals" = y_vals))
}

plot_all_F_approximations <- function(all_protein_F_approximations, colors = NULL, onlyGeo = FALSE, xli = NULL, yli = NULL){
  
  n_proteins = all_protein_F_approximations
  m = length(all_protein_F_approximations[[1]]$F$F_app_list)

  if(is.null(xli) || is.null(yli)){
    obj = get_x_and_y_vals(all_protein_F_approximations,m)
    x_vals = obj$x_vals
    y_vals = obj$y_vals
    
    xli = c(min(x_vals),max(x_vals))
    yli = c(min(y_vals),max(y_vals))
  }  

  if(is.null(colors)) colors = rep("black", length(n_proteins))
  
  plot2d_dist_approx(all_protein_F_approximations[[1]]$F$F_app_list, xli = xli, yli = yli, col = colors[1], name = all_protein_F_approximations[[1]]$name,onlyGeo = onlyGeo)
  for(i in 2:length(all_protein_F_approximations)){
    
    print(i)
    plot2d_dist_approx(all_protein_F_approximations[[i]]$F$F_app_list, add = TRUE, col = colors[i], name = all_protein_F_approximations[[i]]$name, onlyGeo = onlyGeo)
  }

}

getManhattanProjection2d <- function(all_protein_F_approximations){
  prot_number = length(all_protein_F_approximations)
  
  m = length(all_protein_F_approximations[[1]]$F$F_app_list)
  
  df = data.frame(matrix(0,nrow = m*prot_number, ncol = 3))
  colnames(df) = c("name", "x", "y")
  
  for(i in 1:prot_number) {
    print(all_protein_F_approximations[[i]]$name)
    name = all_protein_F_approximations[[i]]$name
    
    ind = (i-1)*m
    for(j in 1:m){
      df[ind+j,1] = name
      df[ind+j,2] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
      df[ind+j,3] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][2]-  all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
    }
    
  }
  
  return(df)
}

# length(all_protein_F_approximations[[1]]$F$F_app_list[[1]])

getManhattanProjection <- function(all_protein_F_approximations){
  prot_number = length(all_protein_F_approximations)
  
  m = length(all_protein_F_approximations[[1]]$F$F_app_list)
  q = length(all_protein_F_approximations[[1]]$F$F_app_list[[1]])
  
  df = data.frame(matrix(0,nrow = m*prot_number, ncol = q+1))
  colnames(df) = c("name", colnames(df)[2:q+1])
  
  for(i in 1:prot_number) {
    print(all_protein_F_approximations[[i]]$name)
    name = all_protein_F_approximations[[i]]$name
    
    ind = (i-1)*m
    for(j in 1:m){
      df[ind+j,1] = name
      
      df[ind+j,1+1] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
      for(k in 2:q){
        df[ind+j,k+1] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][k] - df[ind+j,q]
      }

    }
    
  }
  
  return(df)
}


quickRepSampling <- function(OutputPath, distance_path, fName, pos = TRUE, n = 100, m = 22, q = 21){
  # actually q should be possible to be larger than n
  #
  #-------------------------------------------------------
  
  all_protein_F_approximations = getAll_protein_F_approximations(OutputPath = OutputPath, m = m, q = q, n = n, pos = TRUE)
  
  df = getManhattanProjection(all_protein_F_approximations)
  distances = getDistanceMatrixFromProjection(df)
  
  geos = getGeometricCenters(df)
  
  geoDistances = getDistanceMatrixFromProjection(geos)
  
  # can be loaded with loadDistanceMatrix1()
  write.csv(geoDistances,file = paste(distance_path,"/",fName, ".csv",sep =""))
}


quickRepSampling(OutputPath = "/home/willy/Schreibtisch/106Test/Output/", distance_path = "/home/willy/Schreibtisch/106Test/Output/", fName = "distTest")

library(rbenchmark)
bench = benchmark("quickRepSampling_n_100_m_22_q_2" = quickRepSampling(OutputPath = "/home/willy/Schreibtisch/106Test/Output/",
                                                                distance_path = "/home/willy/Schreibtisch/106Test/Output/",
                                                                fName = "distTest",
                                                                n = 100,
                                                                m = 22,
                                                                q = 2),
          "quickRepSampling_n_100_m_22_q_20" = quickRepSampling(OutputPath = "/home/willy/Schreibtisch/106Test/Output/",
                                                                 distance_path = "/home/willy/Schreibtisch/106Test/Output/",
                                                                 fName = "distTest",
                                                                 n = 100,
                                                                 m = 22,
                                                                 q = 20),
          "quickRepSampling_n_100_m_100_q_20" = quickRepSampling(OutputPath = "/home/willy/Schreibtisch/106Test/Output/",
                                                                distance_path = "/home/willy/Schreibtisch/106Test/Output/",
                                                                fName = "distTest",
                                                                n = 100,
                                                                m = 100,
                                                                q = 20),
          replications = 2,
          columns = c("test", "replications", "elapsed",
                      "relative", "user.self", "sys.self"))

#---------------------------------------------------------------------------------------
#                                 test replications elapsed relative user.self sys.self
# 1  quickRepSampling_n_100_m_22_q_2              2  11.573    1.000    11.411    0.168
# 2 quickRepSampling_n_100_m_22_q_20              2  15.741    1.360    15.550    0.217
# 3 quickRepSampling_n_100_m_100_q_20             2  69.508    6.006    65.665    3.787
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# CompareProteinsCall
#---------------------------------------------------------------------------------------
compareProteinsCall = "/home/willy/PredictingProteinInteractions/MetricGeometry/ComparingProteins/./CompareIsosurfaces.R"
compareProteinsParameterFile = "/home/willy/PredictingProteinInteractions/MetricGeometry/Benchmark/Parameter"
system2(command = compareProteinsCall, args = compareProteinsParameterFile)



library(rbenchmark)
benchmark("CompareProteins" = system2(command = compareProteinsCall, args = compareProteinsParameterFile),
          replications = 1,
          columns = c("test", "replications", "elapsed",
                      "relative", "user.self", "sys.self"))

#---------------------------------------------------------------------------------------
n = 100
m = 100
q = 2
times = 2
OutputPath = "/home/willy/Schreibtisch/106Test/Output/"
all_protein_F_approximations = getAll_protein_F_approximations(OutputPath = OutputPath, m = m, q = q, n = n, pos = TRUE)

functionals = c(getFunctionalProteins(),"000_Trx")
colors = rep("blue", length(all_protein_F_approximations))
for(i in 1:length(all_protein_F_approximations)){
  if(all_protein_F_approximations[[i]]$name %in% functionals) colors[i] = "red"
}

# plot_all_F_approximations(all_protein_F_approximations, colors = colors, onlyGeo = TRUE, xli = c(16, 22), yli = c(1,6))


plot_all_F_approximations(all_protein_F_approximations, colors = colors, onlyGeo = TRUE)



# df = getManhattanProjection2d(all_protein_F_approximations)
df = getManhattanProjection(all_protein_F_approximations)
distances = getDistanceMatrixFromProjection(df)

geos = getGeometricCenters(df)

geos[7,]

xl = geos[7,2] + c(-0.1,1) 
yl = geos[7,3] + c(-1,1) 
unique(geos[intersect(intersect(intersect(which(geos[,2] < xl[2]), which(geos[,2] > xl[1])),which(geos[,3] < yl[2])), which(geos[,3] > yl[1])),1])

geoDistances = getDistanceMatrixFromProjection(geos)

geoDistances[7,]

Emd_distances = emd_parallel(DE_parallel = distances,m = m)


emdRef = readDistanceMatrix3(file = "/home/willy/Schreibtisch/106Test/Output/")
cor(c(Emd_distances[-106,-106]),c(emdRef))


closest = rep(0,106)
for(i in 1:length(closest)){
  closest[i] = length(which(rownames(geoDistances)[which.minn(geoDistances[7,],n = i)] %in% functionals) == TRUE)
}

emdRef = readDistanceMatrix2()
closest2 = rep(0,105)
for(i in 1:length(closest2)){
  closest2[i] = length(which(rownames(emdRef)[which.minn(emdRef[7,],n = i)] %in% functionals) == TRUE)
}
plot(closest2, type = "l", main = "closest neihgbours that are functional", xlab ="NN", ylab = "functionals found", col = "blue")
points(closest, type = "l", main = "closest neihgbours that are functional", xlab ="NN", ylab = "functionals found", col = "red")
legend(x = 35, y = 10,legend = c("RepSamp", "QuickRepSamp"),lty = c(1,1), col = c("blue", "red"))











