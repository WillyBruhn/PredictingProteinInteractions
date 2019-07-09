#-------------------------------------------------------------------------------------------------

plot_vol <- function(box_size, center){
  library(plot3D)
  
  x = prot_file$V6
  y = prot_file$V7
  z = prot_file$V8
  
  par(mfrow=c(1,1))
  points3D(x,y,z)
  
  # b = box_size/2.0
  
  # border3D(x0 = c(center[1] - b), 
  #          y0 = c(center[2] - b),
  #          z0 = c(center[3] - b),
  #          x1 = c(center[1] + b),
  #          y1 = c(center[2] + b),
  #          z1 = c(center[3] + b), add = TRUE, col = "black")
  
}


my_createisosurf <- function(dx_data, dx_meta, level, eps){
  
  # dx_file = paste(prot_name, "/", prot_name, "_pot.dx", sep = "")
  # 
  # dxData <- read.csv(file=paste(folder,"/",dx_file ,sep=""), sep=' ', skip= 11, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
  # 
  # dxData <- head(dxData,-5)
  # # sometimes 5 sometimes 10!!!
  # 
  # 
  # dx_meta = read.csv(file=paste(folder,"/",dx_file ,sep=""), sep=' ', skip= 0, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
  # 
  gridCount = as.numeric(dx_meta$V6[5])
  origin = as.numeric(c(dx_meta$V2[6], dx_meta$V3[6], dx_meta$V4[6]))
  dx = as.numeric(c(dx_meta$V2[7], dx_meta$V3[7], dx_meta$V4[7]))
  dy = as.numeric(c(dx_meta$V2[8], dx_meta$V3[8], dx_meta$V4[8]))
  dz = as.numeric(c(dx_meta$V2[9], dx_meta$V3[9], dx_meta$V4[9]))
  
  lx = gridCount*dx[1]
  ly = gridCount*dy[2]
  lz = gridCount*dz[3]
  
  
  
  
  v1 = as.numeric(dx_data$V1)
  v2 = as.numeric(dx_data$V2)
  v3 = as.numeric(dx_data$V3)
  
  merged = rep(0, 3*length(v1))
  
  for(i in 1:length(v1)){
    merged[1 + ((i-1)*3)] = v1[i]
    merged[1+1+((i-1)*3)] = v2[i]
    merged[1+2+((i-1)*3)] = v3[i]
  }
  
  # gridCount = 129
  # origin = c(-2.079555e+01, -5.045662e+01, -4.726355e+01)
  # dx = c(8.244773e-01, 0.000000e+00, 0.000000e+00)
  # dy = c(0.000000e+00, 6.956269e-01, 0.000000e+00)
  # dz = c(0.000000e+00, 0.000000e+00, 7.680477e-01)
  # 
  # lx = gridCount*dx[1]
  # ly = gridCount*dy[2]
  # lz = gridCount*dz[3]
  
  xVals = c()
  yVals = c()
  zVals = c()
  
  for(i in 1:gridCount){
    for(j in 1:gridCount){
      for(k in 1:gridCount){
        
        ind = ((i-1)*gridCount+(j-1))*gridCount+k
        # print(paste("testing ind ", ind ))
        if(abs(merged[ind] - level) < eps){
          # print(paste("found voxel with level = ", level, sep = ""))
          xVals = c(xVals, i*dx[1] + origin[1])
          yVals = c(yVals, j*dy[2] + origin[2])
          zVals = c(zVals, k*dz[3] + origin[3])
        }
      }
    }
  }
  
  df = data.frame(xVals, yVals, zVals)
  names(df) = c("x","y","z")
  return(df)
  
}

less_resolution <- function(n, points){
  less_res <- points[seq(1, nrow(points), n),]
  return(less_res)
}
plot_with_less_res <- function(theta, phi, n,points, add_flag = FALSE, col = "blue"){
  less_res = less_resolution(n,points)
  points3D(less_res[,1], less_res[,2], less_res[,3], add = add_flag, col = col, theta = theta, phi = phi, cex = 0.01)
}

# # path_to_centerSelect = "/home/sysgen/Documents/LWB/centerSelect/"
# folder = "/home/willy/RedoxChallenges/Redox_old/Output/"
# prot_name="000_Ars"
# 
# dx_file = paste(prot_name, "/", prot_name, "_pot.dx", sep = "")

convertToMatrix <- function(iso){
  X = matrix(rep(0,3*nrow(iso)), nrow = nrow(iso))
  
  for(i in 1:nrow(iso)){
    X[i,1] = as.numeric(iso[i,1])
    X[i,2] = as.numeric(iso[i,2])
    X[i,3] = as.numeric(iso[i,3])
  }
  
  return(X)
}

getPtsFileName <- function(folder,prot_name,positive = TRUE){
  ext = "neg"
  
  if(positive == TRUE){
    ext = "pos"
  }
  
  pts_file = paste(prot_name,"_",ext,".pts",sep="")
  
  return(paste(folder,"/",prot_name,"/",pts_file ,sep=""))
}

readPtsFile <- function(prot_name){
  prot <- read.csv(file=prot_name, sep=';', skip= 0, header=T ,stringsAsFactors=FALSE,  check.names = FALSE, quote ="", dec = ",")
  # ?read.csv
  return(prot)
}


getIsoSurfaceObj <- function(folder,prot_name,eps = 0.1, onlyEveryNth = 1){
  
  iso_pos_name = getPtsFileName(folder,prot_name)
  iso_neg_name = getPtsFileName(folder,prot_name,FALSE)
  
  iso_pos = NULL
  iso_neg = NULL
  
  # 10.1.19
  readFromDx = TRUE
  
  if(readFromDx == TRUE || (!file.exists(iso_pos_name) || !file.exists(iso_neg_name))){
    myPrint("reading from dx ...")
    
    dx_file = paste(prot_name, "/", prot_name, "_pot.dx", sep = "")
    # dx_file = paste(prot_name, "_pot.dx", sep = "")
  
    dxData <- read.csv(file=paste(folder,"/",dx_file ,sep=""), sep=' ', skip= 11, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
  
    dxData <- head(dxData,-5)
    # sometimes 5 sometimes 10!!!
    
    dx_meta = read.csv(file=paste(folder,"/",dx_file ,sep=""), sep=' ', skip= 0, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
    
    
    iso_pos = my_createisosurf(dx_data = dxData, dx_meta, 1.0, eps = eps)
    iso_neg = my_createisosurf(dx_data = dxData, dx_meta, -1.0, eps = eps)
    
    
    write.csv2(iso_pos, file = iso_pos_name, sep = ";", row.names = FALSE)
    write.csv2(iso_neg, file = iso_neg_name, sep = ";", row.names = FALSE)
  } else {
    myPrint("reading from pts ...")
    
    iso_pos = readPtsFile(iso_pos_name)
    iso_neg = readPtsFile(iso_neg_name)
  }

  
  iso_pos = convertToMatrix(iso_pos)
  iso_neg = convertToMatrix(iso_neg)
  
  
  iso_pos = less_resolution(onlyEveryNth,iso_pos)
  iso_neg = less_resolution(onlyEveryNth,iso_neg)
  
  if(FALSE){
  myPrint("Finding nearest neighbours in other ptential ...")
  
  
  NNpositive_name = paste(folder,"/",prot_name,"/",prot_name, "NN_positive.txt", sep = "")
  NNnegative_name = paste(folder,"/",prot_name,"/",prot_name, "NN_negative.txt", sep = "")
  
  NNpositive = NULL
  NNnegative = NULL
  
  if(!file.exists(NNpositive_name) || !file.exists(NNnegative_name)){
    myPrint("Creating NNs ...")
    
    NNpositive = NNInOtherPot(iso_pos,iso_neg)
    NNnegative = NNInOtherPot(iso_neg,iso_pos)
    
    write.table(NNnegative, file = NNnegative_name, row.names = FALSE)
    write.table(NNpositive, file = NNpositive_name, row.names = FALSE)
  } else {

  }
  
    myPrint("reading NNs ...")
    
    NNpositive = read.table(file = NNpositive_name, header = TRUE, dec = ",")
    NNnegative = read.table(file = NNnegative_name, header = TRUE, dec = ",")
  }

  
  # ret = list("name" = prot_name, "positive" = iso_pos, "negative" = iso_neg, "NN_negative" = NNnegative, "NN_positive" = NNpositive)
  
  ret = list("name" = prot_name, "positive" = iso_pos, "negative" = iso_neg)
  
  return(ret)
}


#--------------------------------------------------------------------------------
# iso = my_createisosurf(dx_data = dxData, 1.0, eps = 0.1)
# iso2 = my_createisosurf(dx_data = dxData, -1.0, eps = 0.1)
# 
# iso_less = less_resolution(1,iso)
# iso_less2 = less_resolution(1,iso2)


argmin <- function(v){
  return(which(v == min(v))[1])
}

NNInOtherPot <- function(pot1,pot2){
  positive_minDistances = rep(-1,nrow(pot1))
  for(i in 1:nrow(pot1)){
    distances = rep(100000,nrow(pot2))
    
    for(j in 1:nrow(pot2)){
      distances[j] = euklid_dist(pot1[i,],pot2[j,])
    }
    
    positive_minDistances[i] = argmin(distances)
  }
  
  return(positive_minDistances)
}

# getNearestNeighbourInOtherPotential <- function(pot1, pot2){
#   
#   negative_minDistances = rep(-1,nrow(prot$negative))
#   for(i in 1:nrow(prot$negative)){
#     distances = rep(100000,nrow(prot$positive))
#     
#     for(j in 1:nrow(prot$positive)){
#       distances[j] = euklid_dist(prot$positive)
#     }
#     
#     negative_minDistances[i] = argmin(distances)
#   }
#   
#   ret = list("NN_negative" = negative_minDistances, "NN_positive" = positive_minDistances)
#   
#   return(ret)
# }



plotIsoObj <- function(obj, add = FALSE, secondPlot = FALSE, size = 3.0, extraCol = FALSE){
  if(secondPlot == TRUE){
    open3d()
  }
  
  col = c("red", "blue")
  if(extraCol == TRUE){
    col = c("green", "yellow")
  }

  plot3d(obj$positive, col = col[1], add = add, size = size)
  plot3d(obj$negative, col = col[2], add = TRUE, size = size)
}

plot2IsoObj <- function(X, x_pcs, Y, Y_aligned, y_pcs){
  print("plot2IsoObj ...")
  
  # first close all other windows
  # while (rgl.cur() > 0) { rgl.close() }
  
  # limit = max(max(abs(X$transformed_data)), max(abs(Y$transformed_data))) + 2
  # 
  # limits = c(-limit,limit)
  size = 700
  smallPlotSize = size/2
  offSet = 2000
  # ?par3d()
  # ?open3d()
  open3d()
  plotIsoObj(X)
  
  par3d(windowRect = c(1230,95,1580,445))
  # par3d(windowRect = c(offSet-size,0,offSet-smallPlotSize,smallPlotSize))
  
  open3d()
  plotIsoObj(Y)
  # par3d(windowRect = c(offSet-smallPlotSize,0,offSet,smallPlotSize))
  par3d(windowRect = c(1580,95,1930,445))
  
  # par3d(windowRect = c(offSet-size,offSet-size/2,offSet,size))
  open3d()
  
  # Y_aligned = alignPcaObjs(X = X, Y = Y)
  
  plotIsoObj(X)
  plotIsoObj(Y_aligned,add = TRUE,secondPlot = FALSE,size = 3.0, extraCol = TRUE)
  
  
  
  # par3d(windowRect = c(1580,95,1930,445))
  par3d(windowRect = c(1230,480,1930,1110))
  
}

centerIsoObj <- function(X,center){
  X = translateIsoObj(X,-center)
  return(X)
}

centerIsoObjOnPositive <- function(X){
  X = translateIsoObj(X,-getCentroid(X$positive))
  return(X)
}

centerIsoObjOnNegative <- function(X){
  X = translateIsoObj(X,-getCentroid(X$negative))
  return(X)
}

translateIsoObj <- function(X,v){
  X$positive = translate(X$positive, v)
  X$negative = translate(X$negative, v)
  
  return(X)
}

alignIsoObjWithVectors <- function(X,u,v){
  
  X$positive = alignDataWithVectors(X$positive,u, v)
  X$negative = alignDataWithVectors(X$negative,u, v)
  
  return(X)
}

alignIsoObjWithVectorsInZ0Plane <- function(X,u,v){
  X$positive = alignDataWithVectorsInZ0Plane(X$positive, u, v)
  X$negative = alignDataWithVectorsInZ0Plane(X$negative, u, v)
  
  return(X)
}

alignIsoObjWithFeature <- function(Y,pcs,feature,alignPos){
  feature_center = c()
  if(alignPos == TRUE){
    feature_center = getCentroid(Y$positive[feature,])
  } else  {
    feature_center = getCentroid(Y$negative[feature,])
  }
  
  # center the data on the center of the feature
  Y = translateIsoObj(Y,-feature_center)
  
  # Y_b = Y
  Y_a = alignIsoObjWithVectors(Y,pcs[,3], diag(3)[,3])
  Y = alignIsoObjWithVectorsInZ0Plane(Y_a,pcs[1:2,1], diag(3)[1:2,1])
  
  return(Y)
}

alignIsoObjs <- function(x_pcs,y_pcs, Y ){
  
  Y_a = alignIsoObjWithVectors(Y,x_pcs[,3], y_pcs[,3])
  Y_b = alignIsoObjWithVectorsInZ0Plane(Y_a,x_pcs[1:2,1], y_pcs[1:2,1])
  

  return(Y_b)
}


sampleIsoObj <- function(prot, n, positiveDistribution = c(1:nrow(prot$positive)), negativeDistribution = c(1:nrow(prot$negative))){
  
  # x_indices = sample(c(1:nrow(X)),prob = 1/X_scores,size = sampleSize)
  
  # pos1 = runif(n/2,1,nrow(prot$positive))
  # neg1 = runif(n/2,1,nrow(prot$negative))
  
  pos1 = sample(c(1:nrow(prot$positive)),prob = positiveDistribution, size = n)
  neg1 = sample(c(1:nrow(prot$negative)),prob = negativeDistribution, size = n)
  
  pos2 = prot$NN_negative[neg1,1]
  neg2 = prot$NN_positive[pos1,1]
  
  # print(paste(length(pos1)))
  # print(paste(length(pos2)))
  # print(paste(length(neg1)))
  # print(paste(length(neg2)))
  
  # sample = prot
  # sample$positive = prot$positive[c(pos1,pos2),]
  # sample$negative = prot$negative[c(neg1,neg2),]
  
  
  # ret = list("name" = prot$name, "positive" = iso_pos, "negative" = iso_neg, "NN_negative" = NNnegative, "NN_positive" = NNpositive)
  
  # ret = list("name" = prot$name, "positive" = prot$positive[c(pos1,pos2),], "negative" = prot$negative[c(neg1,neg2),])
  
  ret = list("name" = prot$name, "positive" = prot$positive[round(c(pos1,pos2)),], "negative" = prot$negative[round(c(neg1,neg2)),], "positive_ind" = round(c(pos1,pos2)), "negative_ind" = round(c(neg1,neg2)))
  return(ret)
}


compareProteins <- function(prot1,prot2,sampleSize = 10, repititions = 10000, featureSize = 100, draw = FALSE){
  s = subSamplePseudoPca(prot1, prot2,sampleSize = sampleSize,repitions = repititions,numOfPointsToPlot = featureSize,1.0,FALSE)
  
  return(compare2DataPosAndNeg(s,featureSize,draw = draw,drawAllPoints = TRUE))
}

readInAllProteinsInFolder <- function(folder,proteinNames){
  print("reading in all proteins ..")
  protList = rep(prot1,length(proteinNames))
  for(i in length(proteinNames)){
    print(paste("reading in ", proteinNames[i]))
    prot = getIsoSurfaceObj(folder,proteinNames[i],0.001,1)
    
    protList[i] = prot
  }
  
  return(protList)
}

compare2DataPosAndNeg <- function(comparisonObj, feature_length = 10, draw = TRUE, drawAllPoints = TRUE){
  
  # if(drawAllPoints == FALSE){
  #   X_feature_obj$positive = comparisonObj$X$positive[comparisonObj$X_pos_ordered[c(1:feature_length)],]
  #   X_feature_obj$negative = comparisonObj$X$negative[comparisonObj$X_neg_ordered[c(1:feature_length)],]
  #   
  #   
  #   Y_feature_obj$positive = comparisonObj$Y$positive[comparisonObj$Y_pos_ordered[c(1:feature_length)],]
  #   Y_feature_obj$negative = comparisonObj$Y$negative[comparisonObj$Y_neg_ordered[c(1:feature_length)],]
  # }
  
  X_feature_obj = comparisonObj$X
  Y_feature_obj = comparisonObj$Y
  
  X_feature_obj$positive = comparisonObj$X$positive[comparisonObj$X_pos_ordered[c(1:feature_length)],]
  X_feature_obj$negative = comparisonObj$X$negative[comparisonObj$X_neg_ordered[c(1:feature_length)],]
  
  Y_feature_obj$positive = comparisonObj$Y$positive[comparisonObj$Y_pos_ordered[c(1:feature_length)],]
  Y_feature_obj$negative = comparisonObj$Y$negative[comparisonObj$Y_neg_ordered[c(1:feature_length)],]
  
  X_feature_pos_vec = comparisonObj$X_pos_ordered[c(1:feature_length)]
  X_feature_neg_vec = comparisonObj$X_neg_ordered[c(1:feature_length)]
  
  Y_feature_pos_vec = comparisonObj$Y_pos_ordered[c(1:feature_length)]
  Y_feature_neg_vec = comparisonObj$Y_neg_ordered[c(1:feature_length)]
  
  
  r_positive = compare2Data(X_feature_obj$positive,Y_feature_obj$positive,draw)
  r_negative = compare2Data(X_feature_obj$negative,Y_feature_obj$negative,draw)
  
  
  r = c()
  if(draw == TRUE){
    
    if(r_negative$matching$sum < r_positive$matching$sum){
      
      r_negative_aligned = alignPcaObjs(r_negative$X_pca, r_negative$Y_pca)
      
      # X now centered and rotateted accordingly to the feature of X
      X = alignIsoObjWithFeature(comparisonObj$X, r_positive$X_pca$pcs, X_feature_neg_vec, alignPos = FALSE)
      
      # Y now centered and rotateted accordingly to the feature of Y
      Y = alignIsoObjWithFeature(comparisonObj$Y, r_negative$Y_pca$pcs, Y_feature_neg_vec, alignPos = FALSE)
      Y = alignIsoObjs(r_negative$X_pca$pcs, r_negative$Y_pca$pcs, Y)
      
      
      
      plot2IsoObj(X, r_negative$X_pca$pcs, Y , Y, r_negative$Y_pca$pcs)
      
      drawCoordArrows(center = c(0,0,0), mat = r_negative$X_pca$pcs, col = "purple", name = "")
      drawCoordArrows(center = c(0,0,0), mat = r_negative_aligned$pcs, col = "black", name = "")
      
      
      myPrint("negative ...")
      
      
      
    } else{
      r_positive_aligned = alignPcaObjs(r_positive$X_pca, r_positive$Y_pca)
      
      # X now centered and rotateted accordingly to the feature of X
      X = alignIsoObjWithFeature(comparisonObj$X, r_positive$X_pca$pcs, X_feature_pos_vec, alignPos = TRUE)
      
      # Y now centered and rotateted accordingly to the feature of Y
      Y = alignIsoObjWithFeature(comparisonObj$Y, r_positive$Y_pca$pcs, Y_feature_pos_vec, alignPos = TRUE)
      Y = alignIsoObjs(r_positive$X_pca$pcs, r_positive$Y_pca$pcs, Y)
      
      
      plot2IsoObj(X, r_positive$X_pca$pcs, Y , Y, r_positive$Y_pca$pcs)
      
      drawCoordArrows(center = c(0,0,0), mat = r_positive$X_pca$pcs, col = "purple", name = "")
      drawCoordArrows(center = c(0,0,0), mat = r_positive_aligned$pcs, col = "black", name = "")
      
      myPrint("positive ...")
    }
  }
  
  v = min(r_negative$matching$sum, r_positive$matching$sum) 
  return(v)
}

plotSamplingExample <-function(prot, num = 10){
  samp = sampleIsoObj(prot,num)
  plotIsoObj(prot)
  plotIsoObj(samp,add = TRUE, size = 12.0)
}


subSamplePseudoPca <- function(X,Y,sampleSize,repitions, numOfPointsToPlot, pseudoCount = 1.0, addFeature = TRUE, draw = FALSE){
  # r_base_positive = compare2Data(X$positive,Y$positive,FALSE)
  # r_base_negative = compare2Data(X$negative,Y$negative,FALSE)
  # myPrint("Done with base...")
  #-----------------------------------------
  
  X_scores_positive = rep(pseudoCount,nrow(X$positive))
  x_num_positive = rep(0,nrow(X$positive))
  
  Y_scores_positive = rep(pseudoCount,nrow(Y$positive))
  y_num_positive = rep(0,nrow(Y$positive))
  
  X_scores_negative = rep(pseudoCount,nrow(X$negative))
  x_num_negative = rep(0,nrow(X$negative))
  
  Y_scores_negative = rep(pseudoCount,nrow(Y$negative))
  y_num_negative = rep(0,nrow(Y$negative))
  
  iter = 0
  for(i in 1:repitions){
    iter = iter +1
    myPrint(paste(i/repitions*100, "%"))
    
    # x_indices = sample(c(1:nrow(X)),prob = 1/X_scores,size = sampleSize)
    # X_sample = X[x_indices,]
    # y_indices = sample(c(1:nrow(Y)),prob = 1/Y_scores, size = sampleSize)
    # Y_sample = Y[y_indices,]
    
    # only uniform distribution
    # X_sample = sampleIsoObj(X,sampleSize, 1/X_scores_positive, 1/X_scores_negative)
    # Y_sample = sampleIsoObj(Y,sampleSize, 1/Y_scores_positive, 1/Y_scores_negative)
    
    X_sample = sampleIsoObj(X,sampleSize, seq(1,length(X_scores_positive)), seq(1,length(X_scores_negative)))
    Y_sample = sampleIsoObj(Y,sampleSize, seq(1,length(Y_scores_positive)), seq(1,length(Y_scores_negative)))
    
    
    # r = compare2Data(X_sample,Y_sample,FALSE)
    r_positive = compare2Data(X_sample$positive,Y_sample$positive,FALSE)
    r_negative = compare2Data(X_sample$negative,Y_sample$negative,FALSE)
    
    # vs = r_positive$matching$sum+r_negative$matching$sum+0.1
    
    vs = min(r_positive$matching$sum, r_negative$matching$sum + 0.0001)
    
    # positive
    for(j in 1:length(X_sample$positive_ind)){
      k = X_sample$positive_ind[j]
      
      x_num_positive[k] = x_num_positive[k] + 1
      
      X_scores_positive[k] = (X_scores_positive[k]*(x_num_positive[k]-1) + r_positive$matching$distances[j]/vs)/x_num_positive[k]
    }
    
    for(j in 1:length(Y_sample$positive_ind)){
      k=Y_sample$positive_ind[j]
      
      y_num_positive[k] = y_num_positive[k] + 1
      
      Y_scores_positive[k] = (Y_scores_positive[k]*(y_num_positive[k]-1) + r_positive$matching$inv_distances[j]/vs)/y_num_positive[k]
    }
    
    # negative
    for(j in 1:length(X_sample$negative_ind)){
      k = X_sample$negative_ind[j]
      
      x_num_negative[k] = x_num_negative[k] + 1
      
      X_scores_negative[k] = (X_scores_negative[k]*(x_num_negative[k]-1) + r_negative$matching$distances[j]/vs)/x_num_negative[k]
    }
    
    for(j in 1:length(Y_sample$negative_ind)){
      k=Y_sample$negative_ind[j]
      
      y_num_negative[k] = y_num_negative[k] + 1
      
      Y_scores_negative[k] = (Y_scores_negative[k]*(y_num_negative[k]-1) + r_negative$matching$inv_distances[j]/vs)/y_num_negative[k]
    }
  }
  
  myPrint("myPrinting positive ...")
  myPrint(X_scores_positive[order(X_scores_positive)])
  myPrint(x_num_positive[order(X_scores_positive)])
  
  myPrint(Y_scores_positive[order(Y_scores_positive)])
  myPrint(y_num_positive[order(Y_scores_positive)])
  
  X_feature_positive = order(X_scores_positive)[1:numOfPointsToPlot]
  Y_feature_positive = order(Y_scores_positive)[1:numOfPointsToPlot]
  
  myPrint("myPrinting negative ...")
  myPrint(X_scores_negative[order(X_scores_negative)])
  myPrint(x_num_negative[order(X_scores_negative)])
  
  myPrint(Y_scores_negative[order(Y_scores_negative)])
  myPrint(y_num_negative[order(Y_scores_negative)])
  
  X_feature_negative = order(X_scores_negative)[1:numOfPointsToPlot]
  Y_feature_negative = order(Y_scores_negative)[1:numOfPointsToPlot]
  
  
  if(draw == TRUE){
    plot2ProtFeaturesExtracted(X, Y, X_feature_positive, X_feature_negative, Y_feature_positive, Y_feature_negative,addFeature)
  }
  
  
  feature_distance = (sum(X_feature_negative) + sum(Y_feature_negative) + sum(X_feature_positive) + sum(Y_feature_positive))/(4*numOfPointsToPlot + repitions + sampleSize)
  
  
  return(list("feature_distance" = feature_distance, "sampleSize" = sampleSize, "repitions" = repitions, 
              "featureSize" = numOfPointsToPlot, "pseudoCount" = pseudoCount, 
              "X_pos_ordered" = order(X_scores_positive),
              "X_neg_ordered" = order(X_scores_negative),
              "Y_pos_ordered" = order(Y_scores_positive),
              "Y_neg_ordered" = order(Y_scores_negative),
              "X" = X,
              "Y" = Y))
}

plot2ProtFeaturesExtracted <- function(X, Y, X_feature_positive, X_feature_negative, Y_feature_positive, Y_feature_negative, add = TRUE){
  # while (rgl.cur() > 0) { rgl.close() }
  
  open3d()
  plotIsoObj(X)
  Xs = X
  Xs$positive = X$positive[X_feature_positive,]
  Xs$negative = X$negative[X_feature_negative,]
  plotIsoObj(Xs,add = add, size = 12.0,extraCol = TRUE)
  text3d(0,0,0,text=Xs$name,cex = 1, col = "red")
  
  par3d(windowRect = c(900,95,1400,600))
  
  open3d()
  plotIsoObj(Y)
  Ys = Y
  Ys$positive = Y$positive[Y_feature_positive,]
  Ys$negative = Y$negative[Y_feature_negative,]
  plotIsoObj(Ys,add = add, size = 12.0,extraCol = TRUE)
  text3d(0,0,0,text=Ys$name,cex = 1, col = "red")
  
  
  par3d(windowRect = c(1400,95, 1900,  600))
  
  # par3d(fontname = "Title")
  
  # par3d("windowRect")
}

plot2ProtFeaturesExtractedArbitrary <- function(comparisonObj, length = 10, add = FALSE){
  plot2ProtFeaturesExtracted(comparisonObj$X, comparisonObj$Y, comparisonObj$X_pos_ordered[c(1:length)], comparisonObj$X_neg_ordered[c(1:length)], comparisonObj$Y_pos_ordered[c(1:length)], comparisonObj$Y_neg_ordered[c(1:length)], add = add)
}


plot2Proteins <- function(prot1, prot2){
  plotIsoObj(prot1)
  plotIsoObj(prot2, add = FALSE, secondPlot = TRUE)
}


