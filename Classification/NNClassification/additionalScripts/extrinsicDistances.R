# # https://stackoverflow.com/questions/31820865/error-in-installing-rgl-package
# install.packages("sphereplot")
# install.packages("rgl")

# install.packages("sphereplot")
# library(sphereplot)

# install.packages("plot3D")
# library(plot3D)

# install.packages("clue")
library(clue)

# install.packages("pracma")
library(pracma)

GLOBAL_VERBOSITY = TRUE
myPrint <- function(string){
  if(GLOBAL_VERBOSITY == TRUE){
    print(string)
  }
}


getPolarAngle <- function(u){
  # calculates the polar angle of a given vector
  # in the x,y plane
  
  u_n = u[1]/norm_vec(u)
  dot_prod = u_n
  
  theta = acos(dot_prod)
  
  if(u[2] < 0){
    theta = 2*pi - theta 
  }
  
  return(theta)
}


getRotMatrixToAlign2Vec <- function(u,v){
  
  u_p = getPolarAngle(u)
  v_p = getPolarAngle(v)
  
  theta = u_p - v_p
  
  c_t = cos(theta)
  s_t = sin(theta)
  
  A = matrix(c(c_t, -s_t, s_t, c_t), byrow = TRUE, nrow = 2)
  
  return(A)
}

norm_vec <- function(x) sqrt(sum(x^2))

euklid_dist <- function(a,b){
  return(norm_vec(b-a))
}

getCentroid <- function(X){
  center=c(mean(X[,1]),mean(X[,2]), mean(X[,3]))
  return(center)
}

centerData <- function(X){
  X_c = getCentroid(X)
  
  X_centered = X
  X_centered[,1] = X[,1] - X_c[1]
  X_centered[,2] = X[,2] - X_c[2]
  X_centered[,3] = X[,3] - X_c[3]
  
  return(X_centered)
}

#-----------------------------------------------------------------------------------------------
# 3 D
drawCoordArrows <- function(center, mat, col, name){
  scale = 15
  mat= scale*mat
  
  arrowHead_1 = center + mat[,1]
  arrowHead_2 = center + mat[,2]
  arrowHead_3 = center + mat[,3]
  
  arrow3d(p0 = center, p1 = arrowHead_1, col = col, thickness = 0.1, width = 0.1, type = "rotation")
  arrow3d(p0 = center, p1 = arrowHead_2, col = col, thickness = 0.1, width = 0.1, type = "rotation")
  arrow3d(p0 = center, p1 = arrowHead_3, col = col, thickness = 0.1, width = 0.1, type = "rotation")
  
  ts = 1.3
  t1 = center + mat[,1]*ts
  t2 = center + mat[,2]*ts
  t3 = center + mat[,3]*ts
  
  text3d(t1, text = "Pc1", col = col)
  text3d(t2, text = "Pc2", col = col)
  text3d(t3, text = "Pc3", col = col)
  
  text3d(center, text = name, col, font = 5)
}

translate <- function(X,v){
  
  X[,1] = X[,1] + v[1]
  X[,2] = X[,2] + v[2]
  X[,3] = X[,3] + v[3]
  
  return(X)
}

plotPcaObj <- function(X, col = "red", addFlag = TRUE, limits = c(1,1), s = 3.0){
  plot3d(X$transformed_data, add = addFlag, col = col, xlim = limits, ylim = limits, zlim = limits, main = X$name, size = s)
  drawCoordArrows(getCentroid(X$transformed_data),X$pcs, col = col, name = X$name)
}

plot2Vol <- function(X,Y, Y_aligned,theta=30, phi=30){
  # first close all other windows
  while (rgl.cur() > 0) { rgl.close() }
  
  limit = max(max(abs(X$transformed_data)), max(abs(Y$transformed_data))) + 2
  
  limits = c(-limit,limit)
  size = 700
  smallPlotSize = size/2
  offSet = 2000
  # ?par3d()
  # ?open3d()
  plotPcaObj(X,"red", FALSE, limits)
  par3d(windowRect = c(1230,95,1580,445))
  # par3d(windowRect = c(offSet-size,0,offSet-smallPlotSize,smallPlotSize))
  
  open3d()
  plotPcaObj(Y,"blue", FALSE, limits)
  # par3d(windowRect = c(offSet-smallPlotSize,0,offSet,smallPlotSize))
  par3d(windowRect = c(1580,95,1930,445))

  # par3d(windowRect = c(offSet-size,offSet-size/2,offSet,size))
  open3d()
  
  # Y_aligned = alignPcaObjs(X = X, Y = Y)
  
  plotPcaObj(X,"red", FALSE, limits)
  plotPcaObj(Y_aligned,"blue", TRUE, limits)
  
  # par3d(windowRect = c(1580,95,1930,445))
  par3d(windowRect = c(1230,480,1930,1110))
  
}

plot2VolFeaturesExtracted <- function(X,Y, X_feature, Y_feature){
  # first close all other windows
  while (rgl.cur() > 0) { rgl.close() }
  
  limit = max(max(abs(X$transformed_data)), max(abs(Y$transformed_data))) + 2
  
  limits = c(-limit,limit)
  size = 1400
  smallPlotSize = size/2
  offSet = 2000
  # ?par3d()
  # ?open3d()
  plotPcaObj(X,"red", FALSE, limits)
  plot3d(X$transformed_data[X_feature,], add = TRUE, col = "green", size = 12.0)
  
  
  # par3d(windowRect = c(1230,95,1580,445))
  par3d(windowRect = c(offSet-size,0,offSet-smallPlotSize,smallPlotSize))
  
  open3d()
  plotPcaObj(Y,"blue", FALSE, limits)
  plot3d(Y$transformed_data[Y_feature,], add = TRUE, col = "green", size = 12.0)
  par3d(windowRect = c(offSet-smallPlotSize,0,offSet,smallPlotSize))
  # par3d(windowRect = c(1580,95,1930,445))
  
}

rotateAroundX <- function(X,theta){
  c_t = cos(theta)
  s_t = sin(theta)
  R1 = matrix(c(1, 0, 0,
                0,c_t,-s_t,
                0,s_t,c_t
                ), byrow = TRUE, nrow = 3)
  return(X%*%R1)
}

rotateAroundZ <- function(X,theta){
  c_t = cos(theta)
  s_t = sin(theta)
  R1 = matrix(c(c_t, -s_t, 0,
                s_t,c_t,0,
                0,0,1
  ), byrow = TRUE, nrow = 3)
  return(X%*%R1)
}

isVecSame <- function(x,y){
  return(all(length(x)==length(y)) && all(x==y))
}

rotation = function(x,y){
  
  if(isVecSame(x,y)){
    return(diag(length(x)))
  } 
  
  u=x/sqrt(sum(x^2))
  
  v=y-sum(u*y)*u
  
  # print(v)
  if(norm_vec(v) != 0){
    v = v/norm_vec(v)
  }
  # v=v/sqrt(sum(v^2))
  # print(v)
    
  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
  
  sint=sqrt(1-cost^2);

  # print(u)
  # print(v)
  # print(sint)
    
  diag(length(x)) - u %*% t(u) - v %*% t(v) + 
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
}

# # R = matrix(c(0.8882616, -0.4593379,    0,0.4593379,  0.8882616,    0,0.0000000,0.0000000,1), nrow = 3, byrow = TRUE)
# 
#   
# 
# plot(x = NULL, xlim = c(-5,5), ylim = c(-5,5))
# arrows(0,0,R[1,1],R[2,1])
# arrows(0,0,-R[1,1],-R[2,1],col = "blue")
# 
# arrows(0,0,0,1, col = "red")
# rot = R[,1]%*%rotation(R[,1],-R[,1])
# norm_vec(rot)
# 
# c(0,1,0)%*%rotation(R[,1],-R[,1])
# arrows(0,0,rot[1],rot[2],col = "green")


rotateAroundArbitraryAxis <- function(X,axis,theta){
  # u is a vector
  
  ct = cos(theta)
  ctM = 1-ct
  st = sin(theta)
  
  u = axis/norm_vec(axis)
  
  R = matrix(c(
    ct + u[1]^2*ctM,       u[1]*u[2]*ctM - u[3]*st, u[1]*u[3]*ctM+u[2]*st,
    u[1]*u[2]*ctM+u[3]*st, ct + u[2]^2*ctM,         u[2]*u[3]*ctM-u[1]*st,
    u[3]*u[1]*ctM-u[2]*st, u[3]*u[2]*ctM+u[1]*st,   ct+u[3]^2*ctM
    ), nrow = 3, byrow = TRUE)
  
  return(X%*%R)
}

# R = matrix(c(0.8882616, -0.4593379,    0,0.4593379,  0.8882616,    0,0.0000000,0.0000000,1), nrow = 3, byrow = TRUE)
# plot(x = NULL, xlim = c(-5,5), ylim = c(-5,5))
# arrows(0,0,R[1,1],R[2,1])
# arrows(0,0,-R[1,1],-R[2,1],col = "blue")
# 
# arrows(0,0,0,1, col = "red")
# rot = rotateAroundArbitraryAxis(R[,1],R[,3],pi)
# 
# arrows(0,0,rot[1],rot[2],col = "green")
# 
# 
# for(i in 1:100){
#   print(paste("Test:", i))
#   signs = sign(runif(2,-1.0,1.0))
#   v = c(runif(3,1.0,2.0)*signs[1])
#   u = c(runif(3,1.0,2.0)*signs[2])
#   
#   v = v/norm_vec(v)
#   u = u/norm_vec(u)
#   
#   Rv2u = rotation(v,u)
#   v_t = v %*% Rv2u
#   
#   # print(Rv2u)
#   
#   if(sum(u-v_t) > 0.000001){
#     print("Wrong")
#     print(u-v_t)
# 
#     open3d()
#     arrow3d(p0 = c(0,0,0), p1 = v, col = "red", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = u, col = "blue", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = v_t, col = "green", thickness = 0.1, width = 0.1, type = "rotation")
#     
#     e = diag(3)
#     arrow3d(p0 = c(0,0,0), p1 = e[,1], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = e[,2], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = e[,3], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     return()
#   }
# }

alignDataWithVectors <- function(X,u,v){

  Rv2u = rotation(v,u)
  
  v_t = v %*% Rv2u

  myPrint(paste("u", u[1], u[2], u[3]))
  myPrint(paste("v", v[1], v[2], v[3]))
  
  myPrint(paste("v_t", v_t[1], v_t[2], v_t[3]))
  
  if((v_t[1]*v[1] + v_t[2]*v[2] + v_t[3]*v[3]) == -1){
    myPrint("dot prod is negative")
    Rv2u = rotation(-v,u)
  }
    
  if(sum(u-v_t) > 0.000001){
    print("Wrong")
    print(u-v_t)
    return(FALSE)
  }
  
  X = X %*% Rv2u
  
  return(X)
}

# alignDataWithVectors(X_b,c(0,0,1), c(0,1,0))
# 
# for(i in 1:100){
#   print(paste("Test:", i))
#   signs = sign(runif(2,1.0,1.0))
#   v1 = c(runif(3,1.0,2.0)*signs[1])
#   v2 = c(runif(3,1.0,2.0)*signs[2])
#   
#   v1 = v1/norm_vec(v1)
#   v2 = v2/norm_vec(v2)
#   
#   if(FALSE == alignDataWithVectors(X, v1, v2)){
#     open3d()
#     arrow3d(p0 = c(0,0,0), p1 = v1, col = "red", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = v2, col = "blue", thickness = 0.1, width = 0.1, type = "rotation")
# 
#     e = diag(3)
#     arrow3d(p0 = c(0,0,0), p1 = e[,1], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = e[,2], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = e[,3], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     return()
#   }
# }

alignDataWithVectorsInZ0Plane <- function(X,u,v){
  # u_polar = cart2pol(u)
  # v_polar = cart2pol(v)
  
  u_polar = getPolarAngle(u)
  v_polar = getPolarAngle(v)
  
  myPrint(u_polar)
  myPrint(v_polar)
  
  # phi = -(u_polar[1]-v_polar[1])
  phi = -(u_polar-v_polar)
  
  myPrint(paste("phi", phi))
  # print(paste("u", u[1], u[2], u[3]))
  # print(paste("v", v[1], v[2], v[3]))
  
  # if(v_t[2] < 0){
  #   phi = 2*pi - phi
  # }
  
  v_t = rotateAroundZ(c(v,0),phi)
  
  if(sum(c(u,0)-v_t) > 0.000001){
    print("Wrong")
    print(c(u,0)-v_t)
    return(FALSE)
  }
  
  X = rotateAroundZ(X, phi)
  
  return(X)
}

# compare2Data(X,Y)
# 
# for(i in 1:1000){
#   print(paste("Test:", i))
#   signs = sign(runif(2,-1.0,1.0))
#   v1 = c(runif(2,1.0,2.0)*signs[1])
#   v2 = c(runif(2,1.0,2.0)*signs[2])
#   
#   v1 = v1/norm_vec(v1)
#   v2 = v2/norm_vec(v2)
#   
#   if(FALSE == alignDataWithVectorsInZ0Plane(X, v1, v2)){
#     open3d()
#     arrow3d(p0 = c(0,0,0), p1 = v1, col = "red", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = v2, col = "blue", thickness = 0.1, width = 0.1, type = "rotation")
#     
#     e = diag(3)
#     arrow3d(p0 = c(0,0,0), p1 = e[,1], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = e[,2], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     arrow3d(p0 = c(0,0,0), p1 = e[,3], col = "black", thickness = 0.1, width = 0.1, type = "rotation")
#     return()
#   }
# }



# X_t = alignDataWithVectors(X_b,c(0,0,1), c(1,1,1))
# open3d()
# plot3d(X_t, col = "red")
# plot3d(X_b, add = TRUE, col = "blue")


# Y = rotateAroundZ(X,170)
# Xp = myPCA(X,"X")
# Yp = myPCA(Y,"Y")
# plot2Vol(Xp,Yp)
# 
# print(Xp$transformed_data)
# print(Yp$transformed_data)


projectInPlane0Centered <- function(X,normal){
  # n1*x + n2*y + n3*z + d = 0
  # where x,y,z is a point in the plane
  # that means we can take 0,0,0 since the datais centered
  # => d = 0
  #
  # d_p =  sum(n*p) + d // distance from point to plane
  # p_in_plane = p-d_p*n
  
  X_in_plane = matrix(rep(0,3*nrow(X)), nrow = nrow(X))
  for(i in 1:nrow(X)){
    X_in_plane[i,] = X[i,]-sum(X[i,] %*% normal)*normal
  }

  return(X_in_plane)
}

myPCA <- function(X,name){
  Xc = centerData(X)
  
  # we want to rotate our points in such a way that the 
  # normal-vector points into the z-direction that means (0,0,1)
  # after that we project our points in the plane with z = 0
  m = odregress(Xc[,1:2],Xc[,3])
  m$normal
  
  Xa = alignDataWithVectors(Xc,c(0,0,1), m$normal)
  
  m2 = odregress(Xa[,1],Xa[,2])  
  m2$normal
  
  pcs = diag(3)
  
  pcs[,1] = c(m2$normal[2],-m2$normal[1],0)
  pcs[,2] = c(m2$normal,0)
  pcs[,3] = c(0,0,1)
  
  ret = list("name" = name, "transformed_data" = Xa, "pcs" = pcs)
  return(ret)
}

calcDistanceMatrix <- function(X,Y){
  d_m = matrix(rep(0,nrow(X)*nrow(Y)),nrow = nrow(X))
  
  for(i in 1:nrow(X)){
    for(j in 1:nrow(Y)){
      d = euklid_dist(X[i,], Y[j,])
      d_m[i,j] = d^2
      # d_m[i,j] = euklid_dist(X[i,], Y[j,])
    }
  }
  
  return(d_m)
}

# m = matrix(c(0,1,2,0,0,1,0,0,0), byrow = TRUE, nrow = 3)
# diag(3)
# calcDistanceMatrix(m,m)

getMatching <- function(X,Y){
  d_matrix = calcDistanceMatrix(X,Y)
  
  myPrint(d_matrix)
  
  # print(paste(nrow(X),nrow(Y)))
  
  b = solve_LSAP(d_matrix)

  
  distances = d_matrix[cbind(seq_along(b), b)]
  # o <- order(b)
  
  inv_distances = distances[order(b)]
  
  s = sum(distances)
  
  # print(s)
  
  ret = list("matching" = b, "distances" = distances, "inv_distances" = inv_distances, "sum" = s)
  return(ret)
}

alignPcaObjs <- function(X,Y){
  
  Y_c = centerData(Y$transformed_data)

  Ya = alignDataWithVectors(Y_c,X$pcs[,3], Y$pcs[,3])
  Yb = alignDataWithVectorsInZ0Plane(Ya, X$pcs[1:2,1], Y$pcs[1:2,1])
  
  # Y = myPCA(Yb,"Yaligned")
  
  Y$transformed_data = Yb
  Y$pcs = X$pcs
  Y$name = "Yaligned"
  
  return(Y)
}

compare2Data <- function(X,Y,plot = TRUE){
  Xp = myPCA(X,"X")
  Yp = myPCA(Y,"Y")
  
  Yp_aligned = alignPcaObjs(Xp,Yp)
  
  myPrint("m -------------------------------")
  m = getMatching(Xp$transformed_data,Yp_aligned$transformed_data)
  
  myPrint("m2 -------------------------------")
  Y2 = rotateAroundArbitraryAxis(Yp_aligned$transformed_data, Yp_aligned$pcs[,1],pi)
  m2 = getMatching(Xp$transformed_data,Y2)
  
  myPrint("m3 -------------------------------")
  Y3 = rotateAroundArbitraryAxis(Yp_aligned$transformed_data, Yp_aligned$pcs[,3],pi)
  m3 = getMatching(Xp$transformed_data,Y3)
  
  myPrint("m4 -------------------------------")
  Y4 = rotateAroundArbitraryAxis(Yp_aligned$transformed_data, Yp_aligned$pcs[,2],pi)
  m4 = getMatching(Xp$transformed_data,Y4)
  
  if(m2$sum < m$sum){
    Yp_aligned$transformed_data = Y2
    m = m2
    # print("m2 smaller")
  }
  if(m3$sum < m$sum){
    Yp_aligned$transformed_data = Y3
    m = m3
    # print("m3 smaller")
  }
  if(m4$sum < m$sum){
    Yp_aligned$transformed_data = Y4
    m = m4
    # print("m4 smaller")
  }
  
  myPrint(paste(m$sum, m2$sum, m3$sum, m4$sum))
  myPrint(paste("distance is: ", m$sum))
  
  myPrint(m$matching)
  
  if(plot == TRUE){
    plot2Vol(Xp,Yp,Yp_aligned)
    
    for(i in 1:nrow(Xp$transformed_data)){
      # arrows(Xp$transformed_data[i,1],Xp$transformed_data[i,2], Yp$transformed_data[m$matching[i],1],Yp$transformed_data[m$matching[i],2], col = "blue", length = 0.05)
      
      # start = Xp$transformed_data[i,1],Xp$transformed_data[i,2]
      arrow3d(Xp$transformed_data[i,],Yp_aligned$transformed_data[m$matching[i],], col = "green",thickness = 0.1, width = 0.1, type = "rotation")
      # arrow3d(p0 = center, p1 = arrowHead_1, col = col, thickness = 0.1, width = 0.1, type = "rotation")
      
      # myPrint(Yp_aligned$transformed_data[m$matching[i],])
    }
  }

  
  ret = list("X_pca" = Xp, "Y_pca" = Yp, "matching" = m)
  return(ret)
}

subSamplePseudoPcaOld <- function(X,Y,sampleSize,repitions, numOfPointsToPlot, pseudoCount = 1.0){
  r_base = compare2Data(X,Y,FALSE)
  #-----------------------------------------
  
  
  X_scores = rep(pseudoCount,nrow(X))
  x_num = rep(0,nrow(X))
  
  Y_scores = rep(pseudoCount,nrow(Y))
  y_num = rep(0,nrow(Y))
  
  
  iter = 0
  for(i in 1:repitions){
    iter = iter +1
    print(paste(i/repitions*100, "%"))
    
    # x_indices = sample(c(1:nrow(X)),size = sampleSize)
    # X_sample = X[x_indices,]
    # y_indices = sample(c(1:nrow(Y)), size = sampleSize)
    # Y_sample = Y[y_indices,]
    
    # if(iter > 200){
    x_indices = sample(c(1:nrow(X)),prob = 1/X_scores,size = sampleSize)
    X_sample = X[x_indices,]
    y_indices = sample(c(1:nrow(Y)),prob = 1/Y_scores, size = sampleSize)
    Y_sample = Y[y_indices,]
    # }
    
    
    
    r = compare2Data(X_sample,Y_sample,FALSE)
    
    vs = r$matching$sum+0.1
    
    for(j in 1:length(x_indices)){
      k = x_indices[j]
      
      x_num[k] = x_num[k] + 1
      
      X_scores[k] = (X_scores[k]*(x_num[k]-1) + r$matching$distances[j]/vs)/x_num[k]
    }
    
    for(j in 1:length(y_indices)){
      k=y_indices[j]
      
      y_num[k] = y_num[k] + 1
      
      Y_scores[k] = (Y_scores[k]*(y_num[k]-1) + r$matching$inv_distances[j]/vs)/y_num[k]
    }
  }
  
  print(X_scores[order(X_scores)])
  print(x_num[order(X_scores)])
  
  print(Y_scores[order(Y_scores)])
  print(y_num[order(Y_scores)])
  
  
  # par(mfrow=c(2,2))
  # l = c(-20,20)
  # 
  # xl = c(min(X[,1]),max(Y[,1]))
  # yl = c(min(X[,2]),max(Y[,2]))
  
  # o = max(abs(c(xl,yl)))
  # sl = c(-o,o)
  
  # plot(X, xlim = sl, ylim = sl, col = "red")
  # plot(Y, xlim = sl, ylim = sl, col = "green")
  # 
  # 
  # plot(X[order(X_scores)[1:numOfPointsToPlot],1], X[order(X_scores)[1:numOfPointsToPlot],2], xlim = sl, ylim = sl, col = "red")
  # # points(Y[order(Y_scores)[1:4]], col = "green")
  # plot(Y[order(Y_scores)[1:numOfPointsToPlot],1], Y[order(Y_scores)[1:numOfPointsToPlot],2], xlim = sl, ylim = sl, col = "green")
  
  X_feature = order(X_scores)[1:numOfPointsToPlot]
  
  Y_feature = order(Y_scores)[1:numOfPointsToPlot]
  
  plot2VolFeaturesExtracted(r_base$X_pca,r_base$Y_pca,X_feature,Y_feature)
  
  return(list("X_feature" = X[X_feature,], "Y_feature" = Y[Y_feature,]))
}

generateRandomPoints <- function(n,lim){
  X=matrix(c(runif(n*3,-lim,lim)),nrow = n,byrow = TRUE)
  return(X)
}


plot2dPca <- function(X,col,li, add = FALSE){
  if(add == FALSE){
    plot(X$transformed_data[,1],X$transformed_data[,2], col = col, xlim = li, ylim = li, asp=1, pch = 19)
  } else {
    points(X$transformed_data[,1],X$transformed_data[,2], col = col, xlim = li, ylim = li, asp=1,  pch = 19)
  }


  arrows(0,0,X$pcs[1,1], X$pcs[2,1], col = col, length = 0.1)
  arrows(0,0,X$pcs[1,2], X$pcs[2,2], col = col, length = 0.1)
  
  # ?arrows
  
  sc = 1.3
  text("pc1",x = X$pcs[1,1]*sc, y = X$pcs[2,1]*sc, col = col)
  text("pc2",x = X$pcs[1,2]*sc, y = X$pcs[2,2]*sc, col = col)
}

plot2dOnly <- function(X,Y,Yaligned, YalignedOtherPlot = FALSE, deg = 180){
  l = max(abs(X$transformed_data), abs(Y$transformed_data), abs(Yaligned$transformed_data)) + 2
  li = c(-l,l)
  
  if(YalignedOtherPlot == TRUE){
    par(mfrow=c(2,2))
  } else {
    layout(matrix(c(1,1,1,0,2,2,2,
                    1,1,1,0,2,2,2,
                    1,1,1,0,2,2,2,
                    0,0,0,0,0,0,0,
                    0,3,3,3,3,3,0,
                    0,3,3,3,3,3,0,
                    0,3,3,3,3,3,0,
                    0,3,3,3,3,3,0,
                    0,3,3,3,3,3,0), 3+1+5, 7, byrow = TRUE))
  }

  
  # ?layout
  
  plot2dPca(X,"blue",li)
  
  plot2dPca(Y,"red",li)
  
  plot2dPca(X,"blue",li)
  
  plot2dPca(Yaligned,"red",li, add = TRUE)
  
  if(YalignedOtherPlot != FALSE){
    plot2dPca(X,"blue",li)
    
    YalignedOther = rotateAroundZ(Yaligned$transformed_data, deg)
    Yaligned$transformed_data = YalignedOther
    Yaligned$pcs = rotateAroundZ(Yaligned$pcs, deg)
    plot2dPca(Yaligned,"red",li, add = TRUE)
  }

  
  # points(Y$transformed_data[,1],Y$transformed_data[,2], col = "green")
  
}

plot2dOnlyWithFeatures <- function(X,Y,Yaligned, X_feature, Y_feature){
  plot2dOnly(X,Y,Yaligned)
}
# 
# #-------------------------------------------------------------------------------------------------
# # example 1
# #-------------------------------------------------------------------------------------------------
# feature = matrix(c(1,0,0, 0,0,0, 0,1,0, 1,1,0, 1,2,0, 1,3,0, 0,4,0, 2,4,0), byrow = TRUE, nrow = 8)
# X = matrix(c(1,9,0, 0,9,0, -7,8,0, -4, 8,0), byrow = TRUE, nrow = 4)
# X_b = rbind(X,feature)
# 
# Y = matrix(c(6,-7,0, -2,-5,0, -1,-6,0, -10,-8,0), byrow = TRUE, nrow = 4)
# Y_b = rbind(Y,feature)
# 
# q = compare2Data(X_b,Y_b)
# s = subSamplePseudoPcaOld(X_b,Y_b,sampleSize = 10, repitions = 100, numOfPointsToPlot = 8, pseudoCount = 1.0)
# 
# # plot2dOnly(q$X_pca,q$Y_pca, q$Y_pca)
# plot2dOnlyWithFeatures(q$X_pca,q$Y_pca, q$Y_pca, s$X_feature, s$Y_feature)
# 
# qFeature = compare2Data(s$X_feature,s$Y_feature)


# <----------------------------------
# #-------------------------------------------------------------------------------------------------
# # example 2
# #-------------------------------------------------------------------------------------------------
# feature = matrix(c(1,0,0, 0,0,0, 0,1,0, 1,1,0, 1,2,0, 1,3,0, 0,4,0, 2,4,0), byrow = TRUE, nrow = 8)
# X = matrix(c(1,9,-2, 0,9,5, -7,8,-2, -4, 8,3), byrow = TRUE, nrow = 4)
# X_b = rbind(X,feature)
# 
# Y = matrix(c(6,-7,-2, -2,-5,7, -1,-6,-6, -10,-8,-5), byrow = TRUE, nrow = 4)
# Y_b = rbind(Y,feature)
# 
# q = compare2Data(X_b,Y_b)


# #-------------------------------------------------------------------------------------------------
# # example 3
# #-------------------------------------------------------------------------------------------------
# feature = matrix(c(1,0,0, 0,0,0, 0,1,0, 1,1,0, 1,2,0, 1,3,0, 0,4,0, 2,4,0), byrow = TRUE, nrow = 8)
# feature = rotateAroundX(feature,45)
# X = matrix(c(1,9,0, 0,9,0, -7,8,0, -4, 8,0), byrow = TRUE, nrow = 4)
# X_b = rbind(X,feature)
# 
# Y = matrix(c(6,-7,0, -2,-5,0, -1,-6,0, -10,-8,0), byrow = TRUE, nrow = 4)
# Y_b = rbind(Y,feature)
# 
# compare2Data(X_b,Y_b)
# #-------------------------------------------------------------------------------------------------
# # example 4
# #-------------------------------------------------------------------------------------------------
# feature = matrix(c(1,0,0, 0,0,0, 0,1,0, 1,1,0, 1,2,0, 1,3,0, 0,4,0, 2,4,0), byrow = TRUE, nrow = 8)
# feature = rotateAroundX(feature,45)
# X = matrix(c(1,9,0, 0,9,0, -7,8,0, -4, 8,0), byrow = TRUE, nrow = 4)
# X_b = rbind(X,feature)
# 
# # Y = matrix(c(6,-7,0, -2,-5,0, -1,-6,0, -10,-8,0), byrow = TRUE, nrow = 4)
# Y_b = rotateAroundX(X_b,-50)
# 
# compare2Data(X_b,Y_b)
# #-------------------------------------------------------------------------------------------------
# # example 5
# #-------------------------------------------------------------------------------------------------
# X = matrix(c(-3,0,0, 0,0,0, 3,0,0, 3,4,0), byrow = TRUE, nrow = 4)
# 
# testRotData <- function(X,degree){
#   Y = rotateAroundX(X,degree)
#   
#   ret = compare2Data(X,Y)
# }
# 
# testRotData(X,14)
# 
# #-------------------------------------------------------------------------------------------------
# # example 6
# #-------------------------------------------------------------------------------------------------
# GLOBAL_VERBOSITY = FALSE
# X = generateRandomPoints(100,5)
# testRotData(X,148)
# #-------------------------------------------------------------------------------------------------
# # example 7
# #-------------------------------------------------------------------------------------------------
# GLOBAL_VERBOSITY = FALSE
# feature = matrix(c(1,0,0, 0,0,0, 0,1,0, 1,1,0, 1,2,0, 1,3,0, 0,4,0, 2,4,0), byrow = TRUE, nrow = 8)
# feature = rotateAroundX(feature,45)
# X = matrix(c(1,9,0, 0,9,0, -7,8,0, -4, 8,0), byrow = TRUE, nrow = 4)
# X_b = rbind(X,feature)
# 
# Y = matrix(c(6,-7,0, -2,-5,0, -1,-6,0, -10,-8,0), byrow = TRUE, nrow = 4)
# Y_b = rbind(Y,feature)
# 
# subSamplePseudoPca(X_b,Y_b,8,100, 8, pseudoCount = 1.0)
# 
# #-------------------------------------------------------------------------------------------------
# # example 8
# #-------------------------------------------------------------------------------------------------
# GLOBAL_VERBOSITY = FALSE
# feature = matrix(c(1,0,0, 0,0,0, 0,1,0, 1,1,0, 1,2,0, 1,3,0, 0,4,0, 2,4,0), byrow = TRUE, nrow = 8)
# feature = rotateAroundX(feature,45)
# X = generateRandomPoints(50,8)
# X_b = rbind(X,feature)
# 
# Y = generateRandomPoints(50,8)
# Y_b = rbind(Y,feature)
# 
# subSamplePseudoPca(X_b,Y_b,20,100, 20, pseudoCount = 1.0)

#-------------------------------------------------------------------------------------------------
randomFeatureTest <- function(f,extra){
  lim = c(-9,9)
  plotLim = c(-20,-20)
  
  # feature = matrix(c(1,0, 0,0, 0,1, 1,1, 1,2, 1,3, 0,4, 2,4), byrow = TRUE, nrow = 8)
  feature = matrix(runif(3*f, lim[1], lim[2]), nrow = f, byrow = TRUE)
  
  feature = feature*0.5
  feature[,1] = feature[,1] - 20
  feature[,2] = feature[,2] - 20
  feature[,3] = feature[,3] - 20
  
  
  
  X = matrix(runif(3*extra, lim[1], lim[2]), nrow = extra, byrow = TRUE)
  Y = matrix(runif(3*extra, lim[1], lim[2]), nrow = extra, byrow = TRUE)
  
  X_b = rbind(X,feature)
  Y_b = rbind(Y,feature)
  
  plot(X_b)
  plot(Y_b)
  
  subSamplePseudoPca( Y_b,X_b, 2*f,100,f, pseudoCount = 1)
}

# randomFeatureTest(8,8)









