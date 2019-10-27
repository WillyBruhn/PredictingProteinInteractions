#---------------------------------------------------------------
# examples to compare the geo method with the emd-method
#
#---------------------------------------------------------------

s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s1)

wsPath = "/home/willy/PredictingProteinInteractions/setUp/SourceLoader.R"
# wsPath = as.character(paste(funr::get_script_path(), "/../../setUp/SourceLoader.R", sep = ""))

source(wsPath)
sourceFiles(c("helperFunctions"))
sourceFiles(c("UltraQuickRepeatedSubSampling"))
sourceFiles(c("TriangulateIsoSurface"))
sourceFiles(c("kerasFunctions"))



calcGeo <- function(df_x_y){
  # df_x_y = as.numeric(df_x_y)
  return(c((sum(df_x_y[,1])/length(df_x_y[,1])),sum(df_x_y[,2])/length(df_x_y[,2])))
}

plotExample <- function(df, size = 1){
  xli = as.numeric(c(as.numeric(min(df[,2]))-1,1+as.numeric(max(df[,2]))))
  yli = as.numeric(c(as.numeric(min(df[,3]))-1,1+as.numeric(max(df[,3]))))
  
  plot(df[which(df[,1] == "X"),2:3], pch = 19, xlim = xli, ylim = yli, col = "red", cex = size)
  points(x = df[which(df[,1] == "Y"),2], y = df[which(df[,1] == "Y"),3], col = "blue", pch = 19,  cex = size)
  
  dfX = df[which(df[,1] == "X"),2:3]
  dfX$x = as.numeric(as.character(dfX$x))
  dfX$y = as.numeric(as.character(dfX$y))
  
  geoX = calcGeo(dfX)
  
  dfY = df[which(df[,1] == "Y"),2:3]
  dfY$x = as.numeric(as.character(dfY$x))
  dfY$y = as.numeric(as.character(dfY$y))
  geoY = calcGeo(dfY)
  
  points(x = geoX[1], y = geoX[2], col = "red", pch = 18)
  points(x = geoY[1], y = geoY[2], col = "blue", pch = 18)
  
  return(list("geoX" = geoX, "geoY" = geoY))
}




geoVsEmdXandYdf <- function(df){
  par(mfrow=c(1,2)) 
  
  l = plotExample(df)
  
  dfX = df[which(df$name == "X"),2:3]
  d_all = as.matrix(dist(df[,2:3],method = "manhattan"))
  
  dX = as.matrix(d_all[which(df$name == "X"),which(df$name == "X")])
  dY = d_all[which(df$name == "Y"),which(df$name == "Y")]
  
  d_X_Y = d_all[1:nrow(dX),(nrow(dX)+1):nrow(d_all)]
  
  breaks = seq(from = 0,to = max(dX,d_X_Y)+0.05, by = 0.05)
  
  P = hist(dX[upper.tri(dX)] , breaks = breaks, plot = F)
  Q = hist(d_X_Y , breaks = breaks, plot = F)
  
  # ?hist
  P_normalized = t(as.matrix(P$counts))/sum(P$counts)
  print(P_normalized)
  
  Q_normalized = t(as.matrix(Q$counts))/sum(Q$counts)
  print(Q_normalized)
  
  # ?emd
        
  emdDist = emd2d(P_normalized,Q_normalized)
  
  geoDist = dist(rbind(l$geoX, l$geoY), method = "manhattan")
  
  geoDist2 = geoDist*(max(dX,d_X_Y)*10)
  
  
  print(paste("emd: ",emdDist, ", geo: ", geoDist, ", geo2: ", geoDist2))
  print(paste("Xdistances: "))
  print(dX[upper.tri(dX)])
  print(paste("X_Ydistances: "))
  print(as.vector(d_X_Y))
  
  plot(P, col ="red", main = paste("emd: ",emdDist, "geo: ", geoDist), freq = FALSE)
  plot(Q, col ="blue", add = TRUE, freq = FALSE)
  
  
  return(list("P" = t(as.matrix(P$counts)), "Q"  =t(as.matrix(Q$counts))))
}



#-------------------------------------------------------
# example1

df = data.frame(matrix(0, ncol = 3),stringsAsFactors = FALSE)
colnames(df) = c("name", "x", "y")

df[1,] = c("X", 1, 2)
df[2,] = c("X", 3, 2)
df[3,] = c("Y", 2, 1)
df[4,] = c("Y", 2, 3)

geoVsEmdXandYdf(df)


A <- matrix(1:6 / sum(1:6), 1)
B <- matrix(c(0, 0, 0, 0, 0, 1), 1)

emd2d(A, B)
# if we bring the rows closer, the distance will be reduced
# since mass from the first row has to move down
emd2d(A, B,, 0.1)

#-------------------------------------------------------

#-------------------------------------------------------
# example2

df = data.frame(matrix(0, ncol = 3),stringsAsFactors = FALSE)
colnames(df) = c("name", "x", "y")

df[1,] = c("X", 1, 2)
df[2,] = c("X", 3, 2)
df[3,] = c("Y", 2, 2)
df[4,] = c("Y", 2, 3)

AB = geoVsEmdXandYdf(df)

AB$P
AB$Q

emd2d(AB$P, AB$Q)

A <- matrix(c(0, 0.25, 0, 0, 0, 0.25), 1)
B <- matrix(c(0, 0, 0, 0, 0, 1), 1)

emd2d(A, B)
#-------------------------------------------------------

#-------------------------------------------------------
# example3

df = data.frame(matrix(0, ncol = 3),stringsAsFactors = FALSE)
colnames(df) = c("name", "x", "y")

df[1,] = c("X", 1, 2)
df[2,] = c("X", 3, 3)
df[3,] = c("Y", 2, 2)
df[4,] = c("Y", 2, 3)

geoVsEmdXandYdf(df)


par(mfrow=c(1,1)) 
plotExample(df, size = 3)
points(x = 2, y = 2.5, col = "green", cex = 3, pch = 19)

#-------------------------------------------------------

df = data.frame(matrix(0, ncol = 3),stringsAsFactors = FALSE)
colnames(df) = c("name", "x", "y")

df[1,] = c("X", 1, 2)
df[2,] = c("X", 3, 3)
df[3,] = c("Y", 1, 2+10)
df[4,] = c("Y", 3, 3+10)

geoVsEmdXandYdf(df)


#-------------------------------------------------------
#-------------------------------------------------------
# example4

df = data.frame(matrix(0, ncol = 3),stringsAsFactors = FALSE)
colnames(df) = c("name", "x", "y")

df[1,] = c("X", 1, 2)
df[2,] = c("X", 3, 3)
df[3,] = c("Y", 2, 2)
df[4,] = c("Y", 2, 3)


geoVsEmdXandYdf(df)


#-------------------------------------------------------
#-------------------------------------------------------
# example5

n = 100



dfX = data.frame(matrix(rnorm(n*3), ncol = 3),stringsAsFactors = FALSE)
colnames(dfX) = c("name", "x", "y")
dfX[,1] = rep("X",n) 

dfY = data.frame(matrix(rnorm(n*3,mean = 100,sd = 0.1), ncol = 3),stringsAsFactors = FALSE)
colnames(dfY) = c("name", "x", "y")
dfY[,1] = rep("Y",n) 

df = rbind(dfX,dfY)

geoVsEmdXandYdf(df)

AB = geoVsEmdXandYdf(df)
emd2d(AB$P, AB$Q)

dfXGeo = calcGeo(df_x_y = dfX[,-1])
dfYGeo = calcGeo(df_x_y = dfY[,-1])

manhattan_dist(dfXGeo,dfYGeo)

#-------------------------------------------------------------
projection = read.csv(file = "/home/willy/PredictingProteinInteractions/data/106Test/UltraQuickRepSub/_pos_quickEmd_n_500_m_100_q_1_projection.csv", header = TRUE)
quantiles = projection[,-2]


library(rgl)

labels = read.table("/home/willy/PredictingProteinInteractions/data/106Test/labels.txt", header = TRUE)

geos = read.csv(file = "/home/willy/PredictingProteinInteractions/data/106Test/UltraQuickRepSub/_pos_quickEmd_n_500_m_100_q_1_geometricCenters.csv", header = TRUE)
geos = read.csv(file = "/home/willy/PredictingProteinInteractions/data/106Test/UltraQuickRepSub/_pos_quickEmd_n_100_m_500_q_1_geometricCenters.csv", header = TRUE)

functionals = as.character(labels$name[which(labels$label == "functional")])

functionalInds = which(quantiles[,1] %in% functionals)
notfunctionalInds = c(1:nrow(quantiles))[-functionalInds]

points3d(quantiles[functionalInds,-1], col = "red")
points3d(quantiles[notfunctionalInds,-1], col = "blue")



geos[,-c(1,2,3)]
functionalInds2 = which(geos[,2] %in% functionals)
notfunctionalInds2 = c(1:nrow(geos))[-functionalInds2]
points3d(geos[functionalInds2,-c(1,2,3)],col = "red", size = 15)
points3d(geos[notfunctionalInds2,-c(1,2,3)],col = "blue", size = 15)



quantiles2 = read.csv(file = "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.9_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv", header = TRUE)

functionalInds = which(quantiles2[,1] %in% functionals)
notfunctionalInds = c(1:nrow(quantiles2))[-functionalInds]

points3d(quantiles2[functionalInds,-1], col = "red")
points3d(quantiles2[notfunctionalInds,-1], col = "blue")

#--------------------------------------------------------------------

read










