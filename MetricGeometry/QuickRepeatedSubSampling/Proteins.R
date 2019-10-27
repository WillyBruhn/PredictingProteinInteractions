#!/usr/bin/Rscript

#----------------------------------------------------------------------------------
# Willy Bruhn, 25.10.19
#
# Given the positive and negative potential of the proteins and the active center
# calculate multiple features. These features are the quantiles of the approximation
# of the DE (Distribution of Eccentricities). The features are then used to train
# a neural net.
#
# There are two options:
# 1 training a neural net and evaluating with k-fold-Cross-Validation. Additionally this 
#   model can be exported.
# 2 using a pre-trained model to make predictions on new data.
#
# If you only want to generate the features run with
# onlyGenerateModels = TRUE
#
# Requirements:
# - In the protein-folder the .obj-file has to be present which is delivered by MutComp.
# If the measure should be updated with the information about the active center, then
# the file "active_center.pts" has to be present which is delivered by selectCenter.R.
#
#----------------------------------------------------------------------------------
wsPath = "/home/willy/PredictingProteinInteractions/setUp/SourceLoader.R"
# wsPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/setUp/SourceLoader.R"

# mode = "onlyExperiments2"
# # mode = "onlyGenerateModels

wsPath = as.character(paste(funr::get_script_path(), "/../../setUp/SourceLoader.R", sep = ""))

source(wsPath)
sourceFiles(c("helperFunctions"))
sourceFiles(c("UltraQuickRepeatedSubSampling"))
sourceFiles(c("TriangulateIsoSurface"))
sourceFiles(c("kerasFunctions"))

path106Experiment = getPath("106Experiment")
path120Experiment = getPath("120Experiment")
pathpdbDownloaderExperiment = getPath("pdbDownloaderExperiment")

# pathToExperiment = pathpdbDownloaderExperiment
pathToExperiment = path106Experiment

LABELS = getPath("106ExperimentLabels")
NUMCLASSES = 2

path2Manifold = getPath("Manifold")

# library(keras)
library(readobj)
library(FNN)
library(raster)

library(beepr)

if(is.installed("rgl")) library(rgl)
library(rdist)
library(caret) # F1-score

# hacky way to check if we are on WS
WS_flag = FALSE
if(strsplit(wsPath, "/")[[1]][3] == "sysgen"){
  WS_flag = TRUE
}

SAVE_EXPERIMENTS = TRUE
if(WS_flag == TRUE){
  library(reticulate)
  use_python("/home/sysgen/.pyenv/versions/3.6.3/bin/python3.6", required = TRUE)
  use_virtualenv("/home/sysgen/.pyenv/versions/3.6.3/",required = TRUE)
  SAVE_EXPERIMENTS = TRUE
}

library(keras)
library(doParallel)

# registerDoParallel(detectCores()/2)

numCores = 1
# args = commandArgs(trailingOnly=TRUE)
# if(length(args) > 0){
#   numCores = as.numeric(args[1])
# }

# install.packages("e1071")
library(e1071)

print("done loading ...")


#----------------------------------------------------------------------------------
# interpreting input parameters
library(getopt)
options(warn=-1)
#----------------------------------------------------------------------------------
# Input
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'mode', 'm', 2, "character",
  'onlyGenerateModels'  , 'o', 2, "integer",
  'numClasses', 'n', 2, "integer",
  'pathToExperiment', 'p',2, "character",
  'labels'  , 'l', 2, "character",
  'folds'  , 'f', 2, "integer",
  'cores', 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# mode = "onlyExperiments2"
# mode = "onlyGenerateModels
if ( is.null(opt$onlyGenerateModels    ) ) { opt$onlyGenerateModels    = 0     }
if ( is.null(opt$mode    ) ) { opt$mode    = "training"     }
if(is.null(opt$labels)) {opt$labels    = getPath("106ExperimentLabels")}
if(is.null(opt$numClasses)) {opt$numClasses    = 2}
if(is.null(opt$folds)) {opt$folds    = 10}
if(is.null(opt$pathToExperiment)){opt$pathToExperiment = getPath("106Experiment")}
if ( is.null(opt$cores) ) { 
  opt$cores = Sys.getenv('LSB_MAX_NUM_PROCESSORS')
  if ( is.null(opt$cores) || opt$cores == "") opt$cores = 1
}
if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }

mode = opt$mode
LABELS = opt$labels
numCores = opt$cores
NUMCLASSES = opt$numClasses
pathToExperiment = opt$pathToExperiment
numFolds = opt$folds
mode = opt$mode

print(opt)

print(paste("Using ", numCores, " cores", sep =""))
registerDoParallel(numCores)

quit()
#----------------------------------------------------------------------------------

# install.packages("keras")

# library(keras)

# use_python("/usr/bin/python3.6")
# install_keras()
# 
# library(keras)
# use_python("/usr/bin/python3.7")
# k=backend()
# sess = k$get_session()
# sess$list_devices()
# 
# 
# library(tensorflow)
# install_tensorflow(version = "gpu")

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

localEccentricitiesAndGetQuantiles <- function(model,
                                                n = 0.3,
                                                q = 1)
{
  F_app = getFApproximationsSumMethod_local(mod = model, percent = n, q = q)
  
  quantiles = data.frame(matrix(0,ncol = (q+2)*6+1, nrow = nrow(model$centers)))
  colnames(quantiles) = c("name",
                          as.vector(paste(c("q_pos"),c(1:(q+2)), sep = "")), 
                          as.vector(paste(c("q_neg"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_neg"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_euclid"),c(1:(q+2)), sep = "")), 
                          as.vector(paste(c("q_neg_euclid"),c(1:(q+2)), sep = "")),
                          as.vector(paste(c("q_pos_neg_euclid"),c(1:(q+2)), sep = "")))
  
  quantiles[,1] = rep(model$name, nrow(quantiles))
  for(i in 1:nrow(quantiles)){
    quantiles[i,(1:(q+2))+1] = F_app[[i]]$F_pos_approx
    quantiles[i,(1:(q+2))+1+(q+2)] = F_app[[i]]$F_neg_approx
    quantiles[i,(1:(q+2))+1+(q+2)*2] = F_app[[i]]$F_pos_neg_approx
    
    quantiles[i,(1:(q+2))+1+(q+2)*3] = F_app[[i]]$F_pos_approx_euclid
    quantiles[i,(1:(q+2))+1+(q+2)*4] = F_app[[i]]$F_neg_approx_euclid
    quantiles[i,(1:(q+2))+1+(q+2)*5] = F_app[[i]]$F_pos_neg_approx_euclid
  }
  
  return(quantiles)
}

# localEccentricitiesAndGetQuantiles(models1[[1]])

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


getFApproximationsSumMethod_local <- function(mod, percent = 0.1,q=1){
  #-------------------------------------------------------------------------------
  # for each point in the model take the closest n neighbors and calculate the DE
  #-------------------------------------------------------------------------------
  
  pos_indices = which(mod$posNegVector == TRUE)
  neg_indices = which(mod$posNegVector == FALSE)
  
  n = percent*nrow(mod$d_surface)
  NN_indices = get.knn(mod$d_surface,k = n-1)$nn.index
  
  F_approximations = list()
  # for(i in 1:nrow(mod$centers)){
  F_approximations = lapply(c(1:nrow(mod$centers)), FUN = function(i){
    inds = c(i,NN_indices[i,])
    
    pos_indices_samp = inds[which(mod$posNegVector[inds] == TRUE)]
    neg_indices_samp = inds[which(mod$posNegVector[inds] == FALSE)]
    
    # if(length(pos_indices_samp) < 3 || length(neg_indices_samp) < 3) print(paste("too small",length(pos_indices_samp), length(neg_indices_samp)))
    # 
    #----------------------------------------------------------------------------
    # surface

    F_pos_approx = rep(0,q+2)
    if(length(pos_indices_samp) > 1) {
      d_pos = mod$d_surface[pos_indices_samp,pos_indices_samp]
      F_pos = DistributionOfEccentricities(d_pos, mod$measure[pos_indices_samp])
      F_pos_approx = approximateCDF(F_pos, q)
    }

    F_neg_approx = rep(0,q+2)
    if(length(neg_indices_samp) > 1) {
      d_neg = mod$d_surface[neg_indices_samp,neg_indices_samp]
      F_neg = DistributionOfEccentricities(d_neg, mod$measure[neg_indices_samp])
      F_neg_approx = approximateCDF(F_neg, q)
    }

    merged = c(pos_indices_samp,neg_indices_samp)
    F_pos_neg = DistributionOfEccentricities(mod$d_surface[merged,merged], mod$measure[merged])
    F_pos_neg_approx = approximateCDF(F_pos_neg, q)
    
    #----------------------------------------------------------------------------
    # euclid
    
    F_pos_approx_euclid = rep(0,q+2)
    if(length(pos_indices_samp) > 1) {
      d_pos_euclid = as.matrix(dist(mod$centers[pos_indices_samp,]))
      F_pos_euclid  = DistributionOfEccentricities(d_pos_euclid, mod$measure[pos_indices_samp])
      F_pos_approx_euclid = approximateCDF(F_pos_euclid, q)
    }
    
    F_neg_approx_euclid = rep(0,q+2)
    if(length(neg_indices_samp) > 1) {
      d_neg_euclid = as.matrix(dist(mod$centers[neg_indices_samp,]))
      F_neg_euclid = DistributionOfEccentricities(d_neg_euclid, mod$measure[neg_indices_samp])
      F_neg_approx_euclid = approximateCDF(F_neg_euclid, q)
    }
    
    merged = c(pos_indices_samp,neg_indices_samp)
    F_pos_neg_euclid = DistributionOfEccentricities(as.matrix(dist(mod$centers[merged,])), mod$measure[merged])
    F_pos_neg_approx_euclid = approximateCDF(F_pos_neg_euclid, q)
    
    F_approximations[[i]] = list("F_pos_approx" = F_pos_approx,
                                 "F_neg_approx" = F_neg_approx,
                                 "F_pos_neg_approx" = F_pos_neg_approx,
                                 "F_pos_approx_euclid" = F_pos_approx_euclid,
                                 "F_neg_approx_euclid" = F_neg_approx_euclid,
                                 "F_pos_neg_approx_euclid" = F_pos_neg_approx_euclid)
  })
  # }
    
  return(F_approximations)
}

getFApproximationsSumMethod_local_samllTest <- function(model, percent = 0.3){
  f = getFApproximationsSumMethod_local(model,percent = percent)
  
  m = matrix(0, ncol = 3, nrow = length(f))
  for(i in 1:length(f)){
    m[i,] = f[[i]]$F_pos_neg_approx
  }
  
  return(m)
}


# Rprof(tmp <- tempfile())
# f = getFApproximationsSumMethod_local(models1[[1]],percent = 0.2)
# Rprof()
# summaryRprof(tmp)


# getFApproximationsSumMethod_local_samllTest(models1[[i]], percent = 0.05)

# # small and funny image test
# models1 = models_small
# ms = list()
# 
# for(i in 1:length(models1)){
#   print(i)
#   ms[[i]] = getFApproximationsSumMethod_local_samllTest(models1[[i]], percent = 0.05)
# }
# 
# 
# 
# transformQuants
# 
# 
# for(i in 1:length(ms)){
#   print(i)
#   col = "blue"
#   traf = transformQuants(cbind(rep(0,nrow(ms[[i]])),ms[[i]]),q=1)[,-1]
#   if(models1[[i]]$name %in% functionals) points3d(traf, col = "red")
# 
# }
# 
# i = 2
# models1[[i]]$name %in% functionals
# traf = transformQuants(cbind(rep(0,nrow(ms[[i]])),ms[[i]]),q=1)[,-1]
# points3d(traf, col = "blue", size = 10)
# 
# 
# 
# which(models1$name %in% functionals)
# length(ms)
# for(i in 1:length(ms)){
#   print(i)
#   col = "blue"
# 
#   if(models1[[i]]$name %in% functionals) col = "red"
#   traf = transformQuants(cbind(rep(0,nrow(ms[[i]])),ms[[i]]),q=1)[,-1]
# 
#   while (rgl.cur() > 0) { rgl.close() }
# 
#   points3d(matrix(c(0,0,0,
#                     0,0,1,
#                     0,1,0,
#                     0,1,1,
#                     1,0,0,
#                     1,0,1,
#                     1,1,0,
#                     1,1,1
#                     )*5, ncol = 3, byrow = TRUE), col = "green", size = 10)
# 
#   points3d(traf, col = col)
#   rgl.snapshot(paste("/home/willy/PredictingProteinInteractions/Results/Images/protMeasureProjection/", models1[[i]]$name, ".png",sep = ""))
# }
# 
# 
# 
# mTrx = getFApproximationsSumMethod_local_samllTest(models1[[7]], percent = 0.2)
# m027 = getFApproximationsSumMethod_local_samllTest(models1[[34]], percent = 0.2)
# m001 = getFApproximationsSumMethod_local_samllTest(models1[[8]], percent = 0.2)
# m013 = getFApproximationsSumMethod_local_samllTest(models1[[20]], percent = 0.2)
# m049 = getFApproximationsSumMethod_local_samllTest(models1[[56]], percent = 0.2)
# m046 = getFApproximationsSumMethod_local_samllTest(models1[[53]], percent = 0.2)
# 
# models1[[56]]$name
# 
# points3d(mTrx, col = "green")
# points3d(m027, col = "blue")
# points3d(m001, col = "red")
# points3d(m013, col = "black")
# points3d(m049, col = "brown")
# points3d(m046, col = "yellow")
# 
# points3d(models1[[1]]$centers)
# NN_indices = get.knn(models1[[1]]$d_surface,k = 500)$nn.index
# 
# points3d(models1[[1]]$centers[NN_indices[100,],], col = "green", size = 10)



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
                               alpha = 1,
                               betha = 1,
                               onlyTheseIndices = NULL,
                               recalculate = FALSE){
  
  dirs = list.dirs(path = path, recursive = FALSE, full.names = FALSE)
  if(!is.null(onlyTheseIndices)) dirs = dirs[onlyTheseIndices]
  
  # proteinModels = list()
  proteinModels = foreach(i = 1:length(dirs)) %dopar% {
    # proteinModels[[i]] = 
    getProteinModelStichedSurface(path = path, protName = dirs[i], n_s_euclidean = n_s_euclidean,
                                                       n_s_dijkstra = n_s_dijkstra, plot = FALSE,
                                                       recalculate = recalculate,
                                                       stitchNum =stitchNum,
                                                       measureNearestNeighbors = measureNearestNeighbors,
                                                       alpha = alpha,
                                                       betha = betha)
  }
  
  return(proteinModels)
}


plotProteinModel <- function(fNameOrigingal = NULL,lis, openNew = TRUE, sz = 1, totalMassSize = 30){
  # points3d(points)
  
  # models1[[1]]$atomCoords = atomCoords$atomCoords
  # models1[[1]]$activeCenterInds = atomCoords$activeCenterInds
  
  my_print("plotting surfaces ...")
  if(openNew) {
    rgl.open()
    par3d(windowRect = c(0,0,1000,1000))
  }
  points3d(lis$centers[which(lis$posNegVector == FALSE),], col = "red", size = sz)
  points3d(lis$centers[which(lis$posNegVector == TRUE),], col = "blue", size = sz)
  
  sc = totalMassSize/max(lis$measure)
  
  my_print("plotting measure ...")
  for(i in 1:nrow(lis$centers)){
    # print(lis$measure[i]*sc)
    points3d(x = lis$centers[i,1], y = lis$centers[i,2], z = lis$centers[i,3], col = "green", size = lis$measure[i]*sc)
  }
  
  my_print("plotting atomCoords ...")
  # for(i in 1:nrow(lis$atomCoords)){
    points3d(x = lis$atomCoords[,1], y = lis$atomCoords[,2], z = lis$atomCoords[,3], col = "black", size = sz)
  # }
  
  my_print("plotting active center ...")
  # for(i in lis$activeCenterInds){
  points3d(x = lis$atomCoords[lis$activeCenterInds,1], y = lis$atomCoords[lis$activeCenterInds,2], z = lis$atomCoords[lis$activeCenterInds,3], col = "yellow", size = sz*5)
  # }
  
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
                                          alpha = 1,
                                          betha = 1,
                                          plot = FALSE,
                                          recalculate = FALSE){
  # first stich the model together
  fNameOrigingal = paste(path, "/", protName, "/", protName, ".obj", sep ="")
  fNameStitched = paste(path, "/", protName, "/", protName, "_",stitchNum,"_stitched.obj", sep = "")
  
  fNameModelDownsampled = paste(path, "/", protName, "/", protName, "_sitch_",stitchNum, "_nE_",n_s_euclidean, "_nD_", n_s_dijkstra,
                                "_model_downsampled.rData", sep ="")
  
  fNameModelDownsampled_additional_params = paste(path, "/", protName, "/", protName, 
                                                  "_sitch_",stitchNum, 
                                                  "_nE_",n_s_euclidean, 
                                                  "_nD_", n_s_dijkstra,
                                                  "_NNBoarder_", measureNearestNeighbors,
                                                  "_alpha_", alpha,
                                                  "_betha_",betha,
                                                  "_model_downsampled.rData", sep ="")
  
  if(!file.exists(fNameModelDownsampled_additional_params) || recalculate == TRUE){
    if(!file.exists(fNameModelDownsampled) || recalculate == TRUE){
      if(!file.exists(fNameStitched) || recalculate == TRUE){
        # make watertight
        # obj = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/bathtub/test/bathtub_0110.obj"
        # path2Manifold = "/home/willy/Manifold/build/"
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
      
      
      sampled_indices = myFarthestPointSampling(points, k = n_s_euclidean)
      
      # points3d(points[sampled_indices,], col = "green", size = 20)
      
      d_surface = myShortestDistances(graph,sampled_indices)
      
      fps_surface <- farthest_point_sampling(d_surface)
      sampled_indices2 = fps_surface[1:n_s_dijkstra]
      
   
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
      
      boarderCounts = getBoarderNNCount(lis, neighbors = measureNearestNeighbors)
      
      atomCoords = getACtiveCenterAndAtomCoords(path = path, protName)
      distancesToActiveCenter = getDistancesToActiceCenter(atomCoords = atomCoords$atomCoords,
                                                                        activeCenterInds = atomCoords$activeCenterInds,
                                                                        centers = centers)
      
      measure = getProteinMeasure(distancesToActiveCenter = distancesToActiveCenter,boarderCounts = boarderCounts,alpha = alpha, betha = betha)
      
      # calculate for each starting point the fps-sequence
      fps_sequences = matrix(0,nrow(d_surface),ncol = 100)
      for(i in 1:nrow(d_surface)){
        fps_sequences[i,] = farthest_point_sampling(mat = d_surface,initial_point_index = i,k = 100)
      }
      
      lis = list("centers" = centers,
                 "posNegVector" = posNegVector,
                 "d_surface" = d_surface,
                 "name" = protName,
                 "boarderCounts" = boarderCounts,
                 "measureNearestNeighbors" = measureNearestNeighbors,
                 "atomCoords" = atomCoords$atomCoords,
                 "activeCenterInds" = atomCoords$activeCenterInds,
                 "distancesToActiveCenter" = distancesToActiveCenter,
                 "alpha" = alpha,
                 "betha" = betha,
                 "measure" = measure,
                 "fps_sequences" = fps_sequences)
  
      saveRDS(lis,file = fNameModelDownsampled)
      saveRDS(lis,file = fNameModelDownsampled_additional_params)
    }
    
    
    lis = readRDS(fNameModelDownsampled)
    my_print(fNameModelDownsampled, 1)
    
    if(   is.null(lis$measureNearestNeighbors) 
       || lis$measureNearestNeighbors != measureNearestNeighbors
       || is.null(lis$alpha)
       || lis$alpha != alpha
       || is.null(lis$betha)
       || lis$betha != betha 
       || is.null(lis$atomCoords)
       || is.null(lis$activeCenterInds)
       || is.null(lis$distancesToActiveCenter)){
      my_print(paste("recalculating measure for ", fNameModelDownsampled_additional_params, sep = ""), 1)
      
      if(is.null(lis$boarderCounts) || is.null(lis$measureNearestNeighbors) || lis$measureNearestNeighbors != measureNearestNeighbors){
        lis$boarderCounts = getBoarderNNCount(lis, neighbors = measureNearestNeighbors)
      }

      
      if(is.null(lis$atomCoords) || is.null(lis$activeCenterInds)){
        atomCoords = getACtiveCenterAndAtomCoords(path = path, protName)
        lis$atomCoords = atomCoords$atomCoords
        lis$activeCenterInds = atomCoords$activeCenterInds
      }
      
      if(is.null(lis$distancesToActiveCenter)){
        lis$distancesToActiveCenter = getDistancesToActiceCenter(atomCoords = atomCoords$atomCoords,
                                                             activeCenterInds = atomCoords$activeCenterInds,
                                                             centers = lis$centers)
      }
      
      
      measure = getProteinMeasure(distancesToActiveCenter = lis$distancesToActiveCenter, boarderCounts = lis$boarderCounts,alpha = alpha, betha = betha)
      
      lis = list("centers" = lis$centers,
                 "posNegVector" = lis$posNegVector,
                 "d_surface" = lis$d_surface,
                 "name" = lis$name,
                 "boarderCounts" = lis$boarderCounts,
                 "measureNearestNeighbors" = measureNearestNeighbors,
                 "atomCoords" = lis$atomCoords,
                 "activeCenterInds" = lis$activeCenterInds,
                 "distancesToActiveCenter" = lis$distancesToActiveCenter,
                 "fps_sequences" = lis$fps_sequences,
                 "alpha" = alpha,
                 "betha" = betha,
                 "measure" = measure)
      
      
      saveRDS(lis,file = fNameModelDownsampled_additional_params)
    }
    
  }
    lis = readRDS(fNameModelDownsampled_additional_params)
    my_print(fNameModelDownsampled_additional_params, 1)
  
  if(plot){
    plotProteinModel(fNameOrigingal, lis)
  }
  
  return(lis)
}


get_quantile_name_protein <- function(name, n,m,q, measureNN, alpha, betha, locale){
  paste(name,"_n_",n,"_m_",m,"_q_",q,"_muNN_",measureNN,"_alpha_", alpha, "_betha_", betha, "_loc_", locale, ".csv",sep ="")
}



get_quantiles_protein <- function(path = "/home/willy/Schreibtisch/106Test/Output/", model, n, m, q, measureNN, alpha, betha, recalculate = FALSE, locale = TRUE){
  quantFolder = paste(path = path, "/", model$name, "/Quantiles/", sep ="")
  if(!dir.exists(quantFolder)) dir.create(quantFolder)
  
  quantFileName = paste(quantFolder, get_quantile_name_protein("quantiles", n, m, q, measureNN, alpha = alpha, betha = betha, locale = locale), sep = "")
  
  if(!file.exists(quantFileName) || recalculate == TRUE){
    quantiles = c()
    if(locale){
      quantiles = localEccentricitiesAndGetQuantiles(model = model, n = n, q = q)
    }else{
      quantiles = sampleEccentricitiesAndGetQuantiles(model = model, n = n, m = m, q = q)
    }

    write.csv(quantiles,file = quantFileName, row.names = FALSE)
  }
  
  quantiles = read.csv(quantFileName, header = TRUE, colClasses=c("character",rep("numeric",(q+2)*6)))
  return(quantiles)
}

# install.packages("doParallel")
# library(doParallel)
# registerDoParallel(20)
# 
# 
# foreach(i = 1:1000, .combine=rbind) %dopar%{
#   i^2
# }


get_quantiles_all_proteins <- function(model_vec, path, n, m ,q, recalculate = FALSE, locale = TRUE){
  # quantiles_all = data.frame(matrix(0,ncol = (q+2)*3+1, nrow = m*length(model_vec)))
  # colnames(quantiles_all) = c("name",   as.vector(paste(c("q_pos"),c(1:(q+2)), sep = "")), 
  #                             as.vector(paste(c("q_neg"),c(1:(q+2)), sep = "")),
  #                             as.vector(paste(c("q_pos_neg"),c(1:(q+2)), sep = "")))
  
  num = m
  if(locale) num = nrow(model_vec[[1]]$centers)*length(model_vec)
  
  print(num)
  
  quantiles_all = data.frame(matrix(0,ncol = (q+2)*6+1, nrow = num))
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
  alpha = model_vec[[1]]$alpha
  betha = model_vec[[1]]$betha
  QuantName = paste(QuantFolder,get_quantile_name_protein(name = "All",n = n, m = m, q = q, measureNN = measureNN,alpha = alpha, betha = betha, locale = locale), sep ="")
  
  QuantNameTransformed = paste(QuantFolder,get_quantile_name_protein(name = "All_transformed",n = n, m = m, q = q, measureNN = measureNN,alpha = alpha, betha = betha, locale = locale), sep ="")
  
  
  # if(!file.exists(QuantName) || recalculate == TRUE || !file.exists(QuantNameTransformed)){
    if(!file.exists(QuantName) || recalculate == TRUE){
    if(!dir.exists(QuantFolder)) dir.create(QuantFolder)
    
    # this is a bit hacky, we circumnavigate the combine-methods of foreach
    # in a first round we calculate all the quantiles and in a second loop we collect all the files
    foreach(i = 1:length(model_vec)) %dopar%{
      quant = get_quantiles_protein(path = path, model = model_vec[[i]],n = n, m = m, q = q, recalculate = recalculate, measureNN = measureNN, alpha = alpha, betha = betha, locale = locale)
    }
      
    for(i in 1:length(model_vec)){
    # mat = 
      print(paste(model_vec[[i]]$name, i/length(model_vec)))
      
      quant = get_quantiles_protein(path = path, model = model_vec[[i]],n = n, m = m, q = q, recalculate = FALSE, measureNN = measureNN, alpha = alpha, betha = betha, locale = locale)
  
      start_ind = (i-1)*nrow(quant)+1
      end_ind = start_ind + nrow(quant)-1
  
      quantiles_all[start_ind:end_ind,] = quant
      quantiles_all[start_ind:end_ind,1] = as.character(quant[,1])
    }
    
    write.csv(x = quantiles_all,file = QuantName, row.names = FALSE)
    
    
    # if(!file.exists(QuantNameTransformed)){
    #   quantTransformed = transformQuantsProtein(quantiles_all,q = q)
    #   write.csv(x = quantTransformed,file = QuantNameTransformed, row.names = FALSE)
    # }
  } else {
    quantiles_all = read.csv(file = QuantName, header = TRUE)
  }
  
  return(quantiles_all)
}


# models_small = getAllProteinModels(path = getPath("106Experiment"), n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 10, recalculate = FALSE, alpha = 3,betha = 3)
# quantiles2 = get_quantiles_all_proteins(model_vec = models_small, path = getPath("106Experiment"), n = 0.2, m = 1, q = 1, recalculate = FALSE,locale = TRUE)
# 
# quantiles2


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
  
  if(length(functionalInds) == 0) nonFunctionalInds = c(1:nrow(quantiles))
  
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


getBoarderNNCount <- function(model,neighbors = 20){
  # Count how many points are in the closest neighbors 
  # from the other potential
  #----------------------------------------------------
  
  counts = rep(1,nrow(model$centers))
  for(i in 1:nrow(model$centers)){
    nearestPoints = which.minn(model$d_surface[i,], n = neighbors)
    if(model$posNegVector[i] == FALSE && sum(model$posNegVector[nearestPoints]) > 0){
      counts[i] = sum(model$posNegVector[nearestPoints])
    }
    
    if(model$posNegVector[i] == TRUE && sum(model$posNegVector[nearestPoints]) < neighbors){
      counts[i] = (neighbors - sum(model$posNegVector[nearestPoints]))
    }
  }
  
  return(counts)
}

# fName = "/home/willy/PredictingProteinInteractions/data/120Experiment/Output/1ezk/1ezk_active_center.csv"
# info = file.info(fName)
# info$size
# 
# fName = "/home/willy/PredictingProteinInteractions/data/120Experiment/Output/1aba/1aba_active_center.csv"
# info = file.info(fName)
# info$size
# 
# info
# actCenterMolecules = read.csv2(fName)


getACtiveCenterAndAtomCoords <- function(path, protName){
  # from the pqr-files and after executing
  # centerSelect the exact coordinates of 
  # the active center are extracted
  # returns a list with the coordinates of all
  # atoms, and the indices of the active center
  #------------------------------------------
  filename = paste(path, "/", protName, "/", protName, "_active_center.csv", sep ="")
  info = file.info(filename)
  
  # only then it worked, otherwise we are probably missing a motife still
  if(!is.na(info$size) && info$size > 5){
    print("reading active center ...")
    
    # print(filename)
    
    actCenterMolecules = read.csv2(filename)
    
    filenamePqr = paste(path, "/", protName, "/", protName, "HeadOff.pqr", sep ="")
    # print(filenamePqr)
    
    pqr = read.table(filenamePqr)
    
    
    
    proteinStructure = pqr[,6:8]
    proteinStructure[,1] = proteinStructure[,1] / (max(proteinStructure[,1])- min(proteinStructure[,1]))
    proteinStructure[,2] = proteinStructure[,2] / (max(proteinStructure[,2])- min(proteinStructure[,2]))
    proteinStructure[,3] = proteinStructure[,3] / (max(proteinStructure[,3])- min(proteinStructure[,3]))
    
    geoM = rep(0,3)
    geoM[1] = sum(proteinStructure[,1])/nrow(proteinStructure)
    geoM[2] = sum(proteinStructure[,2])/nrow(proteinStructure)
    geoM[3] = sum(proteinStructure[,3])/nrow(proteinStructure)
    
    proteinStructure[,1] = proteinStructure[,1] - geoM[1]
    proteinStructure[,2] = proteinStructure[,2] - geoM[2]
    proteinStructure[,3] = proteinStructure[,3] - geoM[3]
    
    colnames(proteinStructure) = c("x","y","z")
    
    return(list("atomCoords" = proteinStructure, "activeCenterInds" = actCenterMolecules$start:actCenterMolecules$end))
  } else {
    print("WARNING: Missing active center. Using (0,0,0) as approximation.")
    
    proteinStructure = matrix(0, ncol = 3, nrow = 2)
    return(list("atomCoords" = proteinStructure, "activeCenterInds" = c(1,2)))
  }

}


sampleMultipleTimesFps  <- function(model_vec, quantiles, sampleSize, sampleTimes, m = 100, fps = FALSE){
  # m number of points per model
  #----------------------------------------------
  
  # name in first column
  quants = ncol(quantiles)-1

  subClassNames_train = unique(quantiles[,1])

  relevantQuantiles = as.matrix(quantiles[,(c(1:quants)+1)])

  MatOut = data.frame(matrix(0, nrow = nrow(quantiles)/m*sampleTimes, ncol = ncol(relevantQuantiles)*sampleSize+1))


  for(k in 1:length(subClassNames_train)){
    print(k/length(subClassNames_train))

    startInd = (k-1)*m +1
    endInd = startInd + m -1

    sampled_indices = matrix(0, nrow = sampleTimes, ncol = sampleSize)
    if(fps == TRUE){
      # indices of the rows of which we will take the fps-sequences
      rowinds = sample(c(1:m), sampleTimes, replace = FALSE)
      # print(rowinds)
      sampled_indices = model_vec[[k]]$fps_sequences[rowinds,c(1:sampleSize)]
      
      # print(sampled_indices)
    } else {
      for(i in 1:sampleTimes){
        sampled_indices[i,] = sample(c(1:m), sampleSize, replace = FALSE)
      }
    }

    for(i in 1:sampleTimes){

      
      relevantQuantilesRows = relevantQuantiles[startInd:endInd,]
      
      augmented = relevantQuantilesRows[sampled_indices[i,],]

      rsControl = rowSums(augmented)
      orderingControl = order(rsControl)

      augmented = augmented[orderingControl,]


      startInd2 = (k-1)*sampleTimes+i
      
      augmentedFinal = matrix(as.matrix(augmented), nrow = 1)

      MatOut[startInd2,-1] = augmentedFinal
      # MatOut[startInd2,-1] = matrix(augmented,nrow = 1)
      MatOut[startInd2,1] = as.character(subClassNames_train[k])
    }
  }

  colnames(MatOut)[1] = "name"


  return(MatOut)
}

# mods = getAllProteinModels(path = getPath("106Experiment"),n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000,measureNearestNeighbors = 10,alpha = 0,betha = 0,recalculate = FALSE)
# quantiles2 = get_quantiles_all_proteins(model_vec = mods, path = getPath("106Experiment"), n = 0.1, m = 1000, q = 1, recalculate = FALSE, locale = TRUE)
# 
# ## "reading active center ..."
# #

# sS = 100
# MatOut = sampleMultipleTimesFps(model_vec = mods,quantiles = quantiles2,sampleSize = sS,sampleTimes = 10,m = 1000,fps = TRUE)
# 
# # [,1] [,2] [,3] [,4]
# # [1,]  674  513  200  106
# # [2,]  286  835  748  634
# 
# MatOut[nrow(MatOut), ]
# 
# rowInds = mods[[106]]$fps_sequences[455,1:sS]
# 
# start = 1000*106-1000+1
# subQs = quantiles2[start-1 +rowInds,-1]
# rs = rowSums(subQs)
# ordering = order(rs)
# 
# augmented = subQs[ordering,]
# 
# augmentedFinal = matrix(as.matrix(augmented), nrow = 1)
# 
# 
# MatOutWithoutName = MatOut[nrow(MatOut), -1]
# 
# augmentedFinal
# 
# su = 0
# for(i in 1:ncol(augmentedFinal)){
#   su = su + abs(augmentedFinal[1,i] - MatOutWithoutName[1,i])
# }
# 
# su
# 
# 

# 
# 
# MatOut[1:100,1:5]
# 
# quants = 3*6
# 
# 
# i = 26
# sum(MatOut[i,2:19]) - sum(MatOut[i,20:37])
# 
# ncol(MatOut)
# 
# MatOut[26,]
# 
# nrow(MatOut)
# testFpsSampling(MatOut)
# 
# rS = rowSums(quantiles2[5001:6000,-1])
# ordering = order(rS)
# 
# 
# s = sum(MatOut[26,2:19])
# 
# which(as.vector(rS) ==  s)
# 
# 
# quantiles2[5001:6000,1]
# 
# rS[ordering]
# 
# selectedPoints = mods[[1]]$fps_sequences[1,1:10]
# 
# 
# selections = rep(0,length(selectedPoints))
# for(j in 1:length(selectedPoints)){
#   selections[j] =which(ordering == selectedPoints[j])
# }
# 
# selectedPoints[order(selections)]
# 
# 




testFpsSampling <- function(MatOut){
  half = (ncol(MatOut)-1)/2
  
  
  su = 0
  
  # should be increasingly ordered
  # therefore when the first half is bigger than the second half
  # then something is wrong
  for(i in 1:nrow(MatOut)){
    if(sum(MatOut[i,1:half + 1]) - sum(MatOut[i,(half+1):(2*half) + 1]) > 0 ) {
      print(paste("ERROR", MatOut[i,1], sep = " "))
      print(sum(MatOut[i,1:half + 1]) - sum(MatOut[i,(half+1):(2*half) + 1]))
      # return(FALSE)
      
      su = su+1
    }
  }
  
  print(nrow(MatOut) - su)
  
  return(TRUE)
}


# 
# # "start";"end"
# # 779;823
# 
# 
# MatOut[,1]
# 
# unique(MatOut[,1])
# 
# unique(quantiles2[,1])
# 
# quantiles2[,1] = paste("1_",quantiles2[,1], sep = "")
# 
# getSamplesSurf3(quantiles = quantiles2,model_vec = mods,sampleSize = 10, sampleTimes = 10, numClasses = 2,m = 1000)
# 
# quantiles2[,1]
# 
# # colnames(MatOut)
# 
# 
# prot$fps_sequences
# 
# 
# plotProteinModel(lis = prot,totalMassSize = 1)
# 
# points3d(prot$centers[prot$fps_sequences[1,1:10],],size = 10)
# points3d(prot$centers[prot$fps_sequences[2,],],size = 11, col = "green")
# points3d(prot$centers[prot$fps_sequences[3,],],size = 11, col = "blue")
# points3d(prot$centers[prot$fps_sequences[4,],],size = 11, col = "red")
# points3d(prot$centers[sample(nrow(prot$centers),10),],size = 11, col = "red")
# points3d(prot$centers[sample(nrow(prot$centers),10),],size = 11, col = "blue")


sampleMultipleTimesFpsParallel  <- function(model_vec, quantiles, sampleSize, sampleTimes, m = 100, fps = FALSE){
  # name in first column
  quants = ncol(quantiles)-1
  
  subClassNames_train = unique(quantiles[,1])
  
  relevantQuantiles = as.matrix(quantiles[,(c(1:quants)+1)])
  print("1")
  
  MatOut = data.frame(matrix(0, nrow = nrow(quantiles)/m*sampleTimes, ncol = ncol(relevantQuantiles)*sampleSize+1))
  
  quantTrainMat = t(apply(matrix(c(1:(numPermutations*sampleTimes*length(subClassNames_train))), ncol = 1),1,
                          FUN = function(k){
                            class_ind= ceiling(k/sampleTimes)
                            
                            subClassName_tmp = as.character(subClassNames_train[class_ind])
                            
                            train_indices = ((class_ind-1)*m+1):(class_ind*m)
                            
                            sampled_indices = sample(train_indices, size = sampleSize, replace = FALSE)
                            
                            ordered = relevantQuantiles[sampled_indices,]
                            
                            if(sort == TRUE){
                              
                              
                              ordered = ordered[order(rowSums(ordered)),]
                            }
                            
                            c(subClassName_tmp,as.vector(t(ordered)))
                          }))
  
  quantTrain = data.frame(quantTrainMat, stringsAsFactors=FALSE)
  colnames(quantTrain)[1] = "name"
  
  return(quantTrain)
}





getSamplesSurf3 <- function( quantiles,
                             model_vec,
                             sampleSize, 
                             sampleTimes,
                             numClasses,
                             m,
                             splitPattern = "_",
                             fps = TRUE){
  
  print("starting sampling ...")
  sampledQuantiles = sampleMultipleTimesFps(model_vec = model_vec,
                                            quantiles = quantiles,
                                            sampleSize = sampleSize,
                                            sampleTimes = sampleTimes,
                                            m = m,
                                            fps = fps)
  
  un = unique(sampledQuantiles[,1])
  if(length(un) > 1){
    v = rep(0,length(un))
    for(i in 1:length(un)){
      v[i] = length(which(sampledQuantiles[,1] == un[i]))
    }
    
    print(un)
    
    if(var(v) != 0) return(NULL)
  }

  
  print("generating x and y ...")
  
  y_all_class_names = getClassNamesFromSubClasses(sampledQuantiles[,1],splitPattern = splitPattern)
  
  y = y_all_class_names
  X = as.matrix(sampledQuantiles[,-1])
  
  classLevels = unique(y)
  y = as.numeric(as.factor(y))-1
  y <- to_categorical(y, numClasses)
  
  y_original_names = sampledQuantiles[,1]
  
  return(list("X" = X, "y" = y, "y_original_names" = y_original_names))
}







getDistancesToActiceCenter <- function(atomCoords, activeCenterInds, centers){
  # for each point of the surface calculate the distance
  # to the closest point from the active center
  #-----------------------------------------------
  
  return(as.vector(knnx.dist(atomCoords[activeCenterInds,],centers, k = 1)))
}

getProteinMeasure <- function(distancesToActiveCenter, boarderCounts, alpha = 0, betha = 0){
  # alpha ... sepcifies the importance of the distance to the active site
  # betha ... specifies the importance of a boarder area
  #------------------------------------------------------------------------------

  distNorm = normalize(1/(distancesToActiveCenter))
  boarderNorm = normalize(boarderCounts)
    
  return(normalize(distNorm^alpha* boarderNorm^betha))
}

getQuantilesAlphaBetha <- function(path = path, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, 
                                   measureNearestNeighbors = 10,
                                   recalculate = FALSE,
                                   recalculateQuants = TRUE,
                                   alpha = 2,
                                   betha = 1,
                                   locale = TRUE,
                                   n, m,q){
  models1 = getAllProteinModels(path = path, n_s_euclidean = n_s_euclidean,n_s_dijkstra = n_s_dijkstra,
                                stitchNum = stitchNum, measureNearestNeighbors = measureNearestNeighbors, recalculate = recalculate, alpha = alpha,betha = betha)
  quantiles2 = get_quantiles_all_proteins(model_vec = models1, path = path, n = n, m = m, q = q, recalculate = recalculateQuants, locale = locale)
  return(quantiles2)
}


modelProt1 <- function(TrainTest, epochs = 30, batch_size = 64, weights, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model7 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 300, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 200, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  # history <- model %>% fit(
  #   x_train, y_train, 
  #   epochs = epochs, batch_size = batch_size, 
  #   class_weight = list("1"=106/89,"0"=106/17),
  #   validation_split = 0.05
  # )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  return(model)
}

modelProt2 <- function(TrainTest, epochs = 30, batch_size = 64, sampleSize = NULL, sampleTimes =NULL, q = NULL, weights){
  print("Calling model2 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 2400, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 400, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 200, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  # history <- model %>% fit(
  #   x_train, y_train, 
  #   epochs = epochs, batch_size = batch_size, 
  #   class_weight = list("1"=106/89,"0"=106/17),
  #   validation_split = 0.05
  # )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  
  return(model)
}


modelProt3 <- function(TrainTest, epochs = 30, batch_size = 64, weights, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model3 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 600, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 400, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 200, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 20, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 20, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 20, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  return(model)
}

modelProt4 <- function(TrainTest, epochs = 30, batch_size = 64, weights, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model4 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 1200, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 800, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 400, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 40, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 40, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 40, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  return(model)
}

modelProt5 <- function(TrainTest, epochs = 30, batch_size = 64, weights, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  # for inputs with 900 points
  print("Calling model5 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 1300, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 9000, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 4500, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 450, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 450, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 450, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  # history <- model %>% fit(
  #   x_train, y_train, 
  #   epochs = epochs, batch_size = batch_size, 
  #   class_weight = list("1"=106/89,"0"=106/17),
  #   validation_split = 0.05
  # )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  return(model)
}



# clamp(c(1,2,3,4,2,2), lower = 0, upper = 1)

EPSILON = 0.00000001
recall_keras_metric <- function(y_true, y_pred){
  # """Recall metric.
  # 
  #   Only computes a batch-wise average of recall.
  # 
  #   Computes the recall, a metric for multi-label classification of
  #   how many relevant items are selected.
  # """
  
  # true_positives = sum(round(clamp(y_true * y_pred, 0, 1)))
  # possible_positives = sum(round(clamp(y_true, 0, 1)))
  # recall = true_positives / (possible_positives + EPSILON)
  K <- backend()
  
  true_positives = K$sum(K$round(K$clip(y_true * y_pred, 0, 1)))
  possible_positives = K$sum(K$round(K$clip(y_true, 0, 1)))
  recall = true_positives / (possible_positives + EPSILON)
  
  return(recall)
}

# recall_keras_metric(c(1,0,1,0,1), c(1,1,1,0,0))

precision_keras_metric <- function(y_true, y_pred){
  #   """Precision metric.
  # 
  #         Only computes a batch-wise average of precision.
  # 
  #         Computes the precision, a metric for multi-label classification of
  #         how many selected items are relevant.
  #         """
  # true_positives = sum(round(clamp(y_true * y_pred, 0, 1)))
  # predicted_positives = sum(round(clamp(y_pred, 0, 1)))
  # precision = true_positives / (predicted_positives + EPSILON)
  
  K <- backend()
  
  true_positives = K$sum(K$round(K$clip(y_true * y_pred, 0, 1)))
  predicted_positives = K$sum(K$round(K$clip(y_pred, 0, 1)))
  precision = true_positives / (predicted_positives + EPSILON)
  
  return(precision)
}
# precision_keras_metric(c(1,0,1,0,1), c(1,0,1,0,0))


f1_keras_metric <- function(y_true, y_pred){
  precision = precision_keras_metric(y_true, y_pred)
  recall = recall_keras_metric(y_true, y_pred)
  
  print(precision)
  
  return(2*((precision*recall)/(precision+recall+EPSILON)))
}
# f1_keras_metric(c(1,0,1,0,1), c(1,0,1,0,0))

f1_keras_metric_wrapper <- custom_metric("f1", function(y_true, y_pred) {
  # print("getting called ...")
  return(f1_keras_metric(y_true, y_pred))
})

# f1_keras_metric_wrapper(c(1,0,1,0,1), c(1,0,1,0,0))


modelProt1_f1 <- function(TrainTest, epochs = 30, batch_size = 64, weights, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model1 f1 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 300, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 200, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(),
    metrics = f1_keras_metric_wrapper
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  return(model)
}

conv = FALSE
modelProt_custom <- function(TrainTest, weights, layers = c(10,10), dropOuts = c(0.1,0.1), optimizerFunName = "optimizer_adam",metrics = "f1_keras_metric_wrapper", batch_size = 32, epochs = 20){
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  # define the architecture of the model
  model <- keras_model_sequential()
  if(conv == FALSE){
    # model %>%   layer_dense(units = layers[1],input_shape = c(ncol(x_train))) %>% 
    #   layer_dropout(rate = dropOuts[1]) %>% 
    #   layer_activation("relu")
    # 
    # if(length(layers) > 1){
    #   for(i in 2:length(layers)){
    #     model %>%   layer_dense(units = layers[i])%>%layer_dropout(rate = dropOuts[i]) %>% layer_activation("relu")
    #   }
    # }
    # model %>% layer_dense(units = numClasses, activation = 'softmax')
    
    
    model %>%   layer_dense(units = layers[1],input_shape = c(ncol(x_train))) %>% 
                layer_dropout(rate = dropOuts[1]) %>% 
                layer_activation("relu") %>%
                
                layer_dense(units = layers[2], activation = "relu") %>%
                layer_dropout(rate = dropOuts[2]) %>% 
                layer_activation("relu") %>%
 
                layer_dense(units = layers[3], activation = "relu") %>%
                layer_dropout(rate = dropOuts[3]) %>% 
                layer_activation("relu") %>%     
      
                layer_dense(units = layers[4], activation = "relu") %>%
                layer_dropout(rate = dropOuts[4]) %>% 
                layer_activation("relu") %>%
      
                layer_dense(units = numClasses, activation = 'softmax')
    
    
  }else{
    channels = 6*4
    dropOutsFilters = 0.4
    model %>% 
      layer_reshape(target_shape = c(SAMPLESIZE, ncol(x_train)/(SAMPLESIZE*channels),channels),input_shape = c(ncol(x_train))) %>%
      
      layer_conv_2d(filter = 16, kernel_size = c(3,3), padding = 'same', input_shape = c(SAMPLESIZE, ncol(x_train)/SAMPLESIZE,1)) %>%
      layer_activation("relu") %>%
      layer_average_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
      layer_dropout(dropOutsFilters) %>%

      layer_conv_2d(filter = 16, kernel_size = c(3,3), padding = 'same') %>%
      layer_activation("relu") %>%
      layer_average_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
      layer_dropout(dropOutsFilters) %>%
      
      layer_conv_2d(filter = 16, kernel_size = c(3,3), padding = 'same') %>%
      layer_activation("relu") %>%
      layer_average_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
      layer_dropout(dropOutsFilters) %>%
      
      layer_flatten() %>%
      
      layer_dense(units = 50)%>%
      layer_dropout(0.2) %>%
      layer_activation("relu") %>%
      
      layer_dense(units = 50)%>%
      layer_dropout(0.2) %>%
      layer_activation("relu") %>%
      
      layer_dense(units = numClasses, activation = 'softmax')
    
    
    #   layer_dense(units = layers[1], rate = dropOuts[1]) %>% 
    #   layer_activation("relu")
    #   
    # if(length(layers) > 1){
    #   for(i in 2:length(layers)){
    #     model %>%   layer_dense(units = layers[i])%>%layer_dropout(rate = dropOuts[i]) %>% layer_activation("relu")
    #   }
    # }
    # model %>% layer_dense(units = numClasses, activation = 'softmax')
  }

  
  optimizerFun = optimizer_adam
  print("specifiying optimizer ...")
  if(optimizerFunName == "optimizer_rmsprop") optimizerFun = optimizer_rmsprop
  
  metricsFun = f1_keras_metric_wrapper
  if(metrics == "accuracy") metricsFun = c("accuracy")
  
  # compile the model
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizerFun(),
    metrics = metricsFun
  )
  
  # train the model
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  
  return(model)
}

# 
# TR = readRDS("/home/willy/PredictingProteinInteractions/data/106Test/NNexperiments/Test27sS_10_sT_200_sTt_10_euklid_TRUE_pos_TRUE_neg_TRUE_pos_neg_TRUE_globalToo_TRUE.Rdata")
# 
# model = modelProt_custom(TrainTest = TR,
#                  weights = list("0" = 1, "1" = 1),
#                  layers = c(50,10,10,5),
#                  dropOuts = c(0.4,0.1,0.1,0.1),
#                  optimizerFunName = "optimizer_rmsprop",
#                  metrics = "f1_keras_metric_wrapper")
# 
# model = modelProt_custom(TrainTest = TR,
#                          weights = list("0" = 1, "1" = 1),
#                          layers = c(5,5),
#                          dropOuts = c(0.1,0.1),
#                          optimizerFunName = "optimizer_rmsprop",
#                          metrics = "f1_keras_metric_wrapper")
# 
# 
# lab = read.table(LABELS, header = TRUE)
# folds = createFolds(lab$label, k = 10,list = TRUE)
# 
# y_test_name_inds = folds[[1]]
# 
# test_inds = which(originalNames %in% protNames[y_test_name_inds])
# train_inds = c(1:nrow(TrainTest$X))[-test_inds]
# 
# if(k == 1) train_inds = c(1:nrow(TrainTest$X))
# 
# trainNames = protNames[-y_test_name_inds]
# testNames = protNames[y_test_name_inds]
# 
# Train_X = TrainTest$X[train_inds,]
# Train_y = TrainTest$y[train_inds,]
# 
# Test_X = TrainTest$X[test_inds,]
# Test_y = TrainTest$y[test_inds,]
# 
# Train_X = apply(Train_X, 2, FUN = function(i) as.numeric(as.character(i)))
# Test_X = apply(Test_X, 2, FUN = function(i) as.numeric(as.character(i)))
# 
# Train_y = apply(Train_y, 2, FUN = function(i) as.numeric(as.character(i)))
# Test_y = apply(Test_y, 2, FUN = function(i) as.numeric(as.character(i)))
# 
# 
# if(normalizeInputs){
#   colMins = apply(Train_X,2,min)
#   colMaxs = apply(Train_X,2,max)
#   colRanges = colMaxs - colMins
#   Train_X = sapply(c(1:length(colRanges)), FUN = function(i){ Train_X[,i]/colRanges[i] })
#   Test_X = sapply(c(1:length(colRanges)), FUN = function(i){ Test_X[,i]/colRanges[i] })
# }
# 
# 
# shuf = shuffle(1:nrow(Train_X))
# TrFinal = list("x_train" = Train_X[shuf,], "y_train" = Train_y[shuf,], "x_test" = Test_X, "y_test" = Test_y, "numClasses" = numClasses, "classLevels" = classLevels, "mapping" = mapping)
# 
# classLabels = reverseToCategorical(oneHot = TrFinal$y_train, mapping$name)
# 
# # print(mapping)
# # print(classLabels)
# 
# # classLevels = unique(classLabels)
# weights = rep(0,length(classLevels))
# for(i in 1:length(weights)){
#   weights[i] = 1/(length(which(classLabels == classLevels[i]))/length(classLabels))
# }
# 
# # print(weights)
# weights <- split(weights, mapping$name-1)
# print(weights)
# 
# 
# # return(TrFinal)
# 
# fac = 1
# if(euklid) fac=2
# model = modelFUN(TrainTest = TrFinal,sampleSize = sampleSize,sampleTimes = sampleTimes,q = (q+2)*3*fac,epochs = epochs, batch_size = batch_size,weights = weights)




modelProt2_f1 <- function(TrainTest, epochs = 30, batch_size = 64, weights, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model2 f1 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 300, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 200, activation = 'relu') %>% 
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(),
    metrics = f1_keras_metric_wrapper
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  return(model)
}


modelProt3_f1 <- function(TrainTest, epochs = 30, batch_size = 64, weights, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model3 f1 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 300, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(),
    metrics = f1_keras_metric_wrapper
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  return(model)
}




modelProt4_f1 <- function(TrainTest, epochs = 30, batch_size = 64, weights, sampleSize = NULL, sampleTimes =NULL, q = NULL ){
  print("Calling model4 f1 ...")
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 100, activation = 'relu', input_shape = c(ncol(x_train))) %>% 
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 20, activation = 'relu') %>% 
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 10, activation = 'relu') %>% 
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = numClasses, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(),
    metrics = f1_keras_metric_wrapper
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size, 
    class_weight = weights,
    validation_split = 0.05
  )
  return(model)
}


# ProteinsExperimentKfoldCV( sampleSize = 5,
#                            sampleTimes = 10,
#                            sampleTimes_test = 10,
#                            batch_size = 10,
#                            epochs = 5,
#                            euklid = TRUE,
#                            q = 1,
#                            m = 1000,
#                            numClasses = 2,
#                            fNameTrain = "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                            ExperimentName = "Test999",
#                            modelName = getVarName(modelProt1_f1),
#                            modelFUN = modelProt1_f1,
#                            recalculate = FALSE,
#                            k = 1,
#                            onlySummarizeFolds = FALSE,
#                            normalizeInputs = TRUE,
#                            saveExperiment = TRUE)





convModelProtein <- function(TrainTest, sampleSize, sampleTimes, q, epochs = 30, batch_size = 64, weights){
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_reshape(target_shape = c(sampleSize, ncol(x_train)/sampleSize,1),input_shape = c(ncol(x_train))) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same', input_shape = c(sampleSize, ncol(x_train)/sampleSize,1)) %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 30, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 10, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 10, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 10, kernel_size = c(3,3), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 10, kernel_size = c(4,4), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(3,3),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    
    layer_flatten() %>%
    layer_dense(400) %>%
    layer_dense(200) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    
    layer_dense(100) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(100) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(units = numClasses, activation = 'softmax')
  
  # ?layer_global_average_pooling_1d
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size,
    class_weight = weights,
    validation_split = 0.1
  )
  
  
  return(model)
}

convModelProtein2 <- function(TrainTest, sampleSize, sampleTimes, q, epochs = 30, batch_size = 64, weights){
  
  x_train = TrainTest$x_train
  y_train = TrainTest$y_train
  x_test = TrainTest$x_test
  y_test = TrainTest$y_test
  
  numClasses = TrainTest$numClasses
  
  
  #---------------------------------------------------------
  model <- keras_model_sequential()
  model %>% 
    layer_reshape(target_shape = c(sampleSize, ncol(x_train)/sampleSize,1),input_shape = c(ncol(x_train))) %>%
    
    layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same', input_shape = c(sampleSize, ncol(x_train)/sampleSize,1)) %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_max_pooling_2d(pool_size = c(2,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_average_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_average_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_average_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    layer_conv_2d(filter = 32, kernel_size = c(3,1), padding = 'same') %>%
    layer_activation("relu") %>%
    layer_average_pooling_2d(pool_size = c(3,1),padding = 'same') %>%
    layer_dropout(0.1) %>%
    
    
    layer_flatten() %>%
    layer_dense(400) %>%
    layer_dense(200) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    
    layer_dense(100) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(100) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(10) %>%
    layer_activation("relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(units = numClasses, activation = 'softmax')
  
  # ?layer_global_average_pooling_1d
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size,
    class_weight = weights,
    validation_split = 0.1
  )
  
  
  return(model)
}


getStatsOnAccuracies <- function(path = "/home/willy/PredictingProteinInteractions/Results/TablesProt/"){
  
  files = list.files(path = path, pattern = "Accuracy", full.names = FALSE)
  filesFullNames = list.files(path = path, pattern = "Accuracy", full.names = TRUE)
  
  
  df = data.frame(matrix(0, ncol = 2, nrow = length(files)))
  colnames(df) = c("Test", "Accuracy")
  
  for(i in 1:length(files)){
    s2 = strsplit(files[i],split = ".tex")[[1]][1]
    s3 = strsplit(s2,split = "_")[[1]][2]
    
    df[i,1] = s3
    t = read.table(filesFullNames[i])
    df[i,2] =t[1,1]
  }
  
  df = df[order(df[,2],decreasing = TRUE),]
  return(df)
}

joinStats <- function(path = "/home/willy/PredictingProteinInteractions/Results/TablesProt/"){
  
  files = list.files(path = path, pattern = "stats", full.names = FALSE)
  filesFullNames = list.files(path = path, pattern = "stats", full.names = TRUE)
  
  print(paste("reading ", length(files), " different stats-files ...", sep =""))
  
  # return(files[1])
  
  df = read.csv(filesFullNames[1],header = TRUE)
  
  if(length(files) > 1){
    for(i in 2:length(files)){
      df2 = read.csv(filesFullNames[i],header = TRUE)
      df = rbind(df,df2)
    }
  }
  
  name = rep("", length(files))
  for(i in 1:length(files)){
    name[i] = strsplit(files[i],split = "_")[[1]][1]
  }
  
  df = cbind(name, df)
  
  
  df = df[order(df$f1_score, decreasing = TRUE),]
  
  return(df)
}

#-------------------------------------------------------------------------------------------------------------
# Generate the quantiles
#-------------------------------------------------------------------------------------------------------------
if(mode == "onlyGenerateModels" || mode == "Booth"){
  alphas = c(1)
  bethas = c(3)
  
  alpha_betha_grid = expand.grid(alphas,bethas)
  
  nlocals = c(0.05,0.1,0.2,0.5,0.8)
  # nlocals = c(0.8)
  
  for(alpha in alphas) {
    for(i in 1:length(nlocals)) {
      # tmp = getQuantilesAlphaBetha(alpha = alphas[i],betha = bethas[j], n = 0.2, m = 1,q = 1, locale = TRUE, path = path106Experiment, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 20, recalculate = FALSE,recalculateQuants = TRUE)
      tmp = getQuantilesAlphaBetha(alpha = alpha,betha = alpha, n = nlocals[i], m = 1,q = 1, locale = TRUE, path = pathToExperiment, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 10, recalculate = FALSE,recalculateQuants = FALSE)
    }
  }
}

# tmp = getQuantilesAlphaBetha(alpha = 0,betha = 0, n = 1, m = 1,q = 1, locale = FALSE, path = path106Experiment, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 20, recalculate = FALSE,recalculateQuants = TRUE)
# 
# 
# tmp[1,]
# 
# apply(c(2:ncol(tmp)), 2, FUN = function(i){
#           normalize(tmp)
# })
# 
# # install.packages("BBmisc")
# library(BBmisc)
# normalized = BBmisc::normalize(tmp[,-1],method = "range")
# 
# tNorm = tmp
# tNorm[,-1] = normalized
# 
# 
# functionals = c(getFunctionalProteins(),"000_Trx")
# plot_prot_quants(tNorm,q = 1,functionals = functionals)
# 
# plot_prot_quants(tmp,q = 1,functionals = functionals)
# # 
# # which(is.na(tmp) == TRUE)
# 


#------------------------------------------------------------------------
# Classification with a Neuronal Net model
#------------------------------------------------------------------------
getVarName <- function(v1) {
  deparse(substitute(v1))
}

writeExperimentParametersToFile <- function(pathToStats = "/home/willy/PredictingProteinInteractions/Results/TablesProt/",
                                            sampleSize = 20,
                                            sampleTimes = 10,
                                            sampleTimes_test = 10,
                                            batch_size = 1024,
                                            epochs = 300,
                                            euklid = TRUE,
                                            q = 1,
                                            m = 1000, 
                                            numClasses = 2,
                                            potentials = c("pos","neg","pos_neg"),
                                            fNameTrain = "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_2_betha_1_loc_TRUE.csv",
                                            fNameTest = "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_2_betha_1_loc_TRUE.csv",
                                            ExperimentName = "Test0",
                                            fNameTrain_global = NULL,
                                            fNameTest_global = NULL,
                                            path = "/home/willy/PredictingProteinInteractions/data/106Test/NNexperiments/",
                                            modelName = "NONE",
                                            recalculate = FALSE,
                                            reCalculateTrainTest = FALSE,
                                            labels = "/home/willy/PredictingProteinInteractions/data/labels.txt",
                                            accuracy = 0,
                                            f1_score = 0,
                                            auc = 0
){
  
  
  pos_flag = "pos" %in% potentials
  neg_flag = "neg" %in% potentials
  pos_neg_flag = "pos_neg" %in% potentials
  
  names =  c("accuracy","f1_score","auc","sampleSize", "sampleTimes", "sampleTimes_test", "batch_size",
             "epochs", "euklid","pos_flag", "neg_flag", "pos_neg_flag", "q", "m", "muNN" ,"numClasses", "alpha", "betha", "locPerc",
             "global", "model")
  
  df = data.frame(matrix(0, nrow = 1, ncol = length(names)))
  colnames(df) = names
  
  df$accuracy = accuracy
  df$f1_score = f1_score
  df$auc = auc
  
  df$sampleSize = sampleSize
  df$sampleTimes = sampleTimes
  df$sampleTimes_test = sampleTimes_test
  df$batch_size = batch_size
  df$epochs = epochs
  df$euklid = euklid
  df$q = q
  df$m = m
  df$numClasses = numClasses
  df$pos_flag = pos_flag
  df$neg_flag = neg_flag
  df$pos_neg_flag = pos_neg_flag
  
  alpha = strsplit(fNameTrain, "_")[[1]][11]
  betha = strsplit(fNameTrain, "_")[[1]][13]
  locPerc = strsplit(fNameTrain, "_")[[1]][3]
  muNN = strsplit(fNameTrain, "_")[[1]][9]
  
  df$alpha = as.numeric(alpha)
  df$betha = as.numeric(betha)
  df$locPerc = as.numeric(locPerc)
  df$model = modelName
  df$muNN = as.numeric(muNN)
  
  df$global = FALSE
  if(!is.null(fNameTrain_global))  df$global = TRUE
  
  write.csv(df,file = paste(pathToStats,"/",ExperimentName,"_stats.csv",sep = ""),row.names = FALSE)
  
  return(df)
}

#--------------------
# 
# ProteinsExperiment <- function(sampleSize = 20,
#                                  sampleTimes = 10,
#                                  sampleTimes_test = 10,
#                                  batch_size = 1024,
#                                  epochs = 300,
#                                  euklid = TRUE,
#                                  q = 1,
#                                  m = 1000, 
#                                  numClasses = 2,
#                                  potentials = c("pos","neg","pos_neg"),
#                                  fNameTrain = "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_2_betha_1_loc_TRUE.csv",
#                                  fNameTest = "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_2_betha_1_loc_TRUE.csv",
#                                  ExperimentName = "Test1",
#                                  fNameTrain_global = NULL,
#                                  fNameTest_global = NULL,
#                                  path = "/home/willy/PredictingProteinInteractions/data/106Test/NNexperiments/",
#                                  modelFUN = convModel4,
#                                  modelName,
#                                  recalculate = FALSE,
#                                  reCalculateTrainTest = FALSE,
#                                  labels = "/home/willy/PredictingProteinInteractions/data/labels.txt"){
#   
#   print("------------------------------------------------------")
#   print(paste("Experiment ", ExperimentName))
#   print("------------------------------------------------------")
#   
#   numPermutations = 1
#   
#   pos_flag = "pos" %in% potentials
#   neg_flag = "neg" %in% potentials
#   pos_neg_flag = "pos_neg" %in% potentials
#   
#   global_flag = FALSE
#   if(!is.null(fNameTrain_global) && !is.null(fNameTest_global)) global_flag = TRUE
#   
#   ExperimentFile = paste(path,"/", ExperimentName, "sS_",sampleSize, "_sT_", sampleTimes, "_sTt_", sampleTimes_test, "_euklid_",
#                          euklid,
#                          "_pos_", pos_flag,
#                          "_neg_", neg_flag,
#                          "_pos_neg_", pos_neg_flag,
#                          "_globalToo_", global_flag,
#                          ".Rdata", sep ="")
#   
#   print(ExperimentFile)
#   
#   TrFinal = list()
#   if(!file.exists(ExperimentFile) || recalculate){
#     quantilesTrain = read.csv(file =fNameTrain, header = TRUE)
#     
#     # quantilesTrain[,-1] = apply(quantilesTrain[,-1],2, FUN = function(i) as.numeric(as.character(i)))
#     
#     #-------------------------------------------------------------------------------------------
#     # which features to use
#     #-------------------------------------------------------------------------------------------
#     quants = c(1:(q+2))
#     # keep the name in any case
#     colIndsTokeep = c(1)
#     if(pos_flag) colIndsTokeep = c(colIndsTokeep, quants+1, length(quants)*3+1+quants)
#     if(neg_flag) colIndsTokeep = c(colIndsTokeep, length(quants)*1+1+quants, length(quants)*4+1+quants)
#     if(pos_neg_flag) colIndsTokeep = c(colIndsTokeep, length(quants)*2+1+quants, length(quants)*5+1+quants)
#     colIndsTokeep = sort(colIndsTokeep)
# 
#     #-------------------------------------------------------------------------------------------
#     # classlabels
#     #-------------------------------------------------------------------------------------------
#     protNames = as.character(unique(quantilesTrain[,1]))
#     lab = read.table("/home/willy/PredictingProteinInteractions/data/labels.txt", header = TRUE)
#     
#     # to circumnavigate problems in naming, ("_", is a problem) 
#     mapping = data.frame(matrix(0,ncol = 2, nrow = length(unique(lab$label))))
#     colnames(mapping) = c("originalName", "name")
#     mapping[,1] = sort(unique(lab$label))
#     mapping[,2] = as.numeric(as.factor(sort(unique(lab$label))))
#     
#     protClasses = unlist(lapply(c(1:length(protNames)), FUN = function(i){
#       mapping[which(mapping[,1] == lab$label[which(lab$name == protNames[i])]),2]
#     }))
#     
#     classLabels = unlist(lapply(c(1:nrow(quantilesTrain)), FUN = function(i){
#                                   protClasses[which(protNames == quantilesTrain[i,1])]
#                                 }))
#     
#     quantilesTrain[,1] = unlist(lapply(c(1:length(classLabels)), FUN = function(i){
#       paste(classLabels[i],"_",quantilesTrain[i,1], sep = "")
#     }))
#     
#     # functionals = c(getFunctionalProteins(), "000_Trx")
#     # protNames = unique(quantilesTrain[,1])
#     # functionalInds = which((as.character(as.factor(quantilesTrain[,1])) %in% functionals) == TRUE)
#     # 
#     # 
#     # classLabels = rep("notFunctional", nrow(quantilesTrain))
#     # classLabels[functionalInds] = rep("Functional", length(functionalInds))
#     # 
#     # quantilesTrain[,1] = unlist(lapply(c(1:length(classLabels)), FUN = function(i){
#     #   paste(classLabels[i],"_",quantilesTrain[i,1], sep = "")
#     # }))
#     
#     
#     # quantilesTrain = quantilesTrain[,c(1,8,9,10)]
#     classLevels = mapping$name
#     
#     
#     print("Creating train-set ...")
#     Train = getSamplesSurf2(quantilesTrain,sampleSize = sampleSize,sampleTimes = sampleTimes,euklid = euklid, numPermutations = numPermutations, numClasses = numClasses, m = m,reDo = reCalculateTrainTest)
# 
#     
#     quantilesTest = c()
#     if(fNameTrain != fNameTest){
#       # have to change test-names here too
#       quantilesTest = read.csv(file =fNameTest, header = TRUE)
#       
#       ### TODO! ^^^
#     } else {
#       quantilesTest = quantilesTrain
#     }
#     
#     print("Creating test-set ...")
#     Test = getSamplesSurf2(quantilesTest,sampleSize = sampleSize,sampleTimes = sampleTimes_test,euklid = euklid, numPermutations = numPermutations, numClasses = numClasses, m = m,reDo = reCalculateTrainTest)
#     
#     
#     if(global_flag){
#       #-------------------------------------------------------------------------------------------
#       # classlabels for the global part
#       #-------------------------------------------------------------------------------------------
#       print("adding global information to the model ...")
# 
#       quantilesTrain_global = read.csv(file =fNameTrain_global, header = TRUE)
#       
#       # quantilesTrain_global[,-1] = apply(quantilesTrain_global[,-1],2, FUN = function(i) as.numeric(as.character(i)))
# 
#       protClasses_global = unlist(lapply(c(1:length(protNames)), FUN = function(i){
#         mapping[which(mapping[,1] == lab$label[which(lab$name == protNames[i])]),2]
#       }))
# 
#       classLabels_global = unlist(lapply(c(1:nrow(quantilesTrain_global)),FUN = function(i){
#                                     protClasses_global[which(protNames == quantilesTrain_global[i,1])]
#                                   }))
# 
#       quantilesTrain_global[,1] = unlist(lapply(c(1:length(classLabels_global)), FUN = function(i){
#         paste(classLabels_global[i],"_",quantilesTrain_global[i,1], sep = "")
#       }))
# 
#       print("global Train ...")
#       TrainGlobal = getSamplesSurf2(quantilesTrain_global,sampleSize = 1,sampleTimes = 1,euklid = euklid, numPermutations = 1, numClasses = numClasses, m = 1,reDo = TRUE)
#       # merge the global and the local
#       # the global comes in front of each local-featureset
#       out = lapply(c(1:nrow(Train$X)), FUN = function(i){
#         ind = ceil(i/sampleTimes)
#         c(TrainGlobal$X[ind,],Train$X[i,])
#       })
# 
#       Train$X = matrix(unlist(out), byrow = TRUE, nrow = nrow(Train$X))
# 
#       # #------------------------------------------------------------------------------------
#       # # Test
#       # #------------------------------------------------------------------------------------
#       
#       TestGlobal = c()
#       if(fNameTrain_global != fNameTest_global){
#         # have to change test-names here too
#         quantilesTest_global = read.csv(file =fNameTest_global, header = TRUE)
#         ### TODO! ^^^
#       } else {
#         TestGlobal = TrainGlobal
#       }
#       
#       print("global Test ...")
#        # merge the global and the local
#       # the global comes in front of each local-featureset
#       out = lapply(c(1:nrow(Test$X)), FUN = function(i){
#         ind = ceil(i/sampleTimes_test)
#         c(TestGlobal$X[ind,],Test$X[i,])
#       })
# 
#     Test$X = matrix(unlist(out), byrow = TRUE, nrow = nrow(Test$X))
#     }
#     
#     shuf = shuffle(1:nrow(Train$X))
#     TrFinal = list("x_train" = Train$X[shuf,], "y_train" = Train$y[shuf,], "x_test" = Test$X, "y_test" = Test$y, "numClasses" = numClasses, "classLevels" = classLevels, "mapping" = mapping)
#     
#     saveRDS(TrFinal, ExperimentFile)
#   } else {
#     print("reading from previous experiment ...")
#     TrFinal = readRDS(ExperimentFile)
#   }
#   
#   TrFinal$x_train = apply(TrFinal$x_train,2, FUN = function(i) as.numeric(as.character(i)))
#   TrFinal$x_test = apply(TrFinal$x_test,2, FUN = function(i) as.numeric(as.character(i)))
#   
#   # TrFinal$y_train = apply(TrFinal$y_train,2, FUN = function(i) as.character(i))
#   # TrFinal$y_test = apply(TrFinal$y_test,2, FUN = function(i) as.character(i))
#   # return(TrFinal)
#   
#   fac = 1
#   if(euklid) fac=2
#   model = modelFUN(TrainTest = TrFinal,sampleSize = sampleSize,sampleTimes = sampleTimes,q = (q+2)*3*fac,epochs = epochs, batch_size = batch_size)
#   
#   
#   predictions <- model %>% predict_classes(TrFinal$x_test)
#   
#   # y_origNames = unlist(lapply(c(1:length(Train$y_original_names)), FUN = function(i){
#   #                                                   paste(strsplit(Train$y_original_names[i], "_")[[1]][-1], collapse = "_")
#   # }))
# 
#   
#   pred = predictions+1
#   print(pred)
#   gt = reverseToCategorical(TrFinal$y_test,TrFinal$classLevels)
#   y_test_pred = rep("0",length(pred))
#   su = 0
#   for(i in 1:length(gt)){
#     if(gt[i] == TrFinal$classLevels[pred[i]]) su = su + 1
#     
#     y_test_pred[i] = TrFinal$classLevels[pred[i]]
#   }
#   su/length(gt)
#   y_test_pred
# 
# 
#   confMat = table(factor(as.character(TrFinal$mapping[as.numeric(y_test_pred),1]),
#                          levels=TrFinal$mapping[TrFinal$classLevels,1]),
#                   factor(as.character(TrFinal$mapping[as.numeric(gt),1]),
#                          levels=TrFinal$mapping[TrFinal$classLevels,1]))
#   
#   
#   sum(confMat)
#   confMatNormalized = confMat/colSums(confMat)[col(confMat)]
#   
#   print(confMat)
#   print(confMatNormalized)
#   
#   accuracy = sum(diag(confMat)) / sum(confMat)
#   
#   print("-----------------------------")
#   print(paste("accuracy:", accuracy))
#   print("-----------------------------")
#   
#   # return(list("pred" = as.numeric(y_test_pred), "actual" = as.numeric(gt)))
#   # f1_score = ModelMetrics::f1Score(predicted = as.numeric(y_test_pred), actual = as.numeric(gt))
#   auc = ModelMetrics::auc(predicted = as.numeric(y_test_pred), actual = as.numeric(gt))
#   
#   
#   # f1 from the caret-package
#   pred = as.factor(as.numeric(y_test_pred))
#   act = as.factor(as.numeric(gt))
#   
#   newLevels  = unique(c(as.numeric(y_test_pred),as.numeric(gt)))
#   
#   # return(list("pred" = pred, "act" = act, "newLevels" = newLevels, "gt" = gt, "y_pred" =y_test_pred ))
#   
#   levels(pred) = newLevels
#   levels(act) = newLevels 
#   
#   
#   ?confusionMatrix
#   resu <- confusionMatrix(data = pred, reference = act, mode="prec_recall")
#   f1_score = resu$byClass["F1"]
#   
#   print(paste("f1_score: ", f1_score))
#   
#   write.table(x = signif(accuracy,2),file = paste("/home/willy/PredictingProteinInteractions/Results/TablesProt/Accuracy_", ExperimentName, ".tex", sep = ""),
#               quote = FALSE, col.names = FALSE, row.names = FALSE)
#   
#   print(xtable(x = confMat,caption = "Confusion-matrix 106 Redoxins ",label = "ModelNet10Conf", type = "latex"),
#         file = paste("/home/willy/PredictingProteinInteractions/Results/TablesProt/106TestConf_", ExperimentName, ".tex", sep = ""))
#   
#   
#   print(xtable(x = confMatNormalized,
#                caption = "Confusion-matrix 106 Redoxins (normalized)",
#                label = "106TestConfNormalized",
#                type = "latex"),
#         file = paste("/home/willy/PredictingProteinInteractions/Results/TablesProt/106TestConfNormalized_", ExperimentName, ".tex", sep = ""))
#   
#   
#   writeExperimentParametersToFile(pathToStats = "/home/willy/PredictingProteinInteractions/Results/TablesProt/",
#                                   sampleSize = sampleSize,
#                                   sampleTimes = sampleTimes,
#                                   sampleTimes_test = sampleTimes_test,
#                                   batch_size = batch_size,
#                                   epochs = epochs,
#                                   euklid = euklid,
#                                   q = q,
#                                   m = m, 
#                                   numClasses = numClasses,
#                                   potentials = potentials,
#                                   fNameTrain = fNameTrain,
#                                   fNameTest = fNameTest,
#                                   ExperimentName = ExperimentName,
#                                   fNameTrain_global = fNameTrain_global,
#                                   fNameTest_global = fNameTest_global,
#                                   modelName = modelName,
#                                   accuracy = accuracy,
#                                   f1_score = f1_score,
#                                   auc = auc)
#   
#   stats = joinStats()
#   write.csv(stats, "/home/willy/PredictingProteinInteractions/Results/ProtSummary.csv",row.names = FALSE)
# }
#-----------------

selectFeatures <- function(q, pos_flag, neg_flag, pos_neg_flag, euklid_geo, numOfFeatures = 1){
  #-------------------------------------------------------------------------------------------
  # which features to use
  #-------------------------------------------------------------------------------------------
  # euklid_geo, geo, euklid, both
  # numOfFeatures ... number of matrices that are cbinded
  
  quants = c(1:(q+2))
  # keep the name in any case
  colIndsTokeep = c(1)
  
  colIndsTokeep_tmp = c()
  if(euklid_geo == "both"){
    if(pos_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, quants+1, length(quants)*3+1+quants)
    if(neg_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, length(quants)*1+1+quants, length(quants)*4+1+quants)
    if(pos_neg_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, length(quants)*2+1+quants, length(quants)*5+1+quants)
  } else if(euklid_geo == "euklid"){
    if(pos_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, length(quants)*3+1+quants)
    if(neg_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, length(quants)*4+1+quants)
    if(pos_neg_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, length(quants)*5+1+quants)
  } else if(euklid_geo == "geo"){
    if(pos_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, quants+1)
    if(neg_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, length(quants)*1+1+quants)
    if(pos_neg_flag) colIndsTokeep_tmp = c(colIndsTokeep_tmp, length(quants)*2+1+quants)
  }
    
  for(i in 1:numOfFeatures){
    shift = (length(quants)*6)*(i-1)
    # colIndsTokeep_tmp = colIndsTokeep_tmp + shift
    colIndsTokeep = c(colIndsTokeep, colIndsTokeep_tmp + shift)
  }
  
  colIndsTokeep = sort(colIndsTokeep)
  return(colIndsTokeep)
}

# selectFeatures(q = 20, pos_flag = TRUE, neg_flag = FALSE, pos_neg_flag = FALSE,euklid_geo = "both", numOfFeatures = 4)
# 
# (20+2)*6*4



# example = read.csv("/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_1_betha_1_loc_TRUE.csv", header = TRUE)

# example[1:5,]
# 
# ncol(example)




createModelStatistics <- function(model, TrFinal, expDir, foldNum, testNames, sequentialModel = FALSE){
  
  pred = c()
  if(sequentialModel == TRUE){
    predictions <- model %>% predict_classes(TrFinal$x_test)
    pred = predictions+1
  } else {
    predictions <- model %>% predict(TrFinal$x_test)
    
    for(i in 1:nrow(predictions)){
      pred[i] = which.max(predictions[i,])
    }
  }

  # print(pred)
  gt = reverseToCategorical(TrFinal$y_test,TrFinal$classLevels)
  y_test_pred = rep("0",length(pred))
  su = 0
  for(i in 1:length(gt)){
    if(gt[i] == TrFinal$classLevels[pred[i]]) su = su + 1
    
    y_test_pred[i] = TrFinal$classLevels[pred[i]]
  }
  su/length(gt)
  y_test_pred
  
  # print(y_test_pred)
  # print(gt)
  
  
  confMat = table(factor(as.character(TrFinal$mapping[as.numeric(y_test_pred),1]),
                         levels=TrFinal$mapping[TrFinal$classLevels,1]),
                  factor(as.character(TrFinal$mapping[as.numeric(gt),1]),
                         levels=TrFinal$mapping[TrFinal$classLevels,1]))
  
  
  sum(confMat)
  confMatNormalized = confMat/colSums(confMat)[col(confMat)]
  
  print(confMat)
  print(confMatNormalized)
  
  accuracy = sum(diag(confMat)) / sum(confMat)
  
  print("-----------------------------")
  print(paste("accuracy:", accuracy))
  print("-----------------------------")
  
  # return(list("pred" = as.numeric(y_test_pred), "actual" = as.numeric(gt)))
  # f1_score = ModelMetrics::f1Score(predicted = as.numeric(y_test_pred), actual = as.numeric(gt))
  # auc = ModelMetrics::auc(predicted = as.numeric(y_test_pred), actual = as.numeric(gt))
  
  
  # f1 from the caret-package
  pred = as.factor(as.numeric(y_test_pred))
  act = as.factor(as.numeric(gt))
  newLevels  = unique(c(as.numeric(y_test_pred),as.numeric(gt)))
  levels(pred) = newLevels
  levels(act) = newLevels 
  
  resu <- confusionMatrix(data = pred, reference = act, mode="prec_recall")
  f1_score = resu$byClass["F1"]
  
  print(paste("f1_score: ", f1_score))
  
  write.table(confMat, file = paste(expDir, "/confMat_fold_", foldNum, ".txt", sep = ""))
  write.table(accuracy, file = paste(expDir, "/accuracy_fold_", foldNum, ".txt", sep = ""), row.names = FALSE)
  write.table(f1_score, file = paste(expDir, "/f1_score_", foldNum, ".txt", sep = ""), row.names = FALSE)
  write.table(testNames, file = paste(expDir, "/names_fold_", foldNum, ".txt", sep = ""), row.names = FALSE)
}


addGlobalInformationToModel <- function(protNames,fNameTrain_global, mapping, Train, lab, euklid){
  #-------------------------------------------------------------------------------------------
  # classlabels for the global part
  #-------------------------------------------------------------------------------------------
  print("adding global information to the model ...")
  
  quantilesTrain_global = read.csv(file =fNameTrain_global, header = TRUE)
  
  if(length(which(is.na(quantilesTrain_global) == TRUE)) > 0) {
    print(paste("NAs detected in ", fNameTrain_global, sep = ""))
    return()
  }
  
  # quantilesTrain_global[,-1] = apply(quantilesTrain_global[,-1],2, FUN = function(i) as.numeric(as.character(i)))
  
  protClasses_global = unlist(lapply(c(1:length(protNames)), FUN = function(i){
    mapping[which(mapping[,1] == lab$label[which(lab$name == protNames[i])]),2]
  }))
  
  classLabels_global = unlist(lapply(c(1:nrow(quantilesTrain_global)),FUN = function(i){
    protClasses_global[which(protNames == quantilesTrain_global[i,1])]
  }))
  
  quantilesTrain_global[,1] = unlist(lapply(c(1:length(classLabels_global)), FUN = function(i){
    paste(classLabels_global[i],"_",quantilesTrain_global[i,1], sep = "")
  }))
  
  print("global Train ...")
  TrainGlobal = getSamplesSurf2(quantilesTrain_global,sampleSize = 1,sampleTimes = 1,euklid = euklid, numPermutations = 1, numClasses = numClasses, m = 1,reDo = TRUE)
  # merge the global and the local
  # the global comes in front of each local-featureset
  out = lapply(c(1:nrow(Train$X)), FUN = function(i){
    ind = ceil(i/sampleTimes)
    c(TrainGlobal$X[ind,],Train$X[i,])
  })
  
  Train$X = matrix(unlist(out), byrow = TRUE, nrow = nrow(Train$X))

  return(Train$X)
}

getProtNameFromNameWithClassAsNumber <- function(y_original_names){
  y_original_names_out = rep("", length(y_original_names))
  for(i in 1:length(y_original_names)){
    v = strsplit(y_original_names[i],split = "_")[[1]][-1]
    y_original_names_out[i] = paste(v,collapse = "_")
  }
  
  return(y_original_names_out)
}


createPermutations2 <- function(X, names, y, numPermutations = 1, m = 100){
  #--------------------------------------------------------------------------------
  # for each row in X choose one element with the same name at random
  # this way different representations of the same model are recognized as the same.
  #
  # m ... number of consecutive rows that are from the same protein
  #
  #--------------------------------------------------------------------------------
  
  X_perm = matrix(0,ncol = ncol(X), nrow = numPermutations*nrow(X)/m)
  X_perm_out = matrix(0,ncol = ncol(X), nrow = numPermutations*nrow(X)/m)
  
  inds_perm = rep(0,numPermutations*nrow(X)/m)
  inds_perm_out = rep(0,numPermutations*nrow(X)/m)
  
  names_perm = rep(0,numPermutations*nrow(X)/m)
  y_perm = matrix(0,nrow = numPermutations*nrow(X)/m, ncol = ncol(y))
  
  for(i in 1:(nrow(X)/m)){
    
    # print(i/nrow(X)*m)
    
    start = ceiling(i/m)
    
    start3 = (i-1)*m+1
    end3 = start3+m-1
    
    inds = c(start3:end3)[sample(c(1:m),numPermutations, replace = TRUE)]
    
    # print(paste(start3,end3))
    
    inds2 = c(inds)[shuffle(numPermutations)]
    
    # 
    # print(inds)
    # print(inds2)
    
    startInd = (i-1)*numPermutations+1
    endInd = startInd + numPermutations-1
    
    X_perm[startInd:endInd,] = X[inds,]
    X_perm_out[startInd:endInd,] = X[inds2,]
    
    names_perm[startInd:endInd] = names[inds]
    y_perm[startInd:endInd,] = y[inds,]
  }
  
  # X_perm = X[inds_perm,] 
  # X_perm_out = X[inds_perm_out,] 
  
  return(list("X_perm" = X_perm, "X_perm_out" = X_perm_out, "originalNames" = names_perm, "y" = y_perm))
}


autoEncoder2  <- function(x_train, x_train_permute, epochs = 20, encoderDim = 3, unitNums = c(5,5,5), dropOuts = c(0.1,0.1,0.1), batchSize = 30){
  # set model
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = unitNums[1], activation = "relu", input_shape = ncol(x_train)) %>%
    layer_dropout(dropOuts[1]) %>%
    layer_dense(units = unitNums[2], activation = "relu") %>%
    layer_dropout(dropOuts[2]) %>%
    layer_dense(units = unitNums[3], activation = "relu") %>%
    layer_dropout(dropOuts[3]) %>%
    layer_dense(units = unitNums[4], activation = "relu") %>%
    layer_dropout(dropOuts[4]) %>%
    # layer_dense(units = unitNums[5], activation = "relu") %>%
    # layer_dropout(dropOuts[5]) %>%
    
    layer_dense(units = encoderDim, activation = "relu", name = "bottleneck") %>%
    
    # layer_dense(units = unitNums[5], activation = "relu") %>%
    layer_dense(units = unitNums[4], activation = "relu") %>%
    layer_dense(units = unitNums[3], activation = "relu") %>%
    layer_dense(units = unitNums[2], activation = "relu") %>%
    layer_dense(units = unitNums[1], activation = "relu") %>%
    layer_dense(units = ncol(x_train))
  
  
  
  # view model layers
  summary(model)
  
  
  # compile model
  model %>% compile(
    loss = "mean_squared_error", 
    optimizer = "adam"
  )
  
  # fit model
  model %>% fit(
    x = x_train, 
    y = x_train_permute,
    batch_size = batchSize,
    epochs = epochs,
    verbose = 1
  )
  
  return(model)
}


createPermutations <- function(X, numPermutations = 1, m = 100){
  #--------------------------------------------------------------------------------
  # for each row in X choose one element with the same name at random
  # this way different representations of the same model are recognized as the same.
  #
  # m ... number of consecutive rows that are from the same protein
  #
  #--------------------------------------------------------------------------------
  
  X_perm = matrix(0,ncol = ncol(X), nrow = numPermutations*nrow(X)/m)
  X_perm_out = matrix(0,ncol = ncol(X), nrow = numPermutations*nrow(X)/m)
  
  for(i in 1:nrow(X_perm)){
    
    # print(i/numPermutations)
    
    start = ceiling(i/numPermutations)
    
    start3 = (start-1)*m+1
    end3 = start3+m-1
    
    # print(paste(start3,end3))
    
    ind1 = sample(c(start3:end3),1)
    ind2 = sample(c(start3:end3),1)
    
    X_perm[i,] = X[ind1,]
    X_perm_out[i,] = X[ind2,]
  }
  
  
  return(list("X_perm" = X_perm, "X_perm_out" = X_perm_out))
}

model_built_from_pretrained <- function(GR,modelPreTrained, epochs = 30, batch_size = 64, weights, sampleSize = NULL){
  
  encoder <- keras_model(inputs = modelPreTrained$input, outputs = get_layer(modelPreTrained, "bottleneck")$output)
  # encoder <- keras_model_sequential(inputs = modelPreTrained$input, outputs = get_layer(modelPreTrained, "bottleneck")$output)
  
  # add our custom layers
  predictions <- encoder$output %>% 
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(0.2) %>%
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(0.2) %>% 
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(0.2) %>% 
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(0.2) %>% 
    layer_dense(units = 100, activation = 'relu') %>% 
    layer_dropout(0.2) %>% 
    layer_dense(units = 2, activation = 'softmax')
  
  # this is the model we will train
  fullModel <- keras_model(inputs = encoder$input, outputs = predictions)
  
  # first: train only the top layers (which were randomly initialized)
  freeze_weights(encoder)
  
  fullModel %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  x_train = GR$X_perm
  y_train = GR$y
  
  history <- fullModel %>% fit(
    x_train, y_train, 
    epochs = epochs, batch_size = batch_size,
    validation_split = 0.05)
  
  fullModel %>% summary()
  
  return(fullModel)
}


ProteinsExperimentKfoldCV <- function(sampleSize = 20,
                               sampleTimes = 10,
                               sampleTimes_test = 10,
                               batch_size = 1024,
                               epochs = 300,
                               euklid = "both",
                               q = 1,
                               m = 1000, 
                               numClasses = 2,
                               potentials = c("pos","neg","pos_neg"),
                               fNameTrain = c("/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_2_betha_1_loc_TRUE.csv"),
                               ExperimentName = "Test1",
                               fNameTrain_global = NULL,
                               fNameTest_global = NULL,
                               path = "/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/",
                               modelParameters,
                               recalculate = FALSE,
                               reCalculateTrainTest = FALSE,
                               labels = "/home/willy/PredictingProteinInteractions/data/labels.txt",
                               k = 10,
                               weights_input = NULL,
                               onlySummarizeFolds = FALSE,
                               normalizeInputs = TRUE,
                               saveExperiment = TRUE,
                               splitPattern = "",
                               useColIndsToKeep = TRUE,
                               doParallel = TRUE,
                               sort = TRUE,
                               pre_training = TRUE,
                               numPermutations = 10,
                               pathToModels = "/home/willy/PredictingProteinInteractions/data/106Test/Output/",
                               n_s_euclidean = 1000,
                               n_s_dijkstra = 1000,
                               stitchNum = 2000,
                               fps = TRUE,
                               sampleFunction = 2){
  
  print("------------------------------------------------------")
  print(paste("Experiment ", ExperimentName))
  print("------------------------------------------------------")
  
    
    pos_flag = "pos" %in% potentials
    neg_flag = "neg" %in% potentials
    pos_neg_flag = "pos_neg" %in% potentials
    
    global_flag = FALSE
    if(!is.null(fNameTrain_global) && !is.null(fNameTest_global)) global_flag = TRUE
    
    if(!dir.exists(path)) dir.create(path)
    
    expDir = paste(path,"/", ExperimentName, "/", sep = "")
    if(!dir.exists(expDir)) dir.create(expDir)
    
    ExperimentFile = paste(expDir,"/", ExperimentName,"sS_",sampleSize, "_sT_", sampleTimes, "_sTt_", sampleTimes_test, "_euklid_",
                           euklid,
                           "_pos_", pos_flag,
                           "_neg_", neg_flag,
                           "_pos_neg_", pos_neg_flag,
                           "_globalToo_", global_flag,
                           ".Rdata", sep ="")
    
    print(ExperimentFile)

    lastInd = 19
    df = data.frame(matrix(0,ncol = 2, nrow = lastInd+1 + 2*length(modelParameters$layers)))
    colnames(df) = c("parameter","value")
    df[1,] = c("sampleSize", sampleSize)
    df[2,] = c("sampleTimes", sampleTimes)
    df[3,] = c("sampleTimes_test", sampleTimes_test)
    df[4,] = c("batch_size", modelParameters$batch_size)
    df[5,] = c("epochs", modelParameters$epochs)
    df[6,] = c("euklid", euklid)
    df[7,] = c("q", q)
    df[8,] = c("m", m)
    df[9,] = c("numClasses", numClasses)
    df[10,] = c("pos_flag", pos_flag)
    df[11,] = c("neg_flag", neg_flag)
    df[12,] = c("pos_neg_flag", pos_neg_flag)
    df[13,] = c("fNameTrain_global", fNameTrain_global)
    df[14,] = c("ExperimentName", ExperimentName)
    df[15,] = c("modelName", "customModel")
    df[16,] = c("k", k)
    df[17,] = c("pre_training", pre_training)
    df[18,] = c("fps", fps)
    df[19,] = c("sampleFunction", sampleFunction)
    
    
    for(i in 1:length(fNameTrain)){
      df[lastInd+i,] = c(paste("Ft",i, sep ="_"), fNameTrain[i])
    }
    
    for(i in 1:length(modelParameters$layers)){
      df[lastInd+i + length(fNameTrain),] = c(paste("L",i, sep ="_"), modelParameters$layers[i])
    }
    
    for(i in 1:length(modelParameters$dropOuts)){
      df[lastInd+i+ length(modelParameters$layers) + length(fNameTrain),] = c(paste("d",i, sep ="_"), modelParameters$dropOuts[i])
    }
    

    print(expDir)
    write.table(df, paste(expDir,"/Call.txt", sep =""), row.names = FALSE)
    
    if(!onlySummarizeFolds){
      print("Setting up model ...")
    
    TrainTest = list()
    originalNames = c()
    protNames = c()
    classLevels = c()
    mapping = c()
    lab = read.table(labels, header = TRUE)
    if(!file.exists(ExperimentFile) || recalculate){
      print(paste("reading in ...", fNameTrain[1]))
      quantilesTrain = read.csv(file =fNameTrain[1], header = TRUE)
      
      # fNameTrain is a vector of files that need to be combined
      if(length(fNameTrain) > 1){
        for(i in 2:length(fNameTrain)){
          print(paste("reading in ...", fNameTrain[i]))
          tmp = read.csv(fNameTrain[i], header = TRUE)
          
          quantilesTrain = cbind(quantilesTrain, tmp[,-1])
        }
      }
      
      if(length(which(is.na(quantilesTrain) == TRUE)) > 0) {
        print(paste("NAs detected in ", fNameTrain, sep = ""))
        return()
      }
      
      quantilesTrain_origNames = quantilesTrain[,1]
      
      #-------------------------------------------------------------------------------------------
      # which features to use
      #-------------------------------------------------------------------------------------------
      if(useColIndsToKeep == TRUE){
        colIndsTokeep = selectFeatures(q, pos_flag, neg_flag, pos_neg_flag, euklid, numOfFeatures = length(fNameTrain))
        quantilesTrain = quantilesTrain[,colIndsTokeep]
      }
      
      #-------------------------------------------------------------------------------------------
      # classlabels
      #-------------------------------------------------------------------------------------------
      protNames = as.character(unique(quantilesTrain[,1]))
      
      # to circumnavigate problems in naming, ("_", is a problem) 
      mapping = data.frame(matrix(0,ncol = 2, nrow = length(unique(lab$label))))
      colnames(mapping) = c("originalName", "name")
      mapping[,1] = sort(unique(lab$label))
      mapping[,2] = as.numeric(as.factor(sort(unique(lab$label))))
      
      # remove all proteinNames that are not specified in
      quantilesTrain = quantilesTrain[which(quantilesTrain[,1] %in% lab$name),]
      protNames = as.character(unique(quantilesTrain[,1]))
      
      protClasses = unlist(lapply(c(1:length(protNames)), FUN = function(i){
        mapping[which(mapping[,1] == lab$label[which(lab$name == protNames[i])]),2]
      }))
      
      classLabels = unlist(lapply(c(1:nrow(quantilesTrain)), FUN = function(i){
        protClasses[which(protNames == quantilesTrain[i,1])]
      }))
      
      quantilesTrain[,1] = unlist(lapply(c(1:length(classLabels)), FUN = function(i){
        paste(classLabels[i],"_",quantilesTrain[i,1], sep = "")
      }))
      
      # print(quantilesTrain[1:5,1])
      
      # return(quantilesTrain)
      
      classLevels = mapping$name
      
      TrainTest = c()
      if(sampleFunction == 2){
        TrainTest = getSamplesSurf2(quantilesTrain, sampleSize = sampleSize,sampleTimes = sampleTimes,euklid = euklid, numPermutations = 1, numClasses = numClasses, m = m,reDo = reCalculateTrainTest, splitPattern = splitPattern, sort = TRUE)
      } else if(sampleFunction == 3){
        model_vec = c()
        if(fps == TRUE){
          print(paste("reading models from ", pathToModels, sep = ""))
          model_vec = getAllProteinModels(path = pathToModels,
                                          n_s_euclidean = n_s_euclidean,
                                          n_s_dijkstra = n_s_dijkstra,
                                          stitchNum = stitchNum,
                                          measureNearestNeighbors = 10,
                                          alpha = 0,
                                          betha = 0,
                                          onlyTheseIndices = NULL,
                                          recalculate = FALSE)
        }
        
        TrainTest = getSamplesSurf3(quantiles = quantilesTrain, model_vec = model_vec,sampleSize = sampleSize, sampleTimes = sampleTimes, numClasses = numClasses, m = m, splitPattern = splitPattern, fps = fps)
      }

    
      originalNames = getProtNameFromNameWithClassAsNumber(TrainTest$y_original_names)

      if(global_flag){
        TrainTest$X = addGlobalInformationToModel(protNames = protNames, fNameTrain_global = fNameTrain_global,mapping = mapping,Train = TrainTest,lab = lab,euklid = euklid)
      }
      
        if(saveExperiment == TRUE){
          saveRDS(list("TrainTest" = TrainTest, "originalNames" = originalNames, "protNames" = protNames, "classLevels" = classLevels, "mapping" = mapping), ExperimentFile)
          
        }
      } else {
        print("reading from previous experiment ...")
        TR = readRDS(ExperimentFile)
        
        TrainTest = TR$TrainTest
        originalNames = TR$originalNames
        protNames = TR$protNames
        classLevels = TR$classLevels
        mapping = TR$mapping
      }
    
      # return(list("mapping" = mapping, "originalNames" = originalNames))
    
      #------------------------------------------------------------------------------------------------------------------
      # K-fold-CV
      # create folds with the names
      #------------------------------------------------------------------------------------------------------------------
      folds = createFolds(lab$label, k = k,list = TRUE)
  
      # for(foldInd in 1:length(folds)){
      foreach(foldInd=1:length(folds)) %dopar% {
        
        print(paste("fold", foldInd, "/", k, sep =""))
        
        y_test_name_inds = folds[[foldInd]]
        
        test_inds = which(originalNames %in% protNames[y_test_name_inds])
        train_inds = c(1:nrow(TrainTest$X))[-test_inds]
        
        if(k == 1) train_inds = c(1:nrow(TrainTest$X))
        
        trainNames = protNames[-y_test_name_inds]
        testNames = protNames[y_test_name_inds]
        
        Train_X = TrainTest$X[train_inds,]
        Train_y = TrainTest$y[train_inds,]
        
        Test_X = TrainTest$X[test_inds,]
        Test_y = TrainTest$y[test_inds,]
        
        Train_X = apply(Train_X, 2, FUN = function(i) as.numeric(as.character(i)))
        Test_X = apply(Test_X, 2, FUN = function(i) as.numeric(as.character(i)))

        Train_y = apply(Train_y, 2, FUN = function(i) as.numeric(as.character(i)))
        Test_y = apply(Test_y, 2, FUN = function(i) as.numeric(as.character(i)))

        
        if(normalizeInputs){
          colMins = apply(Train_X,2,min)
          colMaxs = apply(Train_X,2,max)
          colRanges = colMaxs - colMins
          # Train_X = sapply(c(1:length(colRanges)), FUN = function(i){ Train_X[,i]/colRanges[i] })
          Train_X = sapply(c(1:length(colRanges)), FUN = function(i){ (Train_X[,i]- colMins[i])/colRanges[i]})
          Test_X = sapply(c(1:length(colRanges)), FUN = function(i){ (Test_X[,i]- colMins[i])/colRanges[i]})
          # Test_X = sapply(c(1:length(colRanges)), FUN = function(i){ Test_X[,i]/colRanges[i] })
        }
        
        weights = rep(1,length(classLevels))
        
        classLabels = reverseToCategorical(oneHot = Train_y, mapping$name)
        if(modelParameters$metrics == "accuracy"){
          for(i in 1:length(weights)){
            weights[i] = 1/(length(which(classLabels == classLevels[i]))/length(classLabels))
          }
        }
        
        if(!is.null(weights_input) && length(weights_input) == length(classLevels)){
          weights = weights_input * normalize(weights)
        }

        weights <- split(weights, mapping$name-1)
        print(weights)
        
        model <- keras_model_sequential()
        if(pre_training == TRUE){
          # first build an autoEncoder
          print("Building the auto-encoder ...")
          
          GR = createPermutations2(Train_X, names = originalNames,y = Train_y, numPermutations = numPermutations, m = sampleTimes)
          library(permute)
          shuf = shuffle(1:nrow(GR$X_perm))
          GR$X_perm = GR$X_perm[shuf,]
          GR$X_perm_out = GR$X_perm_out[shuf,]
          GR$originalNames = GR$originalNames[shuf]
          GR$y = GR$y[shuf,]
          
          modelPreTrained = autoEncoder2(GR$X_perm, GR$X_perm_out, epochs = 15,encoderDim = 100,dropOuts = c(0.0,0.0,0.0,0.0,0.0), unitNums = c(50,50,50,50),batchSize = modelParameters$batch_size*numPermutations)
          # model %>% save_model_hdf5(paste(outPath,"autoEncodeModel.h5", sep =""))
          # 
          # modelPreTrained <- load_model_hdf5(paste(outPath,"autoEncodeModel.h5", sep =""))
          # modelPreTrained %>% summary()
          
          print("Training the neural-net ...")
          model = model_built_from_pretrained(GR,modelPreTrained,weights = weights, batch_size = modelParameters$batch_size*(numPermutations/sampleTimes) ,epochs = modelParameters$epochs)

          # model %>% save_model_hdf5(paste(expDir,"/my_model.h5", sep = ""))
          
        } else {
          shuf = shuffle(1:nrow(Train_X))
          TrFinal = list("x_train" = Train_X[shuf,], "y_train" = Train_y[shuf,], "x_test" = Test_X, "y_test" = Test_y, "numClasses" = numClasses, "classLevels" = classLevels, "mapping" = mapping)
          # classLabels = reverseToCategorical(oneHot = TrFinal$y_train, mapping$name)
          
          
          model = modelProt_custom(TrainTest = TrFinal,
                                   weights = weights,
                                   layers = modelParameters$layers,
                                   dropOuts = modelParameters$dropOuts,
                                   optimizerFunName = modelParameters$optimizerFunName,
                                   metrics = modelParameters$metrics,
                                   batch_size = modelParameters$batch_size,
                                   epochs = modelParameters$epochs)
          
          # model %>% save_model_hdf5(paste(expDir,"/my_model.h5", sep = ""))
        }
        
        
        # return(model)

        # model %>% save_model_hdf5(paste(expDir,"/my_model.h5", sep = ""))
        # model <- load_model_hdf5(paste(expDir,"my_model.h5", sep =""))
        
        
        TrFinal = list("x_train" = Train_X, "y_train" = Train_y, "x_test" = Test_X, "y_test" = Test_y, "numClasses" = numClasses, "classLevels" = classLevels, "mapping" = mapping)
        
        createModelStatistics(model, TrFinal, expDir, foldInd, testNames = testNames)
      }
  }
  
  # now read in all confusion-matrices and average over the performance
  conf_all = read.table(paste(expDir, "/confMat_fold_", 1, ".txt", sep = ""))
  f1_all = read.table(paste(expDir,"/f1_score_",1,".txt", sep =""), header = TRUE)

  if(k > 1){
    for(i in 2:k){
      conf_i = read.table(paste(expDir, "/confMat_fold_", i, ".txt", sep = ""))
      conf_all = conf_all + conf_i
      
      
      f1_tmp = read.table(paste(expDir,"/f1_score_",i,".txt", sep =""), header = TRUE)
      
      if(is.na(f1_tmp)) f1_tmp = 0
      
      f1_all = f1_all + f1_tmp
    }
  }  

  
  conf_all = as.matrix(conf_all)
  
  accuracy = sum(diag(conf_all)) / sum(conf_all)
  accuracy
  
  
  confMatNormalized = conf_all/colSums(conf_all)[col(conf_all)]
  
  write.table(x = signif(accuracy,2),file = paste(expDir,"/Accuracy.tex", sep = ""),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  print(xtable(x = conf_all,caption = paste("Confusion-matrix from the 106-Redoxins-test-set.",sep =""),label = "Redoxins106TestConf", type = "latex"),
        file = paste(expDir,"Confusion.tex", sep = ""))
  
  print(xtable(x = confMatNormalized,caption = paste("Confusion-matrix from the 106-Redoxins-test-set.",sep =""),label = "Redoxins106TestConfNormalized", type = "latex"),
        file = paste(expDir,"ConfusionNormalized.tex", sep = ""))
  
  
  
  f1_all = f1_all/k
  
  write.table(x = signif(f1_all,2),file = paste(expDir,"/F1.tex", sep = ""),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}

# f1_all = read.table("/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test1/f1_score_1.txt", header = TRUE)
# f1_all
# 
# ProteinsExperimentKfoldCV( sampleSize = 20,
#                            sampleTimes = 200,
#                            sampleTimes_test = 10,
#                            batch_size = 32,
#                            epochs = 5,
#                            euklid = TRUE,
#                            q = 1,
#                            m = 1000,
#                            numClasses = 2,
#                            fNameTrain = "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                            fNameTrain_global = "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_1_m_1_q_1_muNN_10_alpha_0_betha_0_loc_FALSE.csv",
#                            ExperimentName = "Test1",
#                            modelName = getVarName(modelProt1),
#                            modelFUN = modelProt1,
#                            recalculate = FALSE,
#                            k = 10,
#                            onlySummarizeFolds = FALSE,
#                            normalizeInputs = TRUE,
#                            labels = LABELS,
#                            path = "/home/willy/PredictingProteinInteractions/data/120Experiment/NNexperimentsKfoldCV/",
#                            splitPattern = "_")
# 
# unique(tp[,1])
# 
# tp[2001,1]
# 
# read.table(LABELS, header = TRUE)


# modelParameters = list("layers" = c(50,10,10,10), "dropOuts" = c(0.5,0.5,0.5,0.5), "metrics" = "accuracy", "optimizerFunName" = "optimizer_adam", "batch_size" = 32, "epochs" = 20)
# modelParameters = list("layers" = c(500,100,100,10,10), "dropOuts" = c(0.5,0.4,0.4,0.1,0.1), "metrics" = "accuracy", "optimizerFunName" = "optimizer_adam", "batch_size" = 32, "epochs" = 20)
# 
# modelParameters = list("layers" = c(500,100,100,10,10), "dropOuts" = c(0.5,0.4,0.4,0.1,0.1), "metrics" = "accuracy", "optimizerFunName" = "optimizer_adam", "batch_size" = 32, "epochs" = 20)
# ProteinsExperimentKfoldCV( sampleSize = 5,
#                            sampleTimes = 200,
#                            sampleTimes_test = 10,
#                            batch_size = 32,
#                            epochs = 30,
#                            euklid = TRUE,
#                            q = 1,
#                            m = 1000,
#                            numClasses = 2,
#                            fNameTrain = c("/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                           "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                           "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.5_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                           "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv"),
#                            ExperimentName = "Test2",
#                            modelParameters = modelParameters,
#                            recalculate = FALSE,
#                            k = 1,
#                            onlySummarizeFolds = FALSE,
#                            normalizeInputs = TRUE,
#                            saveExperiment = TRUE,
#                            useColIndsToKeep = FALSE)


getExperimentSummary <- function(expDir = "/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/"){
  
  ExperimentFolders = list.dirs(expDir, recursive = FALSE)
  
  if(length(ExperimentFolders) == 0) return(NULL)
  # print(ExperimentFolders)
  
  
  df_summary = data.frame(matrix(0, nrow = length(ExperimentFolders), ncol = 31))
  colnames(df_summary) = c("accuracy",
                           "f1",
                           "alpha",
                           "betha",
                           "q",
                           
                           "sampleSize",
                           "sampleTimes",
                           "epochs",
                           "batch_size",
                           
                           "k",
                           "euklid",
                           "pos_flag",
                           "neg_flag",
                           "pos_neg_flag",
                           "modelName",
                           "ExperimentName",

                           "Ft_1",
                           "Ft_2",
                           "Ft_3",
                           "Ft_4",
                           "Ft_5",
                           
                           "L_1",
                           "L_2",
                           "L_3",
                           "L_4",
                           "L_5",
                           
                           "d_1",
                           "d_2",
                           "d_3",
                           "d_4",
                           "d_5")
  
  
  
  for(i in 1 :length(ExperimentFolders)){
    callFName = paste(ExperimentFolders[i], "/Call.txt", sep  = "")
    if(file.exists(callFName)){
      df = read.table(paste(ExperimentFolders[i], "/Call.txt", sep  = ""), header = TRUE)
      
      fNameTrain = as.character(df[which(df$parameter == "Ft_1"),2])
      alpha = as.numeric(strsplit(fNameTrain, "_")[[1]][11])
      betha = as.numeric(strsplit(fNameTrain, "_")[[1]][13])
      # n_local = as.numeric(strsplit(fNameTrain, "_")[[1]][3])
      q = as.numeric(strsplit(fNameTrain, "_")[[1]][7])
  
      accFile = paste(ExperimentFolders[i],"/Accuracy.tex", sep = "")
      Acc = NA
      if(file.exists(accFile)) Acc = as.numeric(read.table(accFile, header = FALSE))
     
      f1File = paste(ExperimentFolders[i],"/F1.tex", sep = "")
      F1 = NA
      if(file.exists(f1File)) F1 = as.numeric(read.table(f1File, header = FALSE))
  
      df_summary[i,1] = Acc
      df_summary[i,2] = F1
      df_summary[i,3] = alpha
      df_summary[i,4] = betha
      df_summary[i,5] = q
  
      df_summary[i,6] = as.numeric(as.character(df[which(df$parameter == "sampleSize"),2]))
      df_summary[i,7] = as.numeric(as.character(df[which(df$parameter == "sampleTimes"),2]))
      df_summary[i,8] = as.numeric(as.character(df[which(df$parameter == "epochs"),2]))
      df_summary[i,9] = as.numeric(as.character(df[which(df$parameter == "batch_size"),2]))
      
      df_summary[i,10] = as.numeric(as.character(df[which(df$parameter == "k"),2]))
      df_summary[i,11] = as.character(as.character(df[which(df$parameter == "euklid"),2]))
      df_summary[i,12] = as.character(df[which(df$parameter == "pos_flag"),2])
      df_summary[i,13] = as.character(df[which(df$parameter == "neg_flag"),2])
      df_summary[i,14] = as.character(df[which(df$parameter == "pos_neg_flag"),2])
      df_summary[i,15] = as.character(df[which(df$parameter == "modelName"),2])
      df_summary[i,16] = as.character(df[which(df$parameter == "ExperimentName"),2])
      
      lastInd = 16
     
      for(j in 1:5){
        ind =which(df$parameter == paste("Ft", j, sep = "_"))
        
        if(length(ind) == 0 ){
          df_summary[i,lastInd+j] = -1
        } else {
          df_summary[i,lastInd+j] = as.character(df[ind,2])
        }
      }
  
      for(j in 1:5){
        ind = which(df$parameter == paste("L", j, sep = "_"))
        
        if(length(ind) == 0 ){
          df_summary[i,lastInd+j+5] = -1
        } else {
          df_summary[i,lastInd+j+5] = as.numeric(as.character(df[ind,2]))
        }
  
      }
      
      for(j in 1:5){
        ind = which(df$parameter == paste("d", j, sep = "_"))
        
        if(length(ind) == 0 ){
          df_summary[i,lastInd+j+5+5] = -1
        } else {
          df_summary[i,lastInd+j+5+5] = as.numeric(as.character(df[ind,2]))
        }
      }
  
    }
  }
  
  df_summary = df_summary[order(df_summary$f1, decreasing = TRUE),]
  write.csv(df_summary, file = paste(expDir, "/summary.csv", sep =""), row.names = FALSE)
  
  return(df_summary)
}

# suma = getExperimentSummary()



checkIfParametersAreSame <- function(parameters1, parameters2){
  
  relevantParameters = intersect(names(parameters2),names(parameters1))
  
  
  for(i in 1:length(relevantParameters)){
    ind1 = which(names(parameters1) == relevantParameters[i])
    ind2 = which(names(parameters2) == relevantParameters[i])
    if(unlist(parameters1[ind1]) != unlist(parameters2[ind2])) return(FALSE)
  }
  
  return(TRUE)
}

getTrainFName <- function(path, name = "All", n = 0.2, m = 1, q = 1, muNN = 10, alpha = 0, betha = 0, local = TRUE){
  paste(path, "/",name, "_n_", n, "_m_", m, "_q_", q, "_muNN_", muNN, "_alpha_", alpha,"_betha_", betha, "_loc_", local,".csv", sep ="")
}


ExperimentWrapper <- function(parameters, pathKfold, labels, recalculateNAs){
  #--------------------------------------------------------------------------------
  # first check if the experiment with the needed parameters has already been done.
  # If not check if all the necessary files are already generated.
  # Then do the experiment.
  # parameters ... a list with the parameters
  #--------------------------------------------------------------------------------
  
  df_summary = getExperimentSummary(pathKfold)
  
  
  ExperimentName = "EMPTY"
  if(!is.null(df_summary)){
    # check if there is another experiment with the exact same parameters
    for(j in 1:nrow(df_summary)){
      dfList = as.list(df_summary[j,])
      flag = checkIfParametersAreSame(parameters,dfList)
      if(flag == TRUE){
        if(is.na(df_summary$accuracy[j]) && recalculateNAs == TRUE) {
          print(paste("found experiment with same parameters in ", df_summary$ExperimentName[j], ", but the Test was not finished. Will redo the Test.", sep = ""))
          ExperimentName = df_summary$ExperimentName[j]
          
          break
        }
        print(paste("found experiment with same parameters in ", df_summary$ExperimentName[j], " at position ", j," (l. ", j+1,") of ", nrow(df_summary) ,". Skipping this Experiment.", sep = ""))
        return(NULL)
      }
    }
  }
  
  if(ExperimentName == "EMPTY"){
    ExperimentName = getNextExperimentName()
  }
  
  print(ExperimentName) 
  
  
  # if we came until here, then there is no previous experiment with the same parameters
  
  # check if these files exist. If not skip this experiment
  fNameTrain = getTrainFName(path = parameters$path,
                             name = "All",
                             n = parameters$n_local,
                             m = 1,
                             q = parameters$q_local,
                             muNN = 10,
                             alpha = parameters$alpha_local,
                             betha = parameters$betha_local,
                             local = TRUE)
  
  if(!file.exists(fNameTrain)){
    tmp = getQuantilesAlphaBetha(alpha = parameters$alpha_local,
                                 betha = parameters$betha_local,
                                 n = parameters$n_local,
                                 m = 1,
                                 q = parameters$q_local,
                                 locale = TRUE,
                                 path = pathToExperiment,
                                 n_s_euclidean = 1000,
                                 n_s_dijkstra = 1000,
                                 stitchNum = 2000,
                                 measureNearestNeighbors = 10,
                                 recalculate = FALSE,
                                 recalculateQuants = FALSE)
  }
  
  potentials = c()
  if(parameters$pos_flag == TRUE) potentials = c(potentials, "pos")
  if(parameters$neg_flag == TRUE) potentials = c(potentials, "neg")
  if(parameters$pos_neg_flag == TRUE) potentials = c(potentials, "pos_neg")
  
  
  print("starting Experiment ...")
  ProteinsExperimentKfoldCV( sampleSize = parameters$sampleSize,
                             sampleTimes = parameters$sampleTimes,
                             sampleTimes_test = parameters$sampleTimes_test,
                             batch_size = parameters$batch_size,
                             epochs = parameters$epochs,
                             euklid = parameters$euklid,
                             q = parameters$q_local,
                             m = 1000,
                             numClasses = NUMCLASSES,
                             fNameTrain = fNameTrain,
                             fNameTrain_global = fNameTrain_global,
                             ExperimentName = ExperimentName,
                             modelName = parameters$modelName,
                             modelFUN = parameters$modelFun,
                             recalculate = TRUE,
                             potentials = potentials,
                             k = parameters$k,
                             onlySummarizeFolds = FALSE,
                             normalizeInputs = TRUE,
                             saveExperiment = SAVE_EXPERIMENTS,
                             path = pathKfold,
                             labels = labels,
                             reCalculateTrainTest = FALSE)
  
  beep(1)
}


gridSearch <- function( alphas_local = c(0,1,2,3,5),
                        bethas_local = c(0,1,2,3,5),
                        alphas_global = c(0,1,2,3,5),
                        bethas_global = c(0,1,2,3,5),
                        # models = c(modelProt3_f1, modelProt3_f1, modelProt1)
                        sampleSizes = c(5,10,20),
                        sampleTimes = c(200, 400),
                        batch_sizes = c(32,16),
                        epochs = c(20),
                        nlocals = c(0.2, 0.5, 0.1),
                        q_locals = c(1),
                        q_globals = c(1),
                        euklid_val = TRUE,
                        recalculateNAs = TRUE,
                        k = 10){
  for(nloc in nlocals){
    for(alpha_loc in alphas_local){
      for(betha_loc in bethas_local){
        for(alpha_global in alphas_global){
          for(betha_global in bethas_global){
            for(sampleSize in sampleSizes){
              for(sampleTime in sampleTimes){
                for(batch_size in batch_sizes){
                  for(epoch in epochs){
                    for(q_local in q_locals){
                      for(q_global in q_globals){
                        parameters = list("alpha_local" = alpha_loc,
                                          "betha_local" = betha_loc,
                                          "alpha_global" = alpha_global,
                                          "betha_global" = betha_global,
                                          "sampleSize" = sampleSize,
                                          "sampleTimes" = sampleTime,
                                          "sampleTimes_test" = 200,
                                          "batch_size" = batch_size,
                                          "epochs" = epoch,
                                          "euklid" = euklid_val,
                                          "q_local" = q_local,
                                          "q_global" = q_global,
                                          "k" = k,
                                          "pos_flag" = TRUE,
                                          "neg_flag" = TRUE,
                                          "pos_neg_flag" = TRUE,
                                          "modelName" = getVarName(modelProt3_f1),
                                          "modelFun" = modelProt3_f1,
                                          "n_local" = nloc,
                                          "path" = paste(p2,"/Quantiles/", sep = ""))
                        
                        
                        ExperimentWrapper(parameters, NNexperimentsKfoldDir, labels = LABELS, recalculateNAs)
                        
                        
                        if(WS_flag == TRUE) system("/home/sysgen/Documents/LWB/Uploader/Uploader.sh")
                        
                      } 
                    } 
                  } 
                } 
              }
            }
          } 
        } 
      }
    }
  }
}


generateRandomModelParameters <- function(maxLayerNum = 5,
                                          layerSizes = c(5,10,50,20,30,100,500),
                                          dropOuts = c(0.1,0.2,0.7,0.05,0.5,0.3,0.2,0.4),
                                          optimizerFunNames = c("optimizer_adam", "optimizer_rmsprop"),
                                          batch_sizes = c(32,16),
                                          epochs = c(20)){
  
  modelSize = sample(maxLayerNum, 1)
  # layers = rep(0,length(modelSize))
  # layers[1] = chooseOneRandom(layerSizes)
  # 
  # print(layers[1])
  # print(modelSize)
  # if(modelSize > 1){
  #   for(i in 2:length(layers)){
  #     l = chooseOneRandom(layerSizes)
  #     print(l)
  #     print("layer")
  #     print(layers[i-1])
  #     while(l > layers[i-1]){
  #       l = chooseOneRandom(layerSizes)
  #     }
  #     layers[i] = l
  #   }
  # }
  
  layers = sort(sample(layerSizes, modelSize, replace = TRUE),decreasing = TRUE)
  dropOuts = sort(sample(dropOuts, size = modelSize, replace = TRUE), decreasing = TRUE)
  opti = sample(optimizerFunNames, 1)
  
  
  modelParameters = list("layers" = layers, 
                         "dropOuts" = dropOuts, 
                         "metrics" = "accuracy",
                         "optimizerFunName" = opti,
                         "batch_size" = batch_sizes[sample(1:length(batch_sizes), 1)],
                         "epochs" = epochs[sample(1:length(epochs), 1)])
  
  
  return(modelParameters)
}



randomSearch <- function( NNexperimentsKfoldDir = "/home/willy/PredictingProteinInteractions//data/106Test/NNexperimentsKfoldCV/",
                          MAXit = 1,
                          sampleSizes = c(5,10,20),
                          sampleTimes = c(200, 400),
                          alphas = c(3),
                          bethas = c(3),
                          qs = c(1),
                          nlocals = c(0.1,0.2,0.5,0.8),
                          muNN = 10,
                          euklid_val = TRUE,
                          recalculateNAs = TRUE,
                          k = 10,
                          layerSizes = c(5,10,50,20,30,100,500),
                          dropOuts = c(0.1,0.2,0.7,0.05,0.5),
                          metrics = c("accuracy"),
                          optimizerFunNames = c("optimizer_adam", "optimizer_rmsprop"),
                          batch_sizes = c(32,16),
                          epochs = c(20),
                          labels = LABELS){
  
  print("Starting random search ...")
  
  iteration = 0
  while(MAXit > iteration){
    iteration = iteration + 1
    
    df_summary = getExperimentSummary(NNexperimentsKfoldDir)
    startNum = nrow(df_summary)
    
    if(is.null(startNum)) startNum = 0
    
    testName = paste("Test", startNum + 1, sep = "")
    
    modelParameters = generateRandomModelParameters()
    
    sampleTimes = sampleTimes[sample(1:length(sampleTimes), 1)]
    sampleSize = sampleSizes[sample(1:length(sampleSizes), 1)]
    
    
    alpha = alphas[sample(1:length(alphas), 1)]
    betha = bethas[sample(1:length(bethas), 1)]
    q = qs[sample(1:length(qs), 1)]
    
    fNameTrain = c()
    flag = TRUE
    for(i in 1:length(nlocals)){
      p2 = strsplit(pathToExperiment, "/Output/")[[1]][1]
      
      fNameTmp = getTrainFName(path = paste(p2, "/Quantiles/", sep = ""),name = "All",n = nlocals[i],m = 1,q = q,muNN = muNN,alpha = alpha,betha = betha,local = TRUE)
      if(!file.exists(fNameTmp)) {
        print(paste("generating missing file ...", fNameTmp))
        # tmp = getQuantilesAlphaBetha(alpha = alpha,betha = betha, n = nlocals[i], m = 1,q = q, locale = TRUE, path = pathToExperiment, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = muNN, recalculate = FALSE,recalculateQuants = TRUE)

        flag = FALSE
      }
      
      fNameTrain = c(fNameTrain, fNameTmp)
    }
    
    if(flag == TRUE) {
      ProteinsExperimentKfoldCV( sampleSize = sampleSize,
                                 sampleTimes = sampleTimes,
                                 sampleTimes_test = sampleTimes,
                                 batch_size = 32,
                                 epochs = 30,
                                 euklid = TRUE,
                                 q = 1,
                                 m = 1000,
                                 numClasses = 2,
                                 fNameTrain = fNameTrain,
                                 ExperimentName = testName,
                                 modelParameters = modelParameters,
                                 recalculate = FALSE,
                                 k = k,
                                 onlySummarizeFolds = FALSE,
                                 normalizeInputs = TRUE,
                                 saveExperiment = FALSE,
                                 useColIndsToKeep = FALSE,
                                 labels = labels,
                                 path = NNexperimentsKfoldDir)
      
      df_summary = getExperimentSummary(NNexperimentsKfoldDir)
      beep(5)
      if(WS_flag == TRUE) system("/home/sysgen/Documents/LWB/Uploader/Uploader.sh")
    }
    

  }
  

}


#------------------------------------------------------------------------
# p2 = strsplit(pathToExperiment, "/Output/")[[1]][1]
# 
# print(p2)
# NNexperimentsKfoldDir = paste(p2, "/NNexperimentsKfoldCV/", sep = "")
# 
# randomSearch(NNexperimentsKfoldDir = NNexperimentsKfoldDir,
#              k = 10,
#              nlocals = c(0.1,0.2,0.5,0.8),
#              alphas = c(3),
#              bethas = c(3),
#              MAXit = 1)



# df_summary = getExperimentSummary()
# 
# modelParameters = list("layers" = c(50,10,10,10), "dropOuts" = c(0.1,0.1,0.1,0.1), "metrics" = "accuracy", "optimizerFunName" = "optimizer_adam", "batch_size" = 32, "epochs" = 20)
# ProteinsExperimentKfoldCV( sampleSize = 5,
#                            sampleTimes = 400,
#                            sampleTimes_test = 10,
#                            batch_size = 32,
#                            epochs = 30,
#                            euklid = TRUE,
#                            q = 1,
#                            m = 1000,
#                            numClasses = 2,
#                            fNameTrain = c("/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                           "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                           "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.5_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                           "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv"),
#                            ExperimentName = "Test10",
#                            modelParameters = modelParameters,
#                            recalculate = FALSE,
#                            k = 10,
#                            onlySummarizeFolds = FALSE,
#                            normalizeInputs = TRUE,
#                            saveExperiment = FALSE,
#                            useColIndsToKeep = FALSE)
# 
# 
# df_summary = getExperimentSummary()

chooseOneRandom <- function(v){
  return(v[sample(1:length(v), 1)])
}

randomSearch2 <- function(modelParameters, NNexperimentsKfoldDir, ExperimentName, quantilePath,
                          alphas = c(0,1,2,3),
                          bethas = c(0,1,2,3),
                          nlocals = c(0.1,0.2,0.3,0.5,0.8),
                          qs = c(1),
                          sampleSizes = c(5,10,20,40),
                          sampleTimes = c(200, 400, 800),
                          Trials = 1){
  
  for(trial in 1:Trials){
    
    df_summary = getExperimentSummary(expDir = NNexperimentsKfoldDir)
    
    v = as.numeric(unlist(strsplit(df_summary$ExperimentName, split = "Test")))
    nextV = max(v[!is.na(v)]) + 1
    ExperimentName = paste("Test", nextV, sep = "")
    
    
    # choose numbter of feature-matrices at random
    fNameTrainNum = sample(5, size = 1)
    
    nlocals_samp = nlocals[sample(c(1:length(nlocals)), size = fNameTrainNum)]
    
    fNameTrain = c()
    for(i in 1:length(nlocals_samp)){
      alpha = chooseOneRandom(alphas)
      betha = chooseOneRandom(bethas)
      q = chooseOneRandom(qs)
      fName = getTrainFName(path = quantilePath, name = "All", n = nlocals_samp[i], m = 1, q = q,muNN = 10,alpha = alpha, betha = betha,local = TRUE)
      fNameTrain = c(fNameTrain, fName)
      
      if(!file.exists(fName)){
        
        registerDoParallel(16)
        tmp = getQuantilesAlphaBetha(alpha = alpha,betha = betha, n = nlocals_samp[i], m = 1,q = q, locale = TRUE, path = pathToExperiment, n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000, measureNearestNeighbors = 10, recalculate = FALSE,recalculateQuants = FALSE)
        registerDoParallel(numCores)
      }
    }
    
    modelParameters = generateRandomModelParameters()

    ProteinsExperimentKfoldCV( sampleSize = chooseOneRandom(sampleSizes),
                               sampleTimes = chooseOneRandom(sampleTimes),
                               sampleTimes_test = 10,
                               batch_size = 32,
                               epochs = 30,
                               euklid = TRUE,
                               q = q,
                               m = 1000,
                               numClasses = 2,
                               fNameTrain = fNameTrain,
                               ExperimentName = ExperimentName,
                               modelParameters = modelParameters,
                               recalculate = FALSE,
                               k = 10,
                               onlySummarizeFolds = FALSE,
                               normalizeInputs = TRUE,
                               saveExperiment = FALSE,
                               useColIndsToKeep = FALSE,
                               path = NNexperimentsKfoldDir,
                               labels = LABELS)

    df_summary = getExperimentSummary(expDir = NNexperimentsKfoldDir)
    
    if(WS_flag == TRUE) system("/home/sysgen/Documents/LWB/Uploader/Uploader.sh")
  }
}

# mode = "onlyExperiments2"
if(mode == "onlyExperiments2"){
  p2 = strsplit(pathToExperiment, "/Output/")[[1]][1]
  print(p2)
  NNexperimentsKfoldDir = paste(p2, "/NNexperimentsKfoldCV/", sep = "")
  
  df_summary = getExperimentSummary(expDir = NNexperimentsKfoldDir)
  
  
  # 0.05,0.1,0.2,0.3,0.5,0.8)
  
  
  # "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
  # "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
  # "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.05_m_1_q_1_muNN_10_alpha_0_betha_0_loc_TRUE.csv",
  # "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.5_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
  # "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv"

  
  conv = FALSE
  SAMPLESIZE = 20
  modelParameters = list("layers" = c(100,100,100,50,30), "dropOuts" = c(0.2,0.2,0.2,0.2,0.2), "metrics" = "accuracy", "optimizerFunName" = "optimizer_adam", "batch_size" = 16, "epochs" = 20)
  ProteinsExperimentKfoldCV( sampleSize = SAMPLESIZE,
                             sampleTimes = 100,
                             sampleTimes_test = 10,
                             batch_size = 32,
                             epochs = 30,
                             euklid = "both",
                             potentials = c("pos", "neg", "pos_neg"),
                             q = 1,
                             m = 1000,
                             numClasses = 2,
                             fNameTrain = c(  "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.05_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
                                              "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
                                              "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
                                              "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.5_m_1_q_1_muNN_10_alpha_1_betha_1_loc_TRUE.csv",
                                              "/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_1_betha_1_loc_TRUE.csv"),
                                            ExperimentName = "Test120",
                             modelParameters = modelParameters,
                             recalculate = FALSE,
                             k = 10,
                             onlySummarizeFolds = FALSE,
                             normalizeInputs = TRUE,
                             saveExperiment = TRUE,
                             useColIndsToKeep = TRUE,
                             path = NNexperimentsKfoldDir,
                             labels = getPath("106ExperimentLabels"),
                             pre_training = TRUE,
                             numPermutations = 400,
                             fps = TRUE,
                             sampleFunction = 3
                             )

  
  df_summary = getExperimentSummary(expDir = NNexperimentsKfoldDir)
}



if(mode == "onlyExperiments3"){
  p2 = strsplit(pathToExperiment, "/Output/")[[1]][1]
  print(p2)
  NNexperimentsKfoldDir = paste(p2, "/NNexperimentsKfoldCV/", sep = "")
  
  df_summary = getExperimentSummary(expDir = NNexperimentsKfoldDir)
  
  
  # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
  # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
  # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.05_m_1_q_1_muNN_10_alpha_0_betha_0_loc_TRUE.csv",
  # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.5_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
  # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv"
  
  
  conv = FALSE
  SAMPLESIZE = 20
  modelParameters = list("layers" = c(500,100,50,30), "dropOuts" = c(0.2,0.1,0.1,0.1), "metrics" = "accuracy", "optimizerFunName" = "optimizer_adam", "batch_size" = 32, "epochs" = 20)
  ProteinsExperimentKfoldCV( sampleSize = SAMPLESIZE,
                             sampleTimes = 400,
                             sampleTimes_test = 10,
                             batch_size = 32,
                             epochs = 30,
                             euklid = "both",
                             potentials = c("pos", "neg", "pos_neg"),
                             q = 1,
                             m = 1000,
                             numClasses = 2,
                             fNameTrain = c(    "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
                                                "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
                                                "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.05_m_1_q_1_muNN_10_alpha_0_betha_0_loc_TRUE.csv",
                                                "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.5_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
                                                "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv"),
                             ExperimentName = "Test102",
                             modelParameters = modelParameters,
                             recalculate = FALSE,
                             k = 10,
                             onlySummarizeFolds = FALSE,
                             normalizeInputs = TRUE,
                             saveExperiment = TRUE,
                             useColIndsToKeep = TRUE,
                             path = NNexperimentsKfoldDir,
                             labels = getPath("120ExperimentLabels"),
                             pre_training = FALSE,
                             numPermutations = 200,
                             fps = FALSE,
                             sampleFunction = 2
  )
  
  
  df_summary = getExperimentSummary(expDir = NNexperimentsKfoldDir)
}


if(mode == "randomSearch"){
  p2 = strsplit(pathToExperiment, "/Output/")[[1]][1]
  print(p2)
  NNexperimentsKfoldDir = paste(p2, "/NNexperimentsKfoldCV/", sep = "")


  quantilePath = paste(p2, "/Quantiles/", sep = "")
  randomSearch2(modelParameters = modelParameters,
               NNexperimentsKfoldDir = NNexperimentsKfoldDir,
               ExperimentName = ExperimentName,
               quantilePath = quantilePath,
               nlocals = c(0.05,0.1,0.2,0.3,0.8,0.5),
               Trials = 100000)
}


# fNames = list.files("/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/106Test/Quantiles/", recursive = FALSE, pattern = "TRUE")
# 
# nlocVals = c()
# for(i in 1:length(fNames)){
#   nlocVal = as.numeric(strsplit(fNames[i], "_")[[1]][3])
#   nlocVals = c(nlocVals, nlocVal)
# }
# 
# 6*9
# length(nlocVals)

#------------------------------------------------------------------------



# 
# sampleTimes_test = 10
# 
# if(mode == "onlyExperiments" || mode == "both"){
#   
#   p2 = strsplit(pathToExperiment, "/Output/")[[1]][1]
#   
#   
#   NNexperimentsKfoldDir = paste(p2, "/NNexperimentsKfoldCV/", sep = "")
#   if(!dir.exists(NNexperimentsKfoldDir)) dir.create(NNexperimentsKfoldDir)
#   
#   tests = list.dirs(NNexperimentsKfoldDir,recursive = FALSE,full.names = FALSE)
#   
#   print(paste("starting Experiments in", NNexperimentsKfoldDir, "...", sep = " "))
#   
#   EXPERIMENTCOUNT = 0
#   if(length(tests) != 0) {
#     testInds = as.numeric(unlist(strsplit(tests, "Test")))
#     testInds = testInds[!is.na(testInds)]
#     
#     EXPERIMENTCOUNT = max(testInds)
#   }
#   
#   print(paste("starting with ",EXPERIMENTCOUNT, sep = ""))
# 
#   getNextExperimentName <- function(){
#     EXPERIMENTCOUNT <<- EXPERIMENTCOUNT+1
#     
#     return(paste("Test", EXPERIMENTCOUNT, sep = ""))
#   }
#   
#   #------------------------------------------------------------------------
# 
#   alphas_local = c(0,1,2,3,5)
#   bethas_local = c(0,1,2,3,5)
#   
#   alphas_global = c(0,1,2,3,5)
#   bethas_global = c(0,1,2,3,5)
#   
#   # models = c(modelProt3_f1, modelProt3_f1, modelProt1)
#   
#   sampleSizes = c(5,10,20)
#   sampleTimes = c(200, 400)
#   batch_sizes = c(32,16)
#   epochs = c(20)
#   
#   nlocals = c(0.2, 0.5, 0.1)
#   
#   q_locals = c(1)
#   q_globals = c(1)
#   
#   euklid_val = TRUE
#   
#   recalculateNAs = TRUE
#   
#   k = 10
#   
#   
#  
#   
# 
# 
# 
#   summary = getExperimentSummary(NNexperimentsKfoldDir)
# }
# 
# 
# beep(5)

# # stats = joinStats()
# # write.csv(stats, "/home/willy/PredictingProteinInteractions/Results/ProtSummary.csv",row.names = FALSE)
# 
# 
# #-------------------------------------------------------------------------------------------------------------
# # experimental
# #-------------------------------------------------------------------------------------------------------------
# 
# # quit()
# 
# 
# quantiles08 = getQuantilesAlphaBetha(alpha = 1,betha = 1,path = pathToExperiment,
#                                     n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000,
#                                     measureNearestNeighbors = 10, recalculate = FALSE,recalculateQuants = FALSE,
#                                     n = 0.8, m = 1,q = 1, locale = TRUE)
# 
# quantiles02 = getQuantilesAlphaBetha(alpha = 3,betha = 3,path = pathToExperiment,
#                                      n_s_euclidean = 1000,n_s_dijkstra = 1000,stitchNum = 2000,
#                                      measureNearestNeighbors = 10, recalculate = FALSE,recalculateQuants = FALSE,
#                                      n = 0.2, m = 1,q = 1, locale = TRUE)
# 
# plotProteinModel(fNameOrigingal = "/home/willy/PredictingProteinInteractions/data/106Test/Output/000_Trx/000_Trx.obj", lis = models_small[[7]])
# 
# functionals = c(getFunctionalProteins(), "000_Trx")
# plot_prot_quants(quantiles1, q=1, functionals, plotMode = "Booth", withEuclid = TRUE)
# 
# 
# 
# 
# 
# inds = which(quantiles_trafo[,1] == "000_Trx")
# plot_prot_quants(quantiles_trafo[inds,], q=1, functionals, plotMode = "Booth", withEuclid = TRUE)
# 
# inds = which(quantiles_trafo[,1] == "022")
# plot_prot_quants(quantiles_trafo[inds,], q=1, functionals, plotMode = "Booth", withEuclid = TRUE)
# 
# # install.packages("VoxR")
# library(VoxR)
# 
# voxelizeData  <- function(quantiles, res = 0.01){
#   # only works with q = 1
#   
#   vox_pos = vox(quantiles[,2:4],res = res)
#   vox_neg = vox(quantiles[,5:7],res = res)
#   vox_pos_neg = vox(quantiles[,8:10],res = res)
#   
#   vox_pos_euc = vox(quantiles[,11:13],res = res)
#   vox_neg_euc = vox(quantiles[,14:16],res = res)
#   vox_pos_neg_euc = vox(quantiles[,17:19],res = res)
#   
#   return(list("vox_pos" = vox_pos,
#               "vox_neg" = vox_neg,
#               "vox_pos_neg" = vox_pos_neg,
#               
#               "vox_pos_euc" = vox_pos_euc,
#               "vox_neg_euc" = vox_neg_euc,
#               "vox_pos_neg_euc" = vox_pos_neg_euc))
# }
# 
# plotvoxelizeData <- function(voxelizedData, size = 10, newPlot = TRUE){
#   if(newPlot) {
#     open3d()
#   }
#       
#   cols = c("red", "blue","green", "brown", "pink", "lightblue")
#   plot3d(voxelizedData$vox_pos, size = size, col = cols[1], add = !newPlot)
#   plot3d(voxelizedData$vox_neg, size = size, col = cols[2], add = TRUE)
#   plot3d(voxelizedData$vox_pos_neg, size = size, col = cols[3], add = TRUE)
#   
#   plot3d(voxelizedData$vox_pos_euc, size = size, col = cols[4], add = TRUE)
#   plot3d(voxelizedData$vox_neg_euc, size = size, col = cols[5], add = TRUE)
#   plot3d(voxelizedData$vox_pos_neg_euc, size = size, col = cols[6], add = TRUE)
# }
# 
# getTensor <- function(vdata, m = 1000){
#   ar <- array(0, c(resNum,resNum, resNum))
#   for(i in 1:nrow(vdata)){
#     x = vdata[i,1]*resNum
#     y = vdata[i,2]*resNum
#     z = vdata[i,3]*resNum
#     
#     ar[x,y,z] = vdata[i,4]/1000
#   }
#   return(ar)
# }
# 
# pcaTransform <- function(quantiles1){
#   PCA = prcomp(quantiles1[,-1],retx = TRUE,scale. = FALSE)
#   quantiles1_pca = PCA$x
#   
#   Train_X = quantiles1_pca
#   colMins = apply(Train_X,2,min)
#   colMaxs = apply(Train_X,2,max)
#   colRanges = colMaxs - colMins
#   Train_X = sapply(c(1:length(colRanges)), FUN = function(i){ (Train_X[,i]- colMins[i])/colRanges[i]})
#   apply(Train_X,2,max)
#   apply(Train_X,2,min)
#   
#   quantiles_trafo = quantiles1
#   quantiles_trafo[,-1] = Train_X
#   
#   return(quantiles_trafo)
# }
# 
# 
# get3dRepresentationFromProtein <- function(quantiles, res = 0.1, features = c(1,2)){
#   voxelizedData = voxelizeData(quantiles,res = res)
#   
#   # tensorlist = list()
#   # tensorlist = getTensor(voxelizedData$vox_pos)
#   # tensor_pos_euc = getTensor(voxelizedData$vox_pos_euc)
#   
#   x_train <- array(0, c(1,resNum,resNum, resNum,length(features)))
#   
#   # first <- this protein, 
#   # last <- chanel
#   for(i in 1:length(features)){
#     x_train[1,,,,i] = getTensor(matrix(unlist(voxelizedData[features[i]]), ncol = 4, byrow = FALSE))
#   }
#   
#   return(x_train)
# }
# 
# 
# Times = 20
# resNum = 5
# features = c(1,2,3,4,5,6)
# 
# quantiles_trafo = pcaTransform(quantiles08)
# 
# quantiles_trafo2 = pcaTransform(quantiles02)
# 
# inds = which(quantiles_trafo[,1] == "016")
# voxelizedData2 = voxelizeData(quantiles_trafo[inds,],1/resNum)
# plotvoxelizeData(voxelizedData2,newPlot = TRUE)
# 
# inds = which(quantiles_trafo2[,1] == "016")
# voxelizedData2 = voxelizeData(quantiles_trafo2[inds,],1/resNum)
# plotvoxelizeData(voxelizedData2,newPlot = TRUE)
# 
# 
# 
# protNames = as.character(unique(quantiles_trafo[,1]))
# X <- array(dim = c(length(protNames)*Times,resNum,resNum, resNum,length(features)))
# Y = matrix(0, ncol = 2, nrow = length(protNames)*Times)
# Y_orig_names = rep("", length(protNames)*Times)
# for(i in 1:length(protNames)){
#   print(protNames[i])
#   inds = which(quantiles_trafo[,1] == protNames[i])
#   X1 = get3dRepresentationFromProtein(quantiles_trafo[inds,], res = 1/resNum, features = features)
#   
#   ind_start = (i-1)*Times+1
#   ind_end = ind_start + Times -1
#   
#   for(j in ind_start:ind_end){
#     X[j,,,,] = X1
#     Y[j,] <- to_categorical(as.numeric((protNames[i] %in% functionals)), 2)
#     Y_orig_names[j] =   protNames[i]
#   }
# 
# }
# 
# 
# 
# lab = read.table(LABELS, header = TRUE)
# folds = createFolds(lab$label, k = 10,list = TRUE)
# foldInd = 1
# 
# y_test_name_inds = c(folds[[1]], folds[[2]], folds[[3]])
# 
# test_inds = which(Y_orig_names %in% protNames[y_test_name_inds])
# train_inds = c(1:nrow(X))[-test_inds]
# 
# if(k == 1) train_inds = c(1:nrow(X))
# 
# trainNames = protNames[-y_test_name_inds]
# testNames = protNames[y_test_name_inds]
# 
# x_test = X[test_inds,,,,]
# y_test = Y[test_inds,]
# 
# x_train= X[train_inds,,,,]
# y_train = Y[train_inds,]
# 
# shuf = shuffle(1:nrow(x_train))
# x_train = x_train[shuf,,,,]
# y_train = y_train[shuf,]
# 
# 
# filterSize = c(3,3,3)
# filterNum = 16
# filterDropOut = 0.4
# 
# model <- keras_model_sequential()
# model %>% 
#   layer_conv_3d(filter = filterNum, kernel_size = filterSize, padding = 'same', input_shape = c(resNum,resNum,resNum,length(features))) %>%
#   layer_dropout(filterDropOut) %>%
#   layer_activation("relu") %>%
#   
#   layer_conv_3d(filter = filterNum, kernel_size = filterSize, padding = 'same') %>%
#   layer_dropout(filterDropOut) %>%
#   layer_activation("relu") %>%
#   layer_max_pooling_3d(pool_size = c(2,2,2))  %>%
#   
#   layer_conv_3d(filter = filterNum, kernel_size = filterSize, padding = 'same') %>%
#   layer_dropout(filterDropOut) %>%
#   layer_activation("relu") %>%
#   layer_max_pooling_3d(pool_size = c(2,2,2))  %>%
#   
#   layer_flatten() %>%
#   layer_dropout(0.5) %>%
#   
#   layer_dense(50) %>%
#   layer_activation("relu") %>%
#   layer_dropout(0.5) %>%
#   
#   layer_dense(20) %>%
#   layer_activation("relu") %>%
#   layer_dropout(0.5) %>%
#     
#   layer_dense(10) %>%
#   layer_activation("relu") %>%
#   layer_dropout(0.5) %>%
#     
#   layer_dense(units = 2, activation = 'softmax')
# 
# model %>% compile(
#   loss = 'categorical_crossentropy',
#   optimizer = optimizer_adam(),
#   metrics = c('accuracy')
# )
# 
# summary(model)
# 
# 
# weights = c(1,100000)
# weights_in <- split(weights, c(0,1))
# print(weights_in)
# 
# history <- model %>% fit(
#   x_train, y_train, 
#   epochs = 20, batch_size = 16,
#   weights = list("0"=6,"1"=1),
#   validation_split = 0.1
# )
# 
# predictions <- model %>% predict_classes(x_test)
# 
# sum(predictions)
# 
# length(predictions)
# gt = reverseToCategorical(y_test,c(0,1))
# length(which((predictions == gt) == TRUE))/length(predictions)
# 
# TInds = which(gt == "1")
# length(which(predictions[TInds] == "1"))/length(TInds)
# 
# 
# 
# 
# 
# pred = predictions
# # print(pred)
# gt = reverseToCategorical(y_test,c(0,1))
# y_test_pred = rep("0",length(pred))
# su = 0
# for(i in 1:length(gt)){
#   if(gt[i] == as.character(pred[i])) su = su + 1
#   
#   # y_test_pred[i] = TrFinal$classLevels[pred[i]]
# }
# su/length(gt)
# y_test_pred
# 
# 
# 
# confMat = table(factor(y_test_pred),levels = c("0","1"),
#                 factor(gt), levels = c("0","1"))
# 
# confMat


# #-------------------------------------------------------------------------------------------------------------
# # experimental end
# #-------------------------------------------------------------------------------------------------------------