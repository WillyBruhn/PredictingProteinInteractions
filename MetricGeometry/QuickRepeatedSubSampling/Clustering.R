#!/usr/bin/Rscript
#-----------------------------------------------------------------------------------------------------
# Willy Bruhn
# 22.8.19
#-----------------------------------------------------------------------------------------------------
# Train an auto-encoder that reduces the number of features of the proteins. Then the
# bottle-neck-layer is extracted and used to extract a condensed representation of
# the proteins. With these condensed representations then an agglomerative, bottom-up
# clustering is performed. With this clustering then a dendrogram is obtained.
# 
#   inputPath           ... path to a file with a .Rdata-extension that stores the feature-matrix for
#                           all proteins. Can be obtained by running Proteins.R.
# 
#   dendrogramName      ... name of the dendrogram-file. Without the (.pdf)-extension and prefix Dendrogram.
# 
#   numPermutations     ... number of permutations that are created for each representation. For
#                           each protein m rows from the feature-matrix are combined. The order of
#                           this m rows is randomly permutated and numPermutations different
#                           representations are created.
# 
#   m                   ... number of rows from the feature-matrix that are combined for each protein-model.
# 
#   epochs              ... number of epochs to train the autoencoder.
# 
#   batchSize           ... size of a batch. Relevant for the training.
# 
#   l1,l2,l3            ... specifies the encoder-dimensions, that is the size of each layer in the network.
#                           The autoencoder in this implementation has 3 layers.
# 
#   d1,d2,d3            ... in [0,1) specifies the dropout-rates.
#
#-----------------------------------------------------------------------------------------------------
library(permute)

wsPath = "/home/willy/PredictingProteinInteractions/setUp/SourceLoader.R"
# wsPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/setUp/SourceLoader.R"

wsPath = as.character(paste(funr::get_script_path(), "/../../setUp/SourceLoader.R", sep = ""))


source(wsPath)
sourceFiles(c("helperFunctions"))
sourceFiles(c("UltraQuickRepeatedSubSampling"))
sourceFiles(c("TriangulateIsoSurface"))
sourceFiles(c("kerasFunctions"))

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

#-----------------------------------------------------------------------------------------------------
# User-input
#-----------------------------------------------------------------------------------------------------
library(getopt)
options(warn=-1)
#----------------------------------------------------------------------------------
# Input
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'inputPath'  , 'i', 2, "character",
  'numPermutations'   , 'p', 2, "integer",
  'm'   , 'm', 2, "integer",
  'epochs'   , 'e', 2, "integer",
  'batchSize'   , 's', 2, "integer",
  
  'l1'   , 'a', 2, "integer",
  'l2'   , 'b', 2, "integer",
  'l3'   , 'c', 2, "integer",
  
  'd1'   , 'd', 2, "integer",
  'd2'   , 'f', 2, "integer",
  'd3'   , 'g', 2, "integer",
  
  'dendrogramName', 'q', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
if ( is.null(opt$inputPath    ) ) { opt$inputPath    = "/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test118/"}
if ( is.null(opt$numPermutations   ) ) { opt$numPermutations    = 200  }
if ( is.null(opt$m    ) ) { opt$m    = 50  }

if ( is.null(opt$epochs    ) ) { opt$epochs    = 20  }
if ( is.null(opt$batchSize    ) ) { opt$batchSize    = 1024  }
if ( is.null(opt$l1    ) ) { opt$l1    = 100  }
if ( is.null(opt$l2    ) ) { opt$l2    = 100  }
if ( is.null(opt$l3    ) ) { opt$l3    = 100  }

if ( is.null(opt$d1    ) ) { opt$d1    = 0.2  }
if ( is.null(opt$d2    ) ) { opt$d2    = 0.2  }
if ( is.null(opt$d3    ) ) { opt$d3    = 0.2  }

if ( is.null(opt$dendrogramName    ) ) { opt$dendrogramName    = "Autoencoder"  }

if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }

print("--------------------------------------------------------------------------")
print("Using the following parameters")
print(opt)
print("--------------------------------------------------------------------------")
#-----------------------------------------------------------------------------------------------------

getFunctionalProteins <- function(file = "/home/willy/PredictingProteinInteractions/data/106Test/labels.txt"){
  labels = read.table(file = file, header = TRUE)
  return(as.character(labels$name[which(labels$label == "functional")]))
}

getGeos <- function(quantiles){
  
  names = unique(quantiles[,1])
  
  geos = data.frame(matrix(0,ncol = ncol(quantiles), nrow = length(names)))
  colnames(geos) = colnames(quantiles)
  for(i in 1:length(names)){
    print(names[i])
    
    inds = which(quantiles[,1] == names[i])
    geo = colMeans(quantiles[inds,-1])
    # geo = colSums(quantiles[inds,-1])/nrow(quantiles[inds,-1])
    geos[i,-1] = geo
    geos[i,1] = as.character(names[i])
  }
  return(geos)
}

autoEncoder  <- function(x_train, epochs = 20, encoderDim = 3, unitNums = c(5,5,5), dropOuts = c(0.1,0.1,0.1), batchSize = 30){
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
  

  # 
  # # complex but wrong
  # model %>%
  #   layer_dense(units = unitNums[1], activation = "relu", input_shape = ncol(x_train)) %>%
  #   layer_dropout(dropOuts[1])
  # 
  # if(length(unitNums) > 1){
  #   for(i in 2:length(unitNums)){
  #     model %>%
  #       layer_dense(units = unitNums[i], activation = "relu") %>%
  #       layer_dropout(dropOuts[i])
  #   }
  # }
  # 
  # 
  # model %>% layer_dense(units = encoderDim, activation = "relu", name = "bottleneck")
  # 
  # model %>%
  #   layer_dense(units = unitNums[length(unitNums)], activation = "relu") %>%
  #   layer_dropout(dropOuts[length(unitNums)])
  # 
  # if(length(unitNums) > 1){
  #   for(i in (length(unitNums)-1):1){
  #     model %>%
  #       layer_dense(units = unitNums[i], activation = "relu") %>%
  #       layer_dropout(dropOuts[i])
  #   }
  # }
  # 
  # 
  # 
  # model %>% layer_dense(units = ncol(x_train))
  
  
  
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
    y = x_train,
    batch_size = batchSize,
    epochs = epochs,
    verbose = 1
  )
  
  return(model)
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

    
    layer_dense(units = encoderDim, activation = "relu", name = "bottleneck") %>%
    

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
    
    print(i/nrow(X_perm))
    
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

#--------------------------------------------------------------------------
# mat = matrix(c(1:15), ncol = 3)
# mat[c(1,2,1),]

# shuffle(c(100),)

# n = 106
# m = 100
# colN = 100
# 
# X = matrix(rnorm(colN*n*m),ncol = colN, nrow = n*m)
# 
# GR = createPermutations2(X,numPermutations = 10,m = m)
# 
# nrow(X)
# nrow(GR$X_perm)

# n = 10
# m = 2
# X = matrix(c(1:100),ncol = 1, nrow = n*m)
# GR = createPermutations2(X,names = c(1:n*m), y = matrix(rep(0,2*n*m), ncol = 2), numPermutations = 1,m = m)
# 
# GR$X_perm_out



#----------------------------------------------------------------------------------------------------------------

# inputPath = "/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test118/"
# inputPath = "/home/willy/PredictingProteinInteractions/data/120Experiment/NNexperimentsKfoldCV/Test101/"
# inputPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test201/"
inputPath = opt$inputPath

print(paste("reading from ", inputPath, "..."))

TR_fname = list.files(inputPath, pattern = ".Rdata")
TR = readRDS(paste(inputPath,TR_fname, sep = ""))
TrainTest = TR$TrainTest
Train_X = TrainTest$X
Train_X = apply(Train_X, 2, FUN = function(i) as.numeric(as.character(i)))

colMins = apply(Train_X,2,min)
colMaxs = apply(Train_X,2,max)
colRanges = colMaxs - colMins
Train_X = sapply(c(1:length(colRanges)), FUN = function(i){ (Train_X[,i]- colMins[i])/colRanges[i]})

# clear memory inbetween
gc()


print("creating representations ...")
GR = createPermutations2(Train_X, names = TR$originalNames, y = TR$TrainTest$y, numPermutations = opt$numPermutations, m = opt$m)

library(permute)
shuf = shuffle(1:nrow(GR$X_perm))
GR$X_perm = GR$X_perm[shuf,]
GR$X_perm_out = GR$X_perm_out[shuf,]
GR$originalNames = GR$originalNames[shuf]
GR$y = GR$y[shuf,]

suppressPackageStartupMessages(library(keras))
ncol(Train_X)

gc()

model = autoEncoder2(GR$X_perm, GR$X_perm_out, epochs = opt$epochs,encoderDim = 100,dropOuts = c(opt$d1,opt$d2,opt$d3), unitNums = c(opt$l1,opt$l2,opt$l2),batchSize = opt$batchSize)
# model %>% save_model_hdf5(paste(outPath,"autoEncodeModel.h5", sep =""))

# modelPreTrained <- load_model_hdf5(paste(outPath,"autoEncodeModel.h5", sep =""))
# modelPreTrained %>% summary()

# mod <- load_model_hdf5("/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test88/my_model.h5")
# mod %>% summary()


# clustering with the bottleneck
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, GR$X_perm)

length(which(colSums(intermediate_output) == 0))

withNames = data.frame(matrix(0, ncol = ncol(intermediate_output)+1, nrow = nrow(intermediate_output)))
withNames[,-1] = intermediate_output
withNames[,1] = GR$originalNames
# withNames[,1] = d2[,1]

geos = getGeos(withNames)


# library(rgl)
# 
# functionals = getFunctionalProteins()
# functionalInds2 = which(geos[,1] %in% c(functionals, "000_Trx"))
# points3d(geos[functionalInds2, -1], col = "red", size = 10)
# text3d(geos[functionalInds2, -1], texts = geos[functionalInds2, 1])
# 
# nonfunctionalInds2 = c(1:nrow(geos))[-functionalInds2]
# points3d(geos[nonfunctionalInds2,-1], col = "blue", size = 5)  
# text3d(geos[nonfunctionalInds2, -1], texts = geos[nonfunctionalInds2, 1])
# 
# #-------------------------------------------------------------------
# functionalInds = which(TR$originalNames %in% c(functionals, "000_Trx"))
# points3d(intermediate_output[functionalInds,], col = "red")
# 
# nonfunctionalInds = c(1:nrow(intermediate_output))[-functionalInds]
# points3d(intermediate_output[nonfunctionalInds,], col = "blue")
# #-------------------------------------------------------------------

plotBootleNeck <- function(geos, TR, intermediate_output, onlyGeos = FALSE, withNames = TRUE){
  
  functionals = getFunctionalProteins()
  functionalInds2 = which(geos[,1] %in% functionals)
  nonfunctionalInds2 = c(1:nrow(geos))[-functionalInds2]
  
  
  if(ncol(geos) == 3){
    plot(geos[,-1], col = "white")
    
    if(!onlyGeos){
      #-------------------------------------------------------------------
      functionalInds = which(TR$originalNames %in% functionals)
      nonfunctionalInds = c(1:nrow(intermediate_output))[-functionalInds]
      
      
      points(intermediate_output[nonfunctionalInds,], col = "blue", pch = 19, cex = 0.5)

      points(intermediate_output[functionalInds,], col = "red", pch = 19, cex = 0.5)
      

      #-------------------------------------------------------------------
    }
    

    points(geos[nonfunctionalInds2,-1], col = "blue", cex = 1, pch = 19)
    points(geos[functionalInds2, -1], col = "red",cex = 1, pch = 19)
    
    
    if(withNames == TRUE){
      text(geos[nonfunctionalInds2, -1], labels = geos[nonfunctionalInds2, 1], adj = c(0,1))
      text(geos[functionalInds2, -1], labels = geos[functionalInds2, 1], adj = c(0,1))

    }


    
  } else{
    
    points3d(geos[,-1], size = 1)
    
    points3d(geos[nonfunctionalInds2,-1], col = "blue", size = 5)  
    text3d(geos[nonfunctionalInds2, -1], texts = geos[nonfunctionalInds2, 1])
    
    points3d(geos[functionalInds2, -1], col = "red", size = 15)
    text3d(geos[functionalInds2, -1], texts = geos[functionalInds2, 1])
    
    
    
    if(!onlyGeos){
      #-------------------------------------------------------------------
      functionalInds = which(TR$originalNames %in% functionals)
      points3d(intermediate_output[functionalInds,], col = "red")
      
      nonfunctionalInds = c(1:nrow(intermediate_output))[-functionalInds]
      points3d(intermediate_output[nonfunctionalInds,], col = "blue")
      #-------------------------------------------------------------------
    }
  }


}

# citation("ggdendro")
library("cluster")
library("tcltk")
library("ggplot2")
library("ggdendro")
if(!WS_flag) library("plot3D")
library("emdist")

mydendrogramplot2 <- function(outPath, dist, labels,fName){

  # print(labels)
  if(is.null(labels)){
    print("no labels specified ...")
    labels = data.frame(matrix(0, ncol = 2, nrow = nrow(dist)))
    print(nrow(labels))
    print(nrow(dist))
    
    colnames(labels) = c("name","label")
    labels[,2] = rep("label_unknown",nrow(dist))
    labels[,1] = rownames(dist)
  }

  hc2 = hclust(dist(dist), "ave")
  dendr2    <- dendro_data(hc2, type="rectangle") # convert for ggplot
  clust2    <- cutree(hc2,k=2)                    # find 2 clusters
  clust2.df <- data.frame(label=names(clust2), cluster=factor(labels$label))
  
  dendr2[["labels"]] <- merge(dendr2[["labels"]],clust2.df, by="label")
  
  p <- ggplot() + geom_segment(data=segment(dendr2), aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_text(data=label(dendr2), aes(x, y, label=label, hjust=0, color=cluster), 
              size=2) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
  plot(p)
  
  if(fName != ""){
    print(paste(outPath,"/Dendrogram_", fName, ".pdf",sep=""))
    ggsave(filename = paste(outPath,"/Dendrogram_", fName, ".pdf",sep=""),height=8, width = 5)
  }

  
}

getDendrogramFromBootleNeck <- function(geos, labelFname = "/home/willy/PredictingProteinInteractions/data/labels.txt", outPath = "", fName = ""){
  geos = geos[order(geos[,1]),]
  
  geoDistances = as.matrix(dist(x = geos[,-1],upper = TRUE,diag = TRUE, method = "manhattan"))
  rownames(geoDistances) = geos[,1]
  colnames(geoDistances) = geos[,1]
  
  labels = read.table(labelFname, header = TRUE)
  
  labels = labels[which(labels$name %in% geos[,1]),]
  labels = labels[order(labels$name),]
  
  
  print(labels)
  
  num = 1
  geos[num+6,1]
  sort(geoDistances[num+6,], decreasing = TRUE)
  mydendrogramplot2(outPath,geoDistances, labels, fName)
}

if(!WS_flag) library(rgl)
if(!WS_flag) plotBootleNeck(geos,TR, intermediate_output, onlyGeos = FALSE, withNames = FALSE)
# getDendrogramFromBootleNeck(geos, outPath = outPath, fName = "AutoEncoder", labelFname = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/labels.txt")

# getDendrogramFromBootleNeck(geos, outPath = outPath, fName = "AutoEncoder", labelFname = "/home/willy/PredictingProteinInteractions/data/120Experiment/labels120.txt")
getDendrogramFromBootleNeck(geos, outPath = opt$inputPath, fName = opt$dendrogramName)


# ###############################################################################
# # subset
# subset = read.table("//home/sysgen/Documents/LWB/PredictingProteinInteractions/data/106Test/subsetNames.txt", header = FALSE)
# 
# # indices = which(TR$originalNames %in% as.character(subset[,1]))
# # TRSUB = TR
# # TRSUB$TrainTest$X = TR$TrainTest$X[indices,]
# # TRSUB$originalNames = TR$originalNames[indices]
# # TRSUB$TrainTest$y = TR$TrainTest$y[indices,]
# 
# indices2 = which(geos[,1] %in% as.character(subset[,1]))
# 
# length(which( subset[,1] %in% getFunctionalProteins("/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/labels.txt")))
# length(subset[,1])
# 
# 
# labels = read.table("/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/labels.txt", header = TRUE)
# labels
# 
# 
# labels[which(labels$name %in% geos[indices2,1]),]
# 
# 
# getDendrogramFromBootleNeck(geos[indices2,], outPath = outPath, fName = "AutoEncoderSubset", labelFname = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/labels.txt")
# 

