#!/usr/bin/Rscript
#--------------------------------------------------
# Willy Bruhn
# 22.8.19
#
#
#--------------------------------------------------

#--------------------------------------------------------------------------------------------------
library(permute)

wsPath = "/home/willy/PredictingProteinInteractions/setUp/SourceLoader.R"
# wsPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/setUp/SourceLoader.R"

# mode = "onlyExperiments2"
# mode = "onlyGenerateModels"

# mode ="randomSearch"

# mode = "hlat"

wsPath = as.character(paste(funr::get_script_path(), "/../../setUp/SourceLoader.R", sep = ""))


source(wsPath)
sourceFiles(c("helperFunctions"))
sourceFiles(c("UltraQuickRepeatedSubSampling"))
sourceFiles(c("TriangulateIsoSurface"))
sourceFiles(c("kerasFunctions"))

# .restar


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
    
    print(i/numPermutations)
    
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

# n = 106
# m = 100
# colN = 1000
# 
# X = matrix(rnorm(colN*n*m),ncol = colN, nrow = n*m)
# 
# GR = createPermutations(X,numPermutations = 100,m = m)

#----------------------------------------------------------------------------------------------------------------


outPath = "/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test87/"
# outPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test201/"

TR_fname = list.files(outPath, pattern = ".Rdata")
TR = readRDS(paste(outPath,TR_fname, sep = ""))
TrainTest = TR$TrainTest
Train_X = TrainTest$X
Train_X = apply(Train_X, 2, FUN = function(i) as.numeric(as.character(i)))

colMins = apply(Train_X,2,min)
colMaxs = apply(Train_X,2,max)
colRanges = colMaxs - colMins
Train_X = sapply(c(1:length(colRanges)), FUN = function(i){ (Train_X[,i]- colMins[i])/colRanges[i]})


nrow(Train_X)


GR = createPermutations(Train_X, numPermutations = 100, m = 100)

GR$X_perm[1:10,1:10]
GR$X_perm_out[1:5,1:10]


library(permute)
shuf = shuffle(1:nrow(GR$X_perm))
GR$X_perm = GR$X_perm[shuf,]
GR$X_perm_out = GR$X_perm_out[shuf,]
GR$originalNames = TR$originalNames[shuf]
GR$y = TR$TrainTest$y[shuf,]



suppressPackageStartupMessages(library(keras))
ncol(Train_X)

gc()

model = autoEncoder2(GR$X_perm, GR$X_perm_out, epochs = 15,encoderDim = 3,dropOuts = c(0.0,0.0,0.0,0.0,0.0), unitNums = c(500,300,100,100),batchSize = 256)
model %>% save_model_hdf5(paste(outPath,"autoEncodeModel.h5", sep =""))

modelPreTrained <- load_model_hdf5(paste(outPath,"autoEncodeModel.h5", sep =""))
modelPreTrained %>% summary()

mod <- load_model_hdf5("/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test88/my_model.h5")
mod %>% summary()


model_built_from_pretrained <- function(GR, epochs = 30, batch_size = 64, weights, sampleSize = NULL){
 
  encoder <- keras_model(inputs = modelPreTrained$input, outputs = get_layer(modelPreTrained, "bottleneck")$output)
  
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
   
  return(fullModel)
}



model_built_from_pretrained(GR,weights = list("1" = 6, "0" = 1))


# clustering with the bootleneck
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, GR$X_perm)

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

plotBootleNeck <- function(geos, TR, intermediate_output, onlyGeos = FALSE){
  functionals = getFunctionalProteins()
  functionalInds2 = which(geos[,1] %in% functionals)
  points3d(geos[functionalInds2, -1], col = "red", size = 15)
  text3d(geos[functionalInds2, -1], texts = geos[functionalInds2, 1])
  
  nonfunctionalInds2 = c(1:nrow(geos))[-functionalInds2]
  points3d(geos[nonfunctionalInds2,-1], col = "blue", size = 5)  
  text3d(geos[nonfunctionalInds2, -1], texts = geos[nonfunctionalInds2, 1])
  
  if(!onlyGeos){
    #-------------------------------------------------------------------
    functionalInds = which(TR$originalNames %in% functionals)
    points3d(intermediate_output[functionalInds,], col = "red")
    
    nonfunctionalInds = c(1:nrow(intermediate_output))[-functionalInds]
    points3d(intermediate_output[nonfunctionalInds,], col = "blue")
    #-------------------------------------------------------------------
  }

}

library("cluster")
library("tcltk")
library("ggplot2")
library("ggdendro")
library("plot3D")
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

library(rgl)

if(!WS_flag) plotBootleNeck(geos,TR, intermediate_output, onlyGeos = FALSE)
# getDendrogramFromBootleNeck(geos, outPath = outPath, fName = "AutoEncoder", labelFname = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/labels.txt")
getDendrogramFromBootleNeck(geos, outPath = outPath, fName = "AutoEncoder")


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
