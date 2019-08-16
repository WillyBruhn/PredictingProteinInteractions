#!/usr/bin/Rscript
# options <- commandArgs(trailingOnly = TRUE)



#-------------------------------------------------------
#
# install packages
#
#-------------------------------------------------------

print("installing packages ...")

is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}

loadPck <- function(pck){
  if(!is.installed(pck)){install.packages(pck, repos='http://cran.us.r-project.org')}
}

loadPck("RMySQL")
loadPck("plot3D")
loadPck("funr")
loadPck("emdist")
loadPck("ggplot2")
loadPck("ggdendro")
loadPck("tcltk")
loadPck("getopt")
loadPck("rfUtilities")
loadPck("rbenchmark")
loadPck("raster")
loadPck("foreach")
loadPck("FNN")
loadPck("RSNNS")
loadPck("rdist")
loadPck("doMC")
loadPck("misc3d")
loadPck("geometry")
loadPck("xtable")
loadPck("spatstat")
loadPck("igraph")
loadPck("gplots")
loadPck("lubridate")
loadPck("pracma")
loadPck("rgl")
loadPck("keras")
loadPck("readobj")
loadPck("permute")
loadPck("DescTools")
loadPck("doBy")
loadPck("caret")
loadPck("beepr")


# # installation on WS:
# # inside the virtual-env-directories you have to call
# #
# # ./pip3.6 install tensorflow
# #
# # ./pip3.6 install keras
# #
# # and optionally 
# #
# # ./pip3.6 install tensorflow-gpu
# #
# 
# install.packages("devtools")
# devtools::install_github("rstudio/reticulate")
# library(reticulate)
# use_python("/home/sysgen/.pyenv/versions/3.6.3/bin/python3.6", required = TRUE)
# use_virtualenv("/home/sysgen/.pyenv/versions/3.6.3/",required = TRUE)
# # use_virtualenv("~/.virtualenvs/r-reticulate/",required = TRUE)
# 
# reticulate::py_discover_config(required = TRUE)
# py_config()
# reticulate::py_config()
# 
# install.packages("tensorflow")
# library(tensorflow)
# tensorflow::install_tensorflow(envname = "/home/sysgen/.pyenv/versions/3.6.3/")
# 
# devtools::install_github("rstudio/keras")
# library(keras)
# install_keras(method = c("virtualenv"),tensorflow = "/home/sysgen/.pyenv/versions/3.6.3/")
# # install_keras()
# 
# minimalExampleKeras()


#-----------------------------------------------------------------
## /.pyenv/versions/3.6.3/bin> ldd python3.6
# .pyenv/versions/3.6.3/bin> ldd python3.6
# linux-vdso.so.1 (0x00007fff033fd000)
# libpthread.so.0 => /lib64/libpthread.so.0 (0x00007fc736b92000)
# libdl.so.2 => /lib64/libdl.so.2 (0x00007fc73698e000)
# libutil.so.1 => /lib64/libutil.so.1 (0x00007fc73678b000)
# libm.so.6 => /lib64/libm.so.6 (0x00007fc73648e000)
# libc.so.6 => /lib64/libc.so.6 (0x00007fc7360e9000)
# /lib64/ld-linux-x86-64.so.2 (0x00007fc736daf000)
#-----------------------------------------------------------------

#-------------------------------------------------------------
# python:         /usr/bin/python
# libpython:      /usr/lib64/libpython2.7.so.1.0
# pythonhome:     /usr:/usr
# version:        2.7.13 (default, Jan 03 2017, 17:41:54) [GCC]
# numpy:          /usr/lib64/python2.7/site-packages/numpy
# numpy_version:  1.8.0
# 
# python versions found: 
#   /home/sysgen/.virtualenvs/r-reticulate/bin/python
# /usr/bin/python
# /usr/bin/python3



# use_python("/usr/bin/python3", required = TRUE)
## python:         /usr/bin/python3
## libpython:      /usr/lib64/libpython3.4m.so.1.0
## pythonhome:     /usr:/usr
## version:        3.4.6 (default, Mar 22 2017, 12:26:13) [GCC]
## numpy:           [NOT FOUND]
## 
## NOTE: Python version was forced by use_python function

# use_python("/home/sysgen/.virtualenvs/r-reticulate/bin/python3.4", required = TRUE)
# python:         /home/sysgen/.virtualenvs/r-reticulate/bin/python3.4
# libpython:      /usr/lib64/libpython3.4m.so.1.0
# pythonhome:     /usr:/usr
# virtualenv:     /home/sysgen/.virtualenvs/r-reticulate/bin/activate_this.py
# version:        3.4.6 (default, Mar 22 2017, 12:26:13) [GCC]
# numpy:           [NOT FOUND]
# 
# NOTE: Python version was forced by use_python function
#-------------------------------------------------------------


# install.packages("keras")
# library(keras)

# install.packages("lime")
# library(lime)
# 
# # install.packages("tidyquant")
# library(tidyquant)
# 
# # install.packages("rsample")
# # library(rsample)
# 
# # install.packages("recipes")
# library(recipes)
# 
# # install.packages("yardstick")
# library(yardstick)
# 
# # install.packages("corrr")
# # library(corrr)
# 
# library(permute)
# library(xtable)
# install.packages("DescTools")

# for yupiter ...
# https://www.datacamp.com/community/blog/jupyter-notebook-r?utm_source=adwords_ppc&utm_campaignid=898687156&utm_adgroupid=48947256715&utm_device=c&utm_keyword=&utm_matchtype=b&utm_network=g&utm_adpostion=1t1&utm_creative=255798340456&utm_targetid=dsa-473406574235&utm_loc_interest_ms=&utm_loc_physical_ms=1004446&gclid=CjwKCAjw-ITqBRB7EiwAZ1c5U2ue8Ma0thAn3YtgCI-B15RYUtSIZSxCqGx4dPPT8BoKx7gYbnHrchoCL4AQAvD_BwE


#-------------------------------------------------------
#
# sourcable files
#
#-------------------------------------------------------
# PredictingProteinInteractionsFolder = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/"
# PredictingProteinInteractionsFolder = options[1]

path = funr::get_script_path()
vec = strsplit(path, "/setUp")

PredictingProteinInteractionsFolder =  paste(vec[[1]][1],"/", sep ="")
print(PredictingProteinInteractionsFolder)


addSource <- function(table, name, fName){
  sources2 = data.frame(matrix(0,ncol = 2, nrow = 1))
  colnames(sources2) = c("name", "path")
  
  sources2[1, ] = c(name,paste(PredictingProteinInteractionsFolder,fName, sep = ""))
  
  table = rbind(table, sources2)
  return(table)
}

print("creating sources-table ...")
sources = data.frame(matrix(0,ncol = 2, nrow = 0))
colnames(sources) = c("name", "path")

sources = addSource(sources, "helperFunctions", "/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R")
sources = addSource(sources, "UltraQuickRepeatedSubSampling", "/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R")
sources = addSource(sources, "TriangulateIsoSurface", "/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R")
sources = addSource(sources, "QuickRepeatedSubSampling", "/MetricGeometry/QuickRepeatedSubSampling/QuickRepeatedSubSampling.R")
sources = addSource(sources, "BoostedKNN", "/Classification/NNClassification/optimizeDifferentModels/BoostedKNN.R")

sources = addSource(sources, "extrinsicDistances", "/Classification/NNClassification/additionalScripts/extrinsicDistances.R")
sources = addSource(sources, "isoFaces", "/Classification/NNClassification/additionalScripts/isoFaces.R")
sources = addSource(sources, "kerasFunctions", "/Classification/NNClassification/kerasFunctions.R")

write.table(sources, paste(PredictingProteinInteractionsFolder, "setUp/SourcableFiles.txt", sep = ""), row.names = FALSE)


#------------------------------------------------------------------------------------------------------------------------
# other paths and data-sets
#------------------------------------------------------------------------------------------------------------------------

print("creating paths-table ...")
paths = data.frame(matrix(0,ncol = 2, nrow = 0))
colnames(paths) = c("name", "path")

paths = addSource(paths, "Manifold", "/Manifold/build/")
paths = addSource(paths, "ModelNet10", "/data/ModelNet10/ModelNet10/")
paths = addSource(paths, "120Experiment", "/data/120Experiment/Output/")
paths = addSource(paths, "106Experiment", "/data/106Test/Output/")
paths = addSource(paths, "pdbDownloaderExperiment", "/data/pdbDownloaderExperiment/Output/")
paths = addSource(paths, "animals", "/home/willy/PredictingProteinInteractions/data/animals/models/")
paths = addSource(paths, "106ExperimentLabels", "./data/106Test/labels.txt")

write.table(paths, paste(PredictingProteinInteractionsFolder, "setUp/Paths.txt", sep = ""), row.names = FALSE)




# #------------------------------------------------------------------------------------------------------------------------
# # other paths and data-sets
# #------------------------------------------------------------------------------------------------------------------------


minimalExampleKeras <- function(){
  #loading the keras inbuilt mnist dataset
  data<-dataset_mnist()
  
  
  #separating train and test file
  train_x<-data$train$x
  train_y<-data$train$y
  test_x<-data$test$x
  test_y<-data$test$y
  
  rm(data)
  
  
  
  # converting a 2D array into a 1D array for feeding into the MLP and normalising the matrix
  train_x <- array(train_x, dim = c(dim(train_x)[1], prod(dim(train_x)[-1]))) / 255
  test_x <- array(test_x, dim = c(dim(test_x)[1], prod(dim(test_x)[-1]))) / 255
  
  #converting the target variable to once hot encoded vectors using keras inbuilt function
  train_y<-to_categorical(train_y,10)
  test_y<-to_categorical(test_y,10)
  
  #defining a keras sequential model
  model <- keras_model_sequential()
  
  #defining the model with 1 input layer[784 neurons], 1 hidden layer[784 neurons] with dropout rate 0.4 and 1 output layer[10 neurons]
  #i.e number of digits from 0 to 9
  
  model %>%
    layer_dense(units = 784, input_shape = 784) %>%
    layer_dropout(rate=0.4)%>%
    layer_activation(activation = 'relu') %>%
    layer_dense(units = 10) %>%
    layer_activation(activation = 'softmax')
  
  #compiling the defined model with metric = accuracy and optimiser as adam.
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = 'adam',
    metrics = c('accuracy')
  )
  
  #fitting the model on the training dataset
  model %>% fit(train_x, train_y, epochs = 2, batch_size = 128)
  
  #Evaluating model on the cross validation dataset
  loss_and_metrics <- model %>% evaluate(test_x, test_y, batch_size = 128)
  
  return(TRUE)
}

