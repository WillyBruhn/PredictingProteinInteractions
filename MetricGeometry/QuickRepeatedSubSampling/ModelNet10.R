#---------------------------------------------------------------------
# Willy Bruhn, 19.7.19
#
# Read in the model-net10 data-set and create the usual projections.
#
#---------------------------------------------------------------------

s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/helperFunctions.R"
source(s1)

s2 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s2)

library(rgl)

getModel10Net <- function(fName, plot = FALSE){
  
  tab = read.table(fName,skip=2, col.names = c("x","y", "z", "trash"), fill = TRUE)
  tab_xyz = tab[is.na(tab$trash),1:3]
  
  if(plot == TRUE) points3d(tab_xyz)
  return(tab_xyz)
}

getDataSet <- function(datasetPath){
  allFiles = list.files(datasetPath, recursive = TRUE, full.names = TRUE)
  
  dataSet = data.frame(matrix(0,nrow = length(allFiles), ncol = 3))
  colnames(dataSet) = c("model", "trainTest", "file")

  for(i in 1:length(allFiles)){
    vec = strsplit(allFiles[i],split = "/")[[1]]
    vecSmall = vec[(length(vec)-1):length(vec)]
    
    trainTest = vecSmall[1]
    modelName = strsplit(vecSmall[2],split = ".off")[[1]]
    
    dataSet[i,] = c(modelName, trainTest, allFiles[i])
  }  

  return(dataSet)
}

generateF_approximations_3dModel <- function(model_points, n = 100, m = 10, q = 2){
  
  pos13_F_list = list()
  pos13_F_approx_list = list()
  
  for(i in 1:m){
    pos13_F_list[[i]] = samplePointsAndCalculateCDFofEc(all_pts = model_points, n = n,plot = FALSE)
    pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
  }
  
  return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
}

getAllModel_F_approximations <- function(model_vec, n = 100, m = 50, q = 2){
  
  distributions_lists = list()
  
  for(i in 1:length(model_vec)){
    print(i/length(model_vec))
    F_app = generateF_approximations_3dModel(model_vec[[i]]$vert,n = n,q = q, m = m)
    distributions_lists[[i]] =  list("name" = model_vec[[i]]$name,"F" = F_app)
  }
  
  return(distributions_lists)
}

getAllModels <- function(dataSet, maxNum = NULL){
  all_models = list()
  
  if(is.null(maxNum)) maxNum = nrow(dataSet)
  
  for(i in 1:maxNum){
    print(paste(i/maxNum))
    m1 = getModel10Net(dataSet$file[i], FALSE)
    n1 = dataSet$model[i]
    all_models[[i]] = list("vert" = m1, "name" = n1)
  }
  
  return(all_models)
}


#------------------------------------------------------------------------
n = 30
m = 100
q = 10

pathToProjection = "/home/willy/PredictingProteinInteractions/data/ModelNet10/projections/"

datasetPath = "/home/willy/PredictingProteinInteractions/data/ModelNet10/ModelNet10/"



# get all the file names and information if it belongs to train or test
dataSet = getDataSet(datasetPath)


dataSetTrain = dataSet[which(dataSet[,2] == "train"),]

all_models = getAllModels(dataSetTrain)

distributions_all = getAllModel_F_approximations(model_vec = all_models,n = n,m = m,q = q)

projection = getManhattanProjection(distributions_all)

writeProjectionToFile(proj = projection,path = pathToProjection, n = n, m = m, q = q, fName = "proj")



