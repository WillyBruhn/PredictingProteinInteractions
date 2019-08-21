#--------------------------------------------------
# Willy Bruhn
# 22.8.19
#
# Using SOMs to cluster the features.
#
#--------------------------------------------------

install.packages("kohonen")
library(kohonen)

quantiles_01 = read.csv("/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv", header = TRUE)
quantiles_08 = read.csv("/home/willy/PredictingProteinInteractions/data/106Test/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv", header = TRUE)

quantiles_01[1:5,]

?kohonen


data(degelder)
mydata <- list(patterns = degelder$patterns,CellVol = log(degelder$properties[,"cell.vol"]))
## custom distance function
require(Rcpp)
sourceCpp(system.file("Distances", "wcc.cpp", package = "kohonen"))
set.seed(7)
powsom <- supersom(data = mydata, grid = somgrid(6, 4, "hexagonal"),dist.fcts = c("WCCd", "sumofsquares"),keep.data = TRUE)
summary(powsom)


mydata


quantiles_01_list = as.list(quantiles_01[,-1])

length(quantiles_01_list)


powsom <- supersom(data = as.list(quantiles_01[,-1]), grid = somgrid(10, 10, "hexagonal"),
                                dist.fcts = c("manhattan"),keep.data = TRUE)


powsom$unit.classif







?supersom


  
data(yeast)
yeast.supersom <- supersom(yeast, somgrid(6, 6, "hexagonal"),
                           whatmap = c("alpha", "cdc15", "cdc28", "elu"),
                           maxNA.fraction = .5)

plot(yeast.supersom, "changes")

length(yeast)


mygrid <- somgrid(5, 5, "hexagonal")
fakesom <- list(grid = mygrid)
class(fakesom) <- "kohonen"

par(mfrow = c(2,1))
dists <- unit.distances(mygrid)
plot(fakesom, type="property", property = dists[1,],
     main="Distances to unit 1", zlim=c(0,6),
     palette = rainbow, ncolors = 7)

dists <- unit.distances(mygrid, toroidal=TRUE)
plot(fakesom, type="property", property = dists[1,],
     main="Distances to unit 1 (toroidal)", zlim=c(0,6),
     palette = rainbow, ncolors = 7)

#--------------------------------------------------------------------------------------------------

getGeos <- function(quantiles){
  
  names = unique(quantiles[,1])
  
  geos = data.frame(matrix(0,ncol = ncol(quantiles), nrow = length(names)))
  colnames(geos) = colnames(quantiles)
  for(i in 1:length(names)){
    print(names[i])
    
    inds = which(quantiles[,1] == names[i])
    geo = colMeans(quantiles[inds,-1])
    geos[i,-1] = geo
    geos[i,1] = as.character(names[i])
  }
  return(geos)
}

geos = getGeos(quantiles_08)

library(rgl)
functionals = c(getFunctionalProteins(), "000_Trx")
plot_prot_quants(quantiles = geos,q = 1,functionals = functionals,withEuclid = TRUE)




suppressPackageStartupMessages(library(keras))

# set training data
x_train <- as.matrix(quantiles_08[,-1])


array_reshape(x = quantiles_08[,-1],dim = c(106, nrow(quantiles_08)*(ncol(quantiles_08)-1)/106))

quantiles_08_reshaped = matrix(unlist(quantiles_08_num),nrow = 106, ncol = nrow(quantiles_08)*ncol(quantiles_08[,-1])/106, byrow = TRUE)

quantiles_08_num = lapply(quantiles_08[,-1], FUN = function(i) as.numeric(i))









# set model
model <- keras_model_sequential()
model %>%
  layer_dense(units = 6, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 2, activation = "tanh", name = "bottleneck") %>%
  layer_dense(units = 6, activation = "tanh") %>%
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
  y = x_train,
  batch_size = 30,
  epochs = 2,
  verbose = 1
)

# evaluate the performance of the model
mse.ae2 <- evaluate(model, x_train, x_train)
mse.ae2


# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, x_train)





model = load_model_hdf5("/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test51/my_model.h5")

TR = readRDS("/home/willy/PredictingProteinInteractions/data/106Test/NNexperimentsKfoldCV/Test51/Test51sS_5_sT_400_sTt_10_euklid_TRUE_pos_TRUE_neg_TRUE_pos_neg_TRUE_globalToo_FALSE.Rdata")


TrainTest = TR$TrainTest
Train_X = TrainTest$X
Train_X = apply(Train_X, 2, FUN = function(i) as.numeric(as.character(i)))


model


intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "dense_1")$output)
intermediate_output <- predict(intermediate_layer_model, Train_X)



withNames = data.frame(matrix(0, ncol = 4, nrow = nrow(intermediate_output)))
# withNames = cbind(TR$originalNames,intermediate_output)

withNames[,-1] = intermediate_output
withNames[,1] = TR$originalNames
geos = getGeos(withNames)

geos


functionalInds = which(TR$originalNames %in% c(functionals, "000_Trx"))
points3d(intermediate_output[functionalInds,], col = "red")

nonfunctionalInds = c(1:nrow(intermediate_output))[-functionalInds]
points3d(intermediate_output[nonfunctionalInds,], col = "blue")

functionalInds2 = which(geos[,1] %in% c(functionals, "000_Trx"))
points3d(geos[functionalInds2, -1], col = "red", size = 10)
text3d(geos[functionalInds2, -1], texts = geos[functionalInds2, 1])

nonfunctionalInds2 = c(1:nrow(geos))[-functionalInds2]
points3d(geos[nonfunctionalInds2,-1], col = "blue")  
text3d(geos[nonfunctionalInds2, -1], texts = geos[nonfunctionalInds2, 1])





