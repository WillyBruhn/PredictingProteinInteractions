#----------------------------------------------------------
# Using an autoencoder to build a permutation-invariant nn
#----------------------------------------------------------


suppressPackageStartupMessages(library(keras))



weighted_mse <- function(y_true, y_pred, weights){
  K        <- backend()
  weights  <- K$variable(weights)
  # calculate the metric
  loss <- K$sum(weights * (K$pow(y_true - y_pred, 2))) 
  loss
}

metric_weighted_mse <- custom_metric("weighted_mse", function(y_true, y_pred) {
  weighted_mse(y_true, y_pred, weights)
})

model <- model %>% compile(
  loss = 'mse', 
  optimizer = 'rmsprop',
  metrics = metric_weighted_mse)


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

  K <- backend()
  
  true_positives = K$sum(K$round(K$clip(y_true * y_pred, 0, 1)))
  predicted_positives = K$sum(K$round(K$clip(y_pred, 0, 1)))
  precision = true_positives / (predicted_positives + EPSILON)
  
  return(precision)
}
# precision_keras_metric(c(1,0,1,0,1), c(1,0,1,0,0))


permLoss <- function(y_true, y_pred){
  K <- backend()

  y_true_sorted = sort(y_true)
  return(0)
}
# f1_keras_metric(c(1,0,1,0,1), c(1,0,1,0,0))

customPermutationLoss <- custom_metric("permLoss", function(y_true, y_pred) {
  # print("getting called ...")
  return(permLoss(y_true, y_pred))
})




autoEncoderPermutation  <- function(x_train, x_train_perm, epochs = 20, encoderDim = 3, unitNums = c(5,5,5), dropOuts = c(0.1,0.1,0.1)){
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = unitNums[1], activation = "relu", input_shape = ncol(x_train)) %>%
    layer_dense(units = unitNums[2], activation = "relu") %>%
    layer_dense(units = unitNums[3], activation = "relu") %>%
    layer_dense(units = unitNums[4], activation = "relu") %>%

    layer_dense(units = encoderDim, activation = "relu", name = "bottleneck") %>%
    

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
  
  # model %>% compile(
  #   loss = 'mse', 
  #   optimizer = 'rmsprop',
  #   metrics = customPermutationLoss)
  
  # fit model
  model %>% fit(
    x = x_train,
    y = x_train_perm,
    batch_size = 30,
    epochs = epochs,
    verbose = 1
  )
  
  return(model)
}

createPermutations <- function(X, numPermutations = 1, blockSize = 1){
  # blockSize ... size of one block, that means numbers that should not be permutated
  # it has to hold ncol(X) = i * blockSize, where i is a natural number
  
  if(ncol(X) %% blockSize != 0) return(FALSE)
  
  X_perm = matrix(0,ncol = ncol(X), nrow = numPermutations*nrow(X))
  X_notPerm = matrix(0,ncol = ncol(X), nrow = numPermutations*nrow(X))
  
  blockNum = ncol(X)/blockSize
  for(i in 1:numPermutations){
    for(j in 1:nrow(X)){
      ind = (i-1)*nrow(X)+j
      
      blockOrder = shuffle(blockNum)
      
      for(k in 1:blockNum){
          colIndStart = (k-1)*blockSize+1
          colIndEnd = colIndStart+blockSize-1
          
          colIndStart2 = (blockOrder[k]-1)*blockSize+1
          colIndEnd2 = colIndStart2+blockSize-1
          
          # print(colIndStart:colIndEnd)
          
          X_perm[ind,colIndStart:colIndEnd] = X[j,colIndStart2:colIndEnd2]
      }
      
      # X_perm[ind,] = sort(X[j,])
      X_notPerm[ind,] = X[j,]
    }
  }
  
  return(list("X" = X_notPerm, "X_perm" = X_perm))
}

customEvaluation <- function(X, Xperm){
  s = 0
  for(i in 1:nrow(X)){
    s = s + sum(abs(sort(X[i,])-sort(Xperm[i,])))
  }
  
  return(s)
}

#--------------------------------------------------

n = 10
# quants = 3*6*5
quants = 1
features = 3
blocks = 1

# blocks ... number of points combined
# quants ... quantiles
# features ... number of features
dim = quants*blocks*features

X = matrix(rnorm(n = n*dim), ncol = dim)

ncol(X)

library(permute)
Tr = createPermutations(X, numPermutations = 10, blockSize = features*quants)


model = autoEncoderPermutation(Tr$X, Tr$X_perm,unitNums = c(100,100,100,100), dropOuts = c(0.0), encoderDim = 3,epochs = 10)


intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, Tr$X)

intermediate_output_perm <- predict(intermediate_layer_model, Tr$X_perm)

sum(abs(intermediate_output- intermediate_output_perm))


X_output <- predict(model, Tr$X)
customEvaluation(Tr$X, X_output)

customEvaluation(intermediate_output,intermediate_output_perm)



intermediate_output[1:5,]
intermediate_output_perm[1:5,]
Tr$X[1:5,]
X_output[1:5,]
