# s1 = paste(funr::get_script_path(),"/additionalScripts/TriangulateIsoSurface.R", sep ="")
source("/home/willy/PredictingProteinInteractions/Classification/NNClassification/additionalScripts/TriangulateIsoSurface.R")

# s2 = paste(funr::get_script_path(),"/additionalScripts/TriangulateIsoSurface.R", sep ="")
source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/isoFaces.R")
# source(s2)

source("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/extrinsicDistances.R")







readDistanceMatrix1 <- function(file = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/d_matrix.csv"){
  d = read.csv(file = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/d_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
  d = data.matrix(frame = d)
  
  return(d)
}

readDistanceMatrix2 <- function(file = "/home/willy/RedoxChallenges/MasterThesis/memoliModels/FLB_repeated_sampling/FelixMasterThesis/ListOfEMD_positive_100.csv"){
  # Read the distance from a file formated with three columns
  # ProtA ProtB EMDdistances
  #
  # returns a quadratic matrix
  
  q = read.csv(file)
  names = unique(q$NameforEMDA)
  
  d = matrix(0, length(names), length(names))
  colnames(d) = names
  rownames(d) = names
  for(i in 1:nrow(q)){
    
    d[which(rownames(d) == q$NameforEMDA[i]),which(colnames(d) == q$NameforEMDB[i])] = q$EMD[i]
    d[which(colnames(d) == q$NameforEMDB[i]),which(rownames(d) == q$NameforEMDA[i])] = q$EMD[i]
  }
  return(d)
}



readDistanceMatrix3 <- function(file = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Example/EMD_1000_100_10.0000_1_1.0000_1.0000_1.0000_id_test.csv
"){
  q = read.csv(file, sep = ";")
  names = unique(q$protein_name)
  
  # print(names)
  
  d = matrix(0, length(names), length(names))
  colnames(d) = names
  rownames(d) = names
  for(i in 1:nrow(q)){
    
    d[which(rownames(d) == q$protein_name_target[i]),which(colnames(d) == q$protein_name[i])] = q$emd_distance[i]
    d[which(colnames(d) == q$protein_name[i]),which(rownames(d) == q$protein_name_target[i])] = q$emd_distance[i]
    
  }
  return(d)
}

readLabels <- function(file = "/home/willy/PredictingProteinInteractions/data/106Redoxins/Labels/labels.txt"){
  # reads the labels (functional, not_functional) of the proteins.
  
  t = read.table(file = file, header = TRUE)
  return(t)
}



d = readDistanceMatrix3("/home/willy/Schreibtisch/106Test/RepSubOutput/EMD_100_5000_1.0000_1.0000_0.0000_0.0000_id_test_NNact_0.csv")

d = readDistanceMatrix2("/home/willy/RedoxChallenges/MasterThesis/MACD/Output106_FLB_100/ListOfEMD_positive_100.csv")
d = readDistanceMatrix2("/home/willy/RedoxChallenges/MasterThesis/MACD/Output106_FLB_100/ListOfEMD_negative_100.csv")

# d = readDistanceMatrix3("/home/willy/Schreibtisch/120PdbTrxTest/PPI/RepSubOutput/EMD_100_5000_1.0000_1.0000_0.0000_0.0000_id_test_NNact_0.csv")
labels_train = readLabels("/home/willy/Schreibtisch//106Test/Output/labels.txt")

additional = read.table("/home/willy/Schreibtisch/120PdbTrxTest/PPI/RepSubOutput/predictions.txt")

colnames(additional) = colnames(labels_train)

labels_train2 = rbind(labels_train,additional)



getPredictions <- function(d_train,test_names, labels_train,kNN){
  test_indices = rep(0,length(test_names))
  for(i in 1:length(test_indices)){
    test_indices[i] = which(colnames(d_train) == test_names[i])
  }
  
  lab = rep("NA", ncol(d_train))
  
  for(i in 1:ncol(d_train)){
    ind = which(labels_train$name == names(d_train[,1])[i])
    if(length(ind) > 0 )  lab[i] = as.character(labels_train$label[ind])
  }
  
  predictions = memoliNNclassification(d_train,lab,test_indices,n=kNN)
  df = data.frame(test_names,predictions)
  
  return(df)
}


func = which(labels_train$label == "functional")
func2 = which(labels_train2$label == "functional")

scaled = cmdscale(d)
plot(scaled)
points(scaled[func2,], col = "red")
points(scaled[func,], col = "blue")
points(x = scaled[which(rownames(scaled) == "000_Trx"),1], y = scaled[which(rownames(scaled) == "000_Trx"),2], col = "green")


points(x = scaled[which(rownames(scaled) == "045"),1], y = scaled[which(rownames(scaled) == "045"),2], col = "yellow")
points(x = scaled[which(rownames(scaled) == "097"),1], y = scaled[which(rownames(scaled) == "097"),2], col = "brown")

points(c(1,2))

?points

plot(x = c(1), y = c(2))


getFunctionalProteins()

length(which(scaled[,1] > 1))/nrow(scaled)

func2
val = 2
length(which(scaled[func2,1] > val))/length(which(scaled[,1] > val))


length(which(labels_train2[-c(1:106),2] == "functional"))/length(labels_train2[-c(1:106),2])





LeaveOneOutCV <- function(d,classlabels, kNN = 1, positive.class = "functional", negative.class = "not_functional"){
  labels_un = sort(unique(classlabels))
  random_sample = matrix(0,nrow = n, ncol = length(labels_un))
  colnames(random_sample) = labels_un
  for(i in 1:length(labels_un)){
    random_camel = which(classlabels == labels_un[i])
    # print(random_camel)
    random_sample[,i] = sample(random_camel,n, replace = TRUE)
  }
  
  conf = matrix(0,ncol = length(labels_un), nrow = length(labels_un))
  colnames(conf) = paste(labels_un,"_pred",sep="")
  rownames(conf) = labels_un
  
  # print(conf)
  
  for(i in 1:nrow(d)){
    pred = memoliNNclassification(d,classlabels,i,n=kNN)
    
    # print(pred)
    
    # TP
    if(pred == classlabels[i] && classlabels[i] == positive.class){
      conf[1,1] = conf[1,1] + 1
    }
    
    # TN
    if(pred == classlabels[i] && classlabels[i] == negative.class){
      conf[2,2] = conf[2,2] + 1
    }
    
    # FN
    if(pred != classlabels[i] && classlabels[i] == positive.class){
      conf[1,2] = conf[1,2] + 1
    }
    
    # FP
    if(pred != classlabels[i] && classlabels[i] == negative.class){
      conf[2,1] = conf[2,1] + 1
    }
  }
  
  if(sum(conf) != nrow(d)) return(NULL)
  
  return(conf)
}


optimizeAndPlotKnnLOO <- function(d,labels,Kmax = 25){
  f1_scores = rep(0,Kmax)
  accs = rep(0,Kmax)
  TPR = rep(0,Kmax)
  TNR = rep(0,Kmax)
  
  for(k in c(1:Kmax)){
    conf =  LeaveOneOutCV(d,labels, kNN = k)
    
    F1 = F1_score_confusion(conf)
    
    f1_scores[k] = F1$F1
    accs[k] = F1$ACC
    
    TPR[k] = F1$TPR
    TNR[k] = F1$TNR
  }
  
  plot(y = f1_scores, x = c(1:Kmax), type = "l", col = "blue", ylim = c(0,1), ylab = "", xlab = "k nearest neighbors")
  points(y = accs, x = c(1:Kmax), type = "l", col = "red")
  points(y = TPR, x = c(1:Kmax), type = "l", col = "green")
  points(y = TNR, x = c(1:Kmax), type = "l", col = "darkgreen")
  abline(v = which.max(f1_scores))
  legend(x = 15, y = 0.4, legend=c("F_1", "ACC", "TPR", "TNR"),
         col=c("blue", "red", "green", "darkgreen"), lty=rep(1,4), cex=0.8)
}

optimizeAndPlotKnnLOO(d,labels_train$label)

LeaveOneOutCV(d,labels_train$label, kNN = 7)

optimizeAndPlotKnnBootstrap <- function(d,labels,Kmax = 25, n = 1000){
  f1_scores = rep(0,Kmax)
  accs = rep(0,Kmax)
  TPR = rep(0,Kmax)
  TNR = rep(0,Kmax)
  
  for(k in c(1:Kmax)){
    print(k)
    conf =  memoliNNclassificationErrorEstimate(d, labels, n = n, kNN = k,normalized = FALSE)
    
    F1 = F1_score_confusion(conf)
    
    f1_scores[k] = F1$F1
    accs[k] = F1$ACC
    
    TPR[k] = F1$TPR
    TNR[k] = F1$TNR
  }
  
  plot(y = f1_scores, x = c(1:Kmax), type = "l", col = "blue", ylim = c(0,1), ylab = "", xlab = "k nearest neighbors", main = "bootstrap")
  points(y = accs, x = c(1:Kmax), type = "l", col = "red")
  points(y = TPR, x = c(1:Kmax), type = "l", col = "green")
  points(y = TNR, x = c(1:Kmax), type = "l", col = "darkgreen")
  abline(v = which.max(f1_scores))
  legend(x = 15, y = 0.4, legend=c("F_1", "ACC", "TPR", "TNR"),
         col=c("blue", "red", "green", "darkgreen"), lty=rep(1,4), cex=0.8)
}


optimizeAndPlotKnnBootstrap(d,labels_train$label)

# conf = memoliNNclassificationErrorEstimate(d, labels_train2$label, n = 1000, kNN = 10,normalized = FALSE)
conf = memoliNNclassificationErrorEstimate(d, labels_train$label, n = 1000, kNN = 10,normalized = FALSE)
print(conf)



F1_score_confusion <- function(confusion, positive.class="functional") {
  #                functional_pred not_functional_pred
  # functional               0.923               0.077
  # not_functional           0.133               0.867
  #
  # TP FN
  # FP TN
  
  # check for correct format
  positive.class_pred = paste(positive.class,"_pred", sep = "")
  
  # if(which(colnames(confusion) == positive.class_pred) != 1 || which(rownames(confusion) == positive.class) != 1){
  #   print("Error: confusion-matrix does not have the correct format. Aborting ...")
  #   return(NULL)
  # }
  
  pos_pred_col = which(colnames(confusion) == positive.class_pred)
  pos_ind = which(rownames(confusion) == positive.class)
  
  # transform confusion so that in the top-left corner the True-positives are placed
  # TP FN
  # FP TN
  
  confusion2 = confusion
  confusion[1,] = confusion2[pos_ind,]
  confusion[2,] = confusion2[-pos_ind,]
  
  confusion2 = confusion
  confusion[,1] = confusion2[,pos_pred_col]
  confusion[,2] = confusion2[,-pos_pred_col]
  
  
  
  PPV = confusion[1,1]/(confusion[1,1] + confusion[2,1])
  TPR = confusion[1,1]/(confusion[1,1] + confusion[1,2])
  
  TNR = confusion[2,2]/(confusion[2,2] + confusion[2,1])
  
  F1 = 2*(PPV*TPR)/(PPV+TPR)
  
  ACC = (confusion[1,1] + confusion[2,2]) /(sum(confusion))
  
  # print(ACC)
  
  return(list("F1" = F1, "ACC" = ACC, "TPR" = TPR, "TNR" = TNR))
}

positive.class="functional"
length(which(colnames(conf) == positive.class))
length(which(rownames(conf) == positive.class))

which(rownames(conf) == positive.class)

conf2 = conf

conf2[1,] = conf[2,]
conf2[2,] = conf[1,]

rownames(conf2) = c(rownames(conf)[2], rownames(conf)[1])

conf2



F1_score_confusion(conf, positive.class = "not_functional")



conf2 = conf %*% matrix(c(1,0,0,1), nrow = 2, byrow = TRUE)



labels_train


View(d)

install.packages("kknn")
library(kknn)
train.kknn(data = d)

?kknn


install.packages("KernelKnn")
library(KernelKnn)

labels_train = readLabels("/home/willy/Schreibtisch//106Test/Output/labels.txt")
functionals = which(labels_train[,2] == "functional")
labels_train[,2] = 2
labels_train[functionals,2] = 1

weights <- function(x){
  
  print(x)
}


outProtein = distMat.KernelKnnCV(DIST_mat = d,y = as.vector(labels_train$label), k = 3, folds = 10, Levels = unique(labels_train[,2]), seed_num = 24, weights_function = weights)


calculatePredictionAccuracies(outProtein, labels_train[,2])
print((106-16)/106)





data(ionosphere)

X = ionosphere[, -c(2, ncol(ionosphere))]
y = as.numeric(ionosphere[, ncol(ionosphere)])

dist_obj = dist(X)

dist_mat = as.matrix(dist_obj)

out = distMat.KernelKnnCV(dist_mat, y, k = 5, folds = 2, Levels = unique(y))

out$preds
out$folds$fold_1
out$folds$fold_2


y2 = as.character(y)

y2
conf = memoliNNclassificationErrorEstimate(dist_mat, y2, n = 1000, kNN = 10,normalized = TRUE)
print(conf)


length(which(y == 1))/length(y)




calculatePredictionAccuracies <- function(distMatKernelKnnCV_out, trueLabels){
  
  correct = 0
  
  total = 0
  
  for(i in 1:length(distMatKernelKnnCV_out$folds)){
    print(distMatKernelKnnCV_out$folds[i])
    
    predictions = distMatKernelKnnCV_out$preds[i][[1]]
    
    indices = distMatKernelKnnCV_out$folds[[i]]
    
    for(j in 1:nrow(predictions)){
      # print(which.max(predictions[j,]))
      if(which.max(predictions[j,]) == trueLabels[indices[j]]) correct = correct + 1
      
      
      total = total + 1
    }

  }
  # distMatKernelKnnCV_out
  
  print(total)
  print(correct)
  
  print(paste("Acc: ", correct / total, sep = ""))
}


calculatePredictionAccuracies(out, y)



?train.kknn


train.kknn()







data(miete)
(train.con <- train.kknn(nmqm ~ wfl + bjkat + zh, data = miete, 
                         kmax = 25, kernel = c("rectangular", "triangular", "epanechnikov",
                                               "gaussian", "rank", "optimal")))
plot(train.con)









