
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
                             labels = parameters$labelsPath,
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

# if(mode == "onlyExperiments3"){
#   p2 = strsplit(pathToExperiment, "/Output/")[[1]][1]
#   print(p2)
#   NNexperimentsKfoldDir = paste(p2, "/NNexperimentsKfoldCV/", sep = "")
#   
#   df_summary = getExperimentSummary(expDir = NNexperimentsKfoldDir)
#   
#   
#   # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#   # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#   # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.05_m_1_q_1_muNN_10_alpha_0_betha_0_loc_TRUE.csv",
#   # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.5_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#   # "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv"
#   
#   
#   conv = FALSE
#   SAMPLESIZE = 20
#   modelParameters = list("layers" = c(500,100,50,30), "dropOuts" = c(0.2,0.1,0.1,0.1), "metrics" = "accuracy", "optimizerFunName" = "optimizer_adam", "batch_size" = 32, "epochs" = 20)
#   ProteinsExperimentKfoldCV( sampleSize = SAMPLESIZE,
#                              sampleTimes = 400,
#                              sampleTimes_test = 10,
#                              batch_size = 32,
#                              epochs = 30,
#                              euklid = "both",
#                              potentials = c("pos", "neg", "pos_neg"),
#                              q = 1,
#                              m = 1000,
#                              numClasses = 2,
#                              fNameTrain = c(    "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.1_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                                 "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.2_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                                 "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.05_m_1_q_1_muNN_10_alpha_0_betha_0_loc_TRUE.csv",
#                                                 "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.5_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv",
#                                                 "/home/willy/PredictingProteinInteractions/data/120Experiment/Quantiles/All_n_0.8_m_1_q_1_muNN_10_alpha_3_betha_3_loc_TRUE.csv"),
#                              ExperimentName = "Test102",
#                              modelParameters = modelParameters,
#                              recalculate = FALSE,
#                              k = 10,
#                              onlySummarizeFolds = FALSE,
#                              normalizeInputs = TRUE,
#                              saveExperiment = TRUE,
#                              useColIndsToKeep = TRUE,
#                              path = NNexperimentsKfoldDir,
#                              labels = getPath("120ExperimentLabels"),
#                              pre_training = FALSE,
#                              numPermutations = 200,
#                              fps = FALSE,
#                              sampleFunction = 2
#   )
#   
#   
#   df_summary = getExperimentSummary(expDir = NNexperimentsKfoldDir)
# }


# if(mode == "randomSearch"){
#   p2 = strsplit(pathToExperiment, "/Output/")[[1]][1]
#   print(p2)
#   NNexperimentsKfoldDir = paste(p2, "/NNexperimentsKfoldCV/", sep = "")
# 
# 
#   quantilePath = paste(p2, "/Quantiles/", sep = "")
#   randomSearch2(modelParameters = modelParameters,
#                NNexperimentsKfoldDir = NNexperimentsKfoldDir,
#                ExperimentName = ExperimentName,
#                quantilePath = quantilePath,
#                nlocals = c(0.05,0.1,0.2,0.3,0.8,0.5),
#                Trials = 100000)
# }


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