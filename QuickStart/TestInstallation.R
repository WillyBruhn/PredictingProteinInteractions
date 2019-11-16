#!/usr/bin/Rscript
# author: Willy Bruhn, 15.11.2019
#------------------------------------------------------------
#
# Performs the tasks in Quick-Start. Reports for each step 
# if it executed correctly. Reports an error if it fails.
#
# See https://github.com/WillyBruhn/PredictingProteinInteractions
# section QuickStart for a detailed description.
#
#------------------------------------------------------------

checkIfFileExists <- function(file, errorTopic,PredictingProteinInteractionsFolder){
  if(!file.exists(paste(PredictingProteinInteractionsFolder,"/", file, sep =""))){
    print(paste("ERROR: ", errorTopic, "not executed correctly. Missing file", paste(PredictingProteinInteractionsFolder,"/", file, sep ="")))
    return(FALSE)
  }
  
  return(TRUE)
}

path = funr::get_script_path()
vec = strsplit(path, "/")

PredictingProteinInteractionsFolder =  paste(paste(vec[[1]][1:(length(vec[[1]])-1)], collapse ="/"), "/", sep = "")
print("Testing Installation ...")
print(paste("Checking installation in ", PredictingProteinInteractionsFolder))
print("Checking with QuickStart ...")

# just some flags for testing this script. Make Sure to set all to TRUE when shipping
PPIBASIC = TRUE
TRAINING = TRUE
PREDICTION = TRUE

IGNOREERRORS = FALSE
IGNORESTD = TRUE
INGORESTERR = TRUE

if(PPIBASIC){
  print("---------------------------------------------------------------------------------------------")
  print("Step 1: Testing PredictingProteinInteractions.R ...")
  
  
  # Testing PredictingProteinInteractions
  PPICALL = paste(PredictingProteinInteractionsFolder,"/Classification/./predictingProteinInteractions.R --mode SingleDistance --doClustering TRUE --pdb_folder ",
        PredictingProteinInteractionsFolder,"/QuickStart/SmallExample/pdb_train/ --distances_train ",
        PredictingProteinInteractionsFolder,"/QuickStart/SmallExample/Train/UltraQuickRepSub/ --numberOfPoints 4 --rounds 10 --MutCompOutPut ",
        PredictingProteinInteractionsFolder,"/QuickStart/SmallExample/Train/ --doMutComp TRUE --q_val 1 --labels_train ",
        PredictingProteinInteractionsFolder,"/data/labels106.txt", sep ="")
  
  print(PPICALL)
  system(PPICALL,ignore.stdout = IGNORESTD,ignore.stderr = INGORESTERR)
  
  # MutComp checks
  print("Checking MutComp ...")
  if(checkIfFileExists("/QuickStart/SmallExample/Train/Output/", "MutComp", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
  checkIfFileExists("/QuickStart/SmallExample/Train/Output/overview.png", "Imagemaker", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
  checkIfFileExists("/QuickStart/SmallExample/Train/Output/003/003.obj", "MutComp", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
  checkIfFileExists("/QuickStart/SmallExample/Train/Output/003/003_pot.dx", "MutComp", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder)){
    print("... looks good")
  } else {
    if(!IGNOREERRORS) quit()
  }
  
  # PredictingProteinInteractions checks
  print("Checking PredictingProteinInteractions.R ...")
  if(checkIfFileExists("/QuickStart/SmallExample/Train/UltraQuickRepSub/_neg_quickEmd_n_4_m_10_q_1_geo.csv", "Creating Distances", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
     checkIfFileExists("/QuickStart/SmallExample/Train/UltraQuickRepSub/_pos_quickEmd_n_4_m_10_q_1_geo.csv", "Creating Distances", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
     checkIfFileExists("/QuickStart/SmallExample/Train/UltraQuickRepSub/Dendrogramms/Dendrogram__avg_quickEmd_n_4_m_10_q_1.pdf", "Creating Dendrogram", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder)){
    print("... looks good")
  } else {
    if(!IGNOREERRORS) quit()
  }
}

if(TRAINING){
  print("---------------------------------------------------------------------------------------------")
  print("Step 2: Testing Proteins.R (Training) ...")
  
  # Testing Proteins.R
  TRAININGCALL = paste(PredictingProteinInteractionsFolder,"/MetricGeometry/QuickRepeatedSubSampling/./Proteins.R --pathToExperiment ",
                  PredictingProteinInteractionsFolder,"/QuickStart/SmallExample/Train/Output/ --mode evaluation --outPutFolder Test1 --useSmallExample", sep ="")
  
  print(TRAININGCALL)
  system(TRAININGCALL,ignore.stdout = IGNORESTD,ignore.stderr = INGORESTERR)
  
  # MutComp checks
  print("Checking Proteins.R ...")
  if(checkIfFileExists("/QuickStart/SmallExample/Train/Quantiles/All_n_0.01_m_1_q_1_muNN_10_alpha_1_betha_1_loc_TRUE.csv", "Quantile-Calculation", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
     checkIfFileExists("/QuickStart/SmallExample/Train/NNexperimentsKfoldCV/Test1/Accuracy.tex", "Training", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
     checkIfFileExists("/QuickStart/SmallExample/Train/NNexperimentsKfoldCV/Test1/nnModel.h5", "Model-Export", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder)){
    print("... looks good")
  } else {
    if(!IGNOREERRORS) quit()
  }
  
}

if(PREDICTION){
  print("---------------------------------------------------------------------------------------------")
  print("Step 3: Testing Proteins.R (Prediction) ...")
  
  # Testing Proteins.R
  PREPROCESSCALL = paste(PredictingProteinInteractionsFolder, "Classification/./predictingProteinInteractions.R --mode SingleDistance --doClustering FALSE --pdb_folder ",
                         PredictingProteinInteractionsFolder, "/QuickStart/SmallExample/pdb_predict/ --MutCompOutPut ",
                         PredictingProteinInteractionsFolder,"/QuickStart/SmallExample/Predict/ --doMutComp TRUE --q_val 1", sep ="")
  
  # PredictingProteinInteractionsFolder, "/Classification/./predictingProteinInteractions.R --mode SingleDistance --doClustering FALSE --pdb_folder ",
  # PredictingProteinInteractionsFolder, "/QuickStart/SmallExample/pdb_predict/ --MutCompOutPut ",
  # PredictingProteinInteractionsFolder, "/QuickStart/SmallExample/Predict/ --doMutComp TRUE --q_val 1
  
  
  print(PREPROCESSCALL)
  system(PREPROCESSCALL,ignore.stdout = IGNORESTD,ignore.stderr = INGORESTERR)
  
  # MutComp checks
  print("Checking MutComp ...")
  if(checkIfFileExists("/QuickStart/SmallExample/Predict/Output/", "MutComp", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
     checkIfFileExists("/QuickStart/SmallExample/Predict/Output/overview.png", "Imagemaker", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
     checkIfFileExists("/QuickStart/SmallExample/Predict/Output/013/013.obj", "MutComp", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder) &&
     checkIfFileExists("/QuickStart/SmallExample/Predict/Output/013/013_pot.dx", "MutComp", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder)){
    print("... looks good")
  } else {
    if(!IGNOREERRORS) quit()
  }
  
  PREDICTIONCALL = paste(PredictingProteinInteractionsFolder, "MetricGeometry/QuickRepeatedSubSampling/./Proteins.R --outPutFolder TestPred --pathToExperiment ",
                         PredictingProteinInteractionsFolder,"QuickStart/SmallExample/Predict/Output/ --mode prediction --nnModelFolder ",
                         PredictingProteinInteractionsFolder,"/QuickStart/SmallExample/Train/NNexperimentsKfoldCV/Test1/", sep ="")
  
  print(PREDICTIONCALL)
  system(PREDICTIONCALL)
  
  # MutComp checks
  print("Checking Proteins.R ...")
  if(checkIfFileExists("/QuickStart/SmallExample/Predict/NNexperimentsKfoldCV/TestPred/predictions.txt", "Predictions", PredictingProteinInteractionsFolder = PredictingProteinInteractionsFolder)){
    print("... looks good")
  } else {
    if(!IGNOREERRORS) quit()
  }
  
}

print("---------------------------------------------------------------------------------------------")
print("No errors detected. Installation seems fine.")