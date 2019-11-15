#!/usr/bin/Rscript
# author: Willy Bruhn, 14.11.2019
#------------------------------------------------------------
#
# Sets the paths correctly in the DummyParameters.txt
#
#------------------------------------------------------------

changeParameterTo <- function(MutCompParameters,Parameter,Value){
  MutCompParameters[which(MutCompParameters$Parameter == Parameter),2] = Value
  return(MutCompParameters)
  
}


path = funr::get_script_path()
vec = strsplit(path, "/")

PredictingProteinInteractionsFolder =  paste(paste(vec[[1]][1:(length(vec[[1]])-2)], collapse ="/"), "/", sep = "")
print(PredictingProteinInteractionsFolder)


dummyParametersFile = paste(PredictingProteinInteractionsFolder,"/Classification/predictingProteinInteractionsSettings/MutCompParametersDummy.txt", sep = "")
MutCompParameters = read.csv2(dummyParametersFile, header = FALSE)
colnames(MutCompParameters) = c("Parameter", "Value")

MutCompParameters$Value = as.character(MutCompParameters$Value)


# just for manual reference
parametersToBeAdjusted = c("MutCompApbsReference",
                           "MutCompError",
                           "MutCompInputFolder",
                           "MutCompLog",
                           "MutCompOutputFolder",
                           "OutputFolder",
                           "PathToMutComp",
                           "PathToVMD",
                           "vmdRC")

MutCompParameters = changeParameterTo(MutCompParameters = MutCompParameters,
                                      Parameter = "MutCompApbsReference",
                                      paste(PredictingProteinInteractionsFolder, "/PreProcessingProteins/MutComp/APBS/apbsReference.in", sep = ""))

MutCompParameters = changeParameterTo(MutCompParameters = MutCompParameters,
                                      Parameter = "MutCompError",
                                      paste(PredictingProteinInteractionsFolder, "/QuickStart/SmallExample/mutComp.err", sep = ""))

MutCompParameters = changeParameterTo(MutCompParameters = MutCompParameters,
                                      Parameter = "MutCompInputFolder",
                                      paste(PredictingProteinInteractionsFolder, "/QuickStart/SmallExample/", sep = ""))

MutCompParameters = changeParameterTo(MutCompParameters = MutCompParameters,
                                      Parameter = "MutCompLog",
                                      paste(PredictingProteinInteractionsFolder, "/QuickStart/SmallExample/mutComp.log", sep = ""))

MutCompParameters = changeParameterTo(MutCompParameters = MutCompParameters,
                                      Parameter = "MutCompOutputFolder",
                                      paste(PredictingProteinInteractionsFolder, "/QuickStart/SmallExample/Train/", sep = ""))

MutCompParameters = changeParameterTo(MutCompParameters = MutCompParameters,
                                      Parameter = "OutputFolder",
                                      paste(PredictingProteinInteractionsFolder, "/PreProcessingProteins/MutComp/GUI/Parameters/", sep = ""))

MutCompParameters = changeParameterTo(MutCompParameters = MutCompParameters,
                                      Parameter = "PathToMutComp",
                                      paste(PredictingProteinInteractionsFolder, "/PreProcessingProteins/MutComp/", sep = ""))

# write.csv2(x = MutCompParameters, file = dummyParametersFile, row.names = FALSE, quote = FALSE)
write.table(x = MutCompParameters, file = dummyParametersFile, sep=";",  col.names=FALSE, row.names = FALSE, quote = FALSE)
print("Adjusted MutCompParameters.txt")

