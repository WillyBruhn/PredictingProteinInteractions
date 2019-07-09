#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Missing arg", call.=FALSE)
}

fName1_hist = args[1]
fName2_hist = args[2]

emd_passed = as.numeric(args[3])

# emd_passed = as.numeric(args[6])
# fName1 = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/Q_vector.csv"
# fName2 = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/P_vector.csv"
fName1_hist = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/Q_hist.csv"
fName2_hist = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/P_hist.csv"
# bin_num = 2

h1_should = read.csv(fName1_hist, sep = ";", header = FALSE)
h1_should[,1] <- as.numeric(h1_should[,1])

h2_should = read.csv(fName2_hist, sep = ";", header = FALSE)
h2_should[,1] <- as.numeric(h2_should[,1])


dig <- function(number){
  sprintf(number, fmt = '%#.4f')
}


library("emdist")

h1 = t(as.matrix(h1_should))
h2 = t(as.matrix(h2_should))

# print(h1)
# print(h2)


emd = emd2d(h1,h2)

# emd
# print(emd_passed)
# print(dig(emd))


print(dig(emd))

# flag = (dig(emd) == emd_passed)
# 
# if(flag == FALSE){
#   print(dig(emd))
# } else {
#   print(flag)
# }




