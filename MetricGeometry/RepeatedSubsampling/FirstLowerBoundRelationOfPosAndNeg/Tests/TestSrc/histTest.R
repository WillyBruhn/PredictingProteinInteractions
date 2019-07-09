#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Missing arg", call.=FALSE)
}

fName1 = args[1]
fName2 = args[2]

fName1_hist = args[3]
fName2_hist = args[4]

bin_num = as.numeric(args[5])

# emd_passed = as.numeric(args[6])
# fName1 = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/Q_vector.csv"
# fName2 = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/P_vector.csv"
# fName1_hist = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/Q_hist.csv"
# fName2_hist = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/P_hist.csv"
# bin_num = 2


t = read.csv(fName1, sep = ";", header = FALSE)
t[,1] <- as.numeric(t[,1])

t2 = read.csv(fName2, sep = ";", header = FALSE)
t2[,1] <- as.numeric(t2[,1])

bin_width = (max(t2,t)-min(t2,t))/(bin_num)
breaks = seq(min(t2,t), max(t2,t), bin_width)

#?hist

h1 =hist(t[,1], breaks = breaks, plot = FALSE, right = FALSE)
h2 = hist(t2[,1], breaks = breaks, plot = FALSE, right = FALSE)


h1_should = read.csv(fName1_hist, sep = ";", header = FALSE)
h1_should[,1] <- as.numeric(h1_should[,1])

h2_should = read.csv(fName2_hist, sep = ";", header = FALSE)
h2_should[,1] <- as.numeric(h2_should[,1])



# print("c++ R")

flag = TRUE
for(i in 1:length(h1$counts)){
  if(h1_should[i,1] != h1$counts[i]) flag = FALSE
  
  # print(paste(h1_should[i,1], " ", h1$counts[i]))
}
print(flag)


# library("emdist")
# 
# h1 = t(as.matrix(h1$counts))
# h2 = t(as.matrix(h2$counts))
# emd = emd2d(h1,h2)
# 
# # emd
# print(emd_passed)
# # print(emd)







