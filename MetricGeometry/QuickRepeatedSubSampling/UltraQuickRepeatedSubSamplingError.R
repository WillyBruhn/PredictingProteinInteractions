#-----------------------------------
# willy Bruhn, 16.7.19
#
# Check the error of the approximation
#-------------------------------------

s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s1)

s2 = "/home/willy/PredictingProteinInteractions/Classification/NNClassification/optimizeDifferentModels/BoostedKNN.R"
 

getCorrelationsOfModels <- function(models){
  
  
  mat = matrix(0,nrow = length(models), ncol = length(models))

  print(mat)
  names = rep("", length(models))
  for(i in 1:length(models)){
    names[i] = models[[i]]$name
    
    for(j in 1:length(models)){
      print(paste(i,j))
      mat[i,j] = cor(c(models[[i]]$distances),c(models[[j]]$distances))
      mat[j,i] = mat[i,j]
    }
  }

  # colnames(mat) = names
  rownames(mat) = names
  return(mat)
}

#----------------------------------------
# readDistanceMatrix1
# readDistanceMatrix2

emdRef = readDistanceMatrix2(file = "/home/willy/Schreibtisch/106Test/CompareProteinsOutput/ListofEMD_positive_100.csv")
emdQuick1 = readDistanceMatrix1("/home/willy/Schreibtisch/106Test/UltraQuickRepSubSamp/run1_quickEmd_n_10_m_100_q_2_geo.csv")
emdQuick2 = readDistanceMatrix1("/home/willy/Schreibtisch/106Test/UltraQuickRepSubSamp/run2_quickEmd_n_10_m_100_q_2_geo.csv")
emdQuick3 = readDistanceMatrix1("/home/willy/Schreibtisch/106Test/UltraQuickRepSubSamp/run1_quickEmd_n_100_m_100_q_99_geo.csv")
emdQuick4 = readDistanceMatrix1("/home/willy/Schreibtisch/106Test/UltraQuickRepSubSamp/run1_quickEmd_n_100_m_100_q_2_geo.csv")

models = list()
models[[1]] = list("name" = "CompareProteins", "distances" = emdRef)
models[[2]] = list("name" = "run1_quickEmd_n_10_m_100_q_2_geo", "distances" = emdQuick1[-106,-106])
models[[3]] = list("name" = "run2_quickEmd_n_10_m_100_q_2_geo", "distances" = emdQuick2[-106,-106])
models[[4]] = list("name" = "run1_quickEmd_n_100_m_100_q_99_geo", "distances" = emdQuick3[-106,-106])
models[[5]] = list("name" = "run1_quickEmd_n_100_m_100_q_2_geo", "distances" = emdQuick4[-106,-106])

corMat = getCorrelationsOfModels(models)

library(xtable)
corMatLatex = xtable(corMat, caption = "Correlation of different methods",label = "correlationCompareProteins")
print.xtable(corMatLatex, type="latex", file="/home/willy/PredictingProteinInteractions/Results/QuickRepSampling/Correlation1.tex")

#------------------------------------------------------------------------
# pretty plot of 2 distributions
OutputPath = "/home/willy/Schreibtisch/106Test/Output/"
allApprox = getAll_protein_F_approximations(OutputPath = OutputPath,q = 2,m = 100, n = 100)


small2 = list(allApprox[[1]],allApprox[[5]],allApprox[[6]])
proj = getManhattanProjection(all_protein_F_approximations = small2)

pdf("/home/willy/PredictingProteinInteractions/Results/QuickRepSampling/2dPlots/smallExample1.pdf")
plot_all_F_approximations(all_protein_F_approximations = small2,onlyGeo = FALSE, 
                          col = c("blue","green","red"),
                          xli = c(18,28),
                          yli = c(0,10))
dev.off()
#------------------------------------------------------------------------
# emd vs geo
n1 = 100
n2 = 100
mat1 = matrix(rnorm(n1*2), ncol = 2)
mat2 = matrix(rnorm(n1*2,mean = 5), ncol = 2)

mat_both = rbind(mat1,mat2)
plot(mat_both)

d1 = dist(mat1,method = "manhattan")
d2 = dist(mat2,method = "manhattan")

d1_2 = as.matrix(dist(mat_both, method= "manhattan"))[1:n1,(n1+1):(n1+n2)]

breaks = seq(from = 0,to = max(d1,d1_2)+0.05, by = 0.05) 
P <- t(as.matrix(hist(d1 , breaks = breaks, plot = F)$counts))
Q <- t(as.matrix(hist(d1_2 , breaks = breaks, plot = F)$counts))

hist(P, col ="red")
hist(Q, col ="blue", add = TRUE)

emdDist = emd(P,Q)

geo1 = c(sum(mat1[,1])/n1,sum(mat1[,2])/n1) 
geo2 = c(sum(mat2[,2])/n2,sum(mat2[,2])/n2) 

geoDist = dist(rbind(geo1,geo2), method = "manhattan")


print(paste(emdDist, geoDist))






