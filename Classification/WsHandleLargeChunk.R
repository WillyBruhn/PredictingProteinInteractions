

mode = "SingleDistance" 
n = 5
m = 5
q = 1
doMutComp = FALSE



path = "/home/sysgen/server/projects/md-simulations/human_redoxins/"

dirs = list.dirs(path,recursive = FALSE)
notDone = "/home/sysgen/server/projects/md-simulations/human_redoxins//Trxndc8_sep"

dirs = dirs[-which(dirs == notDone)]


print(paste("processing ", dirs))

dirs = dirs[2]

for(k in 1:length(dirs)){
  
  print(paste("processing ", dirs[k], " ", k/length(dirs)*100, " %", sep = ""))
  
  subDirs = list.dirs(dirs[k], recursive = FALSE)
  
  OutPutFolder =""
  pdbFolder = ""
  for(i in 1:length(subDirs)){
    # print(subDirs[i])
    s = strsplit(subDirs[i],"/")
    for(j in 1:length(s)){
      print(s[[j]][length(s[[j]])])
      
      if(s[[j]][length(s[[j]])] == "Output"){
        OutPutFolder = paste(subDirs[i],"/",sep="")
      }
      
      if(s[[j]][length(s[[j]])] == "pdb"){
        pdbFolder = paste(subDirs[i],"/",sep="")
      }
    }
  }
  
  MutCompOutPut = strsplit(OutPutFolder, "Output")[[1]][1]
  
  # MutCompOutPut = OutPutFolder
  distancesFolder = paste(MutCompOutPut, "/QRsampDistances/", sep ="")
  
  print(MutCompOutPut)
  print(distancesFolder)
  
  exec = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/Classification/./predictingProteinInteractions.R "
  fullCall = paste(exec,
                   " --pdb_folder ",pdbFolder,
                   " --MutCompOutPut ",MutCompOutPut,
                   " --mode ", mode,
                   " --doMutComp ", doMutComp ,
                   " --numberOfPoints ", n,
                   " --rounds ", m,
                   " --doCluster TRUE ",
                   " --distances_train ", distancesFolder,
                   " --q ", q, sep ="")
  
  system(fullCall)
  
}


# subDirs = list.dirs(dirs[19], recursive = FALSE)
# 
# OutPutFolder =""
# pdbFolder = ""
# distancesFolder = ""
# for(i in 1:length(subDirs)){
#   # print(subDirs[i])
#   s = strsplit(subDirs[i],"/")
#   for(j in 1:length(s)){
#     print(s[[j]][length(s[[j]])])
#     
#     if(s[[j]][length(s[[j]])] == "Output"){
#       OutPutFolder = paste(subDirs[i],"/",sep="")
#     }
#     
#     if(s[[j]][length(s[[j]])] == "pdb"){
#       pdbFolder = paste(subDirs[i],"/",sep="")
#     }
#   }
# }
# 
# 
# # pdbFolder = ""
# 
# 
# mode = "SingleDistance" 
# n = 10
# m = 5
# q = 2
# doMutComp = FALSE
# 
# 
# MutCompOutPut = strsplit(OutPutFolder, "Output")[[1]][1]
# distancesFolder = paste(MutCompOutPut, "/", "UQDistances", sep = "")
# exec = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/Classification/./predictingProteinInteractions.R "
# fullCall = paste(exec," --pdb_folder ",pdbFolder," --MutCompOutPut ", MutCompOutPut, " --mode ", mode," --doMutComp ", doMutComp ," --numberOfPoints ", n, " --rounds ", m, " --doCluster TRUE --q ", q, " --labels_train NOLABELS --distances_train ", distancesFolder, sep ="")
# 
# system(fullCall)



# t = read.csv("/home/sysgen/server/projects/md-simulations/human_redoxins//Trx1_sep//UQDistances/_quickEmd_n_10_m_5_q_2_geo.csv", header = TRUE, row.names = 1)
# 
# ?read.csv
# 
# nrow(t)
# ncol(t)

 
# t[1:5,1:5]
# t[1,1]


# m = matrix(seq(1:25), ncol = 5)
# m2 = matrix(5, ncol = 5, nrow = 5)
# 
# (m + m2)/2
# 
# 
# m
# m2
# matrix(mapply(m,m2,FUN = max),ncol = 5,nrow = 5)
# 
# 
# d = read.csv("/home/sysgen/server/projects/md-simulations/human_redoxins//Grx1_sep//QRsampDistances//_pos_quickEmd_n_5_m_5_q_1_geo.csv", row.names = 1)
# 
# nrow(d)
# ncol(d)
# 
# d2 = read.csv("/home/sysgen/server/projects/md-simulations/human_redoxins//Grx1_sep//QRsampDistances//_neg_quickEmd_n_5_m_5_q_1_geo.csv", row.names = 1)
# 
# matrix(mapply(as.matrix(d),as.matrix(d2),FUN = max),ncol = ncol(d),nrow = nrow(d))
# 



