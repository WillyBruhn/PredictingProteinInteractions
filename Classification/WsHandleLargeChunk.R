

mode = "SingleDistance" 
n = 5
m = 5
q = 2
doMutComp = FALSE



path = "/home/sysgen/server/projects/md-simulations/human_redoxins/"

dirs = list.dirs(path,recursive = FALSE)
notDone = "/home/sysgen/server/projects/md-simulations/human_redoxins//Trxndc8_sep"

dirs = dirs[-which(dirs == notDone)]


print(paste("processing ", dirs))

dirs = dirs[1:3]

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
  
  write.table(fullCall,paste(MutCompOutPut,"/fullCall.txt", sep = ""),quote = FALSE,col.names = FALSE, row.names = FALSE)
  
}


# combine the summaries from all folders
dirs

# name + statistics = 9 values
summariesTable = data.frame(matrix(0,ncol = 8+5, nrow = 0))
colnames(summariesTable) = c("name", "n", "m", "q", "method", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "mean", "var")

for(k in 1:length(dirs)){
  summaries = list.files(paste(dirs[k], "/QRsampDistances/Summaries/" , sep =""), full.names = TRUE)

  for(i in 1:length(summaries)){
    
    df = data.frame(matrix(0,ncol = 8+5))
    colnames(df) = c("name", "n", "m", "q", "method", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "mean", "var")
    
    
    temp = strsplit(summaries[i],split = ".txt")[[1]]
    vec = strsplit(temp,split = "_")[[1]]
    v2 = vec[(length(vec)-7):length(vec)]
    
    method = v2[1]
    n = as.numeric(v2[4])
    m = as.numeric(v2[6])
    q = as.numeric(v2[8])
    
    sum = read.table(summaries[i],header = TRUE)
    
    
    s1 = strsplit(temp, "/QRsampDistances")[[1]]
    s2 = strsplit(s1, "/")[[1]]
    name = s2[length(s2)]
    
    df[1,1] = name
    df[1,2] = n
    df[1,3] = m
    df[1,4] = q
    df[1,5] = method
    df[1,6:ncol(df)] = sum
    
    print(df)
    
    summariesTable = rbind(summariesTable,df)
  }
  
}

write.table(summariesTable, paste(path, "/summariesAllTogether.txt", sep = ""), row.names = FALSE)







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
d = read.csv("/home/sysgen/server/projects/md-simulations/human_redoxins//Grx1_sep//QRsampDistances//_pos_quickEmd_n_5_m_5_q_1_geo.csv", row.names = 1)

nrow(d)
ncol(d)

d2 = read.csv("/home/sysgen/server/projects/md-simulations/human_redoxins//Grx1_sep//QRsampDistances//_neg_quickEmd_n_5_m_5_q_1_geo.csv", row.names = 1)

mea = matrix(lapply(as.vector(d),as.vector(d2),FUN = mean),ncol = ncol(d),nrow = nrow(d))
mea[1:5,1:5]


mea = (d+d2)/2



d[1:5,1:5]


d2[1:5,1:5]

mea[1:5,1:5]


