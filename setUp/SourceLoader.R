PPISETUP = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/setUp"

sourceFiles <- function(vector){
  t = read.table(file = paste(PPISETUP,"/SourcableFiles.txt", sep = ""), header = TRUE)
  
  # print(t)
  
  for(i in 1:length(vector)){
    s = as.character(t[which(t[,1] == vector[i]),2])
    print(paste("sourcing ", s))
    source(s)
  }
}

printPathsToSources <- function(vector){
  t = read.table(file = paste(PPISETUP,"/SourcableFiles.txt", sep = ""), header = TRUE)
  
  for(i in 1:length(vector)){
    s = as.character(t[which(t[,1] == vector[i]),2])
    print(s)
  }
}
