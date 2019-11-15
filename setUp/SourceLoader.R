# Don't change this line! 
# Automatically generated!
PPISETUP = "/home/willy/Schreibtisch/PPITEST2/PredictingProteinInteractions/setUp"

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

getPath <- function(name){
  t = read.table(file = paste(PPISETUP,"/Paths.txt", sep = ""), header = TRUE)
  
  s = as.character(t[which(t[,1] == name),2])
  
  return(s)
}


