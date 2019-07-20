library(rfUtilities)
library(RSNNS)

generateFileName <- function(n,m,q,fName,path){
  return(paste(path, "/", fName,"_n_", n, "_m_", m, "_q_", q, ".csv", sep = ""))
}

writeProjectionToFile <- function(proj,n,m,q,path = "/home/willy/PredictingProteinInteractions/data/animals/", fName = "proj"){
  fName_final = generateFileName(n=n,m=m,q=q,fName = fName,path = path)
  print(paste("writing projection to ",fName_final, sep =""))
  
  write.csv(proj, file = fName_final, row.names = FALSE)
}

readProjectionFromFile <- function(n,m,q,path,fName){
  fName_final = generateFileName(n=n,m=m,q=q,fName = fName,path = path)
  print(paste("loading projection from ", fName_final, sep =""))
  
  return(read.csv(file = fName_final))
}

writeQuantilesToFileAnimal <- function(model_vec, path = "/home/willy/PredictingProteinInteractions/data/animals/", fName = "quantiles", n,m,q){
  fileName = generateFileName(n = n,m = m,q = q,fName = fName,path = path)  
  
  out = data.frame(matrix(0,nrow = length(model_vec)*m, ncol = q+2+1))
  colnames(out) = c("name", seq(1:(q+2)))
  
  for(i in 1:(length(model_vec))){
    # print(model_vec)
    ind = (i-1)*m
    
    for(j in 1:length(model_vec[[i]]$F$F_app_list)){
      out[ind+j,2:(q+3)] = model_vec[[i]]$F$F_app_list[[j]]
      out[ind+j,1] = model_vec[[i]]$name
      
      print(j/(length(model_vec[[i]]$F$F_app_list)*length(model_vec)))
    }
  }
  
  write.csv(x = out,file = fileName, row.names = FALSE)
}

readQuantilesFromFile <- function(path = "/home/willy/PredictingProteinInteractions/data/animals/", fName = "quantiles", n,m,q){
  fileName = generateFileName(n = n,m = m,q = q,fName = fName,path = path)  
  
  return(read.csv(file = fileName, header = TRUE))
}


getClassNamesFromSubClasses <- function(subClasses, splitPattern = "-"){
  # gets the className from a subclass
  # e.g.
  # Lion-01 --> Lion
  #-------------------------------------------------------
  
  classNames = rep("", length(subClasses))
  for(i in 1:length(subClasses)){
    classNames[i] = strsplit(as.character(subClasses[i]),split = splitPattern)[[1]][1]
  }
  
  return(classNames)
}
