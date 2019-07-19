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