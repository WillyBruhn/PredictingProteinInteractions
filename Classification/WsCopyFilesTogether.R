extractPdbs <- function(extractInds, protDirs){
  pdbsExtracted = c()
  
  for(i in 1:length(protDirs)){
    list.dirs(paste(protDirs[i],"/Output/", sep = ""), recursive = FALSE)
    files = list.files(paste(protDirs[i],"/Input/pdb/", sep = ""), recursive = FALSE)
    files_full = list.files(paste(protDirs[i],"/Input/pdb/", sep = ""), recursive = FALSE, full.names = TRUE)
    
    files[order(files, decreasing = FALSE)]
    
    for(j in 1:length(files)){
      
      if(length(files[j] > 0)){
        print(files[j])
        s = strsplit(files[j], split = "_")
        
        # print(length(s))
        s2 = s[[1]][length(s[[1]])]
        s3 = strsplit(s2, split = ".pdb")
        num = as.numeric(s3[[1]])
        
        if(num %in% extractInds){
          pdbsExtracted = c(pdbsExtracted, files_full[j])
        }
      }
    }
  }
  
  return(pdbsExtracted)
}

copyExtractedPdbs <- function(toDir = "/home/sysgen/server/projects/md-simulations/human_start2/", pdbsExtracted){
  startDir = toDir
  if(!dir.exists(startDir)) dir.create(startDir)
  
  for(i in 1:length(pdbsExtracted)){
    v = strsplit(pdbsExtracted[i], "/")[[1]]
    name = v[length(v)]
    finalDestination = paste(startDir, "/pdb/", name, sep = "")
    file.copy(from = pdbsExtracted[i], to = finalDestination)
  }
}

#--------------------------------------------------------------------------
humanDir = "/home/sysgen/server/projects/md-simulations/human_redoxins/"
protDirs = list.dirs(humanDir, recursive = FALSE)

pdbsExtracted = extractPdbs(extractInds = c(0), protDirs)
copyExtractedPdbs(toDir = "/home/sysgen/server/projects/md-simulations/human_start/", pdbsExtracted = pdbsExtracted)

pdbsExtracted = extractPdbs(extractInds = c(250), protDirs)
copyExtractedPdbs(toDir = "/home/sysgen/server/projects/md-simulations/human_end/", pdbsExtracted = pdbsExtracted)


# folderNames = c()
# for(i in 1:length(pdbsExtracted)){
#   s = unlist(strsplit(pdbsExtracted[i], "/Input/"))[1]
#   s2 = strsplit(s, "/human_redoxins//")[[1]][2]
#   folderNames = c(folderNames,s2)
# }
# 
# 
# folderNames
# protDirsShort= list.dirs(humanDir, recursive = FALSE, full.names = FALSE)
# protDirsShort
# 
# protDirsShort[which(!(protDirsShort %in% folderNames))]

