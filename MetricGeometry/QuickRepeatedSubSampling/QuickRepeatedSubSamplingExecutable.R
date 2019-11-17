#!/usr/bin/Rscript
#---------------------------------------------------------------------------------
# Calculate the repeatedSubSampling Fast. The points are only sampled m times once for each 
# protein. Then the Distributions are calculated. These distributions are then sampled from.
# Willy Bruhn 12.7.2019
# 
#   ProteinsPath        ... path to all proteins as produced by MutComp
#
#   distance_name       ... folder in which all distance-matrices will be stored
#
#   n                   ... number of points to select (see parameters of RepeatedSampling)
#
#   m                   ... sqrt(number of repeatitions)
#
#   cores               ... number of cores to run on
#
# ---------------------------------------------------------------------------------
#  Load the matrix with readMatrix1().
#---------------------------------------------------------------------------------
library(getopt)

options(warn=-1)

#----------------------------------------------------------------------------------
# Input
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'ProteinsPath'  , 'p', 2, "character",
  'distance_name'   , 'r', 2, "character",
  'n'   , 'n', 2, "integer",
  'm'   , 'm', 2, "integer",
  'cores', 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
if ( is.null(opt$ProteinsPath    ) ) { opt$ProteinsPath    = "/home/willy/Schreibtisch/106Test/Output/"     }
if ( is.null(opt$n    ) ) { opt$n    = 100   }
if ( is.null(opt$m    ) ) { opt$m    = 3  }
if ( is.null(opt$cores    ) ) { 
    opt$cores = Sys.getenv('LSB_MAX_NUM_PROCESSORS')
    if ( is.null(opt$cores)) opt$cores = 6
  }
if ( is.null(opt$outPutRepSubSamp    ) ) { opt$distance_name    = paste("quickEmd_n_",opt$n,"_m_",opt$m,".csv",sep ="")    }
if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }

SourceFunctions<-function(file) {
  MyEnv<-new.env()
  source(file=file,local=MyEnv)
  list2env(Filter(f=is.function,x=as.list(MyEnv)),
           envir=parent.env(environment()))
}

# SourceFunctions("/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling.R")
s1 = paste(funr::get_script_path(),"/QuickRepeatedSubSampling.R",sep ="")
# SourceFunctions(s1)
source(s1)
registerDoMC(cores= opt$cores)


packagesLoadedFrom("QuickRepeatedSubSamplingExecutable.R")

#---------------------------------------------------------------------------------------------------------------------
n = opt$n
m = opt$m
times = 1
ProteinsPath = opt$ProteinsPath
distance_name = opt$distance_name
distance_name = paste("quickEmd_n_",n,"_m_",m,".csv",sep ="")
all_protein_models_with_distances = not_vectorized_get_allModels(OutputPath = ProteinsPath, n = n,m = m, times = times)
DE_parallel = computeAllDistancesParallel(all_protein_models_with_distances)
Emd_distance_matrix = emd_parallel(DE_parallel = DE_parallel,m = m)
write.csv(Emd_distance_matrix,file = paste(ProteinsPath,"/",distance_name,sep =""))
#---------------------------------------------------------------------------------------------------------------------

# source("/home/sysgen/Documents/LWB/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling.R")
# n = 100
# m = 22
# times = 1
# OutputPath = "/home/sysgen/Documents/LWB/PredictingProteinInteractions/data/106Model/Proteins/Output/"
# all_protein_models_with_distances = not_vectorized_get_allModels(OutputPath = OutputPath, n = n,m = m, times = times)
# DE_parallel = computeAllDistancesParallel(all_protein_models_with_distances)
# Emd_distance_matrix = emd_parallel(DE_parallel = DE_parallel,m = m)


