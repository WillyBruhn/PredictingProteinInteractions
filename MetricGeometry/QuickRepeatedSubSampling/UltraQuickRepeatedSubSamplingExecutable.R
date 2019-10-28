#!/usr/bin/Rscript
#---------------------------------------------------------------------------------
# Willy Bruhn 15.7.2019
# Calculate the repeatedSubSampling fast. The points are only sampled m times once for each 
# protein. Then the distributions are calculated. Then from these distributions the quantile-
# approximation is calculated. The DE is then calculated with the manhattan distance between
# the quantiles.
# 
#---------------------------------------------------------------------------------
#
#   ProteinsPath        ... path to all proteins as produced by MutComp
#
#   distance_name       ... folder in which all distance-matrices will be stored
#
#   n                   ... number of points to select (see parameters of RepeatedSampling)
#
#   m                   ... sqrt(number of repeatitions)
#
#   q                   ... number of subdivisions of the integral. Basically a higher q
#                           leads to a more accurate approximation. Currently it has to hold
#                           q < n.
#
#   potential           ... pos/neg
#
#   distance_method     ... geo/emd
#
#
#   plot                ... in case of (q == 2) the approximations are ploted into the 2d-plane.
#                           Else no plot is created with a warn-message.
#
#
#   labels              ... a file containing for each protein a label functional/not_functional
#
#   cores               ... number of cores to run on. Trys to detect the number of cores
#                           automatically. If this fails the number of cores is set to 6
#                           by default.
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
  'distance_path'   , 'P', 2, "character",  
  'plot'   , 'o', 2, "logical",  
  'distance_method'   , 'U', 2, "character",  
  'n'   , 'n', 2, "integer",
  'm'   , 'm', 2, "integer",
  'q'   , 'q', 2, "integer",
  'potential'   , 'S', 2, "character",
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
if ( is.null(opt$ProteinsPath    ) ) { opt$ProteinsPath    = paste(funr::get_script_path(),"/../../data/106Model/Proteins/Output/", sep ="")     }
if ( is.null(opt$distance_path    ) ) { opt$distance_path    = paste(funr::get_script_path(),"/../../data/106Model/Proteins/UltraQuickRepSubSamp/", sep ="")     }
if ( is.null(opt$labels    ) ) { opt$labels    = paste(funr::get_script_path(),"/../../data/labels.txt", sep ="")  }
if ( is.null(opt$distance_method    ) ) { opt$distance_method    = "geo"   }
if ( is.null(opt$n    ) ) { opt$n    = 100   }
if ( is.null(opt$m    ) ) { opt$m    = 100  }
if ( is.null(opt$q    ) ) { opt$q    = 2  }
if ( is.null(opt$plot    ) ) { opt$plot    = TRUE  }
if ( is.null(opt$potential    ) ) { opt$potential    = "pos"  }
if ( is.null(opt$cores    ) ) { 
  opt$cores = Sys.getenv('LSB_MAX_NUM_PROCESSORS')
  if ( is.null(opt$cores)) opt$cores = 10
}

# list.files(opt$distance_path)
if ( is.null(opt$distance_name    ) ) { opt$distance_name    =  "run1"}

if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }

SourceFunctions<-function(file) {
  MyEnv<-new.env()
  source(file=file,local=MyEnv)
  list2env(Filter(f=is.function,x=as.list(MyEnv)),
           envir=parent.env(environment()))
}

# s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
s2 = paste(funr::get_script_path(),"/UltraQuickRepeatedSubSampling.R",sep ="")
print(paste("sourcing ", s2, " ..."))
source(s2)
print(paste(" ... done sourcing ", s2))
registerDoMC(cores= opt$cores)
#---------------------------------------------------------------------------------------------------------------------
# Samples the points,
# calculates the distributions,
# approximates them,
# calculate all pairwise distances,
# use the geometric center for the actual distance between the proteins,
# write the distance-matrix to "distTest"
#---------------------------------------------------------------------------------------------------------------------
if(!file.exists(opt$labels)){
  print(paste("Can't open ", opt$labels, sep = ""))
}else{
  labels = read.table(file = opt$labels, header = TRUE)
}

functionals = labels$name[which(labels$label == "functional")]

distName = paste(opt$distance_name,"_quickEmd_n_",opt$n,"_m_",opt$m,"_q_",opt$q,sep ="")

quickRepSampling(OutputPath = opt$ProteinsPath, 
                 distance_path = opt$distance_path,
                 n = opt$n,
                 m = opt$m,
                 q = opt$q,
                 pos = opt$potential,
                 fName = distName,
                 plot = opt$plot,
                 functionals = functionals,
                 distance_method = opt$distance_method)
