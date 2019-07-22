
#---------------------------------------------------------------------------------------------------------
# Willy Bruhn, 15.7.2019
# A distribution is projected into the manhatten-coordinate-system.
# The DE can then be computed very efficiently. Up to 1000 times faster than CompareProteins
#
#---------------------------------------------------------------------------------------------------------

# s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/QuickRepeatedSubSampling.R"
s1 = paste(funr::get_script_path(),"/QuickRepeatedSubSampling.R",sep ="")
print(paste("sourcing ", s1, " ... ", sep =""))
source(s1)

print(paste("... done sourcing ", s1, sep =""))
library(emdist)

#--------------------------------------------------------------------------------
# projectionMethod
#--------------------------------------------------------------------------------

# pts = read_pts_file(OutputPath = OutputPath,protName = "000_Trx")
# CDF = samplePointsAndCalculateCDFofEc(all_pts = pts,n = 100,plot = FALSE)
# 
# integralOfDistribution(CDF)
# x = integrate(CDF, lower = 0, upper = 100)
# x$value

integrateStepFun <- function(F_, a){

  nF = length(knots(F_))
    
  sup = 0
  ind = 0
  for(i in 1:nF){
    if(knots(F_)[i] <= a){
      sup = knots(F_)[i]
      ind = i
    } else{
      break
    }
  }
  
  su = 0
  
  # print(paste("ind",ind))
  # print(sup)
  
  if(ind > 0){
    if(ind > 1){
      for(i in 1:(ind-1)){
        # print(paste("i",i))
        # print(knots(F_)[i+1]-knots(F_)[i])
        su = su + (knots(F_)[i+1]-knots(F_)[i])*F_(knots(F_)[i])
      }
    }


      su = su + (a-knots(F_)[ind])*F_(knots(F_)[ind])
  }
  
  return(su)
}

# test
# integrateStepFun(CDF,300)
# integrateStepFun(CDF,14)


integrateStepFunInverse <- function(F_,b){
  # F_ ... a function
  # b ... integral_{0}^x (F(t) dt) = b, solve for x and return x
  #-------------------------------------------------------------
  
  val_old = 0
  for(i in seq(knots(F_)[1]-0.1,knots(F_)[length(knots(F_))]+0.1,0.1)){
    # print(i)
    val = integrateStepFun(F_,i)
    # print(val)
    if(val > b){
      return(i)
    }
    val_old = val
  }
  return(NULL)
}

approximateCDF <- function(F_, q, type = 1){
  # F_ ... cumulitive distribution function
  # q ... number of equal heights to approximate F_
  #-------------------------------------------------
  knots = knots(F_)
  maxF = max(knots)
  n_F = length(knots)
  
  # print(knots)
  
  
  quants = quantile(F_, probs = seq(0,1,1/(q+1)), type = 1)
  # points(x = quants, y = rep(0,length(quants)), col = "red")
  
  
  # F_app = rep(0,q)
  # # int = integrateStepFun(F_,maxF)
  # for(i in 1:q){
  #   ind = n_F*i/(q+1)
  # 
  #   if(ind - floor(ind) > 0){
  #     F_app[i] = (knots[n_F*i/(q+1)] + knots[(n_F*i+1)/(q+1)])/2
  #   } else{
  #     F_app[i] = knots[n_F*i/(q+1)]
  #   }
  
  # F_app = rep(0,length(quants))
  # 
  # for(i in 1:length(quants)){
  #   F_app[i] = 
  # }
    
    # print((i/(q+1))*(maxF))
    # F_app[i]  = integrateStepFunInverse(F_,v)
    # F_app[i] = integrate(CDF, lower = 0, upper = (i/q)*(maxF-0.1))$value
  # }
  
  return(quants)  
}

# pts = read_pts_file(OutputPath = OutputPath,protName = "000_Trx")
# CDF = samplePointsAndCalculateCDFofEc(all_pts = pts,n = 100,plot = TRUE)
# 
# # ?quantile
# q = 1
# quants = quantile(CDF, probs = seq(0,1,1/(q+1)))
# points(x = quants, y = rep(0,length(quants)), col = "red")
# 
# CDF_approx = approximateCDF(CDF,q = 10, type = 1)
# plotDistWithApproximation(CDF,CDF_approx)
# 
# DifferenceOfIntegral_F_and_Approx(CDF,CDF_approx)


plotDistWithApproximation <- function(F_,F_app,plotToFile = NULL){
  
  if(!is.null(plotToFile)){
    pdf(file = plotToFile)
  }
  
  n = length(knots(F_))
  q = length(F_app)# 1st and last quantile -2???
  sfun = stepfun(F_app , y = c(0,seq(1/q,1,1/q)))
  
  print(knots(F_)[n/2])
  
  
  plot(F_, main = paste("n = ",n , ", q= ", q, sep = ""), xlim = c(0,knots(F_)[n/2]*2))
  plot(sfun, add = TRUE)
  if(!is.null(plotToFile)){
    dev.off()
  }  
}

plotDistWithApproximationDifferentQ <- function(OutputPath, protName, n = 100, qlist = c(1:9), pos =TRUE, plotToFile = NULL){
  
  pts = read_pts_file(OutputPath = OutputPath,protName = protName,pos = pos)
  pos13_F1 = samplePointsAndCalculateCDFofEc(all_pts = pts,n = n,plot = FALSE)
  
  if(!is.null(plotToFile)){
    pdf(file = plotToFile)
  }
  par(mfrow=c(sqrt(length(qlist)),sqrt(length(qlist))))
  
  for(i in 1:length(qlist)){
    F_app1 = approximateCDF(pos13_F1,qlist[i])
    plotDistWithApproximation(pos13_F1,F_app1)
  }
  
  if(!is.null(plotToFile)){
    dev.off()
  }  
  
}

DifferenceOfIntegral_F_and_Approx <- function(F_, F_app){
  q = length(F_app)
  sfun = stepfun(F_app , y = c(0,seq(1/q,1,1/q)))
  
  DifferenceOfIntegral(F_,sfun)
}

manhattan_dist <- function(rating1, rating2){
  distance <- abs(rating1-rating2)
  distance <- sum(distance)
  return(distance)
}

manhattan_dist_distribution_approx <- function(rating1, rating2){
  return(manhattan_dist(rating1,rating2)/(length(rating1)))
}

plot2d_dist_approx <-function(pos13_F_approx_list, name=NULL, col = "red", xli = NULL, yli = NULL, add = FALSE, onlyGeo = FALSE){
  
  if(add == FALSE){
    x_val = pos13_F_approx_list[[1]][1]
    y_val = pos13_F_approx_list[[1]][2]
    
    if(!is.null(xli) && !is.null(yli)) {
      plot(x = x_val, y = y_val-x_val, xlim = xli, ylim = yli,  col = col)
    } else {
      plot(x = x_val, y = y_val-x_val, col = col)
    }
  }
  
  
  x_vals = rep(0,length(pos13_F_approx_list))
  y_vals = rep(0,length(pos13_F_approx_list))
  
  for(i in 1:length(pos13_F_approx_list)){
    x_val = pos13_F_approx_list[[i]][1]
    y_val = pos13_F_approx_list[[i]][2]
    
    if(onlyGeo == FALSE) points(x = x_val, y = y_val-x_val, col = col, pch = 19)
    
    
    x_vals[i] = x_val
    y_vals[i] = y_val-x_val
  }
  
  # plot the name
  if(!is.null(name)){
    geo = c(sum(x_vals)/length(x_vals), sum(y_vals)/length(y_vals))
    
    # print(geo)
    if(onlyGeo == TRUE) {
      text(x = geo[1], y = geo[2], labels = c(name), cex = 2, col = col)
    } else {
      text(x = geo[1], y = geo[2], labels = c(name), cex = 2, col = "black")
    }
    
  }
}

generateF_approximations_3dModelWithMetric <- function(d, n = 100, m = 10, q = 2, pos =TRUE){
  
  # print(nrow(model_points))
  
  pos13_F_list = list()
  pos13_F_approx_list = list()
  
  for(i in 1:m){
    pos13_F_list[[i]] = sampleDistancesAndCalculateCDFofEcWith(d, n = n,plot = FALSE)
    pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
  }
  
  return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
}


generateF_approximations <- function(OutputPath, protName, n = 100, m = 10, q = 2, pos =TRUE){
  
  pos13 = read_pts_file(OutputPath = OutputPath,protName = protName, pos = pos)
  
  pos13_F_list = list()
  pos13_F_approx_list = list()
  
  for(i in 1:m){
    pos13_F_list[[i]] = samplePointsAndCalculateCDFofEc(all_pts = pos13,n = n,plot = FALSE)
    pos13_F_approx_list[[i]] = approximateCDF(pos13_F_list[[i]],q)
  }
  
  return(list("F_list" = pos13_F_list, "F_app_list" = pos13_F_approx_list))
}

getDistanceMatrixFromProjection <- function(df){
  
  # name, n
  q= ncol(df)-2
  
  distances = as.matrix(dist(x = df[,3:q], diag = TRUE, upper = TRUE),method ="manhattan")
  colnames(distances) = df[,1]
  rownames(distances) = df[,1]
  
  return(distances)
}

getGeometricCenters <- function(df,n_){
  
  names = unique(df[,1])
  
  # name, n
  q = ncol(df) -2
  
  geos = data.frame(matrix(0, nrow = length(names), ncol = ncol(df)))
  
  na = c("",q)
  for(i in 1:q){
    na[i] = paste("q_",i,sep="")
  }
  
  colnames(geos) = c("name","n", na)
  
  for(i in 1:length(names)){

    df2 = subset(df, n == n_ & name == names[i],select = c(3:ncol(df)))
    
    geos[i,1] = names[i]
    geos[i,2] = n_
    
    for(j in 1:q){
      geos[i,j+2] = sum(df2[,j])/length(df2[,j])
    }
  }
  
  return(geos)
}
#-------------------------------------------------------------------
## Example

calcError <- function(all_pts, n = 100, q = 2, plot = TRUE, verbose = TRUE, plotToFile = NULL){
  # q = 2
  # n = 100
  pos13_F1 = samplePointsAndCalculateCDFofEc(all_pts = all_pts,n = n,plot = FALSE)
  F_app1 = approximateCDF(pos13_F1,q)
  plotDistWithApproximation(pos13_F1,F_app1,plotToFile)
  DifferenceOfIntegral_F_and_Approx(pos13_F1,F_app1)
  
  pos13_F2 = samplePointsAndCalculateCDFofEc(all_pts = all_pts,n = n,plot = FALSE)
  F_app2 = approximateCDF(pos13_F2,q)
  plotDistWithApproximation(pos13_F2,F_app2,plotToFile)
  DifferenceOfIntegral_F_and_Approx(pos13_F2,F_app2)
  
  DE_exact = DifferenceOfIntegral(pos13_F1,pos13_F2)
  DE_approx = manhattan_dist_distribution_approx(F_app1,F_app2)
  
  if(verbose){
    print(paste("DE_exact = ", DE_exact, ", DE_approx = ", DE_approx, sep =""))
    print(paste("Error =", abs(DE_exact- DE_approx)))
  }
  
  return(abs(DE_exact- DE_approx))
}


calcErrors <- function(all_pts, times = 100, n = 100, q = 2, plot = FALSE){
  errors = rep(0,times)
  for(i in 1:times){
    errors[i] = calcError(all_pts, n, q, plot = plot, verbose = FALSE)
  }
  return(errors)
}

caclErrorCurve <- function(all_pts, times = 100, n = 100, maxQ = 99, plot = FALSE){
  meanErrors = rep(0,maxQ)
  for(i in 1:maxQ){
    print(i)
    meanErrors[i] = mean(calcErrors(pts, times = times, n = n, q = i))
  }
  
  return(meanErrors)
}

# pts = read_pts_file(OutputPath = OutputPath, protName = "000_Trx")
# 
# calcError(pts, n = 100, q = 1, plot = TRUE, verbose = TRUE, 
#           plotToFile = "/home/willy/PredictingProteinInteractions/Results/QuickRepSampling/Examples/ecdf_example1.pdf")
# 
# 
# plotDistWithApproximationDifferentQ(OutputPath = OutputPath, protName = "000_Trx",n = 100, 
#                                     qlist = c(1,2,3,10),
#                                     plotToFile = "/home/willy/PredictingProteinInteractions/Results/QuickRepSampling/Examples/ecdf_exampleMultipleq1.pdf")
# 
# 
# errors = calcErrors(pts, times = 100, n = 100, q = 1)
# boxplot(errors, ylim = c(0,1))
# 
# errors = calcErrors(pts, times = 100, n = 100, q = 10)
# boxplot(errors, ylim = c(0,1))
# 
# errors = calcErrors(pts, times = 100, n = 100, q = 99)
# boxplot(errors, ylim = c(0,1))
# 
# meanErrorsQ = caclErrorCurve(pts, times = 100, n = 3, maxQ = 2)
# 
# plot(x = c(1:length(meanErrorsQ)), y = meanErrorsQ, type ="l", xlab = "q", ylab = "mean(Error)", main = "abs(DE_exact - DE_approx)")
# 
# 
# 
# calcError(pts,q = 10)

#-------------------------------------------------------------------



getAll_protein_F_approximations <- function(OutputPath,  n = 100, m = 50, q = 2, pos = TRUE){
  print(paste("Loading from ", OutputPath, sep =""))
  protein_names = list.dirs(OutputPath, recursive = FALSE, full.names = FALSE)
  
  distributions_lists = list()
  
  for(i in 1:length(protein_names)){
    F_app = generateF_approximations(OutputPath = OutputPath, protName = protein_names[i], n = n, m = m, q = q, pos =pos)
    distributions_lists[[i]] =  list("name" = protein_names[i],"F" = F_app)
  }
  
  return(distributions_lists)
}

get_x_and_y_vals <- function(all_protein_F_approximations,m){
  x_vals = rep(0,length(all_protein_F_approximations)*m)
  y_vals = rep(0,length(all_protein_F_approximations)*m)
  
  for(i in 1:length(all_protein_F_approximations)){
    for(j in 1:m){
      ind = (i-1)*m
      x_vals[ind+j] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
      y_vals[ind+j] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][2]-all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
    }
  }
  
  return(list("x_vals" = x_vals, "y_vals" = y_vals))
}

plot_all_F_approximations <- function(all_protein_F_approximations, colors = NULL, onlyGeo = FALSE, xli = NULL, yli = NULL){
  
  n_proteins = all_protein_F_approximations
  m = length(all_protein_F_approximations[[1]]$F$F_app_list)
  
  if(is.null(xli) || is.null(yli)){
    obj = get_x_and_y_vals(all_protein_F_approximations,m)
    x_vals = obj$x_vals
    y_vals = obj$y_vals
    
    xli = c(min(x_vals),max(x_vals))
    yli = c(min(y_vals),max(y_vals))
  }  
  
  if(is.null(colors)) colors = rep("black", length(n_proteins))
  
  plot2d_dist_approx(all_protein_F_approximations[[1]]$F$F_app_list, xli = xli, yli = yli, col = colors[1], name = all_protein_F_approximations[[1]]$name,onlyGeo = onlyGeo)
  for(i in 2:length(all_protein_F_approximations)){
    
    # print(i)
    plot2d_dist_approx(all_protein_F_approximations[[i]]$F$F_app_list, add = TRUE, col = colors[i], name = all_protein_F_approximations[[i]]$name, onlyGeo = onlyGeo)
  }
  
}

getManhattanProjection2d <- function(all_protein_F_approximations){
  prot_number = length(all_protein_F_approximations)
  
  m = length(all_protein_F_approximations[[1]]$F$F_app_list)
  
  df = data.frame(matrix(0,nrow = m*prot_number, ncol = 3))
  colnames(df) = c("name", "x", "y")
  
  for(i in 1:prot_number) {
    print(all_protein_F_approximations[[i]]$name)
    name = all_protein_F_approximations[[i]]$name
    
    ind = (i-1)*m
    for(j in 1:m){
      df[ind+j,1] = name
      df[ind+j,2] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
      df[ind+j,3] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][2]-  all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
    }
    
  }
  
  return(df)
}

# length(all_protein_F_approximations[[1]]$F$F_app_list[[1]])

getManhattanProjection <- function(all_protein_F_approximations){
  prot_number = length(all_protein_F_approximations)
  
  m = length(all_protein_F_approximations[[1]]$F$F_app_list)
  q = length(all_protein_F_approximations[[1]]$F$F_app_list[[1]])
  n = length(knots(all_protein_F_approximations[[1]]$F$F_list[[1]]))
  
  df = data.frame(matrix(0,nrow = m*prot_number, ncol = q+1+1))
  
  na = c("",q)
  for(i in 1:q){
    na[i] = paste("q_",i,sep="")
  }
  
  colnames(df) = c("name","n", na)
  
  for(i in 1:prot_number) {
    print(paste(all_protein_F_approximations[[i]]$name, i/prot_number))
    name = all_protein_F_approximations[[i]]$name
    
    ind = (i-1)*m
    for(j in 1:m){
      df[ind+j,1] = name
      df[ind+j,2] = n
      
      df[ind+j,1+1+1] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][1]
      for(k in 2:q){
        df[ind+j,k+1+1] = all_protein_F_approximations[[i]]$F$F_app_list[[j]][k] - df[ind+j,k+1]
      }
      
    }
    
  }
  
  return(df)
}


quickRepSampling <- function(OutputPath, distance_path, fName, pos = TRUE, n = 100, m = 22, q = 21, plot = FALSE, functionals = NULL, functionalMainTarget = "000_Trx", distance_method = "emd"){
  # actually q should be possible to be larger than n
  #
  #-------------------------------------------------------
  if(distance_method == "emd"){
    fName_final = paste(distance_path,"/",fName, "_emd.csv",sep ="")
  } else {
    fName_final = paste(distance_path,"/",fName, "_geo.csv",sep ="")
  }

  fName_final_projection = paste(distance_path,"/",fName, "_projection.csv",sep ="")
  fName_final_geometricCenters = paste(distance_path,"/",fName, "_geometricCenters.csv",sep ="")
  
  if(!file.exists(fName_final) || !file.exists(fName_final_projection)){
    print(paste("Creating ",fName_final))
    if(!dir.exists(distance_path)){
      dir.create(distance_path)
    }
    
    all_protein_F_approximations = getAll_protein_F_approximations(OutputPath = OutputPath, m = m, q = q, n = n, pos = TRUE)
    
    df = getManhattanProjection(all_protein_F_approximations)
    write.csv(df,file = fName_final_projection, row.names = FALSE)

    protein_distances = c()
    if(distance_method == "emd"){
      distances = getDistanceMatrixFromProjection(df)
      Emd_distances = emd_parallel(DE_parallel = distances,m = m)
      
      # make symmetric
      for(i in 1:(nrow(Emd_distances)-1)){
        for(j in (i+1):ncol(Emd_distances)){
          v = Emd_distances[i,j] + Emd_distances[j,i]
          Emd_distances[i,j] = v
          Emd_distances[j,i] = v
        }
      }
      
      protein_distances = Emd_distances
      write.csv(Emd_distances,file = fName_final)
    } else {
      geos = getGeometricCenters(df, n = n)
      write.csv(geos, file = fName_final_geometricCenters)
      
      geoDistances = getDistanceMatrixFromProjection(geos)
      protein_distances = geoDistances
      write.csv(geoDistances,file = fName_final)
    }

    # can be loaded with loadDistanceMatrix1()
    
    if(plot == TRUE && q == 2){
      print("ploting ...")
      colors = rep("black", length(all_protein_F_approximations))
      if(!is.null(functionals)){
        colors = rep("blue", length(all_protein_F_approximations))
        for(i in 1:length(all_protein_F_approximations)){
          if(all_protein_F_approximations[[i]]$name %in% functionals) colors[i] = "red"
        }
      }

      fName_final_projection_plot2d = paste(distance_path,"/",fName, "_projection2d_plot.pdf",sep ="")
      pdf(fName_final_projection_plot2d)
      plot_all_F_approximations(all_protein_F_approximations,
                                colors = colors,
                                onlyGeo = TRUE,
                                xli = NULL,
                                yli = NULL)


      dev.off()
      

    }

    fName_final_projection_pseudoRoc = paste(distance_path,"/",fName, "_pseudoRoc.pdf",sep ="")
    pdf(fName_final_projection_pseudoRoc)
    plotPseudoRoc(geoDistances = protein_distances,ind = which(rownames(protein_distances) == functionalMainTarget),functionals = functionals)
    dev.off()
  } 
  
  
  
  print(fName_final)
  protein_distances = read.csv(fName_final, header = TRUE, row.names = 1)
  
  # print(protein_distances)
  
  # print(nrow(protein_distances))
  # print(ncol(protein_distances))
  
  return(protein_distances)
}


# quickRepSampling(OutputPath = "/home/willy/Schreibtisch/106Test/Output/", distance_path = "/home/willy/Schreibtisch/106Test/Output/", fName = "distTest")

# library(rbenchmark)
# bench = benchmark("quickRepSampling_n_100_m_22_q_2" = quickRepSampling(OutputPath = "/home/willy/Schreibtisch/106Test/Output/",
#                                                                        distance_path = "/home/willy/Schreibtisch/106Test/Output/",
#                                                                        fName = "distTest",
#                                                                        n = 100,
#                                                                        m = 22,
#                                                                        q = 2),
#                   "quickRepSampling_n_100_m_22_q_20" = quickRepSampling(OutputPath = "/home/willy/Schreibtisch/106Test/Output/",
#                                                                         distance_path = "/home/willy/Schreibtisch/106Test/Output/",
#                                                                         fName = "distTest",
#                                                                         n = 100,
#                                                                         m = 22,
#                                                                         q = 20),
#                   "quickRepSampling_n_100_m_100_q_20" = quickRepSampling(OutputPath = "/home/willy/Schreibtisch/106Test/Output/",
#                                                                          distance_path = "/home/willy/Schreibtisch/106Test/Output/",
#                                                                          fName = "distTest",
#                                                                          n = 100,
#                                                                          m = 100,
#                                                                          q = 20),
#                   replications = 2,
#                   columns = c("test", "replications", "elapsed",
#                               "relative", "user.self", "sys.self"))

#---------------------------------------------------------------------------------------
#                                 test replications elapsed relative user.self sys.self
# 1  quickRepSampling_n_100_m_22_q_2              2  11.573    1.000    11.411    0.168
# 2 quickRepSampling_n_100_m_22_q_20              2  15.741    1.360    15.550    0.217
# 3 quickRepSampling_n_100_m_100_q_20             2  69.508    6.006    65.665    3.787
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# CompareProteinsCall
#---------------------------------------------------------------------------------------
# compareProteinsCall = "/home/willy/PredictingProteinInteractions/MetricGeometry/ComparingProteins/./CompareIsosurfaces.R"
# compareProteinsParameterFile = "/home/willy/PredictingProteinInteractions/MetricGeometry/Benchmark/Parameter"
# system2(command = compareProteinsCall, args = compareProteinsParameterFile)
# 
# 
# 
# library(rbenchmark)
# benchmark("CompareProteins" = system2(command = compareProteinsCall, args = compareProteinsParameterFile),
#           replications = 1,
#           columns = c("test", "replications", "elapsed",
#                       "relative", "user.self", "sys.self"))

#---------------------------------------------------------------------------------------

# n = 100
# m = 100
# q = 2
# times = 2
# OutputPath = "/home/willy/Schreibtisch/106Test/Output/"
# all_protein_F_approximations = getAll_protein_F_approximations(OutputPath = OutputPath, m = m, q = q, n = n, pos = TRUE)
# 
# functionals = c(getFunctionalProteins(),"000_Trx")
# colors = rep("blue", length(all_protein_F_approximations))
# for(i in 1:length(all_protein_F_approximations)){
#   if(all_protein_F_approximations[[i]]$name %in% functionals) colors[i] = "red"
# }
# 
# # plot_all_F_approximations(all_protein_F_approximations, colors = colors, onlyGeo = TRUE, xli = c(16, 22), yli = c(1,6))
# 
# 
# plot_all_F_approximations(all_protein_F_approximations, colors = colors, onlyGeo = TRUE)
# 
# 
# 
# # df = getManhattanProjection2d(all_protein_F_approximations)
# df = getManhattanProjection(all_protein_F_approximations)
# distances = getDistanceMatrixFromProjection(df)
# 
# geos = getGeometricCenters(df)
# 
# geos[7,]
# 
# xl = geos[7,2] + c(-0.1,1) 
# yl = geos[7,3] + c(-1,1) 
# unique(geos[intersect(intersect(intersect(which(geos[,2] < xl[2]), which(geos[,2] > xl[1])),which(geos[,3] < yl[2])), which(geos[,3] > yl[1])),1])
# 
# geoDistances = getDistanceMatrixFromProjection(geos)
# 
# geoDistances[7,]
# 
# Emd_distances = emd_parallel(DE_parallel = distances,m = m)
# 
# 
# emdRef = readDistanceMatrix3(file = "/home/willy/Schreibtisch/106Test/Output/")
# cor(c(Emd_distances[-106,-106]),c(emdRef))
# 
# 
# closest = rep(0,106)
# for(i in 1:length(closest)){
#   closest[i] = length(which(rownames(geoDistances)[which.minn(geoDistances[7,],n = i)] %in% functionals) == TRUE)
# }
# 
# emdRef = readDistanceMatrix2()
# closest2 = rep(0,105)
# for(i in 1:length(closest2)){
#   closest2[i] = length(which(rownames(emdRef)[which.minn(emdRef[7,],n = i)] %in% functionals) == TRUE)
# }
# plot(closest2, type = "l", main = "closest neihgbours that are functional", xlab ="NN", ylab = "functionals found", col = "blue")
# points(closest, type = "l", main = "closest neihgbours that are functional", xlab ="NN", ylab = "functionals found", col = "red")
# legend(x = 35, y = 10,legend = c("RepSamp", "QuickRepSamp"),lty = c(1,1), col = c("blue", "red"))
# 




plotPseudoRoc <- function(geoDistances,ind,functionals){
  closest = rep(0,nrow(geoDistances))
  for(i in 1:length(closest)){
    closest[i] = length(which(rownames(geoDistances)[which.minn(geoDistances[ind,],n = i)] %in% functionals) == TRUE)
  }
  plot(closest, type = "l", main = "closest neihgbours that are functional", xlab ="NN", ylab = "functionals found", col = "blue")
}






