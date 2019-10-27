#-------------------------------------------------------------------
#
# Benchmarking of CompareProteins and QuickRepeatedSampling
#
#
#
#-------------------------------------------------------------------

s1 = "/home/willy/PredictingProteinInteractions/MetricGeometry/QuickRepeatedSubSampling/UltraQuickRepeatedSubSampling.R"
source(s1)

distance_path = "/home/willy/PredictingProteinInteractions/MetricGeometry/Benchmark/QuickRep/"
OutputPath = "/home/willy/PredictingProteinInteractions/data/106Test/Output/"

quickRepSampling(OutputPath = "/home/willy/PredictingProteinInteractions/data/106Test/Output/",
                 distance_path = "/home/willy/PredictingProteinInteractions/data/106Test/UltraQuickRepSub/",
                 n = 100,
                 m = 22,
                 q = 1,
                 pos = "pos",
                 fName = "quickRepSampling_n_100_m_22_q_1_geo",
                 plot = TRUE,
                 distance_method = "geo",
                 recalculate = TRUE,
                 plotRoc = FALSE)

quickRepSampling(OutputPath = "/home/willy/PredictingProteinInteractions/data/106Test/Output/",
                 distance_path = "/home/willy/PredictingProteinInteractions/data/106Test/UltraQuickRepSub/",
                 n = 100,
                 m = 500,
                 q = 1,
                 pos = "pos",
                 fName = "quickRepSampling_n_100_m_500_q_1_geo",
                 plot = TRUE,
                 distance_method = "geo",
                 recalculate = TRUE,
                 plotRoc = FALSE)


library(rbenchmark)

bench = benchmark("quickRepSampling_n_100_m_22_q_1_emd" = quickRepSampling(OutputPath = OutputPath, 
                                                                       distance_path = distance_path,
                                                                       n = 100,
                                                                       m = 22,
                                                                       q = 1,
                                                                       pos = "pos",
                                                                       fName = "quickRepSampling_n_100_m_22_q_1_emd",
                                                                       plot = TRUE,
                                                                       distance_method = "emd",
                                                                       recalculate = TRUE,
                                                                       plotRoc = FALSE),
                  "quickRepSampling_n_100_m_22_q_1_geo" = quickRepSampling(OutputPath = OutputPath, 
                                                                           distance_path = distance_path,
                                                                           n = 100,
                                                                           m = 22,
                                                                           q = 1,
                                                                           pos = "pos",
                                                                           fName = "quickRepSampling_n_100_m_22_q_1_geo",
                                                                           plot = TRUE,
                                                                           distance_method = "geo",
                                                                           recalculate = TRUE,
                                                                           plotRoc = FALSE),
                  
                  "quickRepSampling_n_100_m_22_q_20_emd" = quickRepSampling(OutputPath = OutputPath, 
                                                                           distance_path = distance_path,
                                                                           n = 100,
                                                                           m = 22,
                                                                           q = 20,
                                                                           pos = "pos",
                                                                           fName = "quickRepSampling_n_100_m_22_q_20_emd",
                                                                           plot = TRUE,
                                                                           distance_method = "emd",
                                                                           recalculate = TRUE,
                                                                           plotRoc = FALSE),
                  "quickRepSampling_n_100_m_22_q_20_geo" = quickRepSampling(OutputPath = OutputPath, 
                                                                           distance_path = distance_path,
                                                                           n = 100,
                                                                           m = 22,
                                                                           q = 20,
                                                                           pos = "pos",
                                                                           fName = "quickRepSampling_n_100_m_22_q_20_geo",
                                                                           plot = TRUE,
                                                                           distance_method = "geo",
                                                                           recalculate = TRUE,
                                                                           plotRoc = FALSE),
                  
                  "quickRepSampling_n_100_m_500_q_1_geo" = quickRepSampling(OutputPath = OutputPath, 
                                                                             distance_path = distance_path,
                                                                             n = 100,
                                                                             m = 500,
                                                                             q = 1,
                                                                             pos = "pos",
                                                                             fName = "quickRepSampling_n_100_m_500_q_1_geo",
                                                                             plot = TRUE,
                                                                             distance_method = "geo",
                                                                             recalculate = TRUE,
                                                                             plotRoc = FALSE),
                  
                  "quickRepSampling_n_100_m_500_q_20_geo" = quickRepSampling(OutputPath = OutputPath, 
                                                                            distance_path = distance_path,
                                                                            n = 100,
                                                                            m = 500,
                                                                            q = 20,
                                                                            pos = "pos",
                                                                            fName = "quickRepSampling_n_100_m_500_q_20_geo",
                                                                            plot = TRUE,
                                                                            distance_method = "geo",
                                                                            recalculate = TRUE,
                                                                            plotRoc = FALSE),

                  replications = 10,
                  columns = c("test", "replications", "elapsed",
                              "relative", "user.self", "sys.self"))

write.table(bench,file = "/home/willy/PredictingProteinInteractions/MetricGeometry/Benchmark/Benchmark.txt",row.names = FALSE)




#---------------------------------------------------------------------------------------
# test                                replications  elapsed relative user.self  sys.self
# quickRepSampling_n_100_m_22_q_1_emd             1  32.345    4.550    31.611    0.937
# quickRepSampling_n_100_m_22_q_1_geo             1   7.109    1.000     6.951    0.228
# quickRepSampling_n_100_m_22_q_20_emd            1  34.646    4.874    33.967    0.887
# quickRepSampling_n_100_m_22_q_20_geo            1  10.034    1.411     9.906    0.202
# quickRepSampling_n_100_m_500_q_20_geo           1 299.356   42.109   295.184    4.669
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
#                                 test replications elapsed relative user.self sys.self
# 1  quickRepSampling_n_100_m_22_q_2              2  11.573    1.000    11.411    0.168
# 2 quickRepSampling_n_100_m_22_q_20              2  15.741    1.360    15.550    0.217
# 3 quickRepSampling_n_100_m_100_q_20             2  69.508    6.006    65.665    3.787
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# test                                 replications elapsed relative user.self sys.self
# quickRepSampling_n_100_m_22_q_1_emd             1  32.939    4.541    32.290    0.862
# quickRepSampling_n_100_m_22_q_1_geo             1   7.254    1.000     7.108    0.215
# quickRepSampling_n_100_m_22_q_20_emd            1  36.026    4.966    35.309    0.965
# quickRepSampling_n_100_m_22_q_20_geo            1  10.166    1.401     9.996    0.259
# quickRepSampling_n_100_m_500_q_1_geo            1 135.448   18.672   133.379    2.418
# quickRepSampling_n_100_m_500_q_20_geo           1 571.133   78.734   566.451    5.000
# quickRepSampling_n_100_m_500_q_20_geo           1 443.418        1   438.875    5.236
#---------------------------------------------------------------------------------------

bench2 = benchmark("quickRepSampling_n_100_m_500_q_20_geo" = quickRepSampling(OutputPath = OutputPath, 
                                                                             distance_path = distance_path,
                                                                             n = 100,
                                                                             m = 500,
                                                                             q = 20,
                                                                             pos = "pos",
                                                                             fName = "quickRepSampling_n_100_m_500_q_20_geo",
                                                                             plot = TRUE,
                                                                             distance_method = "geo",
                                                                             recalculate = TRUE,
                                                                             plotRoc = FALSE),
                  
                  replications = 1,
                  columns = c("test", "replications", "elapsed",
                              "relative", "user.self", "sys.self"))


bench2

#---------------------------------------------------------------------------------------
# CompareProteinsCall
#---------------------------------------------------------------------------------------
compareProteinsCall = "/home/willy/PredictingProteinInteractions/MetricGeometry/ComparingProteins/./CompareIsosurfaces.R"
compareProteinsParameterFile = "/home/willy/PredictingProteinInteractions/MetricGeometry/Benchmark/Parameter"
system2(command = compareProteinsCall, args = compareProteinsParameterFile)



# before benchmarking "CompareProteins" you have to delete the coresponding files from any previous calculations
library(rbenchmark)
benchmark("CompareProteins" = system2(command = compareProteinsCall, args = compareProteinsParameterFile),
          replications = 1,
          columns = c("test", "replications", "elapsed",
                      "relative", "user.self", "sys.self"))

#---------------------------------------------------------------------------------------

bench  = read.table(file = "/home/willy/PredictingProteinInteractions/MetricGeometry/Benchmark/Benchmark.txt", header = TRUE)

benchFinal = bench[,c("test","elapsed", "relative")]
benchFinal$elapsed = bench$elapsed/bench$replications

df = data.frame(matrix(c("CompareProteins_n_100_m_500", 14400, 0), nrow = 1))
colnames(df) = c("test","elapsed", "relative")
benchFinal = rbind(benchFinal,df)
benchFinal$elapsed = as.numeric(benchFinal$elapsed)
benchFinal$elapsed = round(benchFinal$elapsed)

benchFinal
v = benchFinal$elapsed/benchFinal$elapsed[6]
benchFinal$relative = format(round(v, 2), nsmall = 2)

benchFinal

write.table(benchFinal,file = "/home/willy/PredictingProteinInteractions/MetricGeometry/Benchmark/BenchmarkAll.txt",row.names = FALSE,quote = FALSE)

