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
                  
                  "quickRepSampling_n_100_m_500_q_20_emd" = quickRepSampling(OutputPath = OutputPath, 
                                                                            distance_path = distance_path,
                                                                            n = 100,
                                                                            m = 500,
                                                                            q = 20,
                                                                            pos = "pos",
                                                                            fName = "quickRepSampling_n_100_m_500_q_20_emd",
                                                                            plot = TRUE,
                                                                            distance_method = "emd",
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

                  replications = 1,
                  columns = c("test", "replications", "elapsed",
                              "relative", "user.self", "sys.self"))

write.table(bench,file = "/home/willy/PredictingProteinInteractions/MetricGeometry/Benchmark/Benchmark.txt",row.names = FALSE)
#---------------------------------------------------------------------------------------
#                                 test replications elapsed relative user.self sys.self
# 1  quickRepSampling_n_100_m_22_q_2              2  11.573    1.000    11.411    0.168
# 2 quickRepSampling_n_100_m_22_q_20              2  15.741    1.360    15.550    0.217
# 3 quickRepSampling_n_100_m_100_q_20             2  69.508    6.006    65.665    3.787
#---------------------------------------------------------------------------------------

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