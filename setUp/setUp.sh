
# in the parametersFile this is the value for vmdRC
# sudo find . -name ".vmdrc"

#--------------------------------------------------------------------------------------
# compile Repeated Sub Sampling
cd ../MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/

mkdir cmakeBin
cd cmakeBin/

cmake ..

make main
#--------------------------------------------------------------------------------------
