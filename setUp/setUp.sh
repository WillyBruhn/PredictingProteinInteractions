
# in the parametersFile this is the value for vmdRC
sudo find . -name ".vmdrc"

#--------------------------------------------------------------------------------------
# compile Repeated Sub Sampling
cd /home/willy/PredictingProteinInteractions/MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/

mkdir cmakeBin

cd /home/willy/PredictingProteinInteractions/MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/

cmake ..

cd /home/willy/PredictingProteinInteractions/MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/

make
#--------------------------------------------------------------------------------------
