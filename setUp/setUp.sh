
# in the parametersFile this is the value for vmdRC
# sudo find . -name ".vmdrc"

#--------------------------------------------------------------------------------------
# compile Repeated Sub Sampling
#cd ../MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/
#
#mkdir cmakeBin
#cd cmakeBin/
#
#cmake ..
#
#make main
#--------------------------------------------------------------------------------------


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo $DIR

sed -i "s|PPISETUP =.*|PPISETUP = \"$DIR\"|g" SourceLoader.R


./installPackages.R 

# necessary to stich the objects together
#manifoldExe= $(locate /build/manifold)


#manifoldPath=${manifoldExe%/manifold}

#echo "\$manifoldExe"

#rm pathNames.txt
#echo "name path" >> pathNames.txt
#echo "manifold $manifoldExe" >> pathNames.txt
#echo "k k" >> pathNames.txt
















