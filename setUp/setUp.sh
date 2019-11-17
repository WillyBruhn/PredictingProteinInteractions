
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




DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "Using this directory ..."
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



#--------------------------------------------------------------------------------------

antiPrism=$(command -v off2obj)

if [ -z "$antiPrism" ]; then
    echo "installing antiprism (off2obj)... (don't forget to uncomment me)"
    sudo add-apt-repository ppa:antiprism/ppa
    sudo apt-get update
    sudo apt-get install antiprism
fi

cd ..
if [ ! -d "Manifold" ]; then
    git clone --recursive -j8 git://github.com/hjwdzh/Manifold
    echo "installing Manifold (stitching 3d-objects together)..."
    

    cd Manifold
    mkdir build
    cd build
    cmake ..
    make
    cd ..
fi













