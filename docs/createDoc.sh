#!/bin/sh
#-------------------------------------------------------------------
# 14.4.2019
# This script creates the documentation need for the masterthesis.
# - Additionally the masterthesis is copied from the latex-folder.
# - The documentation with mkdocs is built and opened in firefox.
# 
#-------------------------------------------------------------------
PATH2Predicting=$(pwd)
echo $PATH2Predicting

# copy the masterthesis
cp /home/willy/RedoxChallenges/MasterThesis/WillyMT/main.pdf $PATH2Predicting/MasterThesis.pdf
# get the tree
#tree /home/willy/PredictingProteinInteractions/ -d -I site /home/willy/PredictingProteinInteractions/data/Input data/106Redoxins/Output > /home/willy/PredictingProteinInteractions/docs/tree.txt

tree $PATH2Predicting/ -d -L 3 -I site > $PATH2Predicting/docs/tree.txt



# static documentation. TODO: as pdf
cd /home/willy/PredictingProteinInteractions/
mkdocs build


# interactive documentation
cd $PATH2Predicting

# copy to front to see at github
cp docs/index.md README.md

mkdocs serve &
firefox http://127.0.0.1:8000/ &

