# Predicting Protein Interactions
Predict protein-protein-interactions from raw pdb(protein-databank)-files. Find the thesis in **Masterthesis.pdf**. 

This document gives instructions on how to recreate the data that was used in the thesis. As a first step it only serves as a documentation for me personally. The idea is however to transfer this data directly to the CD that has to be submitted with the thesis. 
(Willy Bruhn, 13.4.2019)


## Folder Hierarchy
```
├── Classification
├── data
│   ├── 106Redoxins
│   │   ├── Input
│   │   ├── Labels
│   │   └── Output
│   └── animals
├── docs
├── MetricGeometry
│   ├── ComparingProteins
│   │   ├── EMDandClustering
│   │   └── LowerBounds
│   ├── Downsampling2Step
│   └── RepeatedSubsampling
│       └── FirstLowerBoundRelationOfPosAndNeg
├── PCL
├── PreProcessingProteins
│   └── MutComp
│       ├── APBS
│       ├── deprecated
│       ├── ExampleHierarchy
│       ├── GUI
│       ├── ImageCreator
│       ├── ImageSimilarity
│       ├── MovieCreator
│       ├── mutationTest
│       ├── Output
│       ├── TclScripts
│       └── Tests
└── Results
    ├── Images
    ├── Scripts
    └── Tables
```

## QuickStart
Before you can run this example make sure that all dependencies are installed. A script that trys to do that can be found in *setUp/setUp.sh*. Install all dependencies by typing:
```
./setUp.sh
```

The folder *QuickStart* contains a fully setup example that demonstrates what this implementation can do.
In *ModelTrain/pdb/* pdb-files are stored to train a model. The file *ModelTrain/labels.txt* contains the labels of the proteins (in this case *functional* and *not_functional*, you can specify any label however). In *Test/pdb/* the pdb-files are stored that predictions should be made for. To build the model type:
```
../Classification/./predictingProteinInteractions.R
```

After the script is executed you should have an optimized model in *ModelTrain/bestModel/bestModel.RData*. This file is in a binary-format and can only be read with the R-function *readDs*. Additionally the computed distances are stored in *ModelTrain/RepeatedSubSampling/*.
Now make predictions:

```
../Classification/./predictingProteinInteractions.R --mode Predict --pdb_folder <path/To/QuickStart>/Test/pdb/ --MutCompOutPut <path/To/QuickStart>/Test/Proteins/
```

The predictions are now stored in *Predictions/predictions.txt*. The prediction-file has three columns:
the name of the protein, probability for *functional* probability for *not_functional*.
The output might look like this:

| "ind" | "functional"      | "not_functional"  |
|-------|-------------------|-------------------|
| "013" | 0.554907677356657 | 0.445092322643343 |
|       |                   |                   |
|       |                   |                   |




## 1. Classification
Contains one script that does a classification given one distance-matrix and labels and other parameters for the model.
```
Input:
distance_f_name("")           ... name of the file with the distances

distance_format(1)
                           = 1 ... a simple distance-matrix
                           = 3 ... a 3 column-data-frame. First column ProtA
                                   second column ProtB, third column distance.

 labels_train                   ... name of the file with the labels
                         
 evaluate(TRUE)        
                           = 1 ... with crossvalidation the model is tested
                           = 2 ... the model is tested on the testset
                           = 0 ... the model makes predictions on new data
                                   to which not all classlabels are given.

 evaluate_n                    ... number of times to repeat the classification.
                                   Higher number leads to more accurate estimation.

 kNN                           ... number of nearest neighbors to use for the majority-
                                   voting.

 test_split                    ... if(evaluate == 2) then this file specifies
                                   the names to be tested with a model built from 
                                   the remaining names.
```

## 2. data
### 106Redoxins
Contains the protein-test-set.

### additionalPDBS_1
Contains 19 Redoxins that are very simialr to Thrx.

### 106RedoxinsWithAdditionalPDBS_1
Contains the output of *106Redoxins* and *additionalPDBS_1* together in one folder. The files have to be in one folder for the calculation of the distances.

#### OutputMerged
Folder with all the pts-files.


*RepeatedSubsampling*
```
../../MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/./main /home/willy/PredictingProteinInteractions/data/106RedoxinsWithAddionalPDBS_1/OutputMerged/ /home/willy/PredictingProteinInteractions/data/106RedoxinsWithAddionalPDBS_1/RepSubOutput/ /home/willy/PredictingProteinInteractions/data/106RedoxinsWithAddionalPDBS_1/OutputMerged/names.txt /home/willy/PredictingProteinInteractions/data/106RedoxinsWithAddionalPDBS_1/OutputMerged/names.txt 1 100 500 1 0.0 0.0 test 0
```

*classification of RepeatedSubsampling*
```
../../Classification/./NNClassification.R --evaluate 0 --distance_format 3 --distance_f_name /home/willy/PredictingProteinInteractions/data/106RedoxinsWithAddionalPDBS_1/RepSubOutput/EMD_100_500_1.0000_1.0000_0.0000_0.0000_id_test.csv --kNN 10 --evaluate_n 1000 --test_split /home/willy/PredictingProteinInteractions/data/106Redoxins/Labels/testNames.txt --labels_train /home/willy/PredictingProteinInteractions/data/106Redoxins/Labels/labels.txt > predictions.txt
```

#### RepSubOutput
Folder with the distances obtained with the repeated-subsampling-method.

### Active Center

### animals
Contains the models of the animals also used in Memolis Paper.

## 3. Metric Geometry
*MetricGeometry/ComparingProteins/* contains the original implementation of Felix Berens cloned from Git-hub ([ComparingProteins](https://github.com/BerensF/ComparingProteins), 04.13.2019). 

### ComparingProteins
See for more details [ComparingProteins](https://github.com/BerensF/ComparingProteins).

Modifications:
* Other shebang for linux **#!/usr/bin/Rscript**

### Downsampling2Step
Contains the implementation of the 2-step-downsampling-procedure.

### RepeatedSubsampling
This folder contains the C++ implementation based on the original implementation of the *FLB* [ComparingProteins](https://github.com/BerensF/ComparingProteins) by Felix Berens. An example call would be:
```
pa="/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/"
./main $pa /home/willy/Schreibtisch/repeatedSampling/ /home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Example/proteinsTestSmall.txt /home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Example/proteinsTestSmall.txt 1 2 5 1 0.0 0.0 test 0
```
. More example-calls can be found in **/home/willy/PredictingProteinInteractions/MetricGeometry/RepeatedSubsampling/FirstLowerBoundRelationOfPosAndNeg/example**.

#### Usage
```
  string path = argv[1];							            // path to the folder with the proteins and pts-files
  string outPath = argv[2];						            // path to the folder where the output of this programm should be stored
  string proteinsToCompareFile_target = argv[3];	// proteins to compare to, usually only Trx or all proteins
  string proteinsToCompareFile = argv[4];			    // all other proteins

  measure = stod(argv[5]);						            // model parameter, measure == 1 equals a uniform-distribution
	number_of_selected_points = atoi(argv[6]);		  // model parameters n
	rounds = atoi(argv[7]);							            // model parameters m

  c1 = stod(argv[8]);		                          // model parameters
  c2 = stod(argv[9]);		                          // model parameters
  c3 = stod(argv[10]);		                        // model parameters

  num_bins = log2(rounds)+1;	                    // model parameters

  emd_list_id = argv[11];	                        // for running multiple times with same parameters
	bool allParameterCombinations = (bool)atoi(argv[12]); // for optimizing the modelparameters
	
	int NNtoActCent = atoi(argv[13]);               // = 0 then all points are taken
```

#### Compilation
Go to **/FirstLowerBoundRelationOfPosAndNeg/cmakeBin/** and type
```
cmake ..
```
and then type 
```
make
```
To run the program type:
```
./main
```


## 4. PCL
Contains the implementation that makes use of the ICP(iterative-closest-point-algorithm) of the *Point-Cloud-Library*. 

#### Setting Up
With the script *\IsoSurfSimilarity\setUp.sh* some dependencies are installed and downloaded.

#### Usage
For this implementation a seperate file called *settings.ini* is used

```
ICP = 0                             ... if 0: a genetic algorithm is performed
                                    ... if 1: the ICP-algorithm is performed
                      
ICPiterations = 100                 ... if (ICP = 1): number of iterations to run ICP

K = 0                               ... (if ICP = 0) In the step of finding the nearest neighbour in the potential of the other protein
                                                     this number specifies how many points should be considered. Given two points a similarity is
                                                     calculated depending on the distance to the active center and the local distribution of     
                                                     positive and negative points. Under the K nearest neighbours to the given point in euclidean 
                                                     space the point is chosen that has the highest similarity.
                                                    

WindowSize = 200                    ... size of the small Window for display

displayWeightTH = 10                ... show only points with a weight higher than "displayWeightTH".

evalFunction = simpleSqrdDistances  ... (if ICP = 0) specifies which evaluation-function for the distance should be taken
                                            if "simpleSqrdDistances" then for each point in X the closest point in Y is taken
                                            and the corresponding distance is squared. This is done for all points in X.
                                            All of these values are summed up. Analogously for Y. The sum of the two sums is
                                            the final value.

genNum = 100                        ... (if ICP = 0) specifies how many generations the genetic algorithm should run for

overWriteExistingInd = 1            ... if 1: overwrites a previous (if existing) solution, if the new one is better
                                    ... if 0: if another solution exists this new solution is discarded

path2Data = /home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/
                                    ... path to where the data is stored

popSize = 25                        ... if (ICP = 0): number of individuals in the population for the
                                        genetic algorithm

prot1FileName = 000_Trx             ... name of the protein to be read in. This would for example be
                                    /home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/000_Trx/000_Trx.pcd
    
                                    If you get an error-message showing something like 
                                    "[pcl::PCDReader::readHeader] Could not find file 
                                    '/home/willy/PredictingProteinInteractions/data/106Redoxins/Output//000_Trx/000_Trx.pcd'.
                                    Error loading object/scene file!" 
                                    please make sure to run ptsToPcd.sh beforehand.


prot2FileName = 103                 ... name of the protein to be read in. This would for example be
                                    /home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/103/103.pcd

range = 0                           ... (if ICP = 0) specifies an upper bound for the entries in the translation-vector of the individuals 
                                          in the genetic algorithm.
                                          
                                          if range = 0: then only the rotations are optimized. (The proteins are centered on the active center.)

visInteractive = 1                  ... if 1: an interactive model is displayed
                                    ... if 2: an interactive model is displayed 
                                              and the model is recalculated
                                    ... else: a model is calculated
```

#### Generating pcd(Point-Cloud-Data)-files
The PCL uses its own file-format.
With the scripts *ptsToPcd.sh* and *ptsToPcl.R* the pcd-files can be generated.


In *ptsToPcd.sh* change the variables *path* and *path2ptsToPcd* accordingly:
```
path="/home/willy/PredictingProteinInteractions/data/106Redoxins/Output/"
path2ptsToPcd="/home/willy/PredictingProteinInteractions/PCL/ptsToPcd/"
```

then run *ptsToPcd.sh* by typing:
```
./ptsToPcd.sh
```




### TODO
* Make IsoSurfSimilarity take one argument specifiying the path to the settings-file.
* Uncomment the lines in ptsToPcl.R that install the necessary packages
* move the sources from ptsToPcl.R to a place in PredictingProteinInteractions
* change *ptsToPcd.sh* such that it accepts one argument specifying the path to the data
instead of the user having to edit *ptsToPcd.sh* itself


## 5. PreprocessingProteins

Given a *pdb-file* this script generates the dx-files that are necessary in the further steps to analyse the similarities of the surfaces of the proteins.

### Mutcomp
([MutComp](https://github.com/WillyBruhn/MutComp)). 

The standard path to the parameters-file is (for the use on the WS):
```
/home/sysgen/Documents/LWB/TCL/MutComp/GUI/Parameters/parameters.txt
```
. To run the script type:
```
<pathToMutComp>./process.sh
```
where *<pathToMutComp>* is the path to the script.
If you want to use a different parameter-file just pass it as an argument to the script
```
<pathToMutComp>./process.sh <path/to/a/different/parametersfile>
```
On my machine this call becomes: 
```
/home/willy/PredictingProteinInteractions/PreProcessingProteins/MutComp/./process.sh /home/willy/PredictingProteinInteractions/PreProcessingProteins/MutComp/GUI/Parameters/parametersForThesis.txt
```

I strongly recommmend placing the parametersFile for each data-set in the same folder. 


### centerSelect
([centerSelect](https://github.com/WillyBruhn/centerSelect)). 

## 6. Results

### Scripts
Contains all the scripts that are necessary to reproduce the results that are presented in **Masterthesis.pdf**.

* Classification.R creates the evaluations presented in **Masterthesis**

### Images
Contains the images that were produced in the thesis with the implementations.

### Tables
Contains the tables that were produced in the thesis with the implementations.

# TODO


## Research
* active center
* use the measure parameter
* additional data tested with the model





