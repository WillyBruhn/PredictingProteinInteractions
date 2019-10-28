# Predicting Protein Interactions
Predict protein-interactions from raw pdb(protein-databank)-files. 
This diagram illustrates the different functionallities of this software-tool and how the different parts relate to each other


![Diagram1](docs/Diagram.png)



## Installation
 A script that trys to install all necessary dependencies can be found in *setUp/setUp.sh*. Install all dependencies by typing:
```
./setUp.sh
```
In order to view this documentation on your local machine type:
 
```
./docs/createDoc.sh
```
This requires *mkdocs* to be installed.


## 1. Predicting Protein Interactions


## 2. Preprocessing

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

## 3. Prediction

## 4. Clustering
Create clusterings like this one from the 106-Redoxin-data-set:

![Diagram1](docs/ClusteringExample.png)



### UltraQuickRepeatedSubSampling
Calculate the repeatedSubSampling fast. The points are only sampled m times once for each protein. Then the distributions are calculated. Then from these distributions the quantile-
approximation is calculated. The DE is then calculated with the manhattan distance between
the quantiles.

#### Usage

```
  ProteinsPath        ... path to all proteins as produced by MutComp

  distance_name       ... folder in which all distance-matrices will be stored

  n                   ... number of points to select (see parameters of RepeatedSampling)

  m                   ... sqrt(number of repeatitions)

  q                   ... number of subdivisions of the integral. Basically a higher q
                          leads to a more accurate approximation. Currently it has to hold
                          q < n.

  potential           ... pos/neg

  distance_method     ... geo/emd


  plot                ... in case of (q == 2) the approximations are ploted into the 2d-plane.
                          Else no plot is created with a warn-message.


  labels              ... a file containing for each protein a label functional/not_functional

  cores               ... number of cores to run on. Trys to detect the number of cores
                          automatically. If this fails the number of cores is set to 6
                          by default.
```

### Clustering.R
Train an auto-encoder that reduces the number of features of the proteins. Then the
bottle-neck-layer is extracted and used to extract a condensed representation of
the proteins. With these condensed representations then an agglomerative, bottom-up
clustering is performed. With this clustering then a dendrogram is obtained.


#### Usage
```
  inputPath           ... path to a file with a .Rdata-extension that stores the feature-matrix for
                          all proteins. Can be obtained by running Proteins.R.

  dendrogramName      ... name of the dendrogram-file. Without the (.pdf)-extension and prefix Dendrogram.

  numPermutations     ... number of permutations that are created for each representation. For
                          each protein m rows from the feature-matrix are combined. The order of
                          this m rows is randomly permutated and numPermutations different
                          representations are created.

  m                   ... number of rows from the feature-matrix that are combined for each protein-model.

  epochs              ... number of epochs to train the autoencoder.

  batchSize           ... size of a batch. Relevant for the training.

  l1,l2,l3            ... specifies the encoder-dimensions, that is the size of each layer in the network.
                          The autoencoder in this implementation has 3 layers.

  d1,d2,d3            ... in [0,1) specifies the dropout-rates.
```




