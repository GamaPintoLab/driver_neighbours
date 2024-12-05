# driver_neighbours
Code for the analysis of neighbours of cancer drivers

Data files necessary to run the analysis are available at: 
https://zenodo.org/records/14284408 


Step 1 - prepare data frames for analysis:
Run driver_neighbours_prep.R

Step 2 - find interactions where neighbour expression is influenced by driver mutation (through analysis of paired tumour-normal samples):
Run driver_neighbours_paired.R
Necessary data files: 
1. driver_neib_pairs.csv; 
2. mutationtab1.RData;
3. exptab1_0.RData; 
4. cancertype1.RData
Note: code is parallelized. Change nclust accordingly to the number of available cores.

Step 3 - compute neighbour driver statistical associations:
Run driver_neighbour_associations_tests.R
Necessary data files: 
1. driver_neib_pairs.csv; 
2. mutationtab3.RData; 
3. exptab3.RData; 
4. cancertype3.RData; 
5. dfrespair2.RData; 
Note: code is parallelized. Change nclust accordingly to the number of available cores. Tests with permuted data take long time to run.

Step 4 - analyse associations by neighbour:
Run driver_neighbour_analysis_by_neighbour.R
Necessary data files: 
1. dfrespair2.RData; 
2. dfreslmd3.RData; 
3. ctdtab3.RData; 
4. ctetab3.RData; 
5. cancertype3.RData; 
6. mutationtab3.RData; 
7. exptab3.RData

Step 5 - compute Cox regression analysis for all neighbours and cancer types:
Run driver_neighbour_survival.R
Necessary data files: 
1. clintab2.RData; 
2. exptab2.RData; 
3. neibrho.RData

Step 6 - combine and correlate driver and survival associations:
Run driver_neighbour_analysis_final.R
Necessary data files: 
1. neibsurv.RData; 
2. neibrho.RData; 
3. neibcoef.RData

driver_neighbour_functions.R contains accessory functions used in the remaining R scripts.

