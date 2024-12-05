# driver_neighbours
Code and data for the analysis of neighbours of cancer drivers

Data files necessary to run the analysis is available at: 
https://zenodo.org/records/14284408 


Step 1 - prepare data frames for analysis:
Run driver_neighbours_prep.R

Step 2 - find interactions where neighbour expression is influenced by driver mutation (through analysis of paired tumour-normal samples)
Run driver_neighbours_paired.R
Necessary data files: driver_neib_pairs.csv; mutationtab1.RData;exptab1_0.RData; cancertype1.RData
Note: code is parallelized. Change nclust accordingly to the number of available cores.

Step 3 - compute neighbour driver statistical associations
Run driver_neighbour_associations_tests.R
Necessary data files: driver_neib_pairs.csv; mutationtab3.RData; exptab3.RData; cancertype3.RData; dfrespair2.RData; 
Note: code is parallelized. Change nclust accordingly to the number of available cores. Tests with permuted data take long time to run.

Step 4 - analyse associations by neighbour
Run driver_neighbour_analysis_by_neighbour.R
Necessary data files: dfrespair2.RData; dfreslmd3.RData; ctdtab3.RData; ctetab3.RData; cancertype3.RData; mutationtab3.RData; exptab3.RData

Step 5 - compute Cox regression analysis for all neighbours and cancer types
Run driver_neighbour_survival.R
Necessary data files: clintab2.RData; exptab2.RData; neibrho.RData

Step 6 - combine and correlate driver and survival associations
Run driver_neighbour_analysis_final.R
Necessary data files: neibsurv.RData; neibrho.RData; neibcoef.RData

driver_neighbour_functions.R contains accessory functions used in the remaining R scripts.

