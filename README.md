# driver_neighbours
Code for the analysis of neighbours of cancer drivers<br>

Data files necessary to run the analysis are available at: <br>
https://zenodo.org/records/14674267 <br>


Step 1 - prepare data frames for analysis:<br>
Run driver_neighbours_prep.R<br>

Step 2 - find interactions where neighbour expression is influenced by driver mutation (through analysis of paired tumour-normal samples):<br>
Run driver_neighbours_paired.R<br>
Necessary data files: <br>
1. driver_neib_pairs.csv; 
2. mutationtab1.RData;
3. exptab1_0.RData; 
4. cancertype1.RData<br>
Note: code is parallelized. Change nclust accordingly to the number of available cores.<br>

Step 3 - compute neighbour driver statistical associations:<br>
Run driver_neighbour_associations_tests.R<br>
Necessary data files: <br>
1. driver_neib_pairs.csv; 
2. mutationtab3.RData; 
3. exptab3.RData; 
4. cancertype3.RData; 
5. dfrespair2.RData <br>
Note: code is parallelized. Change nclust accordingly to the number of available cores. Tests with permuted data take long time to run.<br>

Step 4 - analyse associations by neighbour:<br>
Run driver_neighbour_analysis_by_neighbour.R<br>
Necessary data files: <br>
1. dfrespair2.RData; 
2. dfreslmd3.RData; 
3. ctdtab3.RData; 
4. ctetab3.RData; 
5. cancertype3.RData; 
6. mutationtab3.RData; 
7. exptab3.RData<br>

Step 5 - compute Cox regression analysis for all neighbours and cancer types:<br>
Run driver_neighbour_survival.R<br>
Necessary data files: <br>
1. clintab2.RData; 
2. exptab2.RData; 
3. neibrho.RData<br>

Step 6 - combine and correlate driver and survival associations:<br>
Run driver_neighbour_analysis_final.R<br>
Necessary data files: <br>
1. neibsurv.RData; 
2. neibrho.RData; 
3. neibcoef.RData<br>

Stpe 7 - search neighbours in Open Targetd:<br>
Run opentargets.ipynb<br>
Necessary data files: <br>
1. unzip OpenTargets.zip (the resulting OpenTargets folder should be in the same location as opentargets.ipynb); 
<br>

driver_neighbour_functions.R contains accessory functions used in the remaining R scripts.<br>

