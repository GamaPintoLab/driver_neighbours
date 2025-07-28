# driver_neighbours
Code for the analysis of neighbours of cancer drivers.

Raw data files used in this analysis are available on [Zenodo](https://zenodo.org/records/14674267). Extract data.zip to the root of the code directory.

To reproduce the figures in the manuscript use `bct_wct_figures.R`. All the output files needed to reproduce the figures are already available in the Zenodo repository.

To reproduce the entire workflow, first run the python notebooks 1 to 4 (as indicated below) and then the R code.

## Python code
To reproduce python code, use the `environment.yml` file to create a virtual environment using [conda](https://conda.org/), [pixi](https://pixi.sh/latest/) or other similar tool. 

The work is divided into 5 jupyter notebooks:
1. `get_PPIs.ipynb`: this is the first part of our work. It retrieves and processes protein-protein interactions (PPI) data.
2. `preprocessing.ipynb`: This notebook uses the processed PPI data and raw TCGA data to create files for BCT and WCT analysis.
3. `pathway_analysis.ipynb`: Computation of shortest path lenghts between drivers and neighbour genes. 
4. `bct_analysis.ipynb`: Runs BCT analyses. It also includes the tumour mutational burden calculation. 
5. `opentargets_analysis.ipynb`: Contains the last part of our work. Run it after the R code. Using the Open Targets Platform, it searches for potential therapeutic targets among candidate neighbours.  

## R code
R packages needed:
1. doParallel
2. doSnow
3. ggplot2
4. ggpubr
5. arrow
6. tictioc
7. Coselens<br>


To reproduce the figures in the manuscript:<br>
Run bct_wct_figures.R<br>

To run the tumour-normal paired sample analysis and compare it with the signalling pathway analysis:<br>
Run paired_analysis.R<br>
Note: code is parallelized. Change nclust accordingly to the number of available cores.<br>

To run wct analysis (including analysis of shuffled datasets:<br>
Run wct_computation.R<br>
Note: code is parallelized. Change nclust accordingly to the number of available cores. Tests with permuted data take long time to run.<br>

To analyse bct and wct results by neighbour and by driver:<br>
Run bct_wct_analysis.R<br>

To run the validation of wct interactions using Coselens:<br>
Run coselens_validation.R<br>


bct_wct_functions.R contains accessory functions used in the remaining R scripts.<br>

