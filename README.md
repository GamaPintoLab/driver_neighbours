# driver_neighbours
Code for the analysis of neighbours of cancer drivers<br>

Raw data files used in this analysis are available at: <br>
https://zenodo.org/records/14674267 <br>

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

