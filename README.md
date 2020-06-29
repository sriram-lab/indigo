# indigo
This is a readme for indigo repository!

Currently working on better README with code blocks showing how to use code!

Main file to run INDIGO is called indigo.m (in indigoCode). This script executes the following functions:

1. `[summary, druglist] = dataExplore(testFile)`

Takes in the file of interest and returns a summary of descriptive statistics of the data, including synergy and antagonism counts. Also gets a list of unique drugs. Graphs the data to show you distribution of scores.

2. `[sharedInteractions] = dataCompare(compareFiles)`

Takes in a list of files and determines number of shared interactions. Determines correlation between datasets. Only works for datasets containing pairwise drug interactions.

3. `indigoSummary = indigoRun(testFile, trainingData, valMethod, K, standardize, input_type);`

Main part of the process. Takes in the following:
- file with data you want to test
- Validation method: cross validation, holdout, independent
- the number of subsets/iterations you want to do during validation
- orthology file if test data is not e. coli
- trainingData you want to use. Currently can be train on original INDIGO and nature data
- standardize -converts scores to z-scores

4. `[stats, averages, overview] = analyze(indigoSummary)`

Takes output from INDIGO and analyzes results. Determines rank correlation, accuracy, etc. Plots ROC curves, confusion matrices and scatter plots of results (rank sorted).

All data files should be in indigoData

All results from analyze are saved and stored in indigoResults
