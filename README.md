# INDIGO
This is a readme for indigo repository!

INDIGO is a random forest model that makes drug interaction predictions using drug chemogenomics/transcriptomics and drug interaction data. It runs in MATLAB.

I've modified INDIGO so that it can be trained on data from multiple studies representing multiple bacterial strains/species. 

**Google Drive folder containing Data and Results:** https://drive.google.com/drive/folders/13nP4fNvrQv_s-81Li1-Nx59m8m2U8zRT?usp=sharing

Repository contains following folders:

1. `core`
2. `python-scripts`
3. `tb-code`
4. `testing`
5. `validation`

Code required to run INDIGO model is in `validation`. Explanation of how to use code is below.

The following are different functions that can be used.

1. `[summary, druglist] = dataExplore(filename, dataLookup)`

Takes in the file of interest and returns a summary of descriptive statistics of the data, including synergy and antagonism counts. Also gets a list of unique drugs. Graphs the data to show you distribution of scores.

2. `indigoSummary = indigoRun(testData, trainingData, dataLookup, valMethod, K, standardize, modelType, input_type);`

Main part of the process. 

**Inputs:**

- `testData`: File with test data
- `trainingData`: File(s) with training data
- `dataLookup`: Lookup table file with list of data filenames and information such as corresponding orthologs files, synergy cutoff and antagonism cutoff for classification of interaction scores.
- `valMethod`: Either Kfold cross validation, holdout, or independent
- `K`: The number of K subsets/iterations you want to run for cross validation
- `standardize`: Set to `z_score` to convert scores to z-scores
- `modelType`: Either run the E. coli model/original model or the M. Tb model (2). The E. coli model uses E. coli chemogenommics and is trained on E. coli interaction scores to start. The original model also uses E. coli data, but only uses data from [this original paper on INDIGO](https://www.embopress.org/doi/full/10.15252/msb.20156777). The M. Tb model uses E. coli chemogenomics and M. Tb transcriptomics data and is trained on M. Tb interaction scores to start.
- `input_type`: either 1 or 2, default is 2 to make predictions for drug combinations

**Outputs:**

`indigoSummary`: This saves the input settings and contains the predicted scores, which are further analyzed in `analyze`.

3. `[stats, averages, overview] = analyze(indigoSummary,resultIndex, dataFiles))`

Takes output from indigoRun and analyzes results. Determines rank correlation, accuracy, etc. Plots ROC curves, confusion matrices and scatter plots of results (rank sorted). Saves these results to an excel file.

4. `save(saveIndigo(indigoSummary))`

Saves workspace in .mat file. Names .mat file based on test data, training data, validation method and standardization parameter for that run.
