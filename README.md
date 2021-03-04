# INDIGO
This is a readme for indigo repository!

INDIGO is a random forest model that makes drug interaction predictions using drug chemogenomics/transcriptomics and drug interaction data. It runs in MATLAB.

I've modified INDIGO so that it can be trained on data from multiple studies representing multiple bacterial strains/species.

**Google Drive folder containing Data and Results:** https://drive.google.com/drive/folders/13nP4fNvrQv_s-81Li1-Nx59m8m2U8zRT?usp=sharing

Repository contains following folders:

1. `core`
2. `python-scripts`
4. `testing`
5. `validation`

Code required to run INDIGO model is in `validation`. Explanation of how to use code is below.

The following are different functions that can be used.

1. `[summary, drug_list] = data_explore(filename, data_label, data_lookup, number_ref))`

Takes in the file of interest and returns a summary of descriptive statistics of the data, including synergy and antagonism counts. Also gets a list of unique drugs. Graphs the data to show you distribution of scores.

2. `indigo_summary = indigo_run(test_data, training_data, data_lookup, valmethod, K, standardize, model_type, input_type, scoring)`

Main part of the process.

**Inputs:**

- `test_data`: File with test data
- `training_data`: File(s) with training data
- `data_lookup`: Lookup table file with list of data filenames and information such as corresponding orthologs files, synergy cutoff and antagonism cutoff for classification of interaction scores.
- `valmethod`: Either `Kfold`, `independent`, `holdout`, `Kfold_onself`, `holdout_onself`
- `K`: The number of K subsets/iterations you want to run for cross validation
- `standardize`: Set to `z_score` to convert scores to z-scores
- `model_type`: Either run the E. coli model/original model or the M. Tb model (2). The E. coli model uses E. coli chemogenomics and is trained on E. coli interaction scores to start. The original model also uses E. coli data, but only uses data from [this original paper on INDIGO](https://www.embopress.org/doi/full/10.15252/msb.20156777). The M. Tb model uses E. coli chemogenomics and M. Tb transcriptomics data and is trained on M. Tb interaction scores to start.
- `input_type`: either 1 or 2, default is 2 to make predictions for drug combinations
- `scoring`: Either `''`, `'Bliss'`, or `'Loewe'`

**Outputs:**

`indigo_summary`: This saves the input settings and contains the predicted scores, which are further analyzed in `analyze`.

3.  `[stats,averages,overview] = analyze(indigo_summary, result_index, indigo_data, data_files, prediction_idx)`

Takes output from indigo_run and analyzes results. Determines rank correlation, accuracy, etc. Plots ROC curves, confusion matrices and scatter plots of results (rank sorted). Saves these results to excel files.

4. `save(save_indigo(indigo_summary))`

Saves workspace in .mat file. Names .mat file based on test data, training data, validation method and standardization parameter for that run.