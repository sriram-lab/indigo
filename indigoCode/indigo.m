%make sure you are starting in folder one step above code, data and results
%folders
testFile = 'tb_pairwise.xlsx';
trainingData = 'orig_ecoli.xlsx';
valMethod = 'cv';
K = 5;
standardize = '';
input_type = 2;

[summary,drugList] = dataExplore(testFile);
compareFiles = {testFile,'orig_ecoli.xlsx'};
sharedInteractions = dataCompare(compareFiles);

%run indigo on data
indigoSummary = indigoRun(testFile,trainingData,valMethod,K, standardize,input_type);

%analysis
[stats,averages,overview] = analyze(indigoSummary);
%save results
save(saveIndigo(indigoSummary))