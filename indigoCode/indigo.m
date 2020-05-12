%make sure you are starting in folder one step above code, data and results
%folders
testFile = 'ecoli_yeh.xlsx';
valMethod = 'cv';
K = 5;
orthology = '';
trainingData = 'nature';
standardize = '';

[summary,drugList] = dataExplore(testFile);
files = dataFiles();
compareFiles = [testFile;files.natureData];
sharedInteractions = dataCompare(compareFiles);

%run indigo on data
indigoSummary = indigoRun(testFile,valMethod,K,orthology,trainingData,standardize);

%analysis
[stats,averages,overview] = analyze(indigoSummary);

%save results
save(saveIndigo(indigoSummary))