function resultsFile = saveIndigo(indigoSummary)
%save results
dataName = erase(indigoSummary.testFile,'.xlsx');
if exist('indigoSummary.trainingData','var')
    indigoSummary.trainingData = erase(indigoSummary.trainingData,'.xlsx');
    if length(indigoSummary.trainingData) == 6
        indigoSummary.trainingData = 'nature';
    elseif length(indigoSummary.trainingData) == 2
        indigoSummary.trainingData = 'nature_ecoli';
    end

    if indigoSummary.standardized == 1
        resultsFile = strcat('results/',dataName,sprintf('_%s_%s_z.mat', ...
        indigoSummary.trainingData, indigoSummary.valMethod));
    else
        resultsFile = strcat('results/',dataName,sprintf('_%s_%s.mat', ...
        indigoSummary.trainingData, indigoSummary.valMethod));
    end
else
    if indigoSummary.standardized == 1
        resultsFile = strcat('results/',dataName,sprintf('_%s_z.mat', ...
        indigoSummary.valMethod));
    else
        resultsFile = strcat('results/',dataName,sprintf('_%s.mat', ...
        indigoSummary.valMethod));
    end
end