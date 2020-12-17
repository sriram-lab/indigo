function resultsFile = saveIndigo(indigoSummary)
%save results
dataName = erase(indigoSummary.testData,'.xlsx');

if indigoSummary.standardized == 0
    resultsFile = strcat('results/current/','ecoli-model-final/cv/',dataName,sprintf('_%s.mat', indigoSummary.valMethod));
else
    resultsFile = strcat('results/current/','ecoli-model-final/cv-z/',dataName,sprintf('_%s_z.mat', indigoSummary.valMethod));
end