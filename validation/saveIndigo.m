function resultsFile = saveIndigo(indigoSummary)
%save results
dataName = erase(indigoSummary.testData,'.xlsx');

if indigoSummary.standardized == 0
    resultsFile = strcat('results/','mtb-model/bliss/cv/',dataName,sprintf('_%s.mat', indigoSummary.valMethod));
else
    resultsFile = strcat('results/','mtb-model/bliss/cv-z/',dataName,sprintf('_%s_z.mat', indigoSummary.valMethod));
end