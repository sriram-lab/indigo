function resultsFile = saveIndigo(indigoSummary)
    %{
    DESCRIPTION

    Format name of .mat file to save workspace to.
    Function will be called in save() as follows:
    save(saveIndigo(indigoSummary)
    
    Uses the following from indigoSummary to define filepath and filename:
        - modelType
        - standardized
        - valMethod
        - scoring
  
    I/O
    
    REQUIRED INPUTS: 
        1.  indigoSummary:      Model parameters and results 
    %}

    dataName = strcat(erase(indigoSummary.testData,'.xlsx'),'_',indigoSummary.valMethod);
    if strcmp(indigoSummary.standardize,'z_score')
        dataName = sprintf('%s_z', dataName);
    end
    if strcmp(indigoSummary.scoring,'bliss') || strcmp(indigoSummary.scoring,'loewe')
        resultsFile = strcat(indigoSummary.modelType, '/', indigoSummary.scoring, ...
            '/', indigoSummary.valMethod, '/', dataName, '.mat');
    else
        resultsFile = strcat(indigoSummary.modelType, '/', ...
            indigoSummary.valMethod, '/', dataName, '.mat');
    end
end