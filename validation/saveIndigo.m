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

    filename = strcat(erase(indigoSummary.testData,'.xlsx'),'_',indigoSummary.valMethod);
    if strcmp(indigoSummary.standardize,'z_score')
        filename = sprintf('%s_z', filename);
    end
     
    if strcmp(indigoSummary.scoring,'bliss') || strcmp(indigoSummary.scoring,'loewe')
        resultsFile = strcat('results/v2/', indigoSummary.modelType, '/', ...
            indigoSummary.scoring, '/', indigoSummary.valMethod, '/', filename, '.mat');
    else
        resultsFile = strcat('results/v2/',indigoSummary.modelType, '/', ...
            indigoSummary.valMethod, '/', filename, '.mat');
    end
    
end