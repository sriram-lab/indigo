function results_file = save_indigo(indigo_summary)
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

    filename = strcat(erase(indigo_summary.test_data,'.xlsx'),'_',indigo_summary.valmethod);
    if strcmp(indigo_summary.standardize,'z_score')
        filename = sprintf('%s_z', filename);
    end
     
    if strcmp(indigo_summary.scoring,'bliss') || strcmp(indigo_summary.scoring,'loewe')
        results_file = strcat('results/v2/', indigo_summary.model_type, '/', ...
            indigo_summary.scoring, '/', indigo_summary.valmethod, '/', filename, '.mat');
    else
        results_file = strcat('results/v2/',indigo_summary.model_type, '/', ...
            indigo_summary.valmethod, '/', filename, '.mat');
    end
    
end