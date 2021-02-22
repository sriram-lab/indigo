function results_file = save_indigo(indigo_summary, prediction_idx)
    %{
    DESCRIPTION

    Format name of .mat file to save workspace to.
    Function will be called in save() as follows:
    save(save_indigo(indigo_summary)
    
    Uses the following from indigoSummary to define filepath and filename:
        - model_type
        - standardize
        - val_method
        - scoring
  
    I/O
    
    REQUIRED INPUTS: 
        1.  indigo_summary:      Model parameters and results 
    %}
    if prediction_idx == 1
        directory = 'v3_1';
    elseif prediction_idx == 2
        directory = 'v3_2';
    end

    filename = strcat(erase(indigo_summary.test_data,'.xlsx'),'_',indigo_summary.valmethod);

    if strcmp(indigo_summary.standardize,'z_score')
        filename = sprintf('%s_z', filename);    
    end
    if strcmp(indigo_summary.scoring,'bliss') || strcmp(indigo_summary.scoring,'loewe')
        results_file = sprintf(strcat('results/%s/', indigo_summary.model_type, '/', ...
            indigo_summary.scoring, '/', indigo_summary.valmethod, '/', filename,'.mat'), directory);
    else
        results_file = sprintf(strcat('results/%s/',indigo_summary.model_type, '/', ...
            indigo_summary.valmethod, '/', filename,'.mat'), directory);
    end
end