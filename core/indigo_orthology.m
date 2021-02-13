function [deviations, teststaphdat1] = indigo_orthology(phenotype_labels,ecoli_staph_orth,sigma_delta_input,indigo_model)

    %{
    DESCRIPTION

    This function calculates interaction score deviations due to species
    difference (when orthology information is used)

    STEPS
    1. Get E. coli model and orthologs
    2. Determine predicted variable interactions between species

    Author: Murat Cokol
    Created: October 23, 2018
    Updates: August 27, 2020 (Carolina H. Chung)
             January 23, 2021 (David Chang)


    I/O
    
    REQUIRED INPUTS:
      1. phenotype_labels:    labels for phenotype_data (i.e. genes) 
      2. ecoli_staph_orth:    list of ortholog genes b/t species
      3. sigma_delta_input:   numeric matrix of combination profiles
      4. magenta_model:       trained INDIGO model
    
    OUTPUTS:
        1. deviations:        deviations in predicted scores
        2. teststaphdat1:     sigma delta scores after accounting 
                              for non-orthologs. Can be set to 0 or 2,
                              effect is the same.
    %}
  
    %% GET E. COLI MODEL AND ORTHOLOGS
    % Find non-orthologs
    nic_row = [phenotype_labels; phenotype_labels];
    ix = ismember(nic_row, ecoli_staph_orth); % find the orthologs
    nonorthtop = nic_row(~ix);                % get the non-orthologs
%     teststaphdat = sigma_delta_input; 
%     teststaphdat(~ix,:) = 0;                  % modfiy state of non-orthologs
    % Set the sigma scores to be 2 or 0
    ix2 = ismember(phenotype_labels,nonorthtop);
    ix2 = find(ix2);
    teststaphdat1 = sigma_delta_input;
    teststaphdat1(ix2,:) = 2;   
    teststaphdat1(ix2 + length(phenotype_labels),:) = 0; 
    %% DETERMINE PREDICTED VARIABLE INTERACTIONS B/T SPECIES
    testpredictions_staphchem_ecolixns2 = predict(indigo_model,teststaphdat1');
    testpredictions_staphchem_ecolixns20 = predict(indigo_model,sigma_delta_input');
%     testpredictions_staphchem_ecolixns21 = predict(indigo_model,teststaphdat');
%     deviations = testpredictions_staphchem_ecolixns20(:) - testpredictions_staphchem_ecolixns21(:); % output the deviations for the input drugs
    % output the deviations for the input drugs
    deviations = testpredictions_staphchem_ecolixns20(:) - testpredictions_staphchem_ecolixns2(:); 
end