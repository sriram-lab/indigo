function [test_interactions, testinteractions_scores, ...
    indigo_model, sigma_delta_scores, conditions] = ...
    indigo_predict_tb(indigo_model, testdata, input_type, ...
    annotation_filename, chemogenomics_filename, z, ...
    phenotype_data, phenotype_labels, conditions)

    %{
    DESCRIPTION
    This function generates interaction score predictions based on a
    trained INDIGO model.
     
    STEPS
    1. Input processing
    2. Convert and match interaction labels with chemogenomic data
    3. Filter out interactions without chemogenomic data
    4. Define inputs to INDIGO and generate predictions
    
    Author: Murat Cokol
    Created: October 23, 2018
    Updates: Obctober 21, 2020 (Carolina H. Chung)
             August 27, 2020 (Carolina H. Chung)
             January 18, 2021 (David Chang)
    
    I/O
    
    REQUIRED INPUTS: 
        1. indigo_model:            trained MAGENTA model
        2. test_data:               cell array of drug names or interaction
                                    pair names
           -> note: if drug list is provided, then all pairwise drug
                    combinations are determined
        3. input_type:              specify based on test_data provided
                                    (input_type = 1 if drug list)
                                    (input_type = 2 if interaction pairs)
        4. annotation_filename:     filename for matching drug names to 
                                    chemogenomic condition names
        5. chemogenomics_filename:  filename for chemogenomic data
    OPTIONAL INPUTS: 
        1. z:                   threshold value to define significant 
                                effect on fitness (default: z = 2)
        2. phenotype_data:      binary matrix defining sensitive and 
                                resistant phenotypes
        3. phenotype_labels:    labels (i.e. genes) for phenotype data
        4. conditions:          list of conditions (i.e. treatments)
    
    OUTPUTS:
        1. test_interactions:   interaction pairs used for model training
                                (i.e. chemogenomic data available)
        2. testxnscores:        interaction scores corresponding to 
                                train_interactions
        3. magenta_model:       trained MAGENTA model 
        4. sigma_delta_scores:  numerical matrix of combination profiles
        5. ix:                  logical array indicating which entries in 
                                test_data were used for 
        6. phenotype_labels:    labels (i.e. genes) for phenotype data
    %}

    %% INPUT PROCESSING
    if input_type == 1          % drug list mode
        testdrugs = testdata;
    elseif input_type == 2      % interaction pairs mode
        test_interactions1 = testdata;
    else
        error('incorrect input: input_type is 1 (drug) or 2 (interaction)');
    end
    if ~exist('z','var') || isempty(z)
        z = 2;
    end
    if ~exist('phenotype_data','var') || isempty(phenotype_data)
        [phenotype_data, phenotype_labels, conditions] = process_chemgen_tb(chemogenomics_filename,z);
    end

    %% CONVERT AND MATCH INTERACTION LABELS WITH CHEMOGENOMIC DATA
    txt = readcell(annotation_filename);
    [drugxn_id, chemgen_id] = deal(txt(:,1),txt(:,2));

    % Define interaction pairs
    if input_type == 1      % drug list
        drugnames_cell = testdrugs;
        for i = 1:length(drugxn_id)
            drugnames_cell (ismember(drugnames_cell,drugxn_id(i))) = chemgen_id(i);
        end
        drugxn_id1 = unique([drugxn_id;testdrugs]);
        alldrugcombinations = drugxn_id1(nchoosek(1:length(drugxn_id1),2));
        ix = ismember(alldrugcombinations, testdrugs); ix = any(ix,2);
        drugpairsname_cell = alldrugcombinations(ix,:);
        for i = 1:length(drugxn_id)
            drugpairsname_cell(ismember(drugpairsname_cell,drugxn_id(i))) = chemgen_id(i);
        end
    elseif input_type == 2   % interaction pairs
        drugpairsname_cell = test_interactions1;
        for i = 1:length(drugxn_id)
            drugpairsname_cell(ismember(drugpairsname_cell,drugxn_id(i))) = chemgen_id(i);
        end
    end
    
    %% FILTER OUT INTERACTIONS WITHOUT CHEMOGENOMIC DATA
    ix = ismember(drugpairsname_cell, conditions); 
    ix = (sum(~cellfun(@isempty, drugpairsname_cell), 2)) == sum(ix, 2);
    test_interactions = drugpairsname_cell(ix,:); 
    testdrugs = unique(test_interactions); 
    testdrugs(cellfun(@isempty, testdrugs)) = []; 
    
    %% DEFINE INPUTS FOR INDIGO AND GENERATE PREDICTIONS
    [~, pos] = ismember(testdrugs, conditions);
    testchemgen = phenotype_data(:,pos(logical(pos)));
    [testinteractions_scores, indigo_model, sigma_delta_scores] = ...
        indigo_rf_3_tb_ensemble([], [], [], [], testdrugs, testchemgen, ...
        test_interactions, 2, indigo_model);

 end