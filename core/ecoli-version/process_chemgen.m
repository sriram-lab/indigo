function [phenotype_data, phenotype_labels, predictor_names, conditions] = process_chemgen(fname,z)

    % DESCRIPTION 
    % This function processes chemogenomic data to form a binary matrix
    % indicating conditions where KO strains were sensitive or resistant
    % to treatment. 
    % 
    % STEPS 
    % 1. Input processing
    % 2. Load data
    % 3. Convert gene IDs to standard IDs
    % 4. Transform data
    % 5. Define output phenotype variables
    % 
    % Author:   Murat Cokol
    % Created:  October 23, 2018
    % Updates:  September 17, 2020 (Carolina H. Chung)
    %           August 27, 2020 (Carolina H. Chung)

    % I/O
    %{
    OPTIONAL INPUTS: 
        1. fname:   filename for chemogenomic data
        2. z:       threshold value to define significant effect on fitness
                    (default: z = 2)
    
    OUTPUTS:
        1. phenotype_data:      binary matrix containing 1 where E. coli
                                was sensitive or resistant to a treatment
        2. phenotype_labels:    labels (i.e. genes) for phenotype data
        3. conditions:          list of conditions (i.e. treatments)
    %}

    %% INPUT PROCESSING
    if isempty(fname)
        fname = 'ecoli_phenotype_data_cell.xlsx';
    end
    
    if ~exist('z','var') || isempty(z)
        z = 2;
    end
    %% LOAD DATA
    data = readtable(fname,'VariableNamingRule','preserve');
    phenotype_num = data{1:end,2:end};    % numerical data
    probelist = data.Gene;                % gene names (ECK numbers)
    conditions = data.Properties.VariableNames(2:end)';  %list of conditions

    %% CONVERT GENE IDs TO STANDARD IDs
    load ('ecoli_annotation_data1','genenames_array', 'genenames_bnums')
  
    clear plist
    plist = cell(size(probelist));
    for i = 1:length(probelist)
        tem = regexp(probelist{i},'ECK[\d]*-','split');
        try
            plist(i) = tem(2);
        catch
            plist(i) = tem(1);
        end
    end
    plist = regexprep(plist,'''','');
    [ix, pos] = ismember(upper(plist),upper(genenames_array));
    plist_bnums = plist;
    plist_bnums(ix) = genenames_bnums(pos(ix));
    plist_bnums = matlab.lang.makeUniqueStrings(plist_bnums);
    
    % TRANSFORM DATA
    % Sensitive strains
    sensitive_phenotype_num = phenotype_num < -z;
    idx = sum(sensitive_phenotype_num, 2) == 0;
    sensitive_data = sensitive_phenotype_num(~idx,:);
    sensitive_labels = plist_bnums(~idx);
    
    % Resistant strains
    resistant_phenotype_num = phenotype_num > z; 
    idx = sum(resistant_phenotype_num, 2) == 0; 
    resistant_data = resistant_phenotype_num(~idx,:);
    resistant_labels = plist_bnums(~idx);

    % DEFINE OUTPUT PHENOTYPE VARIABLES
    % Sets up row length of sigma delta matrix.
    phenotype_data = [sensitive_data; resistant_data];
    phenotype_labels = [sensitive_labels; resistant_labels];
    
    predictor_names_se = strcat(sensitive_labels, repmat({'_se'},length(sensitive_labels),1));
    predictor_names_r = strcat(resistant_labels, repmat({'_r'},length(resistant_labels),1));
    predictor_names = [predictor_names_se; predictor_names_r];
    predictor_names = [strcat(predictor_names, repmat({'_si'},length(predictor_names),1)); ...
        strcat(predictor_names, repmat({'_d'},length(predictor_names),1))];
end
