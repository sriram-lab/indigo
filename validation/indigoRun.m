function indigoSummary = indigoRun(testData, trainingData, dataLookup, ... 
    valMethod, K, standardize,modelType,input_type, scoring)
arguments
    testData char
    trainingData = ''
    dataLookup = ''
    valMethod char {mustBeMember(valMethod,{'Kfold', 'holdout', ...
                    'holdout_onself', 'Kfold_onself', ...
                    'independent'})} = 'Kfold'
    K {mustBeInteger} = 5
    standardize char {mustBeMember(standardize,{'','z_score'})} = ''
    modelType char {mustBeMember(modelType,{'ecoli_model','mtb_model','original_model'})} = 'ecoli_model';
    input_type {mustBeInteger} = 2;
    scoring char {mustBeMember(scoring,{'bliss', 'loewe', ''})} = ''
end

    %{
    DESCRIPTION

    This function runs INDIGO and returns a model with prediction results.

    STEPS
    1. Input processing for test and training data.
    2. Format interaction scores and sigma delta scores matrices.
    3. Modify sigma delta scores based on orthology.
    4. Build random forest INDIGO model.
    5. Make predictions.

    Author: David Chang
    Created: January 23, 2021

    I/O
    
    REQUIRED INPUTS:
      1. testData:          List of drugs (1) or drug interactions with 
                            scores (2) 

    OPTIONAL INPUTS:
      2. trainingData:      Drug interaction data, can be from multiple
                            species. If model is E. coli, first file must 
                            be E. coli data. If model is M. tb, first file 
                            must be either m. tb or e. coli. Can be single
                            file string or cell array of files.
      3. valMethod:         Validation method, can be the following options:
                                1. holdout_onself: Only use test data
                                2. cv_onself: Only use test data
                                3. independent: Leave all test data out of
                                training
                                4. cv: K-fold cross validation
      4. K:                 Cross-validation parameter 
                            (default K = 5 for Kfold cross-validation)
      5. standardize:       To get zscores or not
      6. modelType:         E. coli model, M. tb model, or original model. 
                            Switches training data and chemogenomics/transcriptomics 
                            data that is used. Original model means using
                            only dataset E. coli dataset from
                            Chandrasekaran et al., 2016.
      7. input_type:        1 for drug list, 2 for drug combinations
      8. scoring            Interaction scoring method. Can be Bliss or
                            Loewe, default is none ('') which means model
                            uses data with both scoring methods.
    OUTPUTS:
      1. indigoSummary:     INDIGO model, predicted scores and input data
    %}

% Initialize output indigoSummary
indigoSummary = struct;

indigoSummary.modelType = modelType;
if strcmp(modelType, 'ecoli_model') || strcmp(modelType, 'original_model')
    %E. coli chemogenomics data
    annotation_file = 'identifiers_match.xlsx';
    chemogenomics_file = 'ecoli_phenotype_data_cell.xlsx';
   
elseif strcmp(modelType, 'mtb_model')
    % M. tb chemogenomics data, scripts have _tb at end
    annotation_file = 'identifiers_match_tb.xlsx';
    chemogenomics_file = [];
    file = 'averaged_data_new.mat';
    S = load(file);
    [normData,col,row] = tb_preprocessing(S.averagedtbdata, ...
                         S.mtb_expression_database_col_ids, ...
                         S.mtb_expression_database_row_ids);
    [phenotype_data, phenotype_labels] = process_transcriptome_tb(normData,row,col);
end

indigoSummary.scoring = scoring;

indigoSummary.testData = testData;
indigoSummary.trainingData = trainingData;
indigoSummary.dataLookup = dataLookup;

% Read in the test data
% Use readtable to handle missing values as {0x0 char}
data = readtable(testData,"ReadVariableNames",false);
scores = data{:,end};
interactions = data{:,1:end-1};

indigoSummary.standardize = standardize;
if strcmp(standardize,'z_score')
    scores = zscores(scores);
end

%get orthologs for test file
testOrthologs = get_orthologs(testData, modelType, dataLookup);

%{ 
COMBINE TRAINING DATA FROM DIFFERENT DATA FILES TO GET OVERALL SIGMA DELTA 
SCORE MATRIX (X) AND INTERACTION SCORE MATRIX (Y) TO BUILD THE MODEL.
%}
if ~isempty(trainingData)
    indigoSummary.trainingData = trainingData;
    interactions_all = [];
    interaction_scores_all = [];
    sigma_delta_scores_all = [];
    
    if ischar(class(trainingData)) || isstring(class(trainingData))
        %convert to cell, accounts for if trainingData is single file
        trainingData = cellstr(trainingData);
    end
    
    for i = 1:length(trainingData)
        fprintf(sprintf('Adding %s to training data\n',trainingData{i}))
        if i == 1
            %Train first
            if strcmp(modelType, 'ecoli_model')
                [train_interactions, train_scores, labels, indigo_model,...
                 sigma_delta_scores, ~] = indigo_train(trainingData{i},standardize, ...
                 annotation_file,chemogenomics_file);
            elseif strcmp(modelType, 'mtb_model')
                [train_interactions, train_scores, labels, indigo_model,...
                 sigma_delta_scores, ~] = indigo_train_tb(trainingData{i},standardize, ...
                 annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
            end
            trainOrthologs = get_orthologs(trainingData{i},modelType);       
            
            if ~isempty(trainOrthologs)
                [~,sigma_delta_scores] = indigo_orthology(labels, ...
                 trainOrthologs,sigma_delta_scores, indigo_model); 
            end
            
        else
            % Use readtable to handle missing values as {0x0 char}
            data = readtable(trainingData{i}, "ReadVariableNames", false);
            train_scores = data{:,end};
            train_interactions = data{:,1:end-1};
            
            if strcmp(standardize, 'standardize')
                train_scores = zscore(train_scores);
            end
            
            if strcmp(modelType, 'ecoli_model')
                [~,~,~,sigma_delta_scores] = indigo_predict(indigo_model,train_interactions, ...
                  input_type,annotation_file,chemogenomics_file);
            elseif strcmp(modelType,'mtb_model')
                [~,~,~,sigma_delta_scores] = indigo_predict_tb(indigo_model,train_interactions, ...
                 input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
            end
            
            trainOrthologs = get_orthologs(trainingData{i},modelType);  
            
            if ~isempty(trainOrthologs)
                [~,sigma_delta_scores] = indigo_orthology(labels, ...
                 trainOrthologs,sigma_delta_scores, indigo_model); 
            end
            
        end
        
        % Adjusting # of columns of training data matrix so that all the 
        % different training data can be combined before building the model
        if size(interactions_all,2) > size(train_interactions,2)
            num_add = size(interactions_all,2) - size(train_interactions,2);
            empty_array = repmat("",length(train_interactions),num_add);
            train_interactions = [train_interactions empty_array];
        elseif size(interactions_all,2) < size(train_interactions,2)
            num_add = size(train_interactions,2) - size(interactions_all,2);
            empty_array = repmat("",length(interactions_all),num_add);
            interactions_all = [interactions_all empty_array];
        end
        
        interaction_scores_all = [interaction_scores_all; train_scores];
        sigma_delta_scores_all = [sigma_delta_scores_all, sigma_delta_scores]; 
        interactions_all = [interactions_all; train_interactions];
    end
end

%% BUILD MODEL AND MAKE PREDICTIONS WITH VALIDATION METHOD OF YOUR CHOOSING
indigoSummary.valMethod = valMethod;

if strcmp(valMethod,'holdout_onself')
    i = 1;
    % K = 0.2 is usually default --> 20% of data is in test set
    [train,test] = crossvalind('HoldOut',length(interactions), K); 
    Xtrain = interactions(train,:);
    Ytrain = scores(train);
    Xtest = interactions(test,:);
    Ytest = scores(test);
    
    % Plug training data into indigo_train
    writecell([Xtrain,num2cell(Ytrain)],'train.xlsx')
    
    if strcmp(modelType, 'ecoli_model')
        [~, ~, labels, indigo_model,...
         sigma_delta_scores, ~] = indigo_train('train.xlsx', standardize,...
         annotation_file,chemogenomics_file);
    elseif strcmp(modelType, 'mtb_model')
        [~, ~, labels, indigo_model,...
         sigma_delta_scores, ~] = indigo_train_tb('train.xlsx',standardize, ...
         annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
    end
    
    % Adjust sigma delta scores based on presenece of orthologous and 
    % nonorthologous genes
    if ~isempty(testOrthologs)
        [~,sigma_delta_scores] = indigo_orthology(labels, testOrthologs, ...
                                 sigma_delta_scores, indigo_model); 
    end

    % Build the model
    tic
    indigo_model = fitrensemble(single(sigma_delta_scores'), ...
                   single(Xtrain),'Method','Bag');
    toc
    
    
    % Make predictions and store results
    predictStep();   
    
elseif strcmp(valMethod, 'holdout')
    i = 1;
    % Add 80% of test set to training and keep 20% as test
    [train,test] = crossvalind('HoldOut',length(scores),0.2); 
    Xtrain = interactions(train,:);
    Ytrain = scores(train);
    Xtest = interactions(test,:);
    Ytest = scores(test);
    
    if strcmp(modelType, 'ecoli_model')
        [~,~,~,sigma_delta_scores] = indigo_predict(indigo_model,Xtrain, ...
         input_type,annotation_file,chemogenomics_file);
    elseif strcmp(modelType, 'mtb_model')
        [~,~,~,sigma_delta_scores] = indigo_predict_tb(indigo_model,Xtrain, ...
         input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
    end


    if ~isempty(testOrthologs)
        [~,sigma_delta_scores] = indigo_orthology(labels, testOrthologs, ...
                                 sigma_delta_scores, indigo_model); 
    end

    interaction_scores = [interaction_scores_all; Ytrain];
    sigma_delta_scores = [sigma_delta_scores_all, sigma_delta_scores];
   
    tic
    indigo_model = fitrensemble(single(sigma_delta_scores'), ...
                   single(interaction_scores),'Method','Bag');
    toc
    
    predictStep();
    
elseif strcmp(valMethod,'independent')
    i = 1;
    
    % Only need to do this if you have more than one file in training data
    % or if the training data has orthologs
    if length(trainingData) > 1 || ~isempty(trainOrthologs)
        tic
        indigo_model = fitrensemble(single(sigma_delta_scores_all'), ...
                   single(interaction_scores_all),'Method','Bag');
        toc
    end
    
    Xtrain = interactions_all;
    Ytrain = interaction_scores_all;
    Xtest = interactions;   
    Ytest = scores;       

    predictStep();
    
elseif strcmp(valMethod,'Kfold_onself')
    indigoSummary.K = K;
    for i = 1:K
        fprintf('Run %d\n',i)
        idx = crossvalind('Kfold',length(scores),K);
        test = (idx == i);
        train = ~test;
        Xtrain = interactions(train,:);
        Ytrain = scores(train);
        Xtest = interactions(test,:);
        Ytest = scores(test);
        writecell([Xtrain,num2cell(Ytrain)],'train.xlsx')
        
        if strcmp(modelType, 'ecoli_model')
            [~, ~, labels, indigo_model,...
             sigma_delta_scores, ~] = indigo_train('train.xlsx',standardize, ...
             annotation_file,chemogenomics_file);
        elseif strcmp(modelType, 'mtb_model')
            [~, ~, labels, indigo_model,...
             sigma_delta_scores, ~] = indigo_train_tb('train.xlsx',standardize, ...
             annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
        end
        
        if ~isempty(testOrthologs)
            [~,sigma_delta_scores] = indigo_orthology(labels, testOrthologs, ...
                                     sigma_delta_scores, indigo_model); 
        end
        
        tic
        indigo_model = fitrensemble(single(sigma_delta_scores'), ...
                   single(Ytrain),'Method','Bag');
        toc
        
        predictStep();   
    end
    
elseif strcmp(valMethod,'Kfold')
    indigoSummary.K = K;
    for i = 1:K
        fprintf('Run %d\n',i)
        idx = crossvalind('Kfold',length(scores),K);
        test = (idx == i);
        train = ~test;
        Xtrain = interactions(train,:);
        Ytrain = scores(train);
        Xtest = interactions(test,:);
        Ytest = scores(test);
        
        if strcmp(modelType, 'ecoli_model')
            [~,~,~,sigma_delta_scores] = indigo_predict(indigo_model,Xtrain, ...
             input_type,annotation_file,chemogenomics_file);
        elseif strcmp(modelType, 'mtb_model')
            [~,~,~,sigma_delta_scores] = indigo_predict_tb(indigo_model,Xtrain, ...
             input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
        end
        
        if ~isempty(testOrthologs)
            [~,sigma_delta_scores] = indigo_orthology(labels, testOrthologs, ...
                                     sigma_delta_scores, indigo_model); 
        end

        interaction_scores = [interaction_scores_all; Ytrain];
        sigma_delta_scores = [sigma_delta_scores_all, sigma_delta_scores];
        
        tic
        indigo_model = fitrensemble(single(sigma_delta_scores'), ...
                       single(interaction_scores),'Method','Bag');
        toc
        
        predictStep();
    end
    
end

% nested function for predicting scores and storing results
function predictStep()
    
    if strcmp(modelType, 'ecoli_model')
         [~,predicted_scores,~,sigma_delta_input] = indigo_predict(indigo_model, ...
          Xtest,input_type,'identifiers_match.xlsx','ecoli_phenotype_data_cell.xlsx');
    elseif strcmp(modelType, 'mtb_model')
         [~,predicted_scores,~,sigma_delta_input] = indigo_predict_tb(indigo_model,Xtest, ...
          input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
    end
    
    if ~isempty(testOrthologs)
        deviations = indigo_orthology(labels, testOrthologs, ... 
                     sigma_delta_input, indigo_model); 

        predicted_scores = predicted_scores - deviations;
    end

    indigoSummary.model{i} = indigo_model;
    indigoSummary.trainPairs{i} = Xtrain;
    indigoSummary.trainScores{i} = Ytrain;
    indigoSummary.sigma_delta_scores{i} = sigma_delta_scores;
    indigoSummary.testPairs{i} = Xtest;
    indigoSummary.testScores{i} = Ytest;
    indigoSummary.predictedScores{i} = predicted_scores;
    fprintf(sprintf('INDIGO predictions subset %d complete!\n',i))
end

if isfile('train.xlsx')
    delete train.xlsx
end

end