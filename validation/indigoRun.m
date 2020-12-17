function indigoSummary = indigoRun(testData,trainingData,valMethod,K,standardize,modelType,input_type)
arguments
    testData char
    trainingData = ''
    valMethod char {mustBeMember(valMethod,{'holdout_onself', 'cv_onself', 'independent', 'cv'})} = 'holdout_onself'
    K {mustBeInteger} = 5
    standardize char {mustBeMember(standardize,{'','standardized'})}= ''
    modelType {mustBeInteger} = 1;
    input_type {mustBeInteger} = 2;
end

%[indigoSummary] = indigoRun(testData,trainingData,valMethod,K, standardize,input_type)
% This function runs indigo and returns a summary of the results
%required argument:
%testData - can be list of drugs (1) or drug combinations (2)
%optional arguments:
%trainingData - single file or cell array of files
%valMethod - name-value argument, options are holdout_onself, cv_onself, independent, cv
%K - fold parameter for K-fold cross validation, default = 5
%standardize - whether data will be standardized or not
%input type, 1 (list of drugs), 2 (drug combos). For 1, use independent
%validation
%modelType 1 == ecoli, 2 == mtb

%scores = drug interaction scores
%interactions = drug pairs
%sigma_delta_scores = sigma delta scores

%Xtrain, Xtest, Ytrain, Ytest - independent/dependent variables during
%holdout or cross validation 

%Initialize output indigoSummary
indigoSummary = struct;

if modelType == 1
    %ecoli chemogenomics data
    indigoSummary.modelType = 1;
    annotation_file = 'identifiers_match.xlsx';
    chemogenomics_file = 'ecoli_phenotype_data_cell.xlsx';
   
elseif modelType == 2
    %mtb chemogenomics data, scripts have _tb at end
    indigoSummary.modelType = 2;
    annotation_file = 'identifiers_match_tb.xlsx';
    chemogenomics_file = [];
    file = 'averaged_data_new.mat';
    S = load(file);
    [normData,col,row] = tb_preprocessing(S.averagedtbdata, ...
                         S.mtb_expression_database_col_ids, ...
                         S.mtb_expression_database_row_ids);
    [phenotype_data, phenotype_labels] = process_transcriptome_tb(normData,row,col);

end

indigoSummary.testData = testData;

sheet = 1;  %for non-standardized data
indigoSummary.standardized = 0;
if strcmp(standardize,'standardized')
    sheet = 2;  %for z score data
    indigoSummary.standardized = 1;
end

[scores,interactions] = xlsread(testData,sheet);

%get orthologs for test file
[testOrthologs,orthologyFile] = get_orthologs(testData,modelType);

if ~isempty(testOrthologs)
    indigoSummary.orthology = orthologyFile;
end


%setting up training data
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
            if modelType == 1
                [train_interactions, train_scores, labels, indigo_model,...
                 sigma_delta_scores, conditions] = indigo_train(trainingData{i},sheet, ...
                 annotation_file,chemogenomics_file);
            elseif modelType == 2
                [train_interactions, train_scores, labels, indigo_model,...
                 sigma_delta_scores, conditions] = indigo_train_tb(trainingData{i},sheet, ...
                 annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
            end
            trainOrthologs = get_orthologs(trainingData{i},modelType);       
            
            if ~isempty(trainOrthologs)
                [~,sigma_delta_scores] = indigo_orthology(labels, ...
                 trainOrthologs,sigma_delta_scores, indigo_model); 
            end
            
        else
            [train_scores, train_interactions] = xlsread(trainingData{i},sheet);
            
            if modelType == 1
                [~,~,~,sigma_delta_scores] = indigo_predict(indigo_model,train_interactions, ...
                  input_type,annotation_file,chemogenomics_file);
            elseif modelType == 2
                [~,~,~,sigma_delta_scores] = indigo_predict_tb(indigo_model,train_interactions, ...
                 input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
            end
            
            trainOrthologs = get_orthologs(trainingData{i},modelType);  
            
            if ~isempty(trainOrthologs)
                [~,sigma_delta_scores] = indigo_orthology(labels, ...
                 trainOrthologs,sigma_delta_scores, indigo_model); 
            end
            
        end
        
        %Training data --> fitrensemble
        if size(interactions_all,2) > size(train_interactions,2)
            num_add = size(interactions_all,2) - size(train_interactions,2);
            empty_array = repmat("",length(train_interactions),num_add);
            train_interactions = [train_interactions empty_array];
        elseif size(interactions_all,2) < size(train_interactions,2)
            num_add = size(train_interactions,2) - size(interactions_all,2);
            empty_array = repmat("",length(interactions_all),num_add);
            interactions_all = [interactions_all empty_array];
        end
        
        interactions_all = [interactions_all; train_interactions];
        interaction_scores_all = [interaction_scores_all; train_scores];
        sigma_delta_scores_all = [sigma_delta_scores_all, sigma_delta_scores];    
    end
end
         
indigoSummary.valMethod = valMethod;

if strcmp(valMethod,'holdout_onself')
    %Only use data from test file
    %Holdout validation
    i = 1;
    [train,test] = crossvalind('HoldOut',length(interactions),0.2); %20% of data is in test set
    Xtrain = interactions(train,:);
    Ytrain = scores(train);
    Xtest = interactions(test,:);
    Ytest = scores(test);
    %plug training data into indigo train
    writecell([Xtrain,num2cell(Ytrain)],'train.xlsx','Sheet',sheet)
    tic
    
    if modelType == 1
        [train_interactions, train_scores, labels, indigo_model,...
         sigma_delta_scores, conditions] = indigo_train('train.xlsx', sheet,...
         annotation_file,chemogenomics_file);
    elseif modelType == 2
        [train_interactions, train_scores, labels, indigo_model,...
         sigma_delta_scores, conditions] = indigo_train_tb('train.xlsx',sheet, ...
         annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
    end
    
    toc
    %Predict and evaluate
    predictStep();   
    
elseif strcmp(valMethod, 'holdout')
    i = 1;
    %Holdout validation - add 80% of test set to training and keep 20% as test
    [train,test] = crossvalind('HoldOut',length(scores),0.2); %20% of data is in test set
    Xtrain = interactions(train,:);
    Ytrain = scores(train);
    Xtest = interactions(test,:);
    Ytest = scores(test);
    
    if modelType == 1
        [~,~,~,sigma_delta_scores] = indigo_predict(indigo_model,Xtrain, ...
         input_type,annotation_file,chemogenomics_file);
    elseif modelType == 2
        [~,~,~,sigma_delta_scores] = indigo_predict_tb(indigo_model,Xtrain, ...
         input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
    end


    if ~isempty(testOrthologs)
        [~,sigma_delta_scores] = indigo_orthology(labels, testOrthologs, ...
                                 sigma_delta_scores, indigo_model); 
    end

    interaction_scores_cv = [interaction_scores_all; Ytrain];
    sigma_delta_scores_cv = [sigma_delta_scores_all, sigma_delta_scores];
   
    %takes a long time!
    tic
    indigo_model = fitrensemble(single(sigma_delta_scores_cv'), ...
                   single(interaction_scores_cv),'Method','Bag');
    toc
    
    %Predict and evaluate
    predictStep();
    
elseif strcmp(valMethod,'independent')
    i = 1;
    %independent validation
    %No need to do this more than once
    %leave out the entire test set and make predictions on it
    
    %Only need to do this if you have more than one file in training data
    if length(trainingData) > 1
        tic
        indigo_model = fitrensemble(single(sigma_delta_scores_all'), ...
                   single(interaction_scores_all),'Method','Bag');
        toc
    end
    
    Xtrain = interactions_all;
    Ytrain = interaction_scores_all;
    Xtest = interactions;   %What you are plugging in to make predictions
    Ytest = scores;         %What you are comparing to at the end

    %Predict and evaluate
    predictStep();
elseif strcmp(valMethod,'cv_onself')
    indigoSummary.K = K;
    for i = 1:K
        fprintf('Run %d\n',i)
        %Only use data from test file
        idx = crossvalind('Kfold',length(scores),K);
        test = (idx == i);
        train = ~test;
        Xtrain = interactions(train,:);
        Ytrain = scores(train);
        Xtest = interactions(test,:);
        Ytest = scores(test);
        writecell([Xtrain,num2cell(Ytrain)],'train.xlsx','Sheet',sheet)

        if ~isempty(testOrthologs)
        [~,sigma_delta_scores] = indigo_orthology(labels, testOrthologs, ...
                                 sigma_delta_scores, indigo_model); 
        end
        
        if modelType == 1
            [train_interactions, train_scores, labels, indigo_model,...
             sigma_delta_scores, conditions] = indigo_train('train.xlsx',sheet, ...
             annotation_file,chemogenomics_file);
        elseif modelType == 2
            [train_interactions, train_scores, labels, indigo_model,...
             sigma_delta_scores, conditions] = indigo_train_tb('train.xlsx',sheet, ...
             annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
        end
        

        %Predict and evaluate
        predictStep();    
    end
elseif strcmp(valMethod,'cv')
    indigoSummary.K = K;
    for i = 1:K
        fprintf('Run %d\n',i)
        %Only method that requires data from test file to be in training data
        %with other sets of data
        %Split up test group into subgroups
        idx = crossvalind('Kfold',length(scores),K);
        test = (idx == i);
        train = ~test;
        Xtrain = interactions(train,:);
        Ytrain = scores(train);
        Xtest = interactions(test,:);
        Ytest = scores(test);
        
        if modelType == 1
            [~,~,~,sigma_delta_scores] = indigo_predict(indigo_model,Xtrain, ...
             input_type,annotation_file,chemogenomics_file);
        elseif modelType == 2
            [~,~,~,sigma_delta_scores] = indigo_predict_tb(indigo_model,Xtrain, ...
             input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
        end
        
        if ~isempty(testOrthologs)
            [~,sigma_delta_scores] = indigo_orthology(labels, testOrthologs, ...
                                     sigma_delta_scores, indigo_model); 
        end

        interaction_scores_cv = [interaction_scores_all; Ytrain];
        sigma_delta_scores_cv = [sigma_delta_scores_all, sigma_delta_scores];
        
        %takes a long time!
        tic
        indigo_model = fitrensemble(single(sigma_delta_scores_cv'), ...
                       single(interaction_scores_cv),'Method','Bag');
        toc
        
        %Predict and evaluate
        predictStep();
    end
end

%nested function for predicting scores and evaluating results
function predictStep()
    %input type changes whether you are predicting all combinations 
    %of list of drugs (1) or given interactions (2)
    
    if modelType == 1
         [~,predicted_scores,~,sigma_delta_input] = indigo_predict(indigo_model, ...
          Xtest,input_type,'identifiers_match.xlsx','ecoli_phenotype_data_cell.xlsx');
    elseif modelType == 2
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
    indigoSummary.trainScores{:,i} = Ytrain;
    indigoSummary.testPairs{i} = Xtest;
    indigoSummary.testScores{:,i} = Ytest;
    indigoSummary.predictedScores{:,i} = predicted_scores;
    fprintf(sprintf('INDIGO predictions subset %d complete!\n',i))
end

if isfile('train.xlsx')
    delete train.xlsx
end

end