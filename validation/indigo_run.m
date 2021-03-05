function indigo_summary = indigo_run(test_data, training_data, data_lookup, ... 
    valmethod, K, standardize, model_type, input_type, scoring)
    arguments
        test_data char
        training_data = ''
        data_lookup = ''
        valmethod char {mustBeMember(valmethod,{'Kfold', 'holdout', ...
                        'holdout_onself', 'Kfold_onself', ...
                        'independent'})} = 'Kfold'
        K {mustBeInteger} = 5
        standardize char {mustBeMember(standardize,{'','z_score'})} = ''
        model_type char {mustBeMember(model_type,{'original_model','ecoli_model','mtb_model'})} = 'ecoli_model';
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
      1. test_data:          List of drugs (1) or drug interactions with 
                            scores (2) 

    OPTIONAL INPUTS:
      2. training_data:      Drug interaction data, can be from multiple
                            species. If model is E. coli, first file must 
                            be E. coli data. If model is M. tb, first file 
                            must be either m. tb or e. coli. Can be single
                            file string or cell array of files.
      3. valmethod:         Validation method, can be the following options:
                                1. holdout_onself: Only use test data
                                2. cv_onself: Only use test data
                                3. independent: Leave all test data out of
                                training
                                4. cv: K-fold cross validation
      4. K:                 Cross-validation parameter 
                            (default K = 5 for Kfold cross-validation)
      5. standardize:       To get zscores or not
      6. model_type:         E. coli model, M. tb model, or original model. 
                            Switches training data and chemogenomics/transcriptomics 
                            data that is used. Original model means using
                            only dataset E. coli dataset from
                            Chandrasekaran et al., 2016.
      7. input_type:        1 for drug list, 2 for drug combinations
      8. scoring            Interaction scoring method. Can be Bliss or
                            Loewe, default is none ('') which means model
                            uses data with both scoring methods.
    OUTPUTS:
      1. indigo_summary:     INDIGO model, predicted scores and input data
    %}

    % Initialize output indigo_summary
    indigo_summary = struct;

    indigo_summary.model_type = model_type;
    if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
        %E. coli chemogenomics data
        chemogenomics_file = 'ecoli_phenotype_data_cell.xlsx';

    elseif strcmp(model_type, 'mtb_model')
        % M. tb chemogenomics data
        chemogenomics_file = [];
        file = 'averaged_data_new.mat';
        S = load(file);
        [normData,col,row] = tb_preprocessing(S.averagedtbdata, ...
                             S.mtb_expression_database_col_ids, ...
                             S.mtb_expression_database_row_ids);
        [phenotype_data, phenotype_labels] = process_transcriptome_tb(normData,row,col);
    end

    indigo_summary.scoring = scoring;
    indigo_summary.test_data = test_data;
    indigo_summary.training_data = training_data;
    indigo_summary.data_lookup = data_lookup;

    % Read in the test data
    % Use readtable to handle missing values as {0x0 char}
    data = readtable(test_data,"ReadVariableNames",false);
    test_scores = data{:,end};
    test_interactions = data{:,1:end-1};

    indigo_summary.standardize = standardize;
    if strcmp(standardize,'z_score')
        test_scores = zscore(test_scores);
    end

    %get orthologs for test file
    test_orthologs = get_orthologs(test_data, model_type, data_lookup);

    %{ 
    COMBINE TRAINING DATA FROM DIFFERENT DATA FILES TO GET OVERALL SIGMA DELTA 
    SCORE MATRIX (X) AND INTERACTION SCORE MATRIX (Y) TO BUILD THE MODEL.
    %}
    if ~isempty(training_data)
        indigo_summary.training_data = training_data;
        train_interactions_all = [];
        train_scores_all = [];
        train_sigma_delta_scores_all = [];

        if ischar(class(training_data)) || isstring(class(training_data))
            %convert to cell, accounts for if trainingData is single file
            training_data = cellstr(training_data);
        end

        for i = 1:length(training_data)
            fprintf(sprintf('Adding %s to training data\n',training_data{i}))
            if i == 1
                %Train first
                if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
                    annotation_file = strcat(erase(training_data{i},'.xlsx'),'_ecoli_match.xlsx');
                    [train_interactions, train_scores, phenotype_labels, indigo_model,...
                     train_sigma_delta_scores, ~] = indigo_train(training_data{i},standardize, ...
                        annotation_file, chemogenomics_file);
                elseif strcmp(model_type, 'mtb_model')
                    annotation_file = strcat(erase(training_data{i},'.xlsx'),'_mtb_match.xlsx');
                    [train_interactions, train_scores, phenotype_labels, indigo_model,...
                     train_sigma_delta_scores, ~] = indigo_train_tb(training_data{i},standardize, ...
                     annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
                end
                train_orthologs = get_orthologs(training_data{i}, model_type, data_lookup);       

                if ~isempty(train_orthologs)
                    [~,train_sigma_delta_scores] = indigo_orthology(phenotype_labels, ...
                     train_orthologs,train_sigma_delta_scores, indigo_model); 
                end

            else
                % Use readtable to handle missing values as {0x0 char}
                data = readtable(training_data{i}, "ReadVariableNames", false);
                train_scores = data{:,end};
                train_interactions = data{:,1:end-1};

                if strcmp(standardize, 'standardize')
                    train_scores = zscore(train_scores);
                end

                if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
                     annotation_file = strcat(erase(training_data{i},'.xlsx'),'_ecoli_match.xlsx');
                    [~,~,~,train_sigma_delta_scores] = indigo_predict(indigo_model,train_interactions, ...
                      input_type,annotation_file,chemogenomics_file);
                elseif strcmp(model_type,'mtb_model')
                    annotation_file = strcat(erase(training_data{i},'.xlsx'),'_mtb_match.xlsx');
                    [~,~,~,train_sigma_delta_scores] = indigo_predict_tb(indigo_model,train_interactions, ...
                     input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
                end

                train_orthologs = get_orthologs(training_data{i},model_type, data_lookup);  

                if ~isempty(train_orthologs)
                    [~,train_sigma_delta_scores] = indigo_orthology(phenotype_labels, ...
                     train_orthologs,train_sigma_delta_scores, indigo_model); 
                end

            end

            % Adjusting # of columns of training data matrix so that all the 
            % different training data can be combined before building the model
            if size(train_interactions_all,2) > size(train_interactions,2)
                num_add = size(train_interactions_all,2) - size(train_interactions,2);
                empty_array = repmat("",length(train_interactions),num_add);
                train_interactions = [train_interactions empty_array];
            elseif size(train_interactions_all,2) < size(train_interactions,2)
                num_add = size(train_interactions,2) - size(train_interactions_all,2);
                empty_array = repmat("",length(train_interactions_all),num_add);
                train_interactions_all = [train_interactions_all empty_array];
            end

            train_sigma_delta_scores_all = [train_sigma_delta_scores_all, train_sigma_delta_scores]; 
            train_scores_all = [train_scores_all; train_scores];        
            train_interactions_all = [train_interactions_all; train_interactions];
        end
    end
    
    predictor_names = [phenotype_labels; phenotype_labels];
    %% BUILD MODEL AND MAKE PREDICTIONS WITH VALIDATION METHOD OF YOUR CHOOSING
    indigo_summary.valmethod = valmethod;
    indigo_summary.K = K;
    if strcmp(valmethod,'holdout_onself')
        i = 1;
        % K = 0.2 is usually default --> 20% of data is in test set
        [train,test] = crossvalind('HoldOut',length(test_interactions), K); 
        Xtrain = test_interactions(train,:);
        Ytrain = test_scores(train);
        Xtest = test_interactions(test,:);
        Ytest = test_scores(test);

        % Plug training data into indigo_train
        writecell([Xtrain,num2cell(Ytrain)],'train.xlsx')

        if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
            annotation_file = strcat(erase(test_data,'.xlsx'),'_ecoli_match.xlsx');
            [~, ~, phenotype_labels, indigo_model,...
             train_sigma_delta_scores, ~] = indigo_train('train.xlsx', standardize, ...
             annotation_file,chemogenomics_file);
        elseif strcmp(model_type, 'mtb_model')
            annotation_file = strcat(erase(test_data,'.xlsx'),'_mtb_match.xlsx');
            [~, ~, phenotype_labels, indigo_model,...
             train_sigma_delta_scores, ~] = indigo_train_tb('train.xlsx',standardize, ...
             annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
        end

        % Adjust sigma delta scores based on presenece of orthologous and 
        % nonorthologous genes
        if ~isempty(test_orthologs)
            [~,train_sigma_delta_scores] = indigo_orthology(phenotype_labels, test_orthologs, ...
                                     train_sigma_delta_scores, indigo_model); 
        end

        % Build the model
        tic
        indigo_model = fitrensemble(single(train_sigma_delta_scores'), ...
                       single(Ytrain),'Method','Bag','PredictorNames',predictor_names);
        toc


        % Make predictions and store results
        predictStep();   

    elseif strcmp(valmethod, 'holdout')
        i = 1;
        % Add 80% of test set to training and keep 20% as test
        [train,test] = crossvalind('HoldOut',length(test_scores),0.2); 
        Xtrain = test_interactions(train,:);
        Ytrain = test_scores(train);
        Xtest = test_interactions(test,:);
        Ytest = test_scores(test);

        if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
            annotation_file = strcat(erase(test_data,'.xlsx'),'_ecoli_match.xlsx');
            [~,~,~,train_sigma_delta_scores] = indigo_predict(indigo_model,Xtrain, ...
             input_type,annotation_file,chemogenomics_file);
        elseif strcmp(model_type, 'mtb_model')
            annotation_file = strcat(erase(test_data,'.xlsx'),'_mtb_match.xlsx');
            [~,~,~,train_sigma_delta_scores] = indigo_predict_tb(indigo_model,Xtrain, ...
             input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
        end


        if ~isempty(test_orthologs)
            [~,train_sigma_delta_scores] = indigo_orthology(phenotype_labels, test_orthologs, ...
                                     train_sigma_delta_scores, indigo_model); 
        end

        train_sigma_delta_scores_all = [train_sigma_delta_scores_all, train_sigma_delta_scores];
        train_scores_all = [train_scores_all; Ytrain];

        tic
        indigo_model = fitrensemble(single(train_sigma_delta_scores_all'), ...
                                    single(train_scores_all), 'Method','Bag', ...
                                    'PredictorNames',predictor_names);
        toc

        predictStep();

    elseif strcmp(valmethod,'independent')
        i = 1;

        % Only need to do this if you have more than one file in training data
        % or if the training data has orthologs

        if length(training_data) > 1 || ~isempty(train_orthologs)
            tic
            indigo_model = fitrensemble(single(train_sigma_delta_scores_all'), ...
                                        single(train_scores_all),'Method','Bag', ...
                                        'PredictorNames',predictor_names);
            toc
        end

        Xtrain = train_interactions_all;
        Ytrain = train_scores_all;
        Xtest = test_interactions;   
        Ytest = test_scores;       

        predictStep();

    elseif strcmp(valmethod,'Kfold_onself')
        for i = 1:K
            fprintf('Run %d\n',i)
            idx = crossvalind('Kfold',length(test_scores),K);
            test = (idx == i);
            train = ~test;
            Xtrain = test_interactions(train,:);
            Ytrain = test_scores(train);
            Xtest = test_interactions(test,:);
            Ytest = test_scores(test);
            writecell([Xtrain,num2cell(Ytrain)],'train.xlsx')

            if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
                annotation_file = strcat(erase(test_data,'.xlsx'),'_ecoli_match.xlsx');
                [~, ~, phenotype_labels, indigo_model,...
                 train_sigma_delta_scores, ~] = indigo_train('train.xlsx',standardize, ...
                 annotation_file,chemogenomics_file);
            elseif strcmp(model_type, 'mtb_model')
                annotation_file = strcat(erase(test_data,'.xlsx'),'_mtb_match.xlsx');
                [~, ~, phenotype_labels, indigo_model,...
                 train_sigma_delta_scores, ~] = indigo_train_tb('train.xlsx',standardize, ...
                 annotation_file,chemogenomics_file,1,phenotype_data,phenotype_labels,col);
            end

            if ~isempty(test_orthologs)
                [~,train_sigma_delta_scores] = indigo_orthology(phenotype_labels, test_orthologs, ...
                                         train_sigma_delta_scores, indigo_model); 
            end

            tic
            indigo_model = fitrensemble(single(train_sigma_delta_scores_all'), ...
                                        single(Ytrain),'Method','Bag', ...
                                        'PredictorNames',predictor_names);
            toc

            predictStep();   
        end

    elseif strcmp(valmethod,'Kfold')
        indigo_summary.K = K;
        for i = 1:K
            fprintf('Run %d\n',i)
            idx = crossvalind('Kfold',length(test_scores),K);
            test = (idx == i);
            train = ~test;
            Xtrain = test_interactions(train,:);
            Ytrain = test_scores(train);
            Xtest = test_interactions(test,:);
            Ytest = test_scores(test);
            if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
                annotation_file = strcat(erase(test_data,'.xlsx'),'_ecoli_match.xlsx');
                [~,~,~,train_sigma_delta_scores] = indigo_predict(indigo_model,Xtrain, ...
                 input_type,annotation_file,chemogenomics_file);
            elseif strcmp(model_type, 'mtb_model')
                annotation_file = strcat(erase(test_data,'.xlsx'),'_mtb_match.xlsx');
                [~,~,~,train_sigma_delta_scores] = indigo_predict_tb(indigo_model,Xtrain, ...
                 input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
            end

            if ~isempty(test_orthologs)
                [~,train_sigma_delta_scores] = indigo_orthology(phenotype_labels, test_orthologs, ...
                                         train_sigma_delta_scores, indigo_model); 
            end

            train_sigma_delta_scores_all = [train_sigma_delta_scores_all, train_sigma_delta_scores];
            train_scores_all = [train_scores_all; Ytrain];

            tic
            indigo_model = fitrensemble(single(train_sigma_delta_scores_all'), ...
                                        single(train_scores_all),'Method','Bag', ...
                                        'PredictorNames',predictor_names);
            toc

            predictStep();
        end
    end 

    % nested function for predicting scores and storing results
    function predictStep()
        
        indigo_summary.model{i} = indigo_model;
        indigo_summary.train_interactions{i} = Xtrain;
        indigo_summary.train_scores{i} = Ytrain;
        indigo_summary.test_interactions{i} = Xtest;
        indigo_summary.test_scores{i} = Ytest;

        if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
            annotation_file = strcat(erase(test_data,'.xlsx'),'_ecoli_match.xlsx');
        [~,predicted_scores,~,test_sigma_delta_scores] = indigo_predict(indigo_model, ...
            Xtest,input_type,annotation_file,chemogenomics_file);
        elseif strcmp(model_type, 'mtb_model')
            annotation_file = strcat(erase(test_data,'.xlsx'),'_mtb_match.xlsx');
            [~,predicted_scores,~,test_sigma_delta_scores] = indigo_predict_tb(indigo_model,Xtest, ...
                input_type,annotation_file,chemogenomics_file, 1, phenotype_data, phenotype_labels, col);
        end

        if ~isempty(test_orthologs)
            deviations = indigo_orthology(phenotype_labels, test_orthologs, ... 
                         test_sigma_delta_scores, indigo_model); 
            
            % first column is for first deviations, second column is second
            predicted_scores_1 = predicted_scores - deviations(:,1);
            predicted_scores_2 = predicted_scores - deviations(:,2);
            indigo_summary.predicted_scores_1{i} = predicted_scores_1;
            indigo_summary.predicted_scores_2{i} = predicted_scores_2;
        else
            indigo_summary.predicted_scores{i} = predicted_scores;
        end
        
        fprintf(sprintf('INDIGO predictions subset %d complete!\n',i))
    end

    if isfile('train.xlsx')
        delete train.xlsx
    end

end