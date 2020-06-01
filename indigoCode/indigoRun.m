function indigoSummary = indigoRun(testFile,valMethod,K,orthology,trainingData,trainingOrthology, standardize,input_type)
arguments
    testFile char
    valMethod char {mustBeMember(valMethod,{'holdout_onself', 'cv_onself', 'independent', 'cv'})} = 'holdout_onself'
    K {mustBeInteger} = 5
    orthology char = ''
    trainingData char = ''
    trainingOrthology char = ''
    standardize char {mustBeMember(standardize,{'','standardized'})}= ''
    input_type {mustBeInteger} = 2;
end

%[indigoSummary] = indigo_run(testFile,valMethod,K,orthology,trainingData,standardize)
% This function runs indigo and returns a summary of the results
%required argument:
%testFile
%optional arguments:
%valMethod - name-value argument, options are holdout_onself, cv_onself, independent, cv
%K - fold parameter for K-fold cross validation, default = 5
%trainingData - options are 'indigo' or 'nature'
%orthology - whether test data is ecoli or not, if not, specify orthology file
%standardize - whether data will be standardized or not

%scores = drug interaction scores
%interactions = drug pairs
%sigma_delta_scores = sigma delta scores

%Xtrain, Xtest, Ytrain, Ytest - independent/dependent variables during
%holdout or cross validation 


%Initialize output. Wherever indigoSummary appears, data is being stored
%for output
indigoSummary = struct;
indigoSummary.testFile = testFile;
% 
% %Only applies if additional training data is present i.e. indigo or nature
% %Determine which training data to use and train Indigo model
% %Store both experimental interaction scores and sigma delta scores here

sheet = 1;  %for non-standardized data
indigoSummary.standardized = 0;
if strcmp(standardize,'standardized')
    sheet = 2;  %for z score data
    indigoSummary.standardized = 1;
end

%store data from testFile
[scores,interactions] = xlsread(testFile,sheet);

%setting up training data
%can specify a file, indigo or nature
if ~isempty(trainingData)
    indigoSummary.trainingData = trainingData;
    sigma_delta_scores_all = [];
    trainFiles = dataFiles();
    if strcmp(trainingData, 'indigo')
        
        [train_interactions, train_scores, labels, indigo_model,...
         sigma_delta_scores, conditions] = indigo_train(trainFiles.indigoData, ...
         sheet,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
         
         interaction_scores_all = train_scores;
         sigma_delta_scores_all = sigma_delta_scores;
    
    elseif strcmp(trainingData,'nature')
        
        [train_interactions, train_scores, labels, indigo_model,...
         sigma_delta_scores, conditions] = indigo_train(trainFiles.natureData{1},sheet, ...
         'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
         
         interaction_scores_all = train_scores;
         sigma_delta_scores_all = sigma_delta_scores;
         for i = 2:length(trainFiles.natureData)
             [train_scores, train_interactions] = xlsread(trainFiles.natureData{i},sheet);

             [~,~,~,sigma_delta_input] = indigo_predict(indigo_model,train_interactions, ...
              2,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');

             [~,sigma_delta_scores] = indigo_orthology(labels, ...
             trainFiles.natureOrthologs{i-1},sigma_delta_input, indigo_model); 

             interaction_scores_all = [interaction_scores_all; train_scores];
             sigma_delta_scores_all = [sigma_delta_scores_all, sigma_delta_scores];
             
         end
         
    elseif strcmp(trainingData,'nature_ecoli')
        for i = 1:length(trainFiles.natureEcoli)
            [train_interactions, train_scores, labels, indigo_model,...
             sigma_delta_scores, conditions] = indigo_train(trainFiles.natureData{i},sheet, ...
             'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
         
            if strcmp(trainFiles.natureEcoli{i},'ecoli_iAi1.xlsx')
                %use orthology
                [~,~,~,sigma_delta_input] = indigo_predict(indigo_model,train_interactions, ...
                2,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');

             [~,sigma_delta_scores] = indigo_orthology(labels, ...
             trainFiles.natureOrthologs{1},sigma_delta_input, indigo_model); 
            end
            interaction_scores_all = train_scores;
            sigma_delta_scores_all = sigma_delta_scores;
        end
    else
        %if not indigo or entire set of nature data and only 1 file
        [train_interactions, train_scores, labels, indigo_model,...
         sigma_delta_scores, conditions] = indigo_train(trainingData,sheet, ...
         'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
         
         if ~isempty(trainingOrthology)
             [~,training_orth] = xlsread(trainingOrthology);  
             indigoSummary.trainingOrthology = trainingOrthology;
             [~,~,~,sigma_delta_input] = indigo_predict(indigo_model,train_interactions, ...
             2,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
             [~,sigma_delta_scores] = indigo_orthology(labels, ...
             training_orth,sigma_delta_input, indigo_model);
         end
         interaction_scores_all = train_scores;
         sigma_delta_scores_all = sigma_delta_scores;
         
    end  
end

%Check to see if orthology applies (data is not ecoli)
if ~isempty(orthology)
    %data you are testing is not ecoli, will require using orthology
    [~,ecoli_orth] = xlsread(orthology);  
    indigoSummary.orthology = orthology;
end
 
indigoSummary.valMethod = valMethod;
indigoSummary.K = K;
for i = 1:K
    fprintf('Run %d\n',i)
    if strcmp(valMethod,'holdout_onself')
        %Only use data from test file
        %Holdout validation
        [train,test] = crossvalind('HoldOut',length(interactions),0.2); %20% of data is in test set
        Xtrain = interactions(train,:);
        Ytrain = scores(train);
        Xtest = interactions(test,:);
        Ytest = scores(test);
        %plug training data into indigo train
        writecell([Xtrain,num2cell(Ytrain)],'train.xlsx','Sheet',sheet)
        
        [train_interactions, train_scores, labels, indigo_model,...
         sigma_delta_scores, conditions] = indigo_train('train.xlsx', sheet,...
         'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');

        %Predict and evaluate
        predictStep();
    elseif strcmp(valMethod,'cv_onself') 
        %Only use data from test file
        idx = crossvalind('Kfold',length(scores),K);
        test = (idx == i);
        train = ~test;
        Xtrain = interactions(train,:);
        Ytrain = scores(train);
        Xtest = interactions(test,:);
        Ytest = scores(test);
        writecell([Xtrain,num2cell(Ytrain)],'train.xlsx')
        
        [train_interactions, train_scores, labels, indigo_model,...
         sigma_delta_scores, conditions] = indigo_train('train.xlsx',sheet, ...
         'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
     
        %Predict and evaluate
        predictStep();    

    elseif strcmp(valMethod,'cv')
        %Only method that requires data from test file to be in training data
        %with other sets of data
        %Split up test group into subgroups
        idx = crossvalind('Kfold',length(scores),K);
        %Train on either indigo, nature or both with subgroup
        test = (idx == i);
        train = ~test;
        Xtrain = interactions(train,:);
        Ytrain = scores(train);
        Xtest = interactions(test,:);
        Ytest = scores(test);

        [~,~,~,sigma_delta_input] = indigo_predict(indigo_model,Xtrain, ...
             2,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');

        if ~isempty(orthology)
            [~,sigma_delta_scores] = indigo_orthology(labels, ecoli_orth, ...
                                     sigma_delta_input, indigo_model); 
        else
            sigma_delta_scores = sigma_delta_input;     %testFile contains ecoli data
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

if strcmp(valMethod,'independent')
    i = 1;
    %independent validation
    %No need to do this more than once
    %leave out the entire test set and make predictions on it
    %train on either indigo, nature or both
    %need to use fitresnsemble
    tic
    indigo_model = fitrensemble(single(sigma_delta_scores_all'), ...
                   single(interaction_scores_all),'Method','Bag');
    toc

    Xtest = interactions;   %What you are plugging in to make predictions
    Ytest = scores;         %What you are comparing to at the end

    %Predict and evaluate
    predictStep();
end

%nested function for predicting scores and evaluating results
%done to shorten lines of code
function predictStep()
    %input type changes whether you are predicting interactions (2) or
    %combinations of specific drugs (1)
    if ~isempty(orthology)
        [~,predicted_ecoli_scores,~,sigma_delta_input] = indigo_predict(indigo_model, ...
         Xtest,input_type,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');

        deviations = indigo_orthology(labels, ecoli_orth, ... 
        sigma_delta_input, indigo_model); 

        predicted_scores = predicted_ecoli_scores - deviations;
    else
        [~, predicted_scores]  = indigo_predict(indigo_model, ...
         Xtest,input_type,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
    end

    indigoSummary.model{i} = indigo_model;
    indigoSummary.trainPairs{i} = Xtrain;
    indigoSummary.trainScores{:,i} = Ytrain;
    indigoSummary.testPairs{i} = Xtest;
    indigoSummary.testScores{:,i} = Ytest;
    indigoSummary.predictedScores{:,i} = predicted_scores;
end

end