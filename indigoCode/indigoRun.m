function indigoSummary = indigoRun(testFile,trainingData,valMethod,K,standardize,input_type)
arguments
    testFile char
    trainingData = ''
    valMethod char {mustBeMember(valMethod,{'holdout_onself', 'cv_onself', 'independent', 'cv'})} = 'holdout_onself'
    K {mustBeInteger} = 5
    standardize char {mustBeMember(standardize,{'','standardized'})}= ''
    input_type {mustBeInteger} = 2;
end

%[indigoSummary] = indigoRun(testFile,trainingData,valMethod,K, standardize,input_type)
% This function runs indigo and returns a summary of the results
%required argument:
%testFile
%optional arguments:
%trainingData - single file or cell array of files
%valMethod - name-value argument, options are holdout_onself, cv_onself, independent, cv
%K - fold parameter for K-fold cross validation, default = 5
%standardize - whether data will be standardized or not
%input type

%scores = drug interaction scores
%interactions = drug pairs
%sigma_delta_scores = sigma delta scores

%Xtrain, Xtest, Ytrain, Ytest - independent/dependent variables during
%holdout or cross validation 


%Initialize output. Wherever indigoSummary appears, data is being stored
%for output
indigoSummary = struct;
indigoSummary.testFile = testFile;

% %Only applies if additional training data is present i.e. indigo or nature
% %Determine which training data to use and train Indigo model
% %Store both experimental interaction scores and sigma delta scores here

sheet = 1;  %for non-standardized data
indigoSummary.standardized = 0;
if strcmp(standardize,'standardized')
    sheet = 2;  %for z score data
    indigoSummary.standardized = 1;
end


files = cellstr(ls('indigoData'));
%store data from testFile
[scores,interactions] = xlsread(testFile,sheet);
orthology = strcat(erase(testFile,'.xlsx'),'_orthologs.xlsx');
if sum(contains(files,orthology)) ~= 0
    %orthology
    [~,test_orth] = xlsread(orthology);
    indigoSummary.orthology = orthology;
end
 
%setting up training data
%can specify a file, indigo or nature
if ~isempty(trainingData)
    indigoSummary.trainingData = trainingData;
    interaction_scores_all = [];
    sigma_delta_scores_all = [];
    
    if ischar(class(trainingData)) || isstring(class(trainingData))
        %convert to cell
        trainingData = cellstr(trainingData);
    end
    for i = 1:length(trainingData)
        orthology = strcat(erase(trainingData{i},'.xlsx'),'_orthologs.xlsx');
        if i == 1
            %Train first
            [train_interactions, train_scores, labels, indigo_model,...
             sigma_delta_scores, conditions] = indigo_train(trainingData{i},sheet, ...
             'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
         
            if sum(contains(files,orthology)) ~= 0
                [~,train_orth] = xlsread(orthology);
                %orthology
                [~,sigma_delta_scores] = indigo_orthology(labels, ...
                 train_orth,sigma_delta_scores, indigo_model); 
            end             
        else
            [train_scores, train_interactions] = xlsread(trainingData{i},sheet);

            [~,~,~,sigma_delta_scores] = indigo_predict(indigo_model,train_interactions, ...
              2,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
             
             if sum(contains(files,orthology)) ~= 0
                [~,train_orth] = xlsread(orthology);
                %orthology
                [~,sigma_delta_scores] = indigo_orthology(labels, ...
                 train_orth,sigma_delta_scores, indigo_model); 
             end
        end
        %Training data --> fitrensemble
        interaction_scores_all = [interaction_scores_all; train_scores];
        sigma_delta_scores_all = [sigma_delta_scores_all, sigma_delta_scores];    
    end
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
        [~,~,~,sigma_delta_scores] = indigo_predict(indigo_model,Xtrain, ...
             2,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');

        if ~isempty(test_orth)
            [~,sigma_delta_scores] = indigo_orthology(labels, test_orth, ...
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
    [~,predicted_scores,~,sigma_delta_input] = indigo_predict(indigo_model, ...
     Xtest,input_type,'identifiers_match_tb.xlsx','ecoli_phenotype_data_cell.xlsx');
    if ~isempty(test_orth)
        deviations = indigo_orthology(labels, test_orth, ... 
        sigma_delta_input, indigo_model); 

        predicted_scores = predicted_scores - deviations;
    end

    indigoSummary.model{i} = indigo_model;
    indigoSummary.trainPairs{i} = Xtrain;
    indigoSummary.trainScores{:,i} = Ytrain;
    indigoSummary.testPairs{i} = Xtest;
    indigoSummary.testScores{:,i} = Ytest;
    indigoSummary.predictedScores{:,i} = predicted_scores;
end

end