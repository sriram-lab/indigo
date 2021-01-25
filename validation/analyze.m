function [stats,averages,overview] = analyze(indigoSummary)

    %{
    DESCRIPTION

    This function analyzes INDIGO prediction results, calculating 
    statistical measures such as Spearman's correlation coefficient,
    accuracy, precision and recall.

    STEPS
    1. Input processing: Read in experimental and predicted scores.
    2. Statistical analysis
    3. Produce figures such as scatterplots, confusion matrices and ROC
    curves.
    4. Save analysis to excel files and .mat files.

    Author: David Chang
    Created: January 23, 2021

    I/O
    
    REQUIRED INPUTS:
      1. indigoSummary:     INDIGO model, predicted scores, validation method 
    OUTPUTS:
      1. stats:             Statistical measures for each subset of 
                            validation results. Only present when K > 1.           
      2. averages:          Mean values of measures contained in stats.
      3. overview:          Statistical measures when combining validation results
                            from all subsets.
    %}

    dataName = erase(indigoSummary.testData,'.xlsx');
    fprintf(sprintf('Results for %s',dataName))
    stats = struct;
    averages = struct;
    
    if strcmp(indigoSummary.valMethod,'Kfold_onself') || strcmp(indigoSummary.valMethod,'Kfold')
        Ytest_total = [];
        Ypred_total = [];
        drugs_total = [];
        for i = 1:indigoSummary.K
            drugs_total = [drugs_total; indigoSummary.testPairs{i}];
            Ytest_total = [Ytest_total; indigoSummary.testScores{i}];
            Ypred_total = [Ypred_total; indigoSummary.predictedScores{i}];
            Ytrain = indigoSummary.trainScores{i};
            Ytest = indigoSummary.testScores{i};
            Ypred = indigoSummary.predictedScores{i};

            %call function
            getStats(Ytrain,Ytest,Ypred,1);   % 1 means analysis per subset
        end
        Ytrain_total = Ytest_total;
    else
        drugs_total = [indigoSummary.testPairs{1}];
        Ytest_total = [indigoSummary.testScores{1}];
        Ypred_total = [indigoSummary.predictedScores{1}];
        Ytrain_total = [indigoSummary.trainScores{1}];
    end

    %overall - a new struct

    overview.drugInteractions = drugs_total;
    overview.experimentalScores = Ytest_total;
    overview.predictedScores = Ypred_total;
    getStats(Ytrain_total,Ytest_total,Ypred_total,2);   % 2 means overall analysis

    %Tables
    varnames = {'Interactions', ...
                    'R (rank)', ...
                    'P value', ...
                    'Accuracy', ...
                    'Absolute error', ...
                    'Precision (synergy)', ...
                    'Recall (synergy)', ...
                    'Precision (antagonism)', ...
                    'Recall (antagonism)', ...
                    'AUC - ROC (synergy)', ...
                    'AUC - ROC (antagonism)'};
                
    if strcmp(indigoSummary.valMethod,'Kfold_onself') || strcmp(indigoSummary.valMethod,'Kfold')
        resultsTable = struct2table(stats);
        resultsTable.Properties.VariableNames = varnames;
        rownames = cell(1,indigoSummary.K);
        for i = 1:indigoSummary.K
            rownames{i} = sprintf('Subset %d', i);
        end
        resultsTable.Properties.RowNames = rownames;

        fields = fieldnames(stats);
        tableArray = zeros(length(fields),1);
        for i = 1:length(fields)
            averages.(fields{i}) = mean(stats.(fields{i}));
            tableArray(i) = mean(stats.(fields{i}));
        end
        averagesTable = table(tableArray);
        averagesTable.Properties.RowNames = varnames; 
        averagesTable.Properties.VariableNames = {'Value'};
    end

    fields = fieldnames(overview);
    tableArray = zeros(length(fields)-3,1);  %don't include first 3 fields
    for i = 4:length(fields)
        tableArray(i-3) = overview.(fields{i});
    end
    overviewTable = table(tableArray);
    overviewTable.Properties.RowNames = varnames;
    overviewTable.Properties.VariableNames = {'Value'};
    
    %% Display analysis and save to .xlsx files
    
    filename = strcat(erase(indigoSummary.testData,'.xlsx'),'_',indigoSummary.valMethod);
    if strcmp(indigoSummary.standardize,'z_score')
        filename = sprintf('%s_z', filename);
    end
    if strcmp(indigoSummary.scoring,'bliss') || strcmp(indigoSummary.scoring,'loewe')
        resultsFile = strcat(indigoSummary.modelType, '/', indigoSummary.scoring, ...
            '/', indigoSummary.valMethod, '/', filename, '.xlsx');
    else
        resultsFile = strcat(indigoSummary.modelType, '/', ...
            indigoSummary.valMethod, '/', filename, '.xlsx');
    end

    %Final output!

    if strcmp(indigoSummary.valMethod,'cv_onself') || strcmp(indigoSummary.valMethod,'cv')
        resultsTable
        averagesTable
        overviewTable
        writetable(resultsTable,resultsFile,'Sheet','results','WriteRowNames',true)
        writetable(averagesTable,resultsFile,'Sheet','averages','WriteRowNames',true)
        writetable(overviewTable,resultsFile,'Sheet','overview','WriteRowNames',true)
    else
        overviewTable
        writetable(overviewTable,resultsFile,'Sheet','overview','WriteRowNames',true)
    end

    %% Conduct statistical analysis
    function getStats(Ytrain, Ytest, Ypred, mode)
        interactionCount = length(Ytest);
        %correlation 
        [R,P] = corr(Ytest, Ypred,'type','Spearman');

        [synergyCutoff, antagonismCutoff] = cutoffs(indigoSummary.testData);
        %synergy - assign value = -1, but remember, synergy is good!

        Ytest(Ytest <= synergyCutoff) = -1;
        Ypred(Ypred <= synergyCutoff) = -1;
        %neutral - assign value = 0
        Ytest(Ytest > synergyCutoff & Ytest < antagonismCutoff) = 0;
        Ypred(Ypred > synergyCutoff & Ypred < antagonismCutoff) = 0;
        %antagonism - assign value = 1, but remember, antagonism is bad!
        Ytest(Ytest >= antagonismCutoff) = 1;
        Ypred(Ypred >= antagonismCutoff) = 1;
        %Get accuracy, precision and recall
        accuracy = sum(Ypred == Ytest)/length(Ypred);
        absError = mean(abs(Ypred-Ytest));
        precisionSynergy = sum(Ypred == -1 & Ytest == -1)/sum(Ypred == -1);
        recallSynergy = sum(Ypred == -1 & Ytest == -1)/sum(Ytest == -1);
        precisionAntagonism = sum(Ypred == 1 & Ytest == 1)/sum(Ypred == 1);
        recallAntagonism = sum(Ypred == 1 & Ytest == 1)/sum(Ytest == 1);

        %Compare model accuracy to 100 random guesses
    %     accuracyGuess = zeros(100,1);
    %     for k = 1:100
    %         random_pred = Ytrain(randperm(length(Ytrain), length(Ytest))); %shuffle values around
    %         accuracyGuess(k) = sum(random_pred == Ytest)/length(Ytest);        
    %     end
    %     
    %     meanGuessAccuracy = mean(accuracyGuess);
    %     [~,p] = ttest2(accuracy,accuracyGuess);
    %     
    
        %% FIGURES
        % ROC curves
        if mode == 1
            figure(1)
            sgtitle(figure(1),sprintf("ROC Curves for %s using INDIGO", ...
            dataName),'Interpreter','none')
            subplot(ceil(indigoSummary.K/2),floor(indigoSummary.K/2),i)    %specific to subset
    %         subplot(2,1,i) 
        else
            sgtitle(figure(4),sprintf("Overall results for %s using INDIGO", ...
            dataName),'Interpreter','none')
            figure(4)
            subplot(3,1,1)
        end
  
        labels = {};  
        if sum(Ytest == -1) > 0  && (sum(Ytest == 0) > 0 || sum(Ytest == 1) > 0)
            %ROC Curve for synergy
            [X_synergy,Y_synergy,T_synergy,AUC_synergy] = perfcurve(Ytest,Ypred,-1);
            plot(X_synergy,Y_synergy)
            labels{end+1} = 'synergy';
        else
            AUC_synergy = nan;
        end
        hold on
        if sum(Ytest == 1) > 0 && (sum(Ytest == 0) > 0 || sum(Ytest == -1) > 0)
            %ROC Curve for antagonism
            [X_antagonism,Y_antagonism,T_antagonism,AUC_antagonism] = perfcurve(Ytest,Ypred,1);
            plot(X_antagonism,Y_antagonism)
            labels{end+1} = 'antagonism';
        else
            AUC_antagonism = nan;
        end
        hold off
        
        legend(labels,'Location','southeast')
        xlabel('False positive rate') 
        ylabel('True positive rate')
       
        if mode == 1
            title(sprintf('Subset %d',i))
        else
            title('Overall ROC curve')
        end
        
        % Confusion matrices
        if mode == 1
            figure(2)
            sgtitle(figure(2),sprintf("Confusion Matrices for %s using INDIGO", ...
            dataName),'Interpreter','none')
            subplot(ceil(indigoSummary.K/2),floor(indigoSummary.K/2),i)     %specific to subset
    %         subplot(2,1,i)
        else
            figure(4)
            subplot(3,1,2)
        end
        C = confusionmat(Ytest,Ypred);
        tempLabels = {};
        if sum(Ytest == -1) > 0 || sum(Ypred == -1) > 0
            tempLabels{end+1} = 'synergistic';
        end
        if sum(Ytest == 0) > 0 || sum(Ypred == 0) > 0
            tempLabels{end+1} = 'neutral';
        end
        if sum(Ytest == 1) > 0 || sum(Ypred == 1) > 0
            tempLabels{end+1}  = 'antagonistic';
        end
        cats = categorical(tempLabels);
        confusionLabels = reordercats(cats,tempLabels);
        cm = confusionchart(C,confusionLabels);
        cm.RowSummary = 'total-normalized';
        cm.ColumnSummary = 'total-normalized';
        if mode == 1
            cm.Title = sprintf('Subset %d',i);   %specific to subset
        else
            cm.Title = 'Overall Confusion Matrix';
        end   
        
        % Scatterplot with line of best fit
        % Rank sorted - larger rank means bigger value
        [Ytest_sorted,I_test] = sort(Ytest);
        Ytest_rank = zeros(1,length(Ytest));
        Ytest_rank(I_test) = 1:length(Ytest);
        [Ypred_sorted,I_pred] = sort(Ypred);
        Ypred_rank = zeros(1,length(Ytest));
        Ypred_rank(I_pred) = 1:length(Ypred);
        
        if mode == 1
            figure(3)
            sgtitle(figure(3),sprintf("Predicted vs experimental scores (rank sorted) for %s using INDIGO", ...
            dataName),'Interpreter','none')
            subplot(ceil(indigoSummary.K/2),floor(indigoSummary.K/2),i) %specific to subset
    %         subplot(2,1,i)
        else
            figure(4)
            subplot(3,1,3)
        end
        scatter(Ytest_rank,Ypred_rank,24,'filled')
        h = lsline; h.Color = 'red'; 
        h.DisplayName = sprintf('R = %.2f', R);
        if mode == 1
            title(sprintf('Subset %d',i))
        else
            title('Overall scatter plot');
        end
        xlabel('Experiment - Rank sorted')
        ylabel('Prediction - Rank sorted')
        legend(h,'Location','northwest')

        if mode == 1
            stats.interactionCount(i,1) = interactionCount;
            stats.rankCorr(i,1) = R;
            stats.pValue(i,1) = P;
            stats.accuracy(i,1) = accuracy;
            stats.absError(i,1) = absError;
            stats.precisionSynergy(i,1) = precisionSynergy;
            stats.recallSynergy(i,1) = recallSynergy;
            stats.precisionAntagonism(i,1) = precisionAntagonism;
            stats.recallAntagonism(i,1) = recallAntagonism;
    %         stats.meanGuessAccuracy(i,1) = meanGuessAccuracy;
    %         stats.compareGuess(i,1) = p;
            stats.AUC_synergy(i,1) = AUC_synergy;
            stats.AUC_antagonism(i,1) = AUC_antagonism;
        else
            overview.interactionCount = interactionCount;
            overview.rankCorr = R;
            overview.pValue = P;
            overview.accuracy = accuracy;
            overview.absError = absError;
            overview.precisionSynergy = precisionSynergy;
            overview.recallSynergy = recallSynergy;
            overview.precisionAntagonism = precisionAntagonism;
            overview.recallAntagonism = recallAntagonism;
    %         overview.meanGuessAccuracy = meanGuessAccuracy;
    %         overview.compareGuess = p;
            overview.AUC_synergy = AUC_synergy;
            overview.AUC_antagonism = AUC_antagonism;
        end
    end
end














%WORK ON THIS NOW!



   %top predictors of synergy
    %rank from lowest to highest
    %take scores that are below synergy cutoff
%     synergyTest = Ytest_sorted(Ytest_sorted < synergyCutoff);
%     synergyPred = Ypred_sorted(Ypred_sorted < synergyCutoff);
%     %get top scores and drug pairs
%     if length(synergyPred) >= 10
%         
%     else 
%         
%     end
%         
%     
%     topSynergyTestScores = synergyTest(1:10);
%     topSynergyPredScores = synergyPred(1:10);
%     topSynergyTestPairs = indigoSummary.testPairs{i}(I_test(1:10),:);
%     topSynergyPredPairs = indigoSummary.testPairs{i}(I_pred(1:10),:);
    
    
    
    %check whether they are correctly classified
    %top predictors of antagonism
    %take scores that are above antagonism cutoff
    %get first 10 scores and drug pairs
    