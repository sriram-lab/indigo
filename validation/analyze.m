function [stats,averages,overview] = analyze(indigoSummary,resultIndex, dataFiles)

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

    %% Structure naming for output files and figures
    dataName = erase(indigoSummary.testData,'.xlsx');
    
    filename = strcat(erase(indigoSummary.testData,'.xlsx'),'_',indigoSummary.valMethod);
    sheetName = indigoSummary.valMethod;
    if strcmp(indigoSummary.standardize,'z_score')
        filename = sprintf('%s_z', filename);
        sheetName = strcat(indigoSummary.valMethod,' (z)');
    end
     
    if strcmp(indigoSummary.scoring,'bliss') || strcmp(indigoSummary.scoring,'loewe')
        overviewFile = strcat('results/v2/', indigoSummary.modelType, '/', ...
            indigoSummary.scoring, '/','overview_results.xlsx');
        resultsFile = strcat('results/v2/', indigoSummary.modelType, '/', ...
            indigoSummary.scoring, '/', indigoSummary.valMethod, '/', filename);
    else
        overviewFile = strcat('results/v2/',indigoSummary.modelType, '/', ...
            'overview_results.xlsx');
        resultsFile = strcat('results/v2/',indigoSummary.modelType, '/', ...
            indigoSummary.valMethod, '/', filename);
    end
    
    %% Process results
    fprintf(sprintf('Results for %s',dataName))
    stats = struct;
    averages = struct;
    
    %% SUBSET DATA FOR KFOLD CV
    if strcmp(indigoSummary.valMethod,'Kfold_onself') || strcmp(indigoSummary.valMethod,'Kfold')
        drugs_total = [];
        Ytest_total = [];
        Ypred_total = [];
        
        f1 = figure(1);
        f1.Visible = 'off';
        t1 = tiledlayout(f1,3,2);
        title(t1,sprintf('ROC Curves for %s\n', dataName), ...
              'FontWeight','bold', 'Interpreter', 'None')
        
        for i = 1:indigoSummary.K
            drugs_total = [drugs_total; indigoSummary.testPairs{i}];
            Ytest_total = [Ytest_total; indigoSummary.testScores{i}];
            Ypred_total = [Ypred_total; indigoSummary.predictedScores{i}];
            Ytest = indigoSummary.testScores{i};
            Ypred = indigoSummary.predictedScores{i};
            
            %correlation 
            [R,P] = corr(Ytest, Ypred,'type','Spearman');
            
            [Ytest, Ypred] = classifyScores(Ytest,Ypred);   
            
            % 1 means analysis per subset
            getStats(Ytest, Ypred, 1)
            
            nexttile
            getROC(Ytest, Ypred, 1)
        end
        
        figFile = strcat(resultsFile,'_roc');
        saveas(f1,figFile,'fig') 
        saveas(f1,figFile,'png')
        close(f1)
        
        f2 = figure(2);
        f2.Visible = 'off';
        t2 = tiledlayout(f2,3,2);
        title(t2, sprintf('Confusion Matrices for %s\n', dataName), ...
              'FontWeight','bold','Interpreter','None')
        
        for i = 1:indigoSummary.K
            nexttile
            getConfusion(Ytest, Ypred, 1)
        end
        
        figFile = strcat(resultsFile,'_cm');
        saveas(f2,figFile,'fig')
        saveas(f2,figFile,'png')
        close(f2)
        
        f3 = figure(3);
        f3.Visible = 'off';
        t3 = tiledlayout(f3,3,2);
        title(t3, sprintf('Scatter Plots (Ranked Sorted) for %s\n', dataName), ...
              'FontWeight','bold','Interpreter','None')
        
        for i = 1:indigoSummary.K     
            nexttile
            getScatter(Ytest, Ypred, 1)  
        end
        
        figFile = strcat(resultsFile,'_sc');
        saveas(f3,figFile,'fig')
        saveas(f3,figFile,'png')
        close(f3)
        
    else
        drugs_total = [indigoSummary.testPairs{1}];
        Ytest_total = [indigoSummary.testScores{1}];
        Ypred_total = [indigoSummary.predictedScores{1}];
    end

    %% OVERALL DATA - ALL RESULTS COMBINED
    
    overview.drugInteractions = drugs_total;
    overview.experimentalScores = Ytest_total;
    overview.predictedScores = Ypred_total;
    
    [R,P] = corr(Ytest_total, Ypred_total,'type','Spearman');

    [Ytest_total, Ypred_total] = classifyScores(Ytest_total, Ypred_total);   

    % 2 means analysis for all results combined
    getStats(Ytest_total, Ypred_total, 2)
    f4 = figure(4);
    f4.Visible = 'off';
    t4 = tiledlayout(f4,3,1);
    title(t4,sprintf('Overall INDIGO Results for %s\n',dataName), ...
        'FontWeight','bold','Interpreter','none')
    nexttile
    getROC(Ytest_total, Ypred_total, 2)    
    nexttile
    getConfusion(Ytest_total, Ypred_total, 2)
    nexttile
    getScatter(Ytest_total, Ypred_total, 2)  
    figFile = strcat(resultsFile,'_overall');
    saveas(f4,figFile,'fig')
    saveas(f4,figFile,'png')
    close(f4)

    %% FORMAT TABLES
    %Tables
    varnames = {'Interaction Count', ...
                    'R (rank)', ...
                    'P value', ...
                    'Accuracy', ...
                    'Absolute Error', ...
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

    if strcmp(indigoSummary.valMethod,'Kfold_onself') || strcmp(indigoSummary.valMethod,'Kfold')
        resultsTable
        averagesTable
        overviewTable
        writetable(resultsTable,strcat(resultsFile,'.xlsx'),'Sheet','results','WriteRowNames',true)
        writetable(averagesTable,strcat(resultsFile,'.xlsx'),'Sheet','averages','WriteRowNames',true)
        writetable(overviewTable,strcat(resultsFile,'.xlsx'),'Sheet','overview','WriteRowNames',true)  
    else
        overviewTable
        writetable(overviewTable,strcat(resultsFile,'.xlsx'),'Sheet','overview','WriteRowNames',true)
    end
    
    overviewTable = rows2vars(overviewTable,'VariableNamingRule','preserve');
    overviewTable = overviewTable(:,2:end);    
    
    writecell(dataFiles, overviewFile, 'Sheet', sheetName, 'Range', ...
        sprintf('A2:A%d', length(dataFiles)+1))
    writecell([{'Datasets'}, varnames], overviewFile, 'Sheet', sheetName, 'Range', 'A1:L1')
    writetable(overviewTable, overviewFile, 'Sheet', sheetName, 'Range', ...
        sprintf('B%d:L%d',resultIndex+1,resultIndex+1), 'WriteVariableNames',false)

    %% CLASSIFY SCORES
    function [Ytest, Ypred] = classifyScores(Ytest, Ypred)
        [synergyCutoff, antagonismCutoff] = cutoffs(indigoSummary.testData, indigoSummary.dataLookup);
        %synergy - assign value = -1, but remember, synergy is good!
        Ytest(Ytest <= synergyCutoff) = -1;
        Ypred(Ypred <= synergyCutoff) = -1;
        %neutral - assign value = 0
        Ytest(Ytest > synergyCutoff & Ytest < antagonismCutoff) = 0;
        Ypred(Ypred > synergyCutoff & Ypred < antagonismCutoff) = 0;
        %antagonism - assign value = 1, but remember, antagonism is bad!
        Ytest(Ytest >= antagonismCutoff) = 1;
        Ypred(Ypred >= antagonismCutoff) = 1;
    end

    %% GET STATISTICS
    function getStats(Ytest, Ypred, mode)
        interactionCount = length(Ytest);
        %Get accuracy, precision and recall
        accuracy = sum(Ypred == Ytest)/length(Ypred);
        absError = mean(abs(Ypred-Ytest));
        precisionSynergy = sum(Ypred == -1 & Ytest == -1)/sum(Ypred == -1);
        recallSynergy = sum(Ypred == -1 & Ytest == -1)/sum(Ytest == -1);
        precisionAntagonism = sum(Ypred == 1 & Ytest == 1)/sum(Ypred == 1);
        recallAntagonism = sum(Ypred == 1 & Ytest == 1)/sum(Ytest == 1);
        
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
        end
    end
    
    %% ROC CURVES
    function getROC(Ytest,Ypred,mode)
        labels = {};  
        if sum(Ytest == -1) > 0
            %ROC Curve for synergy
            [X_synergy,Y_synergy,T_synergy,AUC_synergy] = perfcurve(Ytest,Ypred,-1);
            plot(X_synergy,Y_synergy)
            labels{end+1} = 'synergy';
        else
            AUC_synergy = nan;
        end
        hold on
        if sum(Ytest == 1) > 0
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
        
        if mode == 1
            stats.AUC_synergy(i,1) = AUC_synergy;
            stats.AUC_antagonism(i,1) = AUC_antagonism;
        else
            overview.AUC_synergy = AUC_synergy;
            overview.AUC_antagonism = AUC_antagonism;
        end
    end

    %% CONFUSION MATRICES
    function getConfusion(Ytest,Ypred,mode)
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
        
    end

    %% SCATTER PLOTS
    function getScatter(Ytest,Ypred,mode)
        % Scatter plot with line of best fit
        % Rank sorted - larger rank means bigger value
        [Ytest_sorted,I_test] = sort(Ytest);
        Ytest_rank = zeros(1,length(Ytest));
        Ytest_rank(I_test) = 1:length(Ytest);
        [Ypred_sorted,I_pred] = sort(Ypred);
        Ypred_rank = zeros(1,length(Ytest));
        Ypred_rank(I_pred) = 1:length(Ypred);
        
        scatter(Ytest_rank,Ypred_rank,24,'filled')
        h = lsline; h.Color = 'red'; 
        h.DisplayName = sprintf('R = %.2f', R);
        if mode == 1
            title(sprintf('Subset %d',i))
        else
            title('Overall scatter plot');
        end
        xlabel('Experimental')
        ylabel('Predicted')
        legend(h,'Location','northwest')   
    end
end
