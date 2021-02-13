function [stats,averages,overview] = analyze(indigo_summary,result_index, data_files)

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
    dataname = erase(indigo_summary.test_data,'.xlsx');
    
    filename = strcat(erase(indigo_summary.test_data,'.xlsx'),'_',indigo_summary.valmethod);
    sheetname = indigo_summary.valmethod;
    if strcmp(indigo_summary.standardize,'z_score')
        filename = sprintf('%s_z', filename);
        sheetname = strcat(indigo_summary.valmethod,' (z)');
    end
     
    if strcmp(indigo_summary.scoring,'bliss') || strcmp(indigo_summary.scoring,'loewe')
        overview_file = strcat('results/v2/', indigo_summary.model_type, '/', ...
            indigo_summary.scoring, '/','overview_results.xlsx');
        results_file = strcat('results/v2/', indigo_summary.model_type, '/', ...
            indigo_summary.scoring, '/', indigo_summary.valmethod, '/', filename);
    else
        overview_file = strcat('results/v2/',indigo_summary.model_type, '/', ...
            'overview_results.xlsx');
        results_file = strcat('results/v2/',indigo_summary.model_type, '/', ...
            indigo_summary.valmethod, '/', filename);
    end
    
    %% Process results
    fprintf(sprintf('Results for %s',dataname))
    stats = struct;
    averages = struct;
    
    %% SUBSET DATA FOR KFOLD CV
    if strcmp(indigo_summary.valmethod,'Kfold_onself') || strcmp(indigo_summary.valmethod,'Kfold')
        drugs_total = [];
        Ytest_total = [];
        Ypred_total = [];
        
        f1 = figure(1);
        f1.Visible = 'off';
        t1 = tiledlayout(f1,3,2);
        title(t1,sprintf('ROC Curves for %s\n', dataname), ...
              'FontWeight','bold', 'Interpreter', 'None')
        
        for i = 1:indigo_summary.K
            drugs_total = [drugs_total; indigo_summary.testPairs{i}];
            Ytest_total = [Ytest_total; indigo_summary.testScores{i}];
            Ypred_total = [Ypred_total; indigo_summary.predictedScores{i}];
            Ytest = indigo_summary.testScores{i};
            Ypred = indigo_summary.predictedScores{i};
            
            %correlation 
            [R,P] = corr(Ytest, Ypred,'type','Spearman');
            
            [Ytest_class, Ypred_class] = classify_scores(Ytest,Ypred);   
            
            % 1 means analysis per subset
            get_stats(Ytest_class, Ypred_class, 1)
            
            nexttile
            get_roc(Ytest_class, Ypred_class, 1)
        end
        
        fig_file = strcat(results_file,'_roc');
        saveas(f1,fig_file,'fig') 
        saveas(f1,fig_file,'png')
        close(f1)
        
        f2 = figure(2);
        f2.Visible = 'off';
        t2 = tiledlayout(f2,3,2);
        title(t2, sprintf('Confusion Matrices for %s\n', dataname), ...
              'FontWeight','bold','Interpreter','None')
        
        for i = 1:indigo_summary.K
            nexttile
            get_confusion(Ytest_class, Ypred_class, 1)
        end
        
        fig_file = strcat(results_file,'_cm');
        saveas(f2,fig_file,'fig')
        saveas(f2,fig_file,'png')
        close(f2)
        
        f3 = figure(3);
        f3.Visible = 'off';
        t3 = tiledlayout(f3,3,2);
        title(t3, sprintf('Scatter Plots (Ranked Sorted) for %s\n', dataname), ...
              'FontWeight','bold','Interpreter','None')
        
        for i = 1:indigo_summary.K     
            nexttile
            get_scatter(Ytest, Ypred, 1)  
        end
        
        fig_file = strcat(results_file,'_sc');
        saveas(f3,fig_file,'fig')
        saveas(f3,fig_file,'png')
        close(f3)
        
    else
        drugs_total = [indigo_summary.testPairs{1}];
        Ytest_total = [indigo_summary.testScores{1}];
        Ypred_total = [indigo_summary.predictedScores{1}];
    end

    %% OVERALL DATA - ALL RESULTS COMBINED
    
    overview.drugInteractions = drugs_total;
    overview.experimentalScores = Ytest_total;
    overview.predictedScores = Ypred_total;
    
    [R,P] = corr(Ytest_total, Ypred_total,'type','Spearman');

    [Ytest_class_total, Ypred_class_total] = classify_scores(Ytest_total, Ypred_total);   

    % 2 means analysis for all results combined
    get_stats(Ytest_class_total, Ypred_class_total, 2)
    f4 = figure(4);
    f4.Visible = 'off';
    t4 = tiledlayout(f4,3,1);
    title(t4,sprintf('Overall INDIGO Results for %s\n',dataname), ...
        'FontWeight','bold','Interpreter','none')
    nexttile
    get_roc(Ytest_class_total, Ypred_class_total, 2)    
    nexttile
    get_confusion(Ytest_class_total, Ypred_class_total, 2)
    nexttile
    get_scatter(Ytest_total, Ypred_total, 2)  
    fig_file = strcat(results_file,'_overall');
    saveas(f4,fig_file,'fig')
    saveas(f4,fig_file,'png')
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
                
    if strcmp(indigo_summary.valmethod,'Kfold_onself') || strcmp(indigo_summary.valmethod,'Kfold')
        results_table = struct2table(stats);
        results_table.Properties.VariableNames = varnames;
        rownames = cell(1,indigo_summary.K);
        for i = 1:indigo_summary.K
            rownames{i} = sprintf('Subset %d', i);
        end
        results_table.Properties.RowNames = rownames;

        fields = fieldnames(stats);
        tableArray = zeros(length(fields),1);
        for i = 1:length(fields)
            averages.(fields{i}) = mean(stats.(fields{i}));
            tableArray(i) = mean(stats.(fields{i}));
        end
        averages_table = table(tableArray);
        averages_table.Properties.RowNames = varnames; 
        averages_table.Properties.VariableNames = {'Value'};
    end

    fields = fieldnames(overview);
    tableArray = zeros(length(fields)-3,1);  %don't include first 3 fields
    for i = 4:length(fields)
        tableArray(i-3) = overview.(fields{i});
    end
    overview_table = table(tableArray);
    overview_table.Properties.RowNames = varnames;
    overview_table.Properties.VariableNames = {'Value'};
    
    %% Display analysis and save to .xlsx files

    if strcmp(indigo_summary.valmethod,'Kfold_onself') || strcmp(indigo_summary.valmethod,'Kfold')
        results_table
        averages_table
        overview_table
        writetable(results_table,strcat(results_file,'.xlsx'),'Sheet','results','WriteRowNames',true)
        writetable(averages_table,strcat(results_file,'.xlsx'),'Sheet','averages','WriteRowNames',true)
        writetable(overview_table,strcat(results_file,'.xlsx'),'Sheet','overview','WriteRowNames',true)  
    else
        overview_table
        writetable(overview_table,strcat(results_file,'.xlsx'),'Sheet','overview','WriteRowNames',true)
    end
    
    overview_table = rows2vars(overview_table,'VariableNamingRule','preserve');
    overview_table = overview_table(:,2:end);    
    
    writecell(data_files, overview_file, 'Sheet', sheetname, 'Range', ...
        sprintf('A2:A%d', length(data_files)+1))
    writecell([{'Datasets'}, varnames], overview_file, 'Sheet', sheetname, 'Range', 'A1:L1')
    writetable(overview_table, overview_file, 'Sheet', sheetname, 'Range', ...
        sprintf('B%d:L%d',result_index+1,result_index+1), 'WriteVariableNames',false)

    %% CLASSIFY SCORES
    function [Ytest_class, Ypred_class] = classify_scores(Ytest, Ypred)
        [synergy_cutoff, antagonism_cutoff] = cutoffs(indigo_summary.test_data, indigo_summary.data_lookup);
        Ytest_class = zeros(length(Ytest),1);
        Ypred_class = zeros(length(Ypred),1);
        %synergy - assign value = -1
        Ytest_class(Ytest <= synergy_cutoff) = -1;
        Ypred_class(Ypred <= synergy_cutoff) = -1;
        %antagonism - assign value = 1
        Ytest_class(Ytest >= antagonism_cutoff) = 1;
        Ypred_class(Ypred >= antagonism_cutoff) = 1;
    end

    %% GET STATISTICS
    function get_stats(Ytest_class, Ypred_class, mode)
        interaction_count = length(Ytest_class);
        %Get accuracy, precision and recall
        accuracy = sum(Ypred_class == Ytest_class)/length(Ypred_class);
        abs_error = mean(abs(Ypred_class - Ytest_class));
        precision_synergy = sum(Ypred_class == -1 & Ytest_class == -1)/sum(Ypred_class == -1);
        recall_synergy = sum(Ypred_class == -1 & Ytest_class == -1)/sum(Ytest_class == -1);
        precision_antagonism = sum(Ypred_class == 1 & Ytest_class == 1)/sum(Ypred_class == 1);
        recall_antagonism = sum(Ypred_class == 1 & Ytest_class == 1)/sum(Ytest_class == 1);
        
        if mode == 1
            stats.interaction_count(i,1) = interaction_count;
            stats.rank_corr(i,1) = R;
            stats.p_value(i,1) = P;
            stats.accuracy(i,1) = accuracy;
            stats.abs_error(i,1) = abs_error;
            stats.precision_synergy(i,1) = precision_synergy;
            stats.recall_synergy(i,1) = recall_synergy;
            stats.precision_antagonism(i,1) = precision_antagonism;
            stats.recall_antagonism(i,1) = recall_antagonism;
        else
            overview.interaction_count = interaction_count;
            overview.rank_corr = R;
            overview.p_value = P;
            overview.accuracy = accuracy;
            overview.abs_error = abs_error;
            overview.precision_synergy = precision_synergy;
            overview.recall_synergy = recall_synergy;
            overview.precision_antagonism = precision_antagonism;
            overview.recall_antagonism = recall_antagonism;
        end
    end
    
    %% ROC CURVES
    function get_roc(Ytest_class,Ypred_class,mode)
        labels = {};  
        if sum(Ytest_class == -1) > 0
            %ROC Curve for synergy
            [X_synergy,Y_synergy,T_synergy,AUC_synergy] = perfcurve(Ytest_class,Ypred_class,-1);
            plot(X_synergy,Y_synergy)
            labels{end+1} = 'synergy';
        else
            AUC_synergy = nan;
        end
        hold on
        if sum(Ytest_class == 1) > 0
            %ROC Curve for antagonism
            [X_antagonism,Y_antagonism,T_antagonism,AUC_antagonism] = perfcurve(Ytest_class,Ypred_class,1);
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
    function get_confusion(Ytest_class,Ypred_class,mode)
        C = confusionmat(Ytest_class,Ypred_class);
        tempLabels = {};
        if sum(Ytest_class == -1) > 0 || sum(Ypred_class == -1) > 0
            tempLabels{end+1} = 'synergistic';
        end
        if sum(Ytest_class == 0) > 0 || sum(Ypred_class == 0) > 0
            tempLabels{end+1} = 'neutral';
        end
        if sum(Ytest_class == 1) > 0 || sum(Ypred_class == 1) > 0
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
    function get_scatter(Ytest,Ypred,mode)
        % Scatter plot with line of best fit
        % Rank sorted - larger rank means bigger value
        [~,I_test] = sort(Ytest);
        Ytest_rank = zeros(1,length(Ytest));
        Ytest_rank(I_test) = 1:length(Ytest);
        [~,I_pred] = sort(Ypred);
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
