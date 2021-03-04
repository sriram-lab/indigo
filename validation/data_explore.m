function [summary,drug_list] = data_explore(filename, data_label, data_lookup, number_ref)

    fprintf(sprintf('Exploring data for %s\n',filename));

    figure
    % filename = strcat("/indigoData/",filename);
%     data_name = erase(filename,'.xlsx');
    summary = struct;


    [synergy,antagonism] = cutoffs(filename, data_lookup);

    data = readtable(filename,"ReadVariableNames",false);
    scores = data{:,end};
    drugs = data{:,1:end-1};

    ecoli_orthologs = get_orthologs(filename, 'ecoli_model', data_lookup);
    mtb_orthologs = get_orthologs(filename, 'mtb_model', data_lookup);

    %Descriptive statistics
    summary.min = min(scores);
    summary.max = max(scores);
    summary.range = range(scores);
    summary.median = median(scores);
    summary.mean = mean(scores);
    summary.std = std(scores);
    summary.interaction_count = length(scores);
    summary.synergy_count = sum((scores <= synergy));
    summary.antagonism_count = sum((scores >= antagonism));
    summary.neutral_count = sum((scores > synergy & scores < antagonism));
    summary.ecoli_ortholog_count = length(ecoli_orthologs);
    summary.mtb_ortholog_count = length(mtb_orthologs);

    %Anderson Darling Normality test, if h = 1 and p < 0.05, then scores is not
    %normal
    [summary.reject_normality,summary.p_normality] = adtest(scores);
    %list of unique drugs used
    drug_list = [];
    for i = 1:size(drugs,2)
        drug_list = [drug_list; drugs(:,i)];
    end
    drug_list(~cellfun('isempty',drug_list));
    drug_list = unique(drug_list)
    summary.uniqueDrugs = length(drug_list);
    summary_table = struct2table(summary);
    summary_table
%     writetable(summary_table,'exploratory_stats.xlsx','WriteVariableNames',false,'Range',sprintf('B%d:O%d',number_ref+1,number_ref+1))
    sgtitle(sprintf('Data exploration of %s', data_label),'Interpreter','none');
    subplot(2,2,1)
    histogram(scores)
    hold on
    norm_scores = zscore(scores);
    histogram(norm_scores)
    legend('Raw Scores','Z-Scores')
    hold off
    ax = gca;
    ax.FontWeight = 'bold';
    ax.FontSize = 12;
    xlabel('Interaction scores','FontWeight',"bold","FontSize",16)
    ylabel('Frequency', "FontWeight","bold",'FontSize',16)
    title('Histogram of Interaction Scores', 'FontWeight', 'bold')
    ax.TitleFontSizeMultiplier = 1.5;
    
    subplot(2,2,2)
    normplot(scores)
    ax = gca;
    ax.FontWeight = 'bold';
    ax.FontSize = 12;
    xlabel('Data','FontWeight',"bold","FontSize",16)
    ylabel('Probability', "FontWeight","bold",'FontSize',16)
    title('Normal Probability Plot', 'FontWeight', 'bold')
    ax.TitleFontSizeMultiplier = 1.5;

    subplot(2,2,3)

    boxplot(scores)
    hold on 
    x = ones(size(scores)).*(1+(rand(size(scores))-0.5)/8);
    %scatter plot grouped by synergy, neutral, antagonism
    grouping = zeros(length(scores),1);
    grouping(scores <= synergy) = -1;
    grouping(scores > synergy & scores < antagonism) = 0;
    grouping(scores >= antagonism) = 1;
    gscatter(x,scores,grouping,[0 0 1; 0 0 0; 1 0 0])
    legend('synergistic','neutral','antagonistic')
%     gscatter(x,scores,grouping,[0 0 1; 1 0 0 ])
%     legend('synergistic','antagonistic')
    ax = gca;
    ax.FontWeight = 'bold';
    ax.FontSize = 12;
    xlabel('Interactions','FontWeight',"bold","FontSize",16)
    ylabel('Scores', "FontWeight","bold",'FontSize',16)
    title('Boxplot of Interaction Scores with Scatter Plot Overlay', 'FontWeight', 'bold')
    ax.TitleFontSizeMultiplier = 1.5;

end
