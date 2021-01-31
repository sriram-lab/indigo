function [summary,drugList] = dataExplore(filename, dataLookup)

    fprintf(sprintf('Exploring data for %s\n',filename));

    figure
    % filename = strcat("/indigoData/",filename);
    dataName = erase(filename,'.xlsx');
    summary = struct;


    [synergy,antagonism] = cutoffs(filename, dataLookup);

    data = readtable(filename,"ReadVariableNames",false);
    scores = data{:,end};
    drugs = data{:,1:end-1};

    ecoli_orthologs = get_orthologs(filename, 'ecoli_model', dataLookup);
    mtb_orthologs = get_orthologs(filename, 'mtb_model', dataLookup);

    %Descriptive statistics
    summary.min = min(scores);
    summary.median = median(scores);
    summary.max = max(scores);
    summary.range = range(scores);
    summary.mean = mean(scores);
    summary.std = std(scores);
    summary.interactionCount = length(scores);
    summary.synergyCount = sum((scores <= synergy));
    summary.antagonismCount = sum((scores >= antagonism));
    summary.neutralCount = sum((scores > synergy & scores < antagonism));
    summary.ecoliOrthologCount = length(ecoli_orthologs);
    summary.mtbOrthologCount = length(mtb_orthologs);

    %Anderson Darling Normality test, if h = 1 and p < 0.05, then scores is not
    %normal
    [summary.rejectNormality,summary.pNormality] = adtest(scores);
    %list of unique drugs used
    drugList = [];
    for i = 1:size(drugs,2)
        drugList = [drugList; drugs(:,i)];
    end
    drugList(~cellfun('isempty',drugList));
    drugList = unique(drugList)
    summary.uniqueDrugs = length(drugList);
    summaryTable = struct2table(summary);
    summaryTable

    sgtitle(sprintf('Data exploration of %s', dataName),'Interpreter','none');
    subplot(2,2,1)
    histogram(scores)
    xlabel('Interaction scores')
    ylabel('Frequency')

    hold on
    normScores = zscore(scores);
    histogram(normScores)
    hold off
    legend('original scores','z-scores')
    subplot(2,2,2)

    normplot(scores)

    subplot(2,2,3)

    boxplot(scores)
    xlabel(dataName,'Interpreter','none')
    ylabel('interaction scores')
    hold on 
    x = ones(length(scores),1);
    %scatter plot grouped by synergy, neutral, antagonism
    grouping = zeros(length(scores),1);
    grouping(scores <= synergy) = -1;
    grouping(scores > synergy & scores < antagonism) = 0;
    grouping(scores >= antagonism) = 1;
    gscatter(x,scores,grouping,[0 1 0; 0 0 1; 1 0 0])
    legend('synergistic','neutral','antagonistic')
    xlabel('interactions')
    ylabel('scores')
end