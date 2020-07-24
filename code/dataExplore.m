function [summary,drugList] = dataExplore(filename)

% filename = strcat("/indigoData/",filename);
dataName = erase(filename,'.xlsx');
summary = struct;
%You should make a matlab live script where you run this on each
%scores set and store it as a notebook
%Descriptive statistics
%Either have this plug in to other function, or just make this a mlx file
%and output the analysis for everything so I can always refer to it
%Check distribution and characterize scores
%Can use to determine whether to use zscore or not
%Need to see whether scores has good spread of synergy and antagonism
%which is also based on cutoffs defined by a given paper
[synergy,antagonism] = cutoffs(filename);
%start of descriptive statistics
[scores,drugs] = xlsread(filename);

%Descriptive statistics
summary.min = min(scores);
summary.median = median(scores);
summary.max = max(scores);
summary.range = range(scores);
summary.mean = mean(scores);
summary.std = std(scores);
summary.interactionCount = length(scores);
summary.synergyCount = sum((scores < synergy));
summary.antagonismCount = sum((scores > antagonism));

%if orthology
%count number of orthologs
files = cellstr(ls('data'));
orthology = strcat(erase(filename,'.xlsx'),'_orthologs.xlsx');
if sum(contains(files,orthology)) ~= 0
    [~,orth] = xlsread(orthology);
    summary.orthologsCount = length(orth);
end

%Anderson Darling Normality test, if h = 1 and p < 0.05, then scores is not
%normal
[summary.rejectNormality,summary.pNormality] = adtest(scores);
summaryTable = struct2table(summary);
summaryTable
%list of drugs used
drugList = unique(drugs);

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

% write norm scores to second sheet
filepath = strcat('data/',filename);
writecell([drugs,num2cell(normScores)],filepath,'Sheet',2)
subplot(2,2,2)

normplot(scores)

subplot(2,2,3)


boxplot(scores)
xlabel(dataName,'Interpreter','none')
ylabel('interaction scores')

subplot(2,2,4)
x = 1:length(scores);
%scatter plot grouped by synergy, neutral, antagonism
grouping = zeros(length(scores),1);
grouping(scores < synergy) = -1;
grouping(scores > synergy & scores < antagonism) = 0;
grouping(scores > antagonism) = 1;
gscatter(x,scores,grouping,[0 1 0; 0 0 1; 1 0 0])
legend('synergistic','neutral','antagonistic')
xlabel('interactions')
ylabel('scores')
end