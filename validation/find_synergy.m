function [synergy_scores_ranked, synergy_interactions_ranked] = find_synergy(filename)

data = readcell(filename);
scores = cell2mat(data(:,end));
interactions = data(:,1:end-1);
synergy = cutoffs(filename);

synergy_interactions = interactions(scores < synergy,:);
synergy_scores = scores(scores < synergy);

[synergy_scores_ranked, idx] = sort(synergy_scores);

synergy_interactions_ranked = synergy_interactions(idx,:);
end