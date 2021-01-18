function [antagonism_scores_ranked, antagonism_interactions_ranked] = find_antagonism(filename)

data = readcell(filename);
scores = cell2mat(data(:,end));
interactions = data(:,1:end-1);
[~,antagonism] = cutoffs(filename);

antagonism_interactions = interactions(scores > antagonism,:);
antagonism_scores = scores(scores > antagonism);

[antagonism_scores_ranked, idx] = sort(antagonism_scores,'descend');

antagonism_interactions_ranked = antagonism_interactions(idx,:);
end