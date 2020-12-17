function [antagonism_scores_ranked, antagonism_interactions_ranked] = find_antagonism(filename)

[scores,interactions] = xlsread(filename);
[~,antagonism] = cutoffs(filename);

antagonism_interactions = interactions(scores > antagonism,:);
antagonism_scores = scores(scores > antagonism);

[antagonism_scores_ranked, idx] = sort(antagonism_scores,'descend');

antagonism_interactions_ranked = antagonism_interactions(idx,:);
end