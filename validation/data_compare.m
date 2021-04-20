function shared_interactions_table = data_compare(files)
arguments
    files cell
end
%files = cell array of interaction data files
%can be as many files as you want
%only applies when comparing drug interaction pairs
%function will do pairwise comparison between all data files and orthology
%common use case may be comparing training data file and test data file
names = erase(files,'.xlsx');

%initialize arrays for storing comparison results
R_values = zeros(length(files));
P_values = zeros(length(files));
interaction_counts = zeros(length(files));
shared_interactions = cell(length(files));

%it doesn't matter whether A or B is bigger file
for i = 1:length(files)
    A_label = sprintf('%s score',names{i}); 
    data = readcell(files{i},'Sheet',2);
    A_scores = cell2mat(data(:,end));
    A_drugs = data(:,1:end-1);
    A_drugs = string(A_drugs); 
    for j = 1:length(files)    
        if j > i    %so that you don't get repeat comparisons
            %remove all rows that have more than 2 drugs
            B_label = sprintf('%s score',names{j});
            data = readcell(files{j},'Sheet',2);
            B_scores = cell2mat(data(:,end));
            B_drugs = data(:,1:end-1);
            B_drugs = string(B_drugs); 
            [Lia,Locb] = ismember(A_drugs, B_drugs, 'rows','legacy');
            A_scores_final = A_scores(Lia);
            Locb = nonzeros(Locb); 
            B_scores_final = B_scores(Locb);
            drugs_shared = A_drugs(Lia,:);
            if ~isempty(drugs_shared)
                interaction_counts(i,j) = length(drugs_shared);
                shared_interactions{i,j} = {drugs_shared,A_scores_final,B_scores_final};
                %Show interactions in a table
                T = table(drugs_shared(:,1),drugs_shared(:,2),A_scores_final,B_scores_final);
                T.Properties.VariableNames = {'Drug 1','Drug 2',A_label,B_label};
                shared_interactions_preview = head(T)
                [R_values(i,j),P_values(i,j)] = corr(A_scores_final, B_scores_final,'type','Spearman');
            else
                %if there are no shared interactions
                fprintf('%s and %s have no shared interactions\n\n',names{i},names{j});
            end
        end
    end
end
%tables
%number of shared drug pairs
interaction_counts = array2table(interaction_counts);
interaction_counts.Properties.VariableNames = names;
interaction_counts.Properties.RowNames = names;

%correlation values and significance for shared interactions
R_values = array2table(R_values);
R_values.Properties.VariableNames = names;
R_values.Properties.RowNames = names;
P_values = array2table(P_values);
P_values.Properties.VariableNames = names;
P_values.Properties.RowNames = names;

%show tables
interaction_counts
R_values
P_values

shared_interactions_table = cell2table(shared_interactions);
shared_interactions_table.Properties.VariableNames = files;
shared_interactions_table.Properties.RowNames = files;
end
