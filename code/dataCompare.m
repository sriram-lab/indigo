function sharedInteractions = dataCompare(files)
arguments
    files cell
end
%files = cell array of interaction data files
%can be as many files as you want
%only applies when comparing drug interaction pairs
%function will do pairwise comparison between all data files and orthology
%common use case may be comparing training data file and test data file
names = erase(files,'.xlsx');
allFiles = cellstr(ls('data'));

%initialize arrays for storing comparison results
R_values = zeros(length(files));
P_values = zeros(length(files));
interactionCounts = zeros(length(files));
sharedInteractions = cell(length(files));
orthologCounts = zeros(length(files));
sharedOrthologs = cell(length(files));
%it doesn't matter whether A or B is bigger file
for i = 1:length(files)
    A_label = sprintf('%s score',names{i});   
    [A_scores, A_drugs] = xlsread(files{i});
    A_drugs = string(A_drugs); 
    for j = 1:length(files)    
        if j > i    %so that you don't get repeat comparisons
            %remove all rows that have more than 2 drugs
            B_label = sprintf('%s score',names{j});
            [B_scores, B_drugs] = xlsread(files{j});
            B_drugs = string(B_drugs);
            [Lia,Locb] = ismember(A_drugs, B_drugs, 'rows','legacy');
            A_scores_final = A_scores(Lia);
            Locb = nonzeros(Locb);
            B_scores_final = B_scores(Locb);
            drugsShared = A_drugs(Lia,:);
            if ~isempty(drugsShared)
                interactionCounts(i,j) = length(drugsShared);
                sharedInteractions{i,j} = {drugsShared,A_scores_final,B_scores_final};
                %Show interactions in a table
                T = table(drugsShared(:,1),drugsShared(:,2),A_scores_final,B_scores_final);
                T.Properties.VariableNames = {'Drug 1','Drug 2',A_label,B_label};
                sharedInteractionsPreview = head(T)
                [R_values(i,j),P_values(i,j)] = corr(A_scores_final, B_scores_final,'type','Spearman');
            else
                %if there are no shared interactions
                fprintf('%s and %s have no shared interactions\n\n',names{i},names{j});
            end
                        
            %check if both files have corresponding orthology files
            orthologyFile_A = strcat(erase(files{i},'.xlsx'),'_orthologs.xlsx');
            orthologyFile_B = strcat(erase(files{j},'.xlsx'),'_orthologs.xlsx');
            if sum(contains(allFiles,orthologyFile_A)) ~= 0 && sum(contains(allFiles,orthologyFile_B)) ~= 0
                [~,A_orthologs] = xlsread(orthologyFile_A);
                [~,B_orthologs] = xlsread(orthologyFile_B);
                Lia = ismember(A_orthologs, B_orthologs);
                orthologsShared = A_orthologs(Lia);
                if ~isempty(orthologsShared)
                    orthologCounts(i,j) = length(orthologsShared);
                    sharedOrthologs{i,j} = orthologsShared;
                end          
            end            
        end
    end
end
%tables
%number of shared drug pairs
interactionCounts = array2table(interactionCounts);
interactionCounts.Properties.VariableNames = names;
interactionCounts.Properties.RowNames = names;

%correlation values and significance for shared interactions
R_values = array2table(R_values);
R_values.Properties.VariableNames = names;
R_values.Properties.RowNames = names;
P_values = array2table(P_values);
P_values.Properties.VariableNames = names;
P_values.Properties.RowNames = names;

%orthology
orthologCounts = array2table(orthologCounts);
orthologCounts.Properties.VariableNames = names;
orthologCounts.Properties.RowNames = names;

%show tables
interactionCounts
R_values
P_values
orthologCounts
end
