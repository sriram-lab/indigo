function sharedOrthologs = orthologCompare(orthologs)
arguments
    orthologs cell
end
%orthologs = cell array of orthology files
names = erase(orthologs,'.xlsx');
orthologCounts = zeros(length(orthologs));
sharedOrthologs = cell(length(orthologs));
for i = 1:length(orthologs)
    [~,A_orthologs] = xlsread(orthologs{i});
    for j = 1:length(orthologs)  
        if j > i
            [~,B_orthologs] = xlsread(orthologs{j});
            Lia = ismember(A_orthologs, B_orthologs);
            orthologsShared = A_orthologs(Lia);
            if ~isempty(orthologsShared)
                orthologCounts(i,j) = length(orthologsShared);
                sharedOrthologs{i,j} = orthologsShared;
            else
                %if there are no shared interactions
                fprintf('%s and %s have no shared interactions\n\n',names{i},names{j});
            end
        else
            continue
        end
    end
end
%number of shared orthologous genes
orthologCounts = array2table(orthologCounts);
orthologCounts.Properties.VariableNames = names;
orthologCounts.Properties.RowNames = names;
orthologCounts
end
