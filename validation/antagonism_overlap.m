%Find overlap with other files
function sharedAntagonism = antagonism_overlap(files)
arguments
    files cell
end

sharedAntagonism = cell(length(files));

for i = 1:length(files)
    [antagonism_scores_A, antagonism_interactions_A] = find_antagonism(files{i});
    antagonism_interactions_A = string(antagonism_interactions_A); 
    
    for j = 1:length(files)  
        if j > i
            [antagonism_scores_B, antagonism_interactions_B] = find_antagonism(files{j});
            antagonism_interactions_B = string(antagonism_interactions_B); 

            [Lia,Locb] = ismember(antagonism_interactions_A, antagonism_interactions_B, 'rows','legacy');
            antagonism_interactions_A_final = antagonism_interactions_A(Lia,:);
            antagonism_scores_A_final = antagonism_scores_A(Lia);


            Locb = nonzeros(Locb);
            antagonism_scores_B_final = antagonism_scores_B(Locb);

            drugsShared = antagonism_interactions_A_final;
           
            %Write this to an excel file!
            if ~isempty(drugsShared)
                drugsShared_fullnames = get_fullnames(drugsShared);
                output = [cellstr(drugsShared),drugsShared_fullnames, num2cell(antagonism_scores_A_final),num2cell(antagonism_scores_B_final)];
                sheetName = sprintf('%s_%s', erase(files{i},'.xlsx'),erase(files{j},'.xlsx'));
                writecell(output,'antagonism_interactions.xlsx','Sheet',sheetName);
                sharedAntagonism{i,j} = output;
            end    
        end
    end
end
end
