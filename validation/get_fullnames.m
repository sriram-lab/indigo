%Get fullnames
function fullnames = get_fullnames(interactions)

identifiers_match = readcell('identifiers_match.xlsx');
drugnames = identifiers_match(:,2);
fullnames = cell(size(interactions));
for i = 1:size(interactions,2)
    [~,Locb] = ismember(interactions(:,i),identifiers_match);
    fullnames(:,i) = drugnames(Locb);
end