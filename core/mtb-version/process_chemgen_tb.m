function [phenotype_data, phenotype_labels, conditions] = process_chemgen_tb(fname,z)
%[phenotype_data, phenotype_labels, conditions] = process_chemgen(filename)
% processes chemogenomic data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
if isempty(fname)
    fname = 'ecoli_phenotype_data_cell.xlsx';
end

data = readtable(fname,'VariableNamingRule','preserve');
phenotype_num = data{1:end,2:end};    % numerical data
probelist = data.Gene;                % gene names (ECK numbers)
conditions = data.Properties.VariableNames(2:end)';  %list of conditions
 
if any(isnan(phenotype_num(:,1)))
    phenotype_num(:,1) = '';
end

if ~exist('z','var') || isempty(z)
    z = 1;
end
%disp(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% converting gene ids to standard ids
% load ecoli_annotation_data1 genenames_array      genenames_array1    genenames_bnums
% 
% clear plist
% plist = cell(size(probelist));
% for i = 1:length(probelist)
%     tem = regexp(probelist{i},'-','split');
%     try
%         plist(i) = tem(2);
%     catch
%         plist(i) = tem(1);
%     end
% end
% plist = regexprep(plist,'''','');
% 
% [ix, pos] =ismember(upper(plist),upper(genenames_array));
% plist_bnums = plist;
% plist_bnums(ix) = genenames_bnums(pos(ix));
plist_bnums = probelist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalizing data - ask Sriram!!!
phenotype_num = quantilenorm(phenotype_num);
%z = 2;
for i = 1:length(conditions)
    te = plist_bnums(phenotype_num(:,i) < -z);
    cell_z1_list_t(1:length(te),i) = te;
    lte(i) = length(te);
end
cell_z1_list_t = regexprep(cell_z1_list_t,'''','');
phenotype_labels = unique(cell_z1_list_t(:));

clear nicholslistix_t
for i = 1:length(conditions)
    nicholslistix_t(:,i) = ismember(phenotype_labels,cell_z1_list_t(:,i));
end
phenotype_data = nicholslistix_t;

end
