function [phenotype_data, phenotype_labels] = process_transcriptome_tb(expression_data,expression_data_rowlabels,conditions)

phenotype_num = expression_data;
plist_bnums = expression_data_rowlabels;

z = 1;
for i = 1:length(conditions)
    te = plist_bnums(phenotype_num(:,i) < -z);  % downregulation
    cell_z1_list_t(1:length(te),i) = te;
    lte(i) = length(te);
end

cell_z1_list_t = regexprep(cell_z1_list_t,'''','');
cell_z1_list_t = regexprep(cell_z1_list_t,'[]','');
%phenotype_labels0 = unique(cell_z1_list_t(:));
phenotype_labels0 = unique(cell_z1_list_t (~cellfun(@isempty, cell_z1_list_t)));


clear nicholslistix_t
for i = 1:length(conditions)
     cell_list_column = cell_z1_list_t(:,i);
     nicholslistix_t(:,i) = ismember(phenotype_labels0,cell_list_column(~cellfun('isempty',cell_list_column)));

%    nicholslistix_t(:,i) = ismember(phenotype_labels0,cell_z1_list_t(:,i));
end
%phenotype_data = nicholslistix_t;

for i = 1:length(conditions)
    te1 = plist_bnums(phenotype_num(:,i) > z);  %upregulation
    cell_z1_list_t1(1:length(te1),i) = te1;
    lte1(i) = length(te1);
end
cell_z1_list_t1 = regexprep(cell_z1_list_t1,'''','');
%phenotype_labels1 = unique(cell_z1_list_t1(:));
phenotype_labels1 = unique(cell_z1_list_t1(~cellfun(@isempty, cell_z1_list_t1)));

clear nicholslistix_t1
for i = 1:length(conditions)
    cell_list_column = cell_z1_list_t1(:,i);
    nicholslistix_t1(:,i) = ismember(phenotype_labels1,cell_list_column(~cellfun('isempty',cell_list_column)));

%    nicholslistix_t1(:,i) = ismember(phenotype_labels1,cell_z1_list_t1(:,i));
end
phenotype_data = [nicholslistix_t;nicholslistix_t1];
phenotype_labels = [phenotype_labels0;phenotype_labels1];

end