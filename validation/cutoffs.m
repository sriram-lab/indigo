function [synergy,antagonism] = cutoffs(filename)
%define cutoffs
data_table = readtable('all_data.xlsx');
idx = find(strcmp(filename, data_table.Filename));
synergy = data_table.Synergy(idx);
antagonism = data_table.Antagonism(idx);
end