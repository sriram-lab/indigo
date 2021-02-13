function [synergy,antagonism] = cutoffs(filename, data_lookup)
%define cutoffs
data_table = readtable(data_lookup);
idx = find(strcmp(filename, data_table.Filename));
synergy = data_table.Synergy(idx);
antagonism = data_table.Antagonism(idx);
end