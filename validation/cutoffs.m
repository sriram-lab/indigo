function [synergy,antagonism] = cutoffs(filename, dataLookup)
%define cutoffs
data_table = readtable(dataLookup);
idx = find(strcmp(filename, data_table.Filename));
synergy = data_table.Synergy(idx);
antagonism = data_table.Antagonism(idx);
end