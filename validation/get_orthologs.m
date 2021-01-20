function [orthologs,orthologyFile] = get_orthologs(filename,modelType)

data_table = readtable('all_data_fu.xlsx');

if modelType == 1
    orthology = data_table.EcoliOrthologs;
elseif modelType == 2
    orthology = data_table.MtbOrthologs;
end


idx = find(strcmp(filename, data_table.Filename));

orthologyFile = orthology{idx};

if ~isempty(orthologyFile)
    orthologs = readcell(orthologyFile);
else
    orthologs = '';
    
end
end