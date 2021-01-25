function orthologs = get_orthologs(filename,modelType)

data_table = readtable('all_data.xlsx');

if strcmp(modelType, 'ecoli_model')
    orthology = data_table.EcoliOrthologs;
elseif strcmp(modelType, 'mtb_model')
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