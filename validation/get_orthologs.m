function orthologs = get_orthologs(filename, model_type, data_lookup)

data_table = readtable(data_lookup);

if strcmp(model_type, 'original_model') || strcmp(model_type, 'ecoli_model')
    orthology = data_table.Ecoli_Orthologs;
elseif strcmp(model_type, 'mtb_model')
    orthology = data_table.Mtb_Orthologs;
end


idx = find(strcmp(filename, data_table.Filename));

orthology_file = orthology{idx};

if ~isempty(orthology_file)
    orthologs = readcell(orthology_file);
else
    orthologs = '';
    
end
end