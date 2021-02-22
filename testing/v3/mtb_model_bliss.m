%% independent

% data_lookup = 'all_data.xlsx';
% indigo_data = readtable(data_lookup);
% indigo_data = indigo_data(strcmp(indigo_data.Scoring, 'Bliss'),:);
% data_files = indigo_data.Filename(~indigo_data.Ecoli_Model_Exclusive);
% 
% valmethod = 'independent';
% K = 1;
% standardize = '';
% model_type = 'mtb_model';
% input_type = 2;
% scoring = 'bliss';
% 
% for i=1:length(data_files)
%     test_data = data_files{i};
%     training_data = data_files([1:i-1 i+1:end]);
%     indigo_summary = indigo_run(test_data,training_data,data_lookup,valmethod,K,standardize,model_type,input_type,scoring);
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,1);
%     save(save_indigo(indigo_summary, 1))
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,2);
%     save(save_indigo(indigo_summary, 2))
% end

%% independent (z)

% valmethod = 'independent';
% K = 1;
% standardize = 'z_score';
% model_type = 'mtb_model';
% input_type = 2;
% scoring = 'bliss';
% 
% for i=1:length(data_files)
%     test_data = data_files{i};
%     training_data = data_files([1:i-1 i+1:end]);
%     indigo_summary = indigo_run(test_data,training_data,data_lookup,valmethod,K,standardize,model_type,input_type,scoring);
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,1);
%     save(save_indigo(indigo_summary, 1))
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,2);
%     save(save_indigo(indigo_summary, 2))
% end

%% Kfold

% valmethod = 'Kfold';
% K = 5;
% standardize = '';
% model_type = 'mtb_model';
% input_type = 2;
% scoring = 'bliss';
% 
% for i=1:length(data_files)
%     test_data = data_files{i};
%     training_data = data_files([1:i-1 i+1:end]);
%     indigo_summary = indigo_run(test_data,training_data,data_lookup,valmethod,K,standardize,model_type,input_type,scoring);
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,1);
%     save(save_indigo(indigo_summary, 1))
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,2);
%     save(save_indigo(indigo_summary, 2))
% end

%% Kfold (z)

valmethod = 'Kfold';
K = 5;
standardize = 'z_score';
model_type = 'mtb_model';
input_type = 2;
scoring = 'bliss';

for i=4:length(data_files)
    test_data = data_files{i};
    training_data = data_files([1:i-1 i+1:end]);
    indigo_summary = indigo_run(test_data,training_data,data_lookup,valmethod,K,standardize,model_type,input_type,scoring);
    [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,1);
    save(save_indigo(indigo_summary, 1))
    [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,2);
    save(save_indigo(indigo_summary, 2))
end
