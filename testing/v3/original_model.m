%% independent 

data_lookup = 'all_data.xlsx';
indigo_data = readtable(data_lookup);
data_files = indigo_data.Filename(~indigo_data.Mtb_Model_Exclusive);

% valmethod = 'independent';
% K = 1;
% standardize = '';
% model_type = 'original_model';
% input_type = 2;
% scoring = '';
% 
% for i=1:length(data_files)
%     if strcmp(data_files{i}, 'ecoli_mg1655_chandrasekaran_2016.xlsx')
%         continue
%     end
%     test_data = data_files{i};
%     training_data = 'ecoli_mg1655_chandrasekaran_2016.xlsx';
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
% model_type = 'original_model';
% input_type = 2;
% scoring = '';
% 
% for i=1:length(data_files)
%     if strcmp(data_files{i}, 'ecoli_mg1655_chandrasekaran_2016.xlsx')
%         continue
%     end
%     test_data = data_files{i};
%     training_data = 'ecoli_mg1655_chandrasekaran_2016.xlsx';
%     indigo_summary = indigo_run(test_data,training_data,data_lookup,valmethod,K,standardize,model_type,input_type,scoring);
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,1);
%     save(save_indigo(indigo_summary, 1))
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,2);
%     save(save_indigo(indigo_summary, 2))
% end

%% Kfold

valmethod = 'Kfold';
K = 5;
standardize = '';
model_type = 'original_model';
input_type = 2;
scoring = '';

test_data = 'ecoli_mc4100_loewe_russ_2018.xlsx';
training_data = 'ecoli_mg1655_chandrasekaran_2016.xlsx';
indigo_summary = indigo_run(test_data,training_data,data_lookup,valmethod,K,standardize,model_type,input_type,scoring);
[stats,averages,overview] = analyze(indigo_summary,12, indigo_data, data_files,1);
save(save_indigo(indigo_summary, 1))
[stats,averages,overview] = analyze(indigo_summary,12, indigo_data, data_files,2);
save(save_indigo(indigo_summary, 2))

% for i=1:length(data_files)
%     if strcmp(data_files{i}, 'ecoli_mg1655_chandrasekaran_2016.xlsx')
%         continue
%     end
%     test_data = data_files{i};
%     training_data = 'ecoli_mg1655_chandrasekaran_2016.xlsx';
%     indigo_summary = indigo_run(test_data,training_data,data_lookup,valmethod,K,standardize,model_type,input_type,scoring);
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,1);
%     save(save_indigo(indigo_summary, 1))
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,2);
%     save(save_indigo(indigo_summary, 2))
% end

%% Kfold (z)

% valmethod = 'Kfold';
% K = 5;
% standardize = 'z_score';
% model_type = 'original_model';
% input_type = 2;
% scoring = '';
% 
% for i=1:length(data_files)
%     if strcmp(data_files{i}, 'ecoli_mg1655_chandrasekaran_2016.xlsx')
%         continue
%     end
%     test_data = data_files{i};
%     training_data = 'ecoli_mg1655_chandrasekaran_2016.xlsx';
%     indigo_summary = indigo_run(test_data,training_data,data_lookup,valmethod,K,standardize,model_type,input_type,scoring);
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,1);
%     save(save_indigo(indigo_summary, 1))
%     [stats,averages,overview] = analyze(indigo_summary,i, indigo_data, data_files,2);
%     save(save_indigo(indigo_summary, 2))
% end