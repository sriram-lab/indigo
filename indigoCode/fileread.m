
files = cellstr(ls('indigoData'));
%data = stl_data.xlsx
%ortholog = stl_orthologs.xlsx
pattern = strcat(erase(file,'.xlsx'),'_orthologs.xlsx');
contains(files,pattern)

%File naming scheme
%user inputs all files that they want to be in training data
%nature data can all have same naming scheme
%program checks to see if there are orthologs