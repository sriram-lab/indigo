function F = dataFiles()
%list all the possible data files here - both train and test data
%also all ortholog files
%wrap all names in double quotes

%put everything in a struct
F = struct;
%need to put all your files in one directory
F.indigoData = 'orig_ecoli.xlsx';

F.natureData = {'ecoli_bw25113.xlsx'
               'ecoli_iAi1.xlsx'
               'stlt2.xlsx'
               'st14028.xlsx'
               'pao1.xlsx'
               'pa14.xlsx'};
        
F.natureOrthologs = {'ecoli_iAi1_orthologs.xlsx'
                    'ecoli_stlt2_orthologs.xlsx'
                    'ecoli_st14028_orthologs.xlsx'
                    'ecoli_pao1_orthologs.xlsx'
                    'ecoli_pa14_orthologs.xlsx'};
                
F.natureEcoli = {'ecoli_bw25113.xlsx'
                 'ecoli_iAi1.xlsx'};
                
F.yehData = 'ecoli_yeh.xlsx';
F.asmData = 'asm_data.xlsx';

F.ecoliBlissData = 'ecoli_bliss.xlsx';
F.ecoliLoeweData = 'ecoli_loewe.xlsx';

end