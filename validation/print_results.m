function [row_loc, column_loc] = print_results(indigo_summary, indigo_data, parameter)
    test_data = indigo_summary.test_data;
    model_type = indigo_summary.model_type;
    valmethod = indigo_summary.valmethod;
    scoring = indigo_summary.scoring;
    standardize = indigo_summary.standardize;
    result_idx = indigo_data.Number_Ref(strcmp(indigo_data.Filename, test_data));
    row_loc = result_idx + 3;   % based on dataset

    %% For R
    if strcmp(parameter,'R')
        if strcmp(model_type, 'original_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'E';
            else
                column_loc = 'G';
            end
        elseif strcmp(model_type, 'original_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'I';
            else
                column_loc = 'K';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'M';
            else
                column_loc = 'O';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'Q';
            else
                column_loc = 'S';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'U';
            else
                column_loc = 'W';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'Y';
            else
                column_loc = 'AA';
            end    
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AC';
            else
                column_loc = 'AE';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AG';
            else
                column_loc = 'AI';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'AK';
            else
                column_loc = 'AM';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'AO';
            else
                column_loc = 'AQ';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'AS';
            else
                column_loc = 'AU';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'AW';
            else
                column_loc = 'AY';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'BA';
            else
                column_loc = 'BC';   
            end   
       elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'BE';
            else
                column_loc = 'BG';   
            end     
        end  

    %% For P
    elseif strcmp(parameter, 'P')    
        if strcmp(model_type, 'original_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'F';
            else
                column_loc = 'H';
            end
        elseif strcmp(model_type, 'original_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'J';
            else
                column_loc = 'L';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'N';
            else
                column_loc = 'P';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'R';
            else
                column_loc = 'T';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'V';
            else
                column_loc = 'X';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'Z';
            else
                column_loc = 'AB';
            end    

        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AD';
            else
                column_loc = 'AF';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AH';
            else
                column_loc = 'AJ';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'AL';
            else
                column_loc = 'AN';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'AP';
            else
                column_loc = 'AR';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'AT';
            else
                column_loc = 'AV';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'AX';
            else
                column_loc = 'AZ';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'BB';
            else
                column_loc = 'BD';   
            end   
       elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'BF';
            else
                column_loc = 'BH';   
            end     
        end  
    % everything else
    else
        if strcmp(model_type, 'original_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'E';
            else
                column_loc = 'F';
            end
        elseif strcmp(model_type, 'original_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'G';
            else
                column_loc = 'H';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'I';
            else
                column_loc = 'J';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'K';
            else
                column_loc = 'L';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'M';
            else
                column_loc = 'N';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'O';
            else
                column_loc = 'P';
            end    
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'Q';
            else
                column_loc = 'R';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'S';
            else
                column_loc = 'T';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'U';
            else
                column_loc = 'V';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'W';
            else
                column_loc = 'X';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'Y';
            else
                column_loc = 'Z';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'AA';
            else
                column_loc = 'AB';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AC';
            else
                column_loc = 'AD';   
            end   
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AE';
            else
                column_loc = 'AF';   
            end     
        end  
    end  
end
