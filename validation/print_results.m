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
                column_loc = 'D';
            else
                column_loc = 'F';
            end
        elseif strcmp(model_type, 'original_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'H';
            else
                column_loc = 'J';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'L';
            else
                column_loc = 'N';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'P';
            else
                column_loc = 'R';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'T';
            else
                column_loc = 'V';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'X';
            else
                column_loc = 'Z';
            end    
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AB';
            else
                column_loc = 'AD';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AF';
            else
                column_loc = 'AH';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'AJ';
            else
                column_loc = 'AL';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'AN';
            else
                column_loc = 'AP';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'AR';
            else
                column_loc = 'AT';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'AV';
            else
                column_loc = 'AX';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AZ';
            else
                column_loc = 'BB';   
            end   
       elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'BD';
            else
                column_loc = 'BF';   
            end     
        end  

    %% For P
    elseif strcmp(parameter, 'P')    
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
    % everything else
    else
        if strcmp(model_type, 'original_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'D';
            else
                column_loc = 'E';
            end
        elseif strcmp(model_type, 'original_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'F';
            else
                column_loc = 'G';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'H';
            else
                column_loc = 'I';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'J';
            else
                column_loc = 'K';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'L';
            else
                column_loc = 'M';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'N';
            else
                column_loc = 'O';
            end    
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'P';
            else
                column_loc = 'Q';
            end
        elseif strcmp(model_type, 'ecoli_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'R';
            else
                column_loc = 'S';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'T';
            else
                column_loc = 'U';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, '')
            if strcmp(standardize, '')
                column_loc = 'V';
            else
                column_loc = 'W';
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'X';
            else
                column_loc = 'Y';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'bliss')
            if strcmp(standardize, '')
                column_loc = 'Z';
            else
                column_loc = 'AA';   
            end
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'independent') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AB';
            else
                column_loc = 'AC';   
            end   
        elseif strcmp(model_type, 'mtb_model') && strcmp(valmethod, 'Kfold') && strcmp(scoring, 'loewe')
            if strcmp(standardize, '')
                column_loc = 'AD';
            else
                column_loc = 'AE';   
            end     
        end  
    end  
end
