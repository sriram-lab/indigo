function [normData,col,row] = tb_preprocessing(averagedtbdata,mtb_expression_database_col_ids,mtb_expression_database_row_ids)
% Variables that specify which region of accession numbers and treatments to look at...
firstGene = 1527;
lastGene = 5494;
firstTreatment = 1;
finalTreatment = 248;

% Create batch matrix for combat normalization
batch = (1:248)';

for i = 1:248
    if batch(i) <= 75
        batch(i) = 1;
    elseif batch(i) <= 177
        batch(i) = 2;
    elseif batch(i) <= 210
        batch(i) = 3;
    elseif batch(i) <= 214
        batch(i) = 4;
    elseif batch(i) <= 234
        batch(i) = 5;
    elseif batch(i) <= 237
        batch(i) = 6;
    elseif batch(i) <= 240
        batch(i) = 7;
    else
        batch(i) = 8;
    end
end

data = averagedtbdata(firstGene:lastGene,firstTreatment:finalTreatment);
data(isnan(data)) = 0;

% If you need to impute, uncomment the commands below %
% imputeData = knnimpute(transformData);
% transformData = imputeData;
% Log transform and combine data (will directly convert in input file later)

%{
regData = data(:,1:255);
logTransform = (data(:, 256:266))/log10(2);
transformData = horzcat(regData, logTransform);

%}
normData= ComBat(data, batch, [], 1);

% If you need to normalize data using quantile norm, uncomment the command below %
% normData = quantilenorm(transformData);

col = mtb_expression_database_col_ids(firstTreatment:finalTreatment,1);
row = mtb_expression_database_row_ids(firstGene:lastGene,1);
end

