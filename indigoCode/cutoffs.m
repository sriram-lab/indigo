function [synergy,antagonism] = cutoffs(filename)
%define cutoffs

%indigo, according to article: https://www.embopress.org/doi/full/10.15252/msb.20156777
%loewe
%distribution is decent, but cutoff for antagonism makes things weird
%when doing zscore
%But makes sense to do zscore as scale is off from other scoressets
%Majority of interactions are antagonistic (see Figure 1)
files = dataFiles();

if strcmp(filename, files.indigoData)
    synergy = -0.25;
    antagonism = 1;
    
elseif ismember(filename,files.natureData)
%nature, according to article: https://www.nature.com/articles/s41586-018-0278-9
%bliss
    synergy = -0.1;
    antagonism = 0.1;
    
elseif strcmp(filename, files.yehData)
%Yeh kishony, according to article https://www.nature.com/articles/ng1755#Fig3
%See figure 3
%bliss I think
    synergy = -0.25;
    antagonism = 0.25;

elseif strcmp(filename, files.asmData)
%asm article, according to article https://aac.asm.org/content/58/8/4573
%loewe
%this one was flipped, synergy > 0.25
%antagonsim < -0.25
%So I multiplied everything by -1 to be consistent with other scoring
%schemes
    synergy = -0.25;
    antagonism = 0.25;
    
elseif strcmp(filename,files.ecoliBlissData) || strcmp(filename,files.ecoliLoeweData)
%same as other Kishony paper, according to article https://www.nature.com/articles/s41564-018-0252-1
    synergy = -0.25;
    antagonism = 0.25;
    
elseif strcmp(filename,files.tbData)
    synergy = -0.1;
    antagonism = 0.1;
    
elseif strcmp(filename,files.saureusData)
    synergy = -0.25;
    antagonism = 1;
end
end
