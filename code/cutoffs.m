function [synergy,antagonism] = cutoffs(filename)
%define cutoffs

%indigo, according to article: https://www.embopress.org/doi/full/10.15252/msb.20156777
%loewe
%distribution is decent, but cutoff for antagonism makes things weird
%when doing zscore
%But makes sense to do zscore as scale is off from other scoressets
%Majority of interactions are antagonistic (see Figure 1)
%update so that it is actual file names

natureData = {'ecoli_bw25113.xlsx'
              'ecoli_iAi1.xlsx'
              'stlt2.xlsx'
              'st14028.xlsx'
              'pao1.xlsx'
              'pa14.xlsx'};

if strcmp(filename, 'orig_ecoli.xlsx')
    synergy = -0.25;
    antagonism = 1;
    
elseif ismember(filename,natureData)
%nature, according to article: https://www.nature.com/articles/s41586-018-0278-9
%bliss
    synergy = -0.1;
    antagonism = 0.1;
    
elseif strcmp(filename, 'ecoli_yeh.xlsx')
%Yeh kishony, according to article https://www.nature.com/articles/ng1755#Fig3
%See figure 3
%bliss I think
    synergy = -0.25;
    antagonism = 0.25;

elseif strcmp(filename, 'asm_data.xlsx')
%asm article, according to article https://aac.asm.org/content/58/8/4573
%loewe
%this one was flipped, synergy > 0.25
%antagonsim < -0.25
%So I multiplied everything by -1 to be consistent with other scoring
%schemes
    synergy = -0.25;
    antagonism = 0.25;
    
elseif strcmp(filename,'ecoli_bliss.xlsx') || strcmp(filename,'ecoli_loewe.xlsx')
%same as other Kishony paper, according to article https://www.nature.com/articles/s41564-018-0252-1
    synergy = -0.25;
    antagonism = 0.25;
    
elseif strcmp(filename,'tb.xlsx')
    %FIC
    synergy = -0.1;
    antagonism = 0.1;
    
elseif strcmp(filename,'staph.xlsx')
    %Loewe
    synergy = -0.25;
    antagonism = 1;
    
elseif strcmp(filename,'acinetobacter.xlsx')
    %Loewe-additivity model and Fractional Inhibitory Concentration (FIC)
    %log-FIC, https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006677#
    synergy = -0.2;
    antagonism = 0.2;
    
elseif strcmp(filename,'tb_pairwise.xlsx')
    %Loewe-additivity model and Fractional Inhibitory Concentration (FIC)
    %IC70 (70% inhibitory concentration)
    %??=?log2(FIC)
    %This score is 0, <0 or >0 for additive, synergistic or antagonistic pairs, respectively.
    %https://www.nature.com/articles/s41598-019-48410-y#Sec12
    synergy = 0;
    antagonism = 0;
end
end
