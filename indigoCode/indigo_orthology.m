%indigo_orthology has two inputs
%[predicted scores, sigma delta scores] = indigo_ortholgy
function [deviations, teststaphdat1] = indigo_orthology(phenotype_labels,ecoli_staph_orth,sigma_delta_input,indigo_model)
% deviations = indigo_orthology(phenotype_labels,ecoli_staph_orth,sigma_delta_input,indigo_model);
% % step 1 : get ecoli model 
% % step 2 - get orthologs
% % step3 - get predicted variable interactions between two species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  modify model inputs based on orthology
% [sx, spos] = sort(indigo_model.importance,'descend'); % sort genes based on their importance
nic_row = [phenotype_labels;phenotype_labels]; % row identifiers in the random forest model
    topgenes = nic_row;
    ix = ismember(topgenes,ecoli_staph_orth);  % find the orthologs
    nonorthtop = nic_row(~ix); ix1 = ~ix;   % get the non-orthologs
    teststaphdat = sigma_delta_input; 
    teststaphdat(ix1,:) = 0; % modfiy state of non orthologs
    % set the sigma scores to be 2 or 0
    ix2 = ismember(phenotype_labels,nonorthtop); %sum(ix2)
    ix2 = find(ix2);
    teststaphdat1 = sigma_delta_input;
    teststaphdat1(ix2,:) = 2;
    teststaphdat1(ix2+length(phenotype_labels),:) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  predict interactions that are variable 
    testpredictions_staphchem_ecolixns2 =  predict(indigo_model,teststaphdat1');
    testpredictions_staphchem_ecolixns20 =  predict(indigo_model,sigma_delta_input');
%     disp(testpredictions_staphchem_ecolixns20);
    testpredictions_staphchem_ecolixns21 =  predict(indigo_model,teststaphdat');
%     deviations = testpredictions_staphchem_ecolixns20(:)  - testpredictions_staphchem_ecolixns21(:); % output the deviations for the input drugs
    deviations = testpredictions_staphchem_ecolixns20(:)  - testpredictions_staphchem_ecolixns2(:); % output the deviations for the input drugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end