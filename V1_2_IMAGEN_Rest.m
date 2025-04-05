%% Analyze externalizing and internalizing factors from resting-state FC data

clear; clc;
addpath('path/to/External_Internal_Factor');
load('N3_1_1_make_Ex_In_factor_FC.mat')

%% Extract FC data for externalizing and internalizing factors
rest_exter_fc = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU2_Rest_ex_pp_only_FC_sum');
rest_inter_fc = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU2_Rest_in_nn_only_FC_sum');

%% Calculate mean FC values
rest_exter_fc_mean = mean(rest_exter_fc, 2);
rest_inter_fc_mean = mean(rest_inter_fc, 2);
rest_sub = TableS3_1_1_make_factor_fc.FU2_rest_sub;

%% Load symptom data
[BL_symptom, FU2_symptom, FU3_symptom] = xic_imagen_symptoms;

%% Match subjects and calculate correlations
[~, ind1, ind2] = intersect(rest_sub, FU2_symptom.FU2_subject);

% Correlate externalizing FC with externalizing behavior
[exter_r(1), exter_r(2), exter_r(3)] = xic_corr_one(rest_exter_fc_mean(ind1), ...
    FU2_symptom.FU2_Exter_beha_residual(ind2));

% Correlate internalizing FC with internalizing behavior
[inter_r(1), inter_r(2), inter_r(3)] = xic_corr_one(rest_inter_fc_mean(ind1), ...
    FU2_symptom.FU2_Inter_beha_residual(ind2));

%% Save results
save('IMAGEN_Rest_FC_Correlations.mat', 'exter_r', 'inter_r')