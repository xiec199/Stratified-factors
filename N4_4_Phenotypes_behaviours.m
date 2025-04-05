clear; clc

% Load data
path = '/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr';
load(fullfile(path, 'A5_1_2_Prepare_phnotype_data.mat'));
load('N3_1_1_make_Ex_In_factor_FC.mat');

% Extract relevant data
BL_FC_ex = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.BL_ex_pp_only_FC_sum');
BL_FC_in = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.BL_in_nn_only_FC_sum');
BL_FC_sub = TableS3_1_1_make_factor_fc.BL_subject_match;
BL_FC_ex_sum = mean(BL_FC_ex, 2);
BL_FC_in_sum = mean(BL_FC_in, 2);

% Load covariate data and match subjects
load('/home1/xic_fdu/Datasets/Match_control_subject_All things scale/Cova_subject_jia.mat');
[sub_cova, ind1, ind2] = intersect(cova_subject, BL_FC_sub);

% Perform regression
[~, ~, bl_fc_ex] = regress(BL_FC_ex_sum(ind2), cova_data(ind1, 1:9));
[~, ~, bl_fc_in] = regress(BL_FC_in_sum(ind2), cova_data(ind1, 1:9));

% Regression for task data
for i = 1:4
    [~, ~, bl_fc_ex_task(:, i)] = regress(BL_FC_ex(ind2, i), cova_data(ind1, 1:9));
    [~, ~, bl_fc_in_task(:, i)] = regress(BL_FC_in(ind2, i), cova_data(ind1, 1:9));
end

% Phenotype correlation
Value_Nogender = table;  
Value_Nogender.name = TableS10.cova_name;

for i = 1:52
    sub = TableS10.cova_subject{i, 1}; 
    dat = TableS10.cova_data{i, 1};
    
    % Match subjects
    [~, ind1, ind2] = intersect(sub_cova, sub);
    Value_Nogender.Sample_N(i, 1) = length(ind1);
    
    % Correlation with FC data
    [Value_Nogender.bl_exter_R(i, 1), Value_Nogender.bl_exter_P(i, 1), Value_Nogender.bl_exter_T(i, 1)] = xic_corr(dat(ind2), bl_fc_ex(ind1));
    [Value_Nogender.bl_inter_R(i, 1), Value_Nogender.bl_inter_P(i, 1), Value_Nogender.bl_inter_T(i, 1)] = xic_corr(dat(ind2), bl_fc_in(ind1));          
end 

% Remove specific rows
Value_Nogender([11 20 23], :) = [];

% FDR correction
Value_Nogender.bl_exter_Pfdr = FDR(Value_Nogender.bl_exter_P);
Value_Nogender.bl_inter_Pfdr = FDR(Value_Nogender.bl_inter_P);

% Correlation between FC and phenotypes
[ex_in_r, ex_in_p] = corr(bl_fc_ex, bl_fc_in);
[ph_ex_in_r, ph_ex_in_p] = corr(Value_Nogender.bl_exter_T, Value_Nogender.bl_inter_T);

% Plot results
figure; 
subplot(1, 2, 1);
plot(Value_Nogender.bl_exter_T, Value_Nogender.bl_inter_T, '.'); 
lsline;
title('Phenotype vs. Behavior Correlation');

subplot(1, 2, 2);
plot(bl_fc_ex, bl_fc_in, '.'); 
lsline;
title('FC Ex vs. In Correlation');

% Save results
save('N4_3_1_Phenotypes_behaviours.mat', 'Value_Nogender');
