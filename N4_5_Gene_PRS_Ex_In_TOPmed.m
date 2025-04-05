clear; clc;

% Load data
load('N3_1_1_make_Ex_In_factor_FC.mat');
BL_FC_ex = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.BL_ex_pp_only_FC_sum');
BL_FC_in = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.BL_in_nn_only_FC_sum');
BL_FC_sub = TableS3_1_1_make_factor_fc.BL_subject_match;
FU_FC_data = {...
    'FU2', cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU2_ex_pp_only_FC_sum'), ...
    'FU2_in', cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU2_in_nn_only_FC_sum'), ...
    'FU3', cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU3_ex_pp_only_FC_sum'), ...
    'FU3_in', cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU3_in_nn_only_FC_sum')};

% Compute FC means
BL_FC_ex_sum = mean(BL_FC_ex, 2);
BL_FC_in_sum = mean(BL_FC_in, 2);
FU_FC_ex_sum = cellfun(@(x) mean(x, 2), FU_FC_data(2:4:8), 'UniformOutput', false);
FU_FC_in_sum = cellfun(@(x) mean(x, 2), FU_FC_data(3:4:8), 'UniformOutput', false);

%% Regress covariates
load('/home1/xic_fdu/Datasets/Match_control_subject_All things scale/Cova_subject_jia.mat');
[sub_cova_bl, ind1, ind2] = intersect(cova_subject, BL_FC_sub);
[~, ~, bl_fc_ex] = regress(BL_FC_ex_sum(ind2), cova_data(ind1, 1:9));
[~, ~, bl_fc_in] = regress(BL_FC_in_sum(ind2), cova_data(ind1, 1:9));

% Perform regressions for FU2 and FU3
FU_data = {'fu2', FU_FC_data{2}, FU_FC_data{4}; 'fu3', FU_FC_data{6}, FU_FC_data{8}};
for i = 1:size(FU_data, 1)
    [sub_cova_fu, ind1, ind2] = intersect(cova_subject, eval([FU_data{i, 1}, '_FC_sub']));
    [~, ~, eval([FU_data{i, 1}, '_fc_ex'])] = regress(eval([FU_data{i, 1}, '_FC_ex_sum'])(ind2), cova_data(ind1, 1:9));
    [~, ~, eval([FU_data{i, 1}, '_fc_in'])] = regress(eval([FU_data{i, 1}, '_FC_in_sum'])(ind2), cova_data(ind1, 1:9));
end

% External factors
load('/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr/A5_1_1_2_Network_FC.mat');
fc_sub = TableS9_2.BL_FC_tale.subject;
fc_mat = TableS9_2.BL_FC_tale.oppo_pp;

% Read ADHD PRS data
prs_path = '/home1/xic_fdu/Datasets/Gwas_Chao/Exteranlising_PRS/ADHD/';
[adhd_prs, prs_sub] = prs_read(prs_path);

% Perform external factor regressions
[sub, ind1, ind2] = intersect(fc_sub, sub_cova_bl);
[~, ~, exter_residual] = regress(bl_fc_ex(ind1), [ones(length(sub_cova_bl), 1), fc_mat(ind2)]);
[~, ~, inter_residual] = regress(bl_fc_in(ind1), [ones(length(sub_cova_bl), 1), fc_mat(ind2)]);

% Correlations with ADHD PRS
[~, r, p, t] = arrayfun(@(x) xic_corr_one(adhd_prs(ind1), eval([x, '_residual'])(ind2)), {'exter', 'inter'}, 'UniformOutput', false);

% Internal factors (PRS)
folder = dir('/home1/xic_fdu/Datasets/Gwas_Chao/Previous_PRS/PRS_Mats/*.mat');
[PRS_inter_bl, PRS_dat] = xic_PRS_fc(folder, sub_cova_bl, bl_fc_in);

% Correlations for internal factors
[~, ind1, ind2] = intersect(sub, PRS_dat.sub{13});
inter_resi = inter_residual(ind1);
inter_fc = bl_fc_in(ind1);
exter_fc = bl_fc_ex(ind1);
prs_dat = PRS_dat.dat{13}(ind2);

% Calculate correlations
[r, p, t] = arrayfun(@(x) xic_corr_one(eval([x, '_fc']), prs_dat), {'inter', 'exter'}, 'UniformOutput', false);

