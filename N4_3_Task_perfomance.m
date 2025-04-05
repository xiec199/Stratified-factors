clear; clc

% Add necessary paths
addpath('/public/home/xic_fdu/Ongoing_project/PET_01_Pattern_SCZ_Hub/customcolormap/');
addpath('/home1/xic_fdu/xic_analysis/');

% Load relevant data
load('/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr/A5_1_1_2_Network_FC.mat');
fc_sub = TableS9_2.BL_FC_tale.subject;
fc_mat = TableS9_2.BL_FC_tale.oppo_pp;

load('N3_1_1_make_Ex_In_factor_FC.mat');
[BL_symptom, FU2_symptom, FU3_symptom] = xic_imagen_symptoms;

% Load and process FC data
BL_FC_ex = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.BL_ex_pp_only_FC_sum');
BL_FC_in = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.BL_in_nn_only_FC_sum');
BL_FC_sub = TableS3_1_1_make_factor_fc.BL_subject_match;
BL_FC_ex_sum = mean(BL_FC_ex, 2);
BL_FC_in_sum = mean(BL_FC_in, 2);

FU2_FC_ex = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU2_ex_pp_only_FC_sum');
FU2_FC_in = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU2_in_nn_only_FC_sum');
FU2_FC_sub = TableS3_1_1_make_factor_fc.FU2_subject_match;
FU2_FC_ex_sum = mean(FU2_FC_ex, 2);
FU2_FC_in_sum = mean(FU2_FC_in, 2);

FU3_FC_ex = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU3_ex_pp_only_FC_sum');
FU3_FC_in = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU3_in_nn_only_FC_sum');
FU3_FC_sub = TableS3_1_1_make_factor_fc.FU3_subject_match;
FU3_FC_ex_sum = mean(FU3_FC_ex, 2);
FU3_FC_in_sum = mean(FU3_FC_in, 2);

%% MID performance
load('/share/home1/ISTBI_data/IMAGEN/IMAGEN_behavior/BL_MID_RT_ACC.mat');
[MID_sub, ind1, ind2] = xic_intersect(ACC_id, BL_FC_sub, fc_sub);

BL_FC_MID = BL_FC_ex_sum(ind2);
BL_ACC_MID = mean(ACC(ind1, :), 2);
BL_ACC_MID_all = ACC(ind1, :);

% Remove outliers
[BL_ACC_MID, ind] = rmoutliers(BL_ACC_MID);
BL_FC_MID = BL_FC_MID(ind == 0);
RT = RT(ind == 0);
BL_MID_SUB = MID_sub(ind == 0);

% Correlation analysis
[r, p, t] = xic_corr(BL_FC_MID, BL_ACC_MID);
[r, p, t] = xic_corr(BL_FC_MID, RT);

[sub_mid, ind1, ind2] = intersect(fc_sub, BL_MID_SUB);
[r, p] = partialcorr(BL_FC_MID(ind2), BL_ACC_MID(ind2), fc_mat(ind1));
xic_r2t(r, length(ind2));

% Regress FC data
[~, ~, bl_fc_residua1] = regress(BL_FC_MID(ind2), [ones(length(ind1), 1), fc_mat(ind1)]);
mid_acc = BL_ACC_MID(ind2);

%% Symptoms (External Behavior)
load('/share/home1/ISTBI_data/IMAGEN/IMAGEN_behavior/BL_MID_RT_ACC.mat');
[MID_sub, ind1, ind2] = intersect(ACC_id, BL_symptom.BL_subject);

BL_sym = BL_symptom.BL_Exter_beha_residual(ind2);
BL_ACC_MID = mean(ACC(ind1, :), 2);
BL_ACC_MID_all = ACC(ind1, :);

% Remove outliers
[BL_ACC_MID, ind] = rmoutliers(BL_ACC_MID);
BL_sym = BL_sym(ind == 0);
RT = RT(ind == 0);
BL_MID_SUB = MID_sub(ind == 0);

[r, p, t] = xic_corr(BL_sym, BL_ACC_MID);

%% SST (Stop Signal Task)
load('/share/home1/ISTBI_data/IMAGEN/IMAGEN_behavior/SSRT_0918_BL.mat');
ssrt_nan = SSRT > 0;
ID = ID(ssrt_nan);
SSRT = SSRT(ssrt_nan);
PerGOsuc = PerGOsuc(ssrt_nan);

[SST_SUB, ind1, ind2] = intersect(ID, BL_FC_sub);

BL_FC_SST = BL_FC_in_sum(ind2);
BL_ACC_SST = PerGOsuc(ind1, :);
[BL_ACC_SST, ind] = rmoutliers(BL_ACC_SST);

% Remove outliers and align data
BL_FC_SST = BL_FC_SST(ind == 0);
SSRT = SSRT(ind == 0);
BL_SST_SUB = SST_SUB(ind == 0);

[r, p, t] = xic_corr(BL_ACC_SST, BL_FC_SST);
[r, p, t] = xic_corr(SSRT, BL_FC_SST);

[sub_sst, ind1, ind2] = intersect(fc_sub, BL_SST_SUB);
[r, p] = partialcorr(BL_FC_SST(ind2), BL_ACC_SST(ind2), fc_mat(ind1));
xic_r2t(r, length(ind2));

% Regress FC data
[~, ~, bl_fc_residua2] = regress(BL_FC_SST(ind2), [ones(length(ind1), 1), fc_mat(ind1)]);
sst_acc = BL_ACC_SST(ind2);

%% Symptoms (Internal Behavior)
load('/share/home1/ISTBI_data/IMAGEN/IMAGEN_behavior/SSRT_0918_BL.mat');
ssrt_nan = SSRT > 0;
ID = ID(ssrt_nan);
SSRT = SSRT(ssrt_nan);
PerGOsuc = PerGOsuc(ssrt_nan);

[SST_SUB, ind1, ind2] = intersect(ID, BL_symptom.BL_subject);

BL_sym = BL_symptom.BL_Inter_beha_residual(ind2);
BL_ACC_SST = PerGOsuc(ind1, :);

% Remove outliers
[BL_ACC_SST, ind] = rmoutliers(BL_ACC_SST);
BL_sym = BL_sym(ind == 0);
SSRT = SSRT(ind == 0);
BL_SST_SUB = SST_SUB(ind == 0);

[r, p, t] = xic_corr(BL_ACC_SST, BL_sym);

%% Combine SST and MID results
[~, ind1, ind2] = intersect(sub_sst, sub_mid);

[r, p] = corr(sst_acc(ind1), bl_fc_residua2(ind1));
[r, p] = corr(mid_acc(ind2), bl_fc_residua1(ind2));
