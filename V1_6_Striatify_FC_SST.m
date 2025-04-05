%% Analysis of Externalizing and Internalizing Factors in STRATIFY Dataset

clear; clc;
load('N2_1_5_make_Ex_In_factor_eachFC.mat')
mid = load('N5_3_2_1_Striatify_FC_MID_v2.mat');
sst = load('N5_3_2_1_Striatify_FC_SST_v2.mat');

%% Extract subject IDs and functional connectivity data
sst_sub = sst.Striatify_FC_SST.subID;
sst_dat_ex = sst.Striatify_FC_SST.ex_pp;
sst_dat_in = sst.Striatify_FC_SST.in_nn;
mid_sub = mid.Striatify_FC_MID.subID;
mid_dat_ex = mid.Striatify_FC_MID.ex_pp;
mid_dat_in = mid.Striatify_FC_MID.in_nn;

%% Load demographic and clinical data
data = readtable('STRATIFY_demographics.csv');
adhd = readtable('STRATIFY-STRATIFY_ADHD-BASIC_DIGEST.csv');
adhd_dat = table2array(adhd(:, 29:end));
adhd_dat_age = table2array(adhd(:, 5))/365;
adhd_dat_sum = nansum(adhd_dat, 2);
adhd_sub = adhd.UserCode;

%% Prepare covariates for regression
% Process site information
sub_site = (string(data.recruitmentSite));
sub_site(contains(sub_site, 'BERLIN')) = 1;
sub_site(contains(sub_site, 'LONDON')) = 2;
sub_site(contains(sub_site, 'SOUTHAMPTON')) = 3;

% Create covariate matrix
sub_gender = data.sex;
sub_cova = [ones(length(sub_gender), 1), contains(sub_gender, 'F'), dummyvar(str2double(sub_site))];
sub_cova(:, end) = [];
sub_id = data.PSC2;

% Process diagnostic group information
sub_type = string(data.patientGroup);
for i = 1:length(sub_type)
    if (sub_type(i)) == '';
        sub_type(i) = 'NaN';
    end
end

% Match ADHD data with subject IDs
[sub_idAll, ind1, ind2] = intersect(adhd_sub, sub_id);
adhd_dat_sum = adhd_dat_sum(ind1);
sub_cova = sub_cova(ind2, :);
sub_type = sub_type(ind2);
sub_age = adhd_dat_age(ind1);

%% Extract SST functional connectivity data
[~, ind1, ind2] = intersect(sst_sub, sub_idAll);
sub_fc = sst_dat_in(ind1);

%% Uncomment for MID analysis instead of SST
% [~, ind1, ind2] = intersect(mid_sub, sub_idAll);
% sub_fc = mid_dat_in(ind1);

%% Compare FC differences between diagnostic groups
sub_fc_type = sub_type(ind2);
sub_age_fc = sub_age(ind2);
sub_cova_fc = sub_cova(ind2, :);
sub_sex = sub_cova(:, 2);

% Regress out covariates
[~, ~, sub_fcR] = regress(sub_fc, sub_cova_fc);

% Identify diagnostic groups
control_id = contains(sub_fc_type, 'Control') + contains(sub_fc_type, 'NaN');
AN_id = contains(sub_fc_type, 'AN');
AUD_id = contains(sub_fc_type, 'AUD');
BN_id = contains(sub_fc_type, 'BN');
MDD_id = contains(sub_fc_type, 'MDD');

% Extract FC values for each group
sub_fc_control = sub_fcR(control_id == 1);
sub_fc_case_in = sub_fcR((AN_id + BN_id + MDD_id) == 1);
sub_fc_AN = sub_fcR(AN_id == 1);
sub_fc_AUD = sub_fcR(AUD_id == 1);
sub_fc_BN = sub_fcR(BN_id == 1);
sub_fc_MDD = sub_fcR(MDD_id == 1);

%% Perform statistical comparisons
% AUD vs Control
[tvalue_AUD, p_AUD, dvalue_AUD] = xic_ttest2(sub_fc_AUD, sub_fc_control);

% All internalizing disorders vs Control
[tvalue_IN, p_IN, dvalue_IN] = xic_ttest2(sub_fc_case_in, sub_fc_control);

% AN vs Control
[tvalue_AN, p_AN, dvalue_AN] = xic_ttest2(sub_fc_AN, sub_fc_control);

% BN vs Control
[tvalue_BN, p_BN, dvalue_BN] = xic_ttest2(sub_fc_BN, sub_fc_control);

% MDD vs Control
[tvalue_MDD, p_MDD, dvalue_MDD] = xic_ttest2(sub_fc_MDD, sub_fc_control);

%% Save results
results = struct();
results.AUD = struct('tvalue', tvalue_AUD, 'pvalue', p_AUD, 'dvalue', dvalue_AUD);
results.Internalizing = struct('tvalue', tvalue_IN, 'pvalue', p_IN, 'dvalue', dvalue_IN);
results.AN = struct('tvalue', tvalue_AN, 'pvalue', p_AN, 'dvalue', dvalue_AN);
results.BN = struct('tvalue', tvalue_BN, 'pvalue', p_BN, 'dvalue', dvalue_BN);
results.MDD = struct('tvalue', tvalue_MDD, 'pvalue', p_MDD, 'dvalue', dvalue_MDD);

save('STRATIFY_FC_Analysis.mat', 'results');