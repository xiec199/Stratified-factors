%% Analysis of Depression Data Using Externalizing and Internalizing Factors

clear; clc;
load('N2_1_5_make_Ex_In_factor_eachFC.mat')
Xinan = load('Xinan_FC_Map_268.mat');

%% Load subject information
case_info = readtable('Patient_info.csv');
con_info = readtable('Con_info.csv');

%% Prepare functional connectivity data
xinanFC = reshape(Xinan.fc_pnc, [], size(Xinan.fc_pnc, 3));
xinanSub = Xinan.fc_sub_name;
exter = exter_pospos_type;
inter = inter_negneg_type;

%% Calculate mean FC for externalizing and internalizing factors
for i = 1:4
    exter_fc_rest(:, i) = nanmean(xinanFC(exter{i}, :))';
    inter_fc_rest(:, i) = nanmean(xinanFC(inter{i}, :))';
end
inter_fc_mean = mean(inter_fc_rest, 2);

%% Compare depression patients vs controls
Control = contains(xinanSub, 'sub-con');
[h, p, c, t] = ttest2(inter_fc_mean(Control == 0, :), inter_fc_mean(Control == 1, :));
Tvalue = t.tstat;
Pvalue = p/2;
Cohen_D = cohen_t2(t.tstat, sum(Control == 0), sum(Control == 1));

%% Calculate mean age across all participants
mean_age = mean([case_info.Age; con_info.Age]);

%% Save results
results = struct('Tvalue', Tvalue, 'Pvalue', Pvalue, 'Cohen_D', Cohen_D, 'mean_age', mean_age);
save('Depression_FC_Analysis.mat', 'results');