%% Neural Factor Analysis for ABCD Dataset
% Load and analyze externalizing and internalizing factors from functional connectivity data

clear; clc;

%% Load necessary data
load('N2_1_5_make_Ex_In_factor_eachFC.mat')
load('ABCD_MID_CONN_v2.mat')
load('ABCD_SST_CONN_v3.mat');
addpath('path/to/task_template')
addpath('path/to/external_internal_factor')

%% Prepare data structures
exter = exter_pospos_type;
inter = inter_negneg_type;
Sub_SSTs = strrep(Sub_SST, 'sub-NDAR', 'NDAR_');

%% Extract functional connectivity data
% Reshape and organize SST and MID task data
sst_Data{1} = reshape(ABCD_SST_StopSuccess, [], size(ABCD_SST_StopSuccess, 3));
sst_Data{2} = reshape(ABCD_SST_StopFailure, [], size(ABCD_SST_StopFailure, 3));
mid_Data{1} = reshape(ABCD_MID_feed, [], size(ABCD_MID_feed, 3));
mid_Data{2} = reshape(ABCD_MID_anti, [], size(ABCD_MID_anti, 3));

%% Calculate mean FC for externalizing and internalizing factors
% Externalizing factors
exter_fc_sst(:,1) = mean(sst_Data{1,1}(exter{1},:))';
exter_fc_sst(:,2) = mean(sst_Data{1,2}(exter{2},:))';
exter_fc_mid(:,1) = mean(mid_Data{1,1}(exter{3},:))';
exter_fc_mid(:,2) = mean(mid_Data{1,2}(exter{4},:))';

% Internalizing factors 
inter_fc_sst(:,1) = mean(sst_Data{1,1}(inter{1},:))';
inter_fc_sst(:,2) = mean(sst_Data{1,2}(inter{2},:))';
inter_fc_mid(:,1) = mean(mid_Data{1,1}(inter{3},:))';
inter_fc_mid(:,2) = mean(mid_Data{1,2}(inter{4},:))';

%% Analyze baseline behavioral data 
load('N5_3_2_ABCD_Tables_v2.mat')
dat_table = abcd_dat_resi{5, 1};
[bl_sub, ind1, ind2, ind3] = xic_intersect(dat_table.sub, Sub_MID, Sub_SSTs);

% Extract behavioral and FC data
cbcl_bl = table2array(dat_table(ind1, 2:end));
exter_bl_mid = mean(exter_fc_mid(ind2, 2:end), 2);
exter_bl_sst = mean(exter_fc_sst(ind3, 2:end), 2);
inter_bl_mid = mean(inter_fc_mid(ind2, 2:end), 2);
inter_bl_sst = mean(inter_fc_sst(ind3, 2:end), 2);

% Correlate behavior with externalizing factors
[r1, p1, t1] = xic_corr_one(cbcl_bl, exter_bl_mid);
[r2, p2, t2] = xic_corr_one(cbcl_bl, exter_bl_sst);
[r3, p3, t3] = xic_corr_one(cbcl_bl, exter_bl_mid + exter_bl_sst);

% Correlate behavior with internalizing factors
[r4, p4, t4] = xic_corr_one(cbcl_bl, inter_bl_mid);
[r5, p5, t5] = xic_corr_one(cbcl_bl, inter_bl_sst);
[r6, p6, t6] = xic_corr_one(cbcl_bl, inter_bl_mid + inter_bl_sst);

% Create externalizing results table
R_exter = table;
R_exter.Name = dat_table.Properties.VariableNames(2:end)';
R_exter.R_mid_exter = r1;
R_exter.P_mid_exter = p1;
R_exter.T_mid_exter = t1;
R_exter.R_sst_exter = r2;
R_exter.P_sst_exter = p2;
R_exter.T_sst_exter = t2;
R_exter.R_sstMID_exter = r3;
R_exter.P_sstMID_exter = p3;
R_exter.T_sstMID_exter = t3;

% Create internalizing results table
R_inter = table;
R_inter.Name = dat_table.Properties.VariableNames(2:end)';
R_inter.R_mid_inter = r4;
R_inter.P_mid_inter = p4;
R_inter.T_mid_inter = t4;
R_inter.R_sst_inter = r5;
R_inter.P_sst_inter = p5;
R_inter.T_sst_inter = t5;
R_inter.R_sstMID_inter = r6;
R_inter.P_sstMID_inter = p6;
R_inter.T_sstMID_inter = t6;

%% Analyze follow-up behavioral data
load('N5_3_2_ABCD_Tables_v2_fu2.mat')
dat_table = abcd_dat_resi{5, 1};
[bl_sub, ind1, ind2, ind3] = xic_intersect(dat_table.sub, Sub_MID, Sub_SSTs);

% Extract behavioral and FC data
cbcl_bl = table2array(dat_table(ind1, 2:end));
exter_bl_mid = mean(exter_fc_mid(ind2, 2:end), 2);
exter_bl_sst = mean(exter_fc_sst(ind3, 2:end), 2);
inter_bl_mid = mean(inter_fc_mid(ind2, 2:end), 2);
inter_bl_sst = mean(inter_fc_sst(ind3, 2:end), 2);

% Correlate behavior with externalizing factors
[r1, p1, t1] = xic_corr_one(cbcl_bl, exter_bl_mid);
[r2, p2, t2] = xic_corr_one(cbcl_bl, exter_bl_sst);
[r3, p3, t3] = xic_corr_one(cbcl_bl, exter_bl_mid + exter_bl_sst);

% Correlate behavior with internalizing factors
[r4, p4, t4] = xic_corr_one(cbcl_bl, inter_bl_mid);
[r5, p5, t5] = xic_corr_one(cbcl_bl, inter_bl_sst);
[r6, p6, t6] = xic_corr_one(cbcl_bl, inter_bl_mid + inter_bl_sst);

% Create follow-up externalizing results table
R_exter_fu = table;
R_exter_fu.Name = dat_table.Properties.VariableNames(2:end)';
R_exter_fu.R_mid_exter = r1;
R_exter_fu.P_mid_exter = p1;
R_exter_fu.T_mid_exter = t1;
R_exter_fu.R_sst_exter = r2;
R_exter_fu.P_sst_exter = p2;
R_exter_fu.T_sst_exter = t2;
R_exter_fu.R_sstMID_exter = r3;
R_exter_fu.P_sstMID_exter = p3;
R_exter_fu.T_sstMID_exter = t3;

% Create follow-up internalizing results table
R_inter_fu = table;
R_inter_fu.Name = dat_table.Properties.VariableNames(2:end)';
R_inter_fu.R_mid_inter = r4;
R_inter_fu.P_mid_inter = p4;
R_inter_fu.T_mid_inter = t4;
R_inter_fu.R_sst_inter = r5;
R_inter_fu.P_sst_inter = p5;
R_inter_fu.T_sst_inter = t5;
R_inter_fu.R_sstMID_inter = r6;
R_inter_fu.P_sstMID_inter = p6;
R_inter_fu.T_sstMID_inter = t6;

%% Save results
save('ABCD_FC_CBCL.mat', 'R_inter_fu', 'R_exter_fu', 'R_inter', 'R_exter')