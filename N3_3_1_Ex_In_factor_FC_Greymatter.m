
clear;clc
bl  = readtable('/public/mig_old_storage/home1/ISTBI_data/IMAGEN/Trajectory/BL/aparc_GrayVol.txt');
fu2  = readtable('/public/mig_old_storage/home1/ISTBI_data/IMAGEN/Trajectory/FU2/aparc_GrayVol.txt');
fu3  = readtable('/public/mig_old_storage/home1/ISTBI_data/IMAGEN/Trajectory/FU3/aparc_GrayVol.txt');

bl_sub = bl.SubID;
bl_mean = mean(table2array(bl(:,2:end)),2);
bl_sub(isnan(bl_mean)) =[];
bl_mean(isnan(bl_mean)) =[];

fu2_sub = fu2.SubID;
fu2_mean = mean(table2array(fu2(:,2:end)),2);
fu2_sub(isnan(fu2_mean)) =[];
fu2_mean(isnan(fu2_mean)) =[];

fu3_sub = fu3.SubID;
fu3_mean = mean(table2array(fu3(:,2:end)),2);
fu3_sub(isnan(fu3_mean)) =[];
fu3_mean(isnan(fu3_mean)) =[];

path = '/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr';
addpath(path);

load('N3_1_1_make_Ex_In_factor_FC.mat')
BL_FC_ex = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.BL_ex_pp_only_FC_sum');
BL_FC_in = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.BL_in_nn_only_FC_sum');
BL_FC_sub = TableS3_1_1_make_factor_fc.BL_subject_match;
BL_FC_ex_sum = mean(BL_FC_ex,2);
BL_FC_in_sum = mean(BL_FC_in,2);

FU2_FC_ex = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU2_ex_pp_only_FC_sum');
FU2_FC_in = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU2_in_nn_only_FC_sum');
FU2_FC_sub = TableS3_1_1_make_factor_fc.FU2_subject_match;
FU2_FC_ex_sum = mean(FU2_FC_ex,2);
FU2_FC_in_sum = mean(FU2_FC_in,2);

FU3_FC_ex = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU3_ex_pp_only_FC_sum');
FU3_FC_in = cell2mat(TableS3_1_1_make_factor_fc.Across_FC.FU3_in_nn_only_FC_sum');
FU3_FC_sub = TableS3_1_1_make_factor_fc.FU3_subject_match;
FU3_FC_ex_sum = mean(FU3_FC_ex,2);
FU3_FC_in_sum = mean(FU3_FC_in,2);

% regres the covariates
load('/home1/xic_fdu/Datasets/Match_control_subject_All things scale/Cova_subject_jia.mat')
[sub_cova_bl,ind1,ind2] = intersect(cova_subject,BL_FC_sub);
[~,~,bl_fc_ex] = regress(BL_FC_ex_sum(ind2),cova_data(ind1,1:9));
[~,~,bl_fc_in] = regress(BL_FC_in_sum(ind2),cova_data(ind1,1:9));

[sub_cova_fu2,ind1,ind2] = intersect(cova_subject,FU2_FC_sub);
[~,~,fu2_fc_ex] = regress(FU2_FC_ex_sum(ind2),cova_data(ind1,1:9));
[~,~,fu2_fc_in] = regress(FU2_FC_in_sum(ind2),cova_data(ind1,1:9));

[sub_cova_fu3,ind1,ind2] = intersect(cova_subject,FU3_FC_sub);
[~,~,fu3_fc_ex] = regress(FU3_FC_ex_sum(ind2),cova_data(ind1,1:9));
[~,~,fu3_fc_in] = regress(FU3_FC_in_sum(ind2),cova_data(ind1,1:9));

%% comapre the changes 
[~,ind1,ind2] = intersect(sub_cova_bl,bl_sub);
[r p] = corr(bl_fc_ex(ind1),bl_mean(ind2))
[r p] = corr(bl_fc_in(ind1),bl_mean(ind2))

[~,ind1,ind2] = intersect(sub_cova_bl,fu2_sub);
[r p] = corr(bl_fc_in(ind1),c(ind2))

[~,ind1,ind2] = intersect(sub_cova_fu3,fu3_sub);
[r p] = corr(fu3_fc_ex(ind1),fu3_mean(ind2))


%% comapre the trajectory data
[fc_bl_fu3,ind1,ind2] = intersect(sub_cova_bl,sub_cova_fu3); 
ex_diff = bl_fc_ex(ind1) - fu3_fc_ex(ind2);
in_diff  = bl_fc_in(ind1) - fu3_fc_in(ind2);


[all_sub,ind1,ind2,ind3] = xic_intersect(fc_bl_fu3,bl_sub,fu3_sub);
ex_diff =  ex_diff(ind1);
in_diff =  in_diff(ind1);

bl_grey = bl_mean(ind2);
fu3_grey = fu3_mean(ind3);

[r p t]  = xic_corr(ex_diff,bl_grey-fu3_grey);
[r p t]  = xic_corr(in_diff,bl_grey-fu3_grey);


grey_diff = bl_grey-fu3_grey;

[bl_ex_sort,ind] = sort(bl_in,'descend');

bl_ex_high  = ind(1:226);
bl_ex_low  = ind((length(ind)-225):end);

ex_high_grey = grey_diff(bl_ex_high);
ex_low_grey = grey_diff(bl_ex_low);

mean(ex_high_grey)
mean(ex_low_grey)




