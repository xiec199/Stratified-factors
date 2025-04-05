
clear;clc
addpath('/public/home/xic_fdu/Ongoing_project/2023_External_Internal_Factor')
load('N3_1_1_make_Ex_In_factor_FC.mat')
[BL_symptom,FU2_symptom,FU3_symptom] = xic_imagen_symptoms;

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

[~,ind1,ind2] = intersect(BL_symptom.BL_subject,BL_FC_sub);

[r p] = corr(BL_symptom.BL_Exter_beha_residual(ind1),BL_FC_in_sum(ind2));

[r p] = corr(BL_symptom.BL_Inter_beha_residual(ind1),BL_FC_ex_sum(ind2));

%% developmental changes

[sub_all,ind1,ind2,ind3] = xic_intersect(BL_FC_sub,FU2_FC_sub,FU3_FC_sub);
BL_ex  = BL_FC_ex_sum(ind1);  BL_in  = BL_FC_in_sum(ind1);
FU2_ex = FU2_FC_ex_sum(ind2); FU2_in = FU2_FC_in_sum(ind2);
FU3_ex = FU3_FC_ex_sum(ind3); FU3_in = FU3_FC_in_sum(ind3);

subplot(1,2,1);
bar([mean(BL_ex),mean(FU2_ex),mean(FU3_ex)])
title('Externalising')
subplot(1,2,2);
bar([mean(BL_in),mean(FU2_in),mean(FU3_in)])
title('Internalising');


for i=1:size(BL_ex,1)
    plot(1,BL_ex(i),'.','MarkerSize',10);hold on
    plot(2,FU2_ex(i),'.','MarkerSize',10);hold on
    plot(3,FU3_ex(i),'.','MarkerSize',10);hold on
    
end
xlim([0 4])

%% Correlation between the BL fc and FU2 and FU3 symptoms
[~,ind1,ind2] = intersect(BL_FC_sub,FU2_symptom.FU2_subject);
[ex_fu2(1), ex_fu2(2),ex_fu2(3)] = xic_corr_one(BL_FC_ex_sum(ind1),FU2_symptom.FU2_Exter_beha_residual(ind2));
[in_fu2(1), in_fu2(2),in_fu2(3)] = xic_corr_one(BL_FC_in_sum(ind1),FU2_symptom.FU2_Inter_beha_residual(ind2));

for i=1:4
    [ex_fu2_task(i,1), ex_fu2_task(i,2),ex_fu2_task(i,3)] = xic_corr_one(BL_FC_ex(ind1,i),FU2_symptom.FU2_Exter_beha_residual(ind2));
    [in_fu2_task(i,1), in_fu2_task(i,2),in_fu2_task(i,3)] = xic_corr_one(BL_FC_in(ind1,i),FU2_symptom.FU2_Inter_beha_residual(ind2));
end

[~,ind1,ind2] = intersect(BL_FC_sub,FU3_symptom.FU3_subject);
[ex_fu3(1), ex_fu3(2),ex_fu3(3)]  = xic_corr_one(BL_FC_ex_sum(ind1),FU3_symptom.FU3_Exter_beha_residual(ind2));
[in_fu3(1), in_fu3(2),in_fu3(3)]  = xic_corr_one(BL_FC_in_sum(ind1),FU3_symptom.FU3_Inter_beha_residual(ind2));

for i=1:4
    [ex_fu3_task(i,1), ex_fu3_task(i,2),ex_fu3_task(i,3)] = xic_corr_one(BL_FC_ex(ind1,i),FU3_symptom.FU3_Exter_beha_residual(ind2));
    [in_fu3_task(i,1), in_fu3_task(i,2),in_fu3_task(i,3)] = xic_corr_one(BL_FC_in(ind1,i),FU3_symptom.FU3_Inter_beha_residual(ind2));
end

%% Correlation between the FU2 fc and FU2 FU3 symptoms
[~,ind1,ind2] = intersect(FU2_FC_sub,FU2_symptom.FU2_subject);
[fu2_ex_fu2(1), fu2_ex_fu2(2),fu2_ex_fu2(3)] = xic_corr_one(FU2_FC_ex_sum(ind1),FU2_symptom.FU2_Exter_beha_residual(ind2));
[fu2_in_fu2(1), fu2_in_fu2(2),fu2_in_fu2(3)] = xic_corr_one(FU2_FC_in_sum(ind1),FU2_symptom.FU2_Inter_beha_residual(ind2));

for i=1:4
    [fu2_ex_fu2_task(i,1), fu2_ex_fu2_task(i,2),fu2_ex_fu2_task(i,3)] = xic_corr_one(FU2_FC_ex(ind1,i),FU2_symptom.FU2_Exter_beha_residual(ind2));
    [fu2_in_fu2_task(i,1), fu2_in_fu2_task(i,2),fu2_in_fu2_task(i,3)] = xic_corr_one(FU2_FC_in(ind1,i),FU2_symptom.FU2_Inter_beha_residual(ind2));
end

[~,ind1,ind2] = intersect(FU2_FC_sub,FU3_symptom.FU3_subject);
[fu2_ex_fu3(1), fu2_ex_fu3(2),fu2_ex_fu3(3)]  = xic_corr(FU2_FC_ex_sum(ind1),FU3_symptom.FU3_Exter_beha_residual(ind2));
[fu2_in_fu3(1), fu2_in_fu3(2),fu2_in_fu3(3)]  = xic_corr(FU2_FC_in_sum(ind1),FU3_symptom.FU3_Inter_beha_residual(ind2));

for i=1:4
    [fu2_ex_fu3_task(i,1), fu2_ex_fu3_task(i,2),fu2_ex_fu3_task(i,3)] = xic_corr_one(FU2_FC_ex(ind1,i),FU3_symptom.FU3_Exter_beha_residual(ind2));
    [fu2_in_fu3_task(i,1), fu2_in_fu3_task(i,2),fu2_in_fu3_task(i,3)] = xic_corr_one(FU2_FC_in(ind1,i),FU3_symptom.FU3_Inter_beha_residual(ind2));
end

%% %% Correlation between the FU3 fc and FU3 symptoms
[~,ind1,ind2] = intersect(FU3_FC_sub,FU3_symptom.FU3_subject);
[fu3_ex_fu3(1), fu3_ex_fu3(2),fu3_ex_fu3(3)] = xic_corr_one(FU3_FC_ex_sum(ind1),FU3_symptom.FU3_Exter_beha_residual(ind2));
[fu3_in_fu3(1), fu3_in_fu3(2),fu3_in_fu3(3)] = xic_corr_one(FU3_FC_in_sum(ind1),FU3_symptom.FU3_Inter_beha_residual(ind2));

for i=1:4
    [fu3_ex_fu3_task(i,1), fu3_ex_fu3_task(i,2),fu3_ex_fu3_task(i,3)] = xic_corr_one(FU3_FC_ex(ind1,i),FU3_symptom.FU3_Exter_beha_residual(ind2));
    [fu3_in_fu3_task(i,1), fu3_in_fu3_task(i,2),fu3_in_fu3_task(i,3)] = xic_corr_one(FU3_FC_in(ind1,i),FU3_symptom.FU3_Inter_beha_residual(ind2));
end



%% the correlation between BL and FU2 brain and behaviour
[~,ind1,ind2,ind3] = xic_intersect(BL_FC_sub,FU2_FC_sub,FU3_FC_sub);
[ex_brain_r,ex_brain_p] = corr([BL_FC_ex_sum(ind1),FU2_FC_ex_sum(ind2),FU3_FC_ex_sum(ind3)]);
[in_brain_r,in_brain_p] = corr([BL_FC_in_sum(ind1),FU2_FC_in_sum(ind2),FU3_FC_in_sum(ind3)]);

[~,ind1,ind2,ind3] = xic_intersect(BL_symptom.BL_subject,FU2_symptom.FU2_subject,FU3_symptom.FU3_subject);
[ex_beha_r,ex_beha_p] = corr([BL_symptom.BL_Exter_beha_residual(ind1),FU2_symptom.FU2_Exter_beha_residual(ind2),FU3_symptom.FU3_Exter_beha_residual(ind3)]);
[in_beha_r,in_beha_p] = corr([BL_symptom.BL_Inter_beha_residual(ind1),FU2_symptom.FU2_Inter_beha_residual(ind2),FU3_symptom.FU3_Inter_beha_residual(ind3)]);

%% Correlation between the BL fc and FU2 and FU3 symptoms controlling the BL symptoms
[~,ind1,ind2,ind3] = xic_intersect(BL_FC_sub,FU2_symptom.FU2_subject,BL_symptom.BL_subject);
[ex_fu2_r_con, ex_fu2_p_con] = partialcorr(BL_FC_ex_sum(ind1),FU2_symptom.FU2_Exter_beha_residual(ind2),BL_symptom.BL_Exter_beha_residual(ind3));
[in_fu2_r_con, in_fu2_p_con] = partialcorr(BL_FC_in_sum(ind1),FU2_symptom.FU2_Inter_beha_residual(ind2),BL_symptom.BL_Inter_beha_residual(ind3));

[~,ind1,ind2,ind3] = xic_intersect(BL_FC_sub,FU3_symptom.FU3_subject,BL_symptom.BL_subject);
[ex_fu3_r_con, ex_fu3_p_con] = partialcorr(BL_FC_ex_sum(ind1),FU3_symptom.FU3_Exter_beha_residual(ind2),BL_symptom.BL_Exter_beha_residual(ind3));
[in_fu3_r_con, in_fu3_p_con]  = partialcorr(BL_FC_in_sum(ind1),FU3_symptom.FU3_Inter_beha_residual(ind2),BL_symptom.BL_Inter_beha_residual(ind3));


