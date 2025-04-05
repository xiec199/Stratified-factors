

clear;clc
addpath('/public/home/xic_fdu/Ongoing_project/PET_01_Pattern_SCZ_Hub/customcolormap/');
addpath('/home1/xic_fdu/xic_analysis/');
addpath('/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr')
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
%%
[all,ind1,ind2] = intersect(BL_FC_sub,BL_symptom.BL_subject);
[r p] = corr(BL_FC_ex_sum(ind1),BL_symptom.BL_Exter_beha_residual(ind2))
[r p] = corr(BL_FC_in_sum(ind1),BL_symptom.BL_Inter_beha_residual(ind2))


%% basic FC states
[ex_bl(1),ex_bl(2),ex_bl(3)] = xic_ttest(BL_FC_ex_sum);
[ex_fu2(1),ex_fu2(2),ex_fu2(3)] = xic_ttest(FU2_FC_ex_sum);
[ex_fu3(1),ex_fu3(2),ex_fu3(3)] = xic_ttest(FU3_FC_ex_sum);

[in_bl(1),in_bl(2),in_bl(3)] = xic_ttest(BL_FC_in_sum);
[in_fu2(1),in_fu2(2),in_fu2(3)] = xic_ttest(FU2_FC_in_sum);
[in_fu3(1),in_fu3(2),in_fu3(3)] = xic_ttest(FU3_FC_in_sum);

%% developmental changes

[sub_all,ind1,ind2,ind3] = xic_intersect(BL_FC_sub,FU2_FC_sub,FU3_FC_sub);
BL_ex  = BL_FC_ex_sum(ind1);  BL_in  = BL_FC_in_sum(ind1);
FU2_ex = FU2_FC_ex_sum(ind2); FU2_in = FU2_FC_in_sum(ind2);
FU3_ex = FU3_FC_ex_sum(ind3); FU3_in = FU3_FC_in_sum(ind3);

exter_data = [BL_ex,FU2_ex,FU3_ex];
inter_data = [BL_in,FU2_in,FU3_in];

[num_people, num_timepoints] = size(exter_data);

save('N3_2_2_Longitudinal_Change.mat','exter_data','inter_data');

figure;
for i = 1:667
    disp(i)
    person_ex = exter_data(i,:);
    person_in = inter_data(i,:);
    
    subplot(2,1,1);
    plot(1:num_timepoints, person_ex, '-o','Color', [[221 62 74]/255,0.2], 'MarkerFaceColor', [221 62 74]/255);    
    hold on; 
    
    subplot(2,1,2);
    plot(1:num_timepoints, person_in, '-o','Color', [[68 122 188]/255,0.2], 'MarkerFaceColor',[68 122 188]/255);    
    hold on; 
   
    coeff_ex = polyfit(1:num_timepoints, person_ex, 1);
    coeff_in = polyfit(1:num_timepoints, person_in, 1);
    slop_ex(i,1) = coeff_ex(1);
    slop_in(i,1) = coeff_in(1);
end

subplot(2,1,1);
xlim([0.5 3.5])
title('Externalising')
xticks([1:3]);
xticklabels({'Age-14','Age-19','Age-23'});

subplot(2,1,2);
xlim([0.5 3.5])
title('Internalising')
xticks([1:3]);
xticklabels({'Age-14','Age-19','Age-23'});

%% changes difference

slop_diff = slop_ex - slop_in; 
[r p] = corr(slop_ex, slop_in);
[h p c, t] = ttest(slop_ex,slop_in);
[h p c, t] = ttest(slop_in);


d = cohen_t(t.tstat,length(slop_ex))


[~,ind1,ind2,ind3] = xic_intersect(sub_all,BL_symptom.BL_subject,FU2_symptom.FU2_subject);
[r p t] = xic_corr(slop_ex(ind1),BL_symptom.BL_Exter_beha_residual(ind2));
[r p t] = xic_corr(slop_in(ind1),BL_symptom.BL_Inter_beha_residual(ind2));

[r p] = corr(slop_diff(ind1),BL_symptom.BL_Exter_beha_residual(ind2));
[r p] = corr(slop_diff(ind1),FU2_symptom.FU2_Exter_beha_residual(ind3));



