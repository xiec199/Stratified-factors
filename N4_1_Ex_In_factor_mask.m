
clear;clc
load('N2_1_5_make_Ex_In_factor_eachFC.mat')
addpath('/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr')
out_path = '/public/home/xic_fdu/Ongoing_project/2023_External_Internal_Factor/exter_inter_mask';

exter_all = unique(cell2mat(exter_pospos_type'));
inter_all = unique(cell2mat(inter_negneg_type'));

%% all condition
exter_mask = zeros(268,268);
exter_mask(exter_all) = 1;
exter_mask = exter_mask + exter_mask';

inter_mask = zeros(268,268);
inter_mask(inter_all) = 1;
inter_mask = inter_mask + inter_mask';

dlmwrite('N4_1_1_exter_all.txt', exter_mask, 'delimiter', '\t');
dlmwrite('N4_1_1_inter_all.txt', inter_mask, 'delimiter', '\t');

exter_roi = sum(exter_mask,2); 
inter_roi = sum(inter_mask,2); 
[r,p] = corr(exter_roi, inter_roi,'type','spearman');

exter_roi(exter_roi<18) = 0;
inter_roi(inter_roi<8) = 0;

xic_shen_template(exter_roi,'exter_all_threshold.nii',out_path);
xic_shen_template(inter_roi,'inter_all_threshold.nii',out_path);
    
[degree_r, degree_p] = corr(sum(exter_mask,2),sum(inter_mask,2));

%% single condition
task = {'SST_stopSuccess','SST_stopFailure','MID_feedhit','MID_antici_hit'};
for i=1:4
    exter_task = zeros(268,268);
    exter_task(exter_pospos_type{i}) = 1;
    exter_task = exter_task + exter_task';
    exter_task_all(:,:,i) = exter_task;
    exter_task_degree(:,i) = sum(exter_task,2);
    
    inter_task = zeros(268,268);
    inter_task(inter_negneg_type{i}) = 1;
    inter_task = inter_task + inter_task';
    inter_task_all(:,:,i) = inter_task;
    inter_task_degree(:,i) = sum(inter_task,2);
    
    exter_out = ['N4_1_1_exter_',task{i},'.txt'];
    inter_out = ['N4_1_1_inter_',task{i},'.txt'];
    
    dlmwrite(exter_out, exter_task, 'delimiter', '\t');
    dlmwrite(inter_out, inter_task, 'delimiter', '\t');
    
    exter_task_nii = ['exter_',task{i},'.nii'];
    inter_task_nii = ['inter_',task{i},'.nii'];
    
    xic_shen_template(sum(exter_task,2),exter_task_nii,out_path);
    xic_shen_template(sum(inter_task,2),inter_task_nii,out_path);
end

[task_r, task_p] = corr([exter_task_degree,inter_task_degree]);
col2 = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'},60);
colormap(col2);

imagesc(task_r)

