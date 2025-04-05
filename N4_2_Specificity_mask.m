
clear;clc
cd('/public/home/xic_fdu/Ongoing_project/2023_External_Internal_Factor')
addpath('/public/home/xic_fdu/Ongoing_project/PET_01_Pattern_SCZ_Hub/customcolormap');
load('N4_1_2_Ex_In_factor_mask_plot.mat')

%% regional level
oppo_path = '/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr';
load(fullfile(oppo_path,'A2_2_1_2_make_mask_across_single_v2.mat'))
share_pp = oppo_pp_dimen{1} + oppo_pp_dimen{2}  + oppo_pp_dimen{3} +oppo_pp_dimen{4};
share_pp_net = share_pp + tril(share_pp,-1)';
share_degree = sum(share_pp_net,2);

load('N2_1_5_make_Ex_In_factor_eachFC.mat')
exter_all = unique(cell2mat(exter_pospos_type'));
inter_all = unique(cell2mat(inter_negneg_type'));

exter_mask = zeros(268,268);
exter_mask(exter_all) = 1;
exter_mask = exter_mask + tril(exter_mask,-1)';
exter_degree = sum(exter_mask,2);

inter_mask = zeros(268,268);
inter_mask(inter_all) = 1;
inter_mask = inter_mask + tril(inter_mask,-1)';
inter_degree = sum(inter_mask,2);

share_weight = share_degree./sum(share_degree)+1;
exter_weight = exter_degree./sum(exter_degree)+1;
inter_weight = inter_degree./sum(inter_degree)+1;

share_speci = zscore(share_weight./exter_weight + share_weight./inter_weight);
exter_speci = zscore(exter_weight./share_weight + exter_weight./inter_weight);
inter_speci = zscore(inter_weight./share_weight + inter_weight./exter_weight);

speci_all = [share_speci,exter_speci,inter_speci];

for i=1:268
    roi = speci_all(i,:);
    roi_dis = abs(min(roi)-0);
    if (min(roi)) < 0 
        roi_nor = roi + roi_dis ;
    else
        roi_nor = roi - roi_dis ;
    end
    roi_all(i,:) = roi_nor; 
end

for i=1:268
    x_roi = roi_all(i,:);
    [~,ind(i,1)] = max(x_roi);
end
[ind1,num]  = sort(ind);
roi_all_reorder = roi_all(num,:);

roi_type_sort = zeros(268,3);
for i=1:3
    roi_type = roi_all_reorder(ind1==i,:);
    
    clear type_dis
    for j=1:length(roi_type)
        type_dis(j,1) = max(roi_type(j,:)) -min(roi_type(j,:));    
    end
    [ind,num1] = sort(type_dis,'descend');
    roi_type_sort(ind1==i,:) = roi_type(num1,:);
end


fig1 = figure(1);
set(fig1, 'Position', [100, 100, 1200, 600]);
set(fig1, 'Color', 'w');

plot(roi_type_sort(:,1),'.','MarkerSize',10,'Color',[93 168 78]/255);hold on
plot(roi_type_sort(:,2),'.','MarkerSize',10,'Color',[225 45 61]/255);hold on
plot(roi_type_sort(:,3),'.','MarkerSize',10,'Color',[57 122 188]/255); hold on

for i=1:268
    lin1 = roi_type_sort(i,:);
    roi_lin = sort(lin1);
    lin_min = roi_lin(2);
    lin_max = roi_lin(3);
    y = [lin_min,lin_max];
    x = [i,i];
    plot(x, y,'color',[180 180 180]/255);
end

ylim([-0.2 14])
xlim([-10 280])

%% Proportion
share_speci_info = roi_all_reorder(ind1==1,:);
exter_speci_info = roi_all_reorder(ind1==2,:);
inter_speci_info = roi_all_reorder(ind1==3,:);

share_roi = num(ind1 == 1);
exter_roi = num(ind1 == 2);
inter_roi = num(ind1 == 3);


[p, tbl, stats] = anova1([share_speci_info(:,1)', exter_speci_info(:,2)', inter_speci_info(:,3)']', ind1, 'off');
disp(tbl);
[h,p,c,t] = ttest2(share_speci_info(:,1),exter_speci_info(:,2));
[h,p,c,t] = ttest2(share_speci_info(:,1),inter_speci_info(:,3));


for i=1:length(share_speci_info)
    roi_share = share_speci_info(i,:);
    [~,share_ind(i,1)] = min(roi_share);
end

for i=1:length(exter_speci_info)
    roi_exter = exter_speci_info(i,:);
    [~,exter_ind(i,1)] = min(roi_exter);
end


for i=1:length(inter_speci_info)
    roi_inter = inter_speci_info(i,:);
    [~,inter_ind(i,1)] = min(roi_inter);
end


fig2 = figure(2);
set(fig2, 'Position', [100, 100, 1200, 600]);
set(fig2, 'Color', 'w');

subplot(1,3,1)
pie([sum(share_ind==2),sum(share_ind==3)]);
title('Share:Ex-In')
subplot(1,3,2)
pie([sum(exter_ind==1),sum(exter_ind==3)]);
title('Exter:Sh-In')
subplot(1,3,3)
pie([sum(inter_ind==1),sum(inter_ind==2)]);
title('Inter:Sh-Ex')

%% ROI mask 
addpath('/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr');

share_speci = roi_all_reorder(ind1==1,1);
share_speci_order = sort(share_speci,'descend');
share_speci_mask = zeros(268,1);
share_speci_mask(num(ind1==1,:)) = share_speci;
share_speci_mask(share_speci_mask<4.0) = 0;

exter_speci = roi_all_reorder(ind1==2,2);
exter_speci_order = sort(exter_speci,'descend');
exter_speci_mask = zeros(268,1);
exter_speci_mask(num(ind1==2,:)) = exter_speci;
exter_speci_mask(exter_speci_mask<2.42) = 0;

inter_speci = roi_all_reorder(ind1==3,3);
inter_speci_order = sort(inter_speci,'descend');
inter_speci_mask = zeros(268,1);
inter_speci_mask(num(ind1==3,:)) = inter_speci;
inter_speci_mask(inter_speci_mask<2.62) = 0;

out_path = '/public/home/xic_fdu/Ongoing_project/2023_External_Internal_Factor/exter_inter_mask';

xic_shen_template(share_speci_mask,'Specific_share_mask_filter.nii',out_path);
xic_shen_template(exter_speci_mask,'Specific_exter_mask_filter.nii',out_path);
xic_shen_template(inter_speci_mask,'Specific_inter_mask_filter.nii',out_path);


%% Network -level 
load('N4_1_2_Ex_In_factor_mask_plot.mat')
mask_ind = tril(ones(10,10))==1;
mask_num = zeros(10,10);
mask_num(mask_ind > 0) = 1:55;

exter_num = exter_mat(mask_ind >0);
inter_num = inter_mat(mask_ind >0);
share_num = oppo_mat(mask_ind >0);

exter_weight = exter_num./sum(exter_num)+1;
inter_weight = inter_num./sum(inter_num)+1;
share_weight = share_num./sum(share_num)+1;

share_speci = zscore(share_weight./exter_weight + share_weight./inter_weight);
exter_speci = zscore(exter_weight./share_weight + exter_weight./inter_weight);
inter_speci = zscore(inter_weight./share_weight + inter_weight./exter_weight);


speci_all = [share_speci,exter_speci,inter_speci];
clear roi_all
for i=1:55
    roi = speci_all(i,:);
    roi_dis = abs(min(roi)-0);
    if (min(roi)) < 0 
        roi_nor = roi + roi_dis ;
    else
        roi_nor = roi - roi_dis ;
    end
    roi_all(i,:) = roi_nor; 
end

clear ind
for i=1:55
    x_roi = roi_all(i,:);
    [~,ind(i,1)] = max(x_roi);
end
[ind1,num]  = sort(ind);
roi_all_reorder = roi_all(num,:);

roi_type_sort = zeros(55,3);
for i=1:3
    roi_type = roi_all_reorder(ind1==i,:);
    
    clear type_dis
    for j=1:length(roi_type)
        type_dis(j,1) = max(roi_type(j,:)) -min(roi_type(j,:));    
    end
    [ind,num1] = sort(type_dis,'descend');
    roi_type_sort(ind1==i,:) = roi_type(num1,:);
end


fig3 = figure(3);
set(fig3, 'Position', [100, 100, 1200, 600]);
set(fig3, 'Color', 'w');

plot(roi_type_sort(:,1),'.','MarkerSize',10,'Color',[93 168 78]/255);hold on
plot(roi_type_sort(:,2),'.','MarkerSize',10,'Color',[225 45 61]/255);hold on
plot(roi_type_sort(:,3),'.','MarkerSize',10,'Color',[57 122 188]/255); hold on

for i=1:55
    lin1 = roi_type_sort(i,:);
    roi_lin = sort(lin1);
    lin_min = roi_lin(2);
    lin_max = roi_lin(3);
    y = [lin_min,lin_max];
    x = [i,i];
    plot(x, y,'color',[180 180 180]/255);
end

ylim([-0.5 8])
xlim([-5 60])

share_speci_info = roi_all_reorder(ind1==1,:);
exter_speci_info = roi_all_reorder(ind1==2,:);
inter_speci_info = roi_all_reorder(ind1==3,:);

share_roi = num(ind1 == 1);
exter_roi = num(ind1 == 2);
inter_roi = num(ind1 == 3);


[p, tbl, stats] = anova1([share_speci_info(:,1)', exter_speci_info(:,2)', inter_speci_info(:,3)']', ind1, 'off');
disp(tbl);
[h,p,c,t] = ttest2(share_speci_info(:,1),exter_speci_info(:,2));


clear share_ind  exter_ind inter_ind
for i=1:length(share_speci_info)
    roi_share = share_speci_info(i,:);
    [~,share_ind(i,1)] = min(roi_share);
end

for i=1:length(exter_speci_info)
    roi_exter = exter_speci_info(i,:);
    [~,exter_ind(i,1)] = min(roi_exter);
end


for i=1:length(inter_speci_info)
    roi_inter = inter_speci_info(i,:);
    [~,inter_ind(i,1)] = min(roi_inter);
end

fig4 = figure(4);
set(fig4, 'Position', [100, 100, 1200, 600]);
set(fig4, 'Color', 'w');

subplot(1,3,1)
pie([sum(share_ind==2),sum(share_ind==3)]);
title('Share:Ex-In')
subplot(1,3,2)
pie([sum(exter_ind==1),sum(exter_ind==3)]);
title('Exter:Sh-In')
subplot(1,3,3)
pie([sum(inter_ind==1),sum(inter_ind==2)]);
title('Inter:Sh-Ex')

share_speci = roi_all_reorder(ind1==1,1);
share_speci_order = sort(share_speci,'descend');
share_speci_mask = zeros(55,1);
share_speci_mask(num(ind1==1,:)) = share_speci;
share_speci_mask_sort = sort(share_speci_mask,'descend');
share_net = zeros(10,10);
share_net(mask_ind >0) = share_speci_mask;
share_net(share_net<2.3) = 0;

exter_speci = roi_all_reorder(ind1==2,2);
exter_speci_order = sort(exter_speci,'descend');
exter_speci_mask = zeros(55,1);
exter_speci_mask(num(ind1==2,:)) = exter_speci;
exter_net = zeros(10,10);
exter_net(mask_ind >0) = exter_speci_mask;
exter_net(exter_net<2.8) = 0;

inter_speci = roi_all_reorder(ind1==3,3);
inter_speci_order = sort(inter_speci,'descend');
inter_speci_mask = zeros(55,1);
inter_speci_mask(num(ind1==3,:)) = inter_speci;
inter_net = zeros(10,10);
inter_net(mask_ind >0) = inter_speci_mask;
inter_net(inter_net<2.1) = 0;

fig5 = figure(5);
set(fig5, 'Position', [100, 100, 1200, 400]);
set(fig5, 'Color', 'w');

red = customcolormap([0 0.5 1], {'#F40000','#F78E8D','#FFFFFF'});   %
blue = customcolormap([0 0.5 1], {'#2165A6','#4F80FF','#FFFFFF'});   %

subplot(1,3,1)
imagesc(share_net);
title('Share')
subplot(1,3,2)
imagesc(exter_net);
title('Exter')
subplot(1,3,3)
imagesc(inter_net);
title('Inter')
colormap(red);
colorbar
