
clear;clc
load('N1_1_1_make_Ex_In_factor_part1.mat')
load('N2_1_2_make_Ex_In_factor_Permute.mat')

condi = [1,2,5,6];
%% calculate the number of edges in external and internal specific brain network

for i=1:4
    % exter and inter
    [exter_pos_pos,exter_pos_neg,exter_neg_neg,exter_neg_pos] = type_fc_exter_condi(TableS1_1,condi(i));               
    [inter_pos_pos,inter_pos_neg,inter_neg_neg,inter_neg_pos] = type_fc_inter_condi(TableS1_1,condi(i));
    [share_pos_pos,share_pos_neg,share_neg_neg,share_neg_pos] = type_fc_share_condi(TableS1_1,condi(i));
         
    % exter and inter specific
    [exter_pospos_type{i}] = setdiff(exter_pos_pos,share_pos_pos);
    [exter_posneg_type{i}] = setdiff(exter_pos_neg,share_pos_neg);
    [exter_negneg_type{i}] = setdiff(exter_neg_neg,share_neg_neg);   
    [exter_negpos_type{i}] = setdiff(exter_neg_pos,share_neg_pos);   
           
    [inter_pospos_type{i}] = setdiff(inter_pos_pos,share_pos_pos);
    [inter_posneg_type{i}] = setdiff(inter_pos_neg,share_pos_neg);
    [inter_negneg_type{i}] = setdiff(inter_neg_neg,share_neg_neg);   
    [inter_negpos_type{i}] = setdiff(inter_neg_pos,share_neg_pos);   
     
    exter_pospos_type_num(i,1) = length(exter_pospos_type{i});
    exter_posneg_type_num(i,1) = length(exter_posneg_type{i});
    exter_negneg_type_num(i,1) = length(exter_negneg_type{i});
    exter_negpos_type_num(i,1) = length(exter_negpos_type{i});

    inter_pospos_type_num(i,1) = length(inter_pospos_type{i});
    inter_posneg_type_num(i,1) = length(inter_posneg_type{i});
    inter_negneg_type_num(i,1) = length(inter_negneg_type{i});
    inter_negpos_type_num(i,1) = length(inter_negpos_type{i});
    
%% Permute calculate the number of edges in external and internal specific brain network
for j=1:1000
    exter_pospos = TableS2_1_perm.Exter_PosPos(j,:);
    exter_posneg = TableS2_1_perm.Exter_PosNeg(j,:);
    exter_negneg = TableS2_1_perm.Exter_NegNeg(j,:);
    exter_negpos = TableS2_1_perm.Exter_NegPos(j,:);
    
    inter_pospos = TableS2_1_perm.Inter_PosPos(j,:);
    inter_posneg = TableS2_1_perm.Inter_PosNeg(j,:);
    inter_negneg = TableS2_1_perm.Inter_NegNeg(j,:);
    inter_negpos = TableS2_1_perm.Inter_NegPos(j,:);

    share_pospos = TableS2_1_perm.Share_PosPos(j,:);
    share_posneg = TableS2_1_perm.Share_PosNeg(j,:);
    share_negneg = TableS2_1_perm.Share_NegNeg(j,:);
    share_negpos = TableS2_1_perm.Share_NegPos(j,:);
  
    % exter
    exter_pos_pos_perm = [exter_pospos{1, condi(i)}];
    exter_pos_neg_perm = [exter_posneg{1, condi(i)}];
    exter_neg_neg_perm = [exter_negneg{1, condi(i)}];
    exter_neg_pos_perm = [exter_negpos{1, condi(i)}];
                     
    % inter
    inter_pos_pos_perm = [inter_pospos{1, condi(i)}];
    inter_pos_neg_perm = [inter_posneg{1, condi(i)}];
    inter_neg_neg_perm = [inter_negneg{1, condi(i)}];
    inter_neg_pos_perm = [inter_negpos{1, condi(i)}];
       
    % share
    share_pos_pos_perm = [share_pospos{1, condi(i)}];
    share_pos_neg_perm = [share_posneg{1, condi(i)}];
    share_neg_neg_perm = [share_negneg{1, condi(i)}];
    share_neg_pos_perm = [share_negpos{1, condi(i)}];
     
    % exter and inter specific
    [exter_pospos_type_perm] = setdiff(exter_pos_pos_perm,share_pos_pos_perm);
    [exter_posneg_type_perm] = setdiff(exter_pos_neg_perm,share_pos_neg_perm);
    [exter_negneg_type_perm] = setdiff(exter_neg_neg_perm,share_neg_neg_perm);   
    [exter_negpos_type_perm] = setdiff(exter_neg_pos_perm,share_neg_pos_perm);   
           
    [inter_pospos_type_perm] = setdiff(inter_pos_pos_perm,share_pos_pos_perm);
    [inter_posneg_type_perm] = setdiff(inter_pos_neg_perm,share_pos_neg_perm);
    [inter_negneg_type_perm] = setdiff(inter_neg_neg_perm,share_neg_neg_perm);   
    [inter_negpos_type_perm] = setdiff(inter_neg_pos_perm,share_neg_pos_perm); 
    
    exter_pospos_perm_num(j,1) = length(exter_pospos_type_perm);
    exter_posneg_perm_num(j,1) = length(exter_posneg_type_perm);
    exter_negneg_perm_num(j,1) = length(exter_negneg_type_perm);
    exter_negpos_perm_num(j,1) = length(exter_negpos_type_perm);
    
    inter_pospos_perm_num(j,1) = length(inter_pospos_type_perm);
    inter_posneg_perm_num(j,1) = length(inter_posneg_type_perm);
    inter_negneg_perm_num(j,1) = length(inter_negneg_type_perm);
    inter_negpos_perm_num(j,1) = length(inter_negpos_type_perm);
end

%% correction for the significant overlap of exter. and inter. FC for each condition

Corr_P_exter_pospos(i,1) = 1- sum(exter_pospos_type_num(i,1)>inter_pospos_perm_num)/1000;
Corr_P_exter_posneg(i,1) = 1- sum(exter_posneg_type_num(i,1)>exter_posneg_perm_num)/1000;
Corr_P_exter_negneg(i,1) = 1- sum(exter_negneg_type_num(i,1)>exter_negneg_perm_num)/1000;
Corr_P_exter_negpos(i,1) = 1- sum(exter_negpos_type_num(i,1)>exter_negpos_perm_num)/1000;

Corr_P_inter_pospos(i,1) = 1- sum(inter_pospos_type_num(i,1)>inter_pospos_perm_num)/1000;
Corr_P_inter_posneg(i,1) = 1- sum(inter_posneg_type_num(i,1)>inter_posneg_perm_num)/1000;
Corr_P_inter_negneg(i,1) = 1- sum(inter_negneg_type_num(i,1)>inter_negneg_perm_num)/1000;
Corr_P_inter_negpos(i,1) = 1- sum(inter_negpos_type_num(i,1)>inter_negpos_perm_num)/1000;
end

save N2_1_5_make_Ex_In_factor_eachFC.mat exter_pospos_type inter_negneg_type

