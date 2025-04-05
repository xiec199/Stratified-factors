
clear;clc
load('N1_1_1_make_Ex_In_factor_part1.mat')
load('N2_1_2_make_Ex_In_factor_Permute.mat')

condi = [1,2,5,6];
%% calculate the number of edges in external and internal specific brain network

% exter and inter
[exter_pos_pos,exter_pos_neg,exter_neg_neg,exter_neg_pos] = type_fc_exter(TableS1_1,condi);               
[inter_pos_pos,inter_pos_neg,inter_neg_neg,inter_neg_pos] = type_fc_inter(TableS1_1,condi);
[share_pos_pos,share_pos_neg,share_neg_neg,share_neg_pos] = type_fc_share(TableS1_1,condi);
         
% exter and inter specific
[exter_pospos_type] = setdiff(exter_pos_pos,share_pos_pos);
[exter_posneg_type] = setdiff(exter_pos_neg,share_pos_neg);
[exter_negneg_type] = setdiff(exter_neg_neg,share_neg_neg);   
[exter_negpos_type] = setdiff(exter_neg_pos,share_neg_pos);   
           
[inter_pospos_type] = setdiff(inter_pos_pos,share_pos_pos);
[inter_posneg_type] = setdiff(inter_pos_neg,share_pos_neg);
[inter_negneg_type] = setdiff(inter_neg_neg,share_neg_neg);   
[inter_negpos_type] = setdiff(inter_neg_pos,share_neg_pos);   
 
exter_pospos_num = length(exter_pospos_type);
exter_posneg_num = length(exter_posneg_type);
exter_negneg_num = length(exter_negneg_type);
exter_negpos_num = length(exter_negpos_type);
    
inter_pospos_num = length(inter_pospos_type);
inter_posneg_num = length(inter_posneg_type);
inter_negneg_num = length(inter_negneg_type);
inter_negpos_num = length(inter_negpos_type);

save N2_1_4_make_Ex_In_factor_sumFC.mat exter_pospos_type inter_negneg_type
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
    exter_pos_pos_perm = [exter_pospos{1, 1},exter_pospos{1, 2},exter_pospos{1, 5},exter_pospos{1, 6}];
    exter_pos_neg_perm = [exter_posneg{1, 1},exter_posneg{1, 2},exter_posneg{1, 5},exter_posneg{1, 6}];
    exter_neg_neg_perm = [exter_negneg{1, 1},exter_negneg{1, 2},exter_negneg{1, 5},exter_negneg{1, 6}];
    exter_neg_pos_perm = [exter_negpos{1, 1},exter_negpos{1, 2},exter_negpos{1, 5},exter_negpos{1, 6}];
                     
    % inter
    inter_pos_pos_perm = [inter_pospos{1, 1},inter_pospos{1, 2},inter_pospos{1, 5},inter_pospos{1, 6}];
    inter_pos_neg_perm = [inter_posneg{1, 1},inter_posneg{1, 2},inter_posneg{1, 5},inter_posneg{1, 6}];
    inter_neg_neg_perm = [inter_negneg{1, 1},inter_negneg{1, 2},inter_negneg{1, 5},inter_negneg{1, 6}];
    inter_neg_pos_perm = [inter_negpos{1, 1},inter_negpos{1, 2},inter_negpos{1, 5},inter_negpos{1, 6}];
       
    % share
    share_pos_pos_perm = [share_pospos{1, 1},share_pospos{1, 2},share_pospos{1, 5},share_pospos{1, 6}];
    share_pos_neg_perm = [share_posneg{1, 1},share_posneg{1, 2},share_posneg{1, 5},share_posneg{1, 6}];
    share_neg_neg_perm = [share_negneg{1, 1},share_negneg{1, 2},share_negneg{1, 5},share_negneg{1, 6}];
    share_neg_pos_perm = [share_negpos{1, 1},share_negpos{1, 2},share_negpos{1, 5},share_negpos{1, 6}];
     
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

Corr_P_exter_pospos = 1- sum((exter_pospos_num)>exter_pospos_perm_num)/1000;
Corr_P_exter_posneg = 1- sum((exter_posneg_num)>exter_posneg_perm_num)/1000;
Corr_P_exter_negneg = 1- sum((exter_negneg_num)>exter_negneg_perm_num)/1000;
Corr_P_exter_negpos = 1- sum((exter_negpos_num)>exter_negpos_perm_num)/1000;

Corr_P_inter_pospos = 1- sum((inter_pospos_num)>inter_pospos_perm_num)/1000;
Corr_P_inter_posneg = 1- sum((inter_posneg_num)>inter_posneg_perm_num)/1000;
Corr_P_inter_negneg = 1- sum((inter_negneg_num)>inter_negneg_perm_num)/1000;
Corr_P_inter_negpos = 1- sum((inter_negpos_num)>inter_negpos_perm_num)/1000;


