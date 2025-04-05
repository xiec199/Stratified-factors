
clear;clc
load('N1_1_1_make_Ex_In_factor_part1.mat')
load('N2_1_2_make_Ex_In_factor_Permute.mat')

%% calculate the number of edges in external and internal specific brain network
for i=1:8
    exter_all = [TableS1_1.exter_pos_pos_all_dimension_id{i}; ...
                 TableS1_1.exter_pos_neg_all_dimension_id{i}; ...
                 TableS1_1.exter_neg_neg_all_dimension_id{i}; ...
                 TableS1_1.exter_neg_pos_all_dimension_id{i}];
                 
    exter_all_uni = unique(exter_all);
    
    inter_all = [TableS1_1.inter_pos_pos_all_dimension_id{i}; ...
                 TableS1_1.inter_pos_neg_all_dimension_id{i}; ...
                 TableS1_1.inter_neg_neg_all_dimension_id{i}; ...
                 TableS1_1.inter_neg_pos_all_dimension_id{i}];
                 
    inter_all_uni = unique(inter_all);   
        
    share_all = [TableS1_1.share_pos_pos_all_dimension_id{i}; ...
                 TableS1_1.share_pos_neg_all_dimension_id{i}; ...
                 TableS1_1.share_neg_neg_all_dimension_id{i}; ...
                 TableS1_1.share_neg_pos_all_dimension_id{i}];
                 
    share_all_uni = unique(share_all);
    
    exter_uni = setdiff(exter_all_uni,share_all_uni);
    inter_uni = setdiff(inter_all_uni,share_all_uni);
    
    exter_uni_num(i,1) = length(exter_uni);
    inter_uni_num(i,1) = length(inter_uni);
    
end

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
    
for i=1:8
    disp([j,i])
    
    exter_all = [exter_pospos{i},exter_posneg{i},exter_negneg{i},exter_negpos{i}];             
    exter_all_uni = unique(exter_all);
    
    inter_all = [inter_pospos{i}, inter_posneg{i}, inter_negneg{i},inter_negpos{i}];             
    inter_all_uni = unique(inter_all);
    
    share_all = [share_pospos{i}, share_posneg{i}, share_negneg{i},share_negpos{i}];             
    share_all_uni = unique(share_all);
    
    exter_uni_perm = setdiff(exter_all_uni,share_all_uni);
    inter_uni_perm = setdiff(inter_all_uni,share_all_uni);
    
    exter_uni_num_perm(j,i) = length(exter_uni_perm);
    inter_uni_num_perm(j,i) = length(inter_uni_perm);
    
end
end

%% correction for the significant overlap of exter. and inter. FC for each condition

for i=1:8
    Corr_P_exter(i,1) = 1- sum(exter_uni_num(i)>exter_uni_num_perm(:,i))/1000;
    Corr_P_inter(i,1) = 1- sum(inter_uni_num(i)>inter_uni_num_perm(:,i))/1000;
end


