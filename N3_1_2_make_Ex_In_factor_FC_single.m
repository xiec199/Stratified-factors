
clear;clc

load('N2_1_5_make_Ex_In_factor_eachFC.mat')

[Matrix_name,BL_matrix_all,BL_subject,BL_subject_all, ...
             FU2_matrix_all,FU2_subject,FU2_subject_all, ...
             FU3_matrix_all,FU3_subject,FU3_subject_all]= read_fc_ex_in;
         
BL_subject_match = BL_subject{1}(BL_subject_all(:,1));
FU2_subject_match = FU2_subject{1}(FU2_subject_all(:,1));
FU3_subject_match = FU3_subject{1}(FU3_subject_all(:,1));

name = {'SST_stopsucc';'SST_stopfail';'MID_feedhit';'MID_anticihit'};


load('/home1/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr/FU2_shen_rest/Shen_Rest_matrix.mat');

fc_fu2_rest = reshape(FU2_matrix_rest,[],size(FU2_matrix_rest,3))';

Across_FC = table;
Across_FC.name = name;

for i=1:4
       fc_bl = BL_matrix_all{i}; 
       fc_fu2 = FU2_matrix_all{i};
       fc_fu3 = FU3_matrix_all{i};
    
       % BL       
       eval('Across_FC.BL_ex_pp{i} = (fc_bl(BL_subject_all(:,i),exter_pospos_type{1,i}));');
       eval('Across_FC.BL_in_nn{i} = (fc_bl(BL_subject_all(:,i),inter_negneg_type{1,i}));');
       
       % FU2
       eval('Across_FC.FU2_ex_pp{i} = (fc_fu2(FU2_subject_all(:,i),exter_pospos_type{1,i}));');
       eval('Across_FC.FU2_in_nn{i} = (fc_fu2(FU2_subject_all(:,i),inter_negneg_type{1,i}));');
       
       % FU2_Rest
       eval('Across_FC.FU2_Rest_ex_pp{i} = (fc_fu2_rest(:,exter_pospos_type{1,i}));');
       eval('Across_FC.FU2_Rest_in_nn{i} = (fc_fu2_rest(:,inter_negneg_type{1,i}));');
       
       % FU3
       eval('Across_FC.FU3_ex_pp{i} = (fc_fu3(FU3_subject_all(:,i),exter_pospos_type{1,i}));');
       eval('Across_FC.FU3_in_nn{i} = (fc_fu3(FU3_subject_all(:,i),inter_negneg_type{1,i}));');
       
end

TableS3_1_1_make_factor_fc.Across_FC = Across_FC;
TableS3_1_1_make_factor_fc.BL_subject_match = BL_subject_match';
TableS3_1_1_make_factor_fc.FU2_subject_match = FU2_subject_match';
TableS3_1_1_make_factor_fc.FU3_subject_match = FU3_subject_match';
TableS3_1_1_make_factor_fc.FU2_rest_sub  = FU2_rest_sub;

save N3_1_2_make_Ex_In_factor_FC_single TableS3_1_1_make_factor_fc




