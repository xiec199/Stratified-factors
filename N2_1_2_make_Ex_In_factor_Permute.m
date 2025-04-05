
%%
clear;clc;close all

model_path  = '/CPM_permute/';

model_name = {'_SST_stopsuc.mat','_SST_stopfai.mat' ...
        'SST_gowrong.mat','feed_miss.mat','feed_hit.mat' ...
        'antici_hit.mat','EFT_neutral.mat','EFT_angry.mat'};
 
name = {'SST:stop success','SST:stop failure','SST:go wrong',...
       'MID:feedmiss','MID:feedhit','MID: antici hit', ...
       'EFT:neutral','EFT:angry'};


number=[3 1 4 7 2 5 6 8];
dis_names  = {'asd','adhd','cd' ,'od','anxiety','dep','eat','speph'};
  
load('/public/home/xic_fdu/Finished_project/2019_NP_factor/IAMGEN_Develop_diagnostic_0814/Trans_diagnostic_p_facotr/A2_0_Mean_brain_behavior.mat')

model_id = {'SST_stop_sucess_CPM_Result_permute','SST_stop_failure_CPM_Result_permute','SST_gowrong_CPM_Result_permute', ...               
            'BL_MID_feed_miss_CPM_Result_permute','MID_feed_hit_CPM_Result_permute','MID_antici_hit_CPM_Result_permute', ...
            'EFT_neutral_CPM_Result_permute','EFT_angry_CPM_Result_permute' };
        
mask_ind = tril(reshape(1:64,8,8),-1);
exter_mask = mask_ind(1:4,1:4); exter_mask = exter_mask(exter_mask~=0);
inter_mask = mask_ind(5:8,5:8); inter_mask = inter_mask(inter_mask~=0);
across_mask = mask_ind(5:8,1:4); across_mask = across_mask(across_mask~=0);

% permute edges

for k = 1:8  
    
    models = dir(fullfile(model_path,['*',model_name{k}]));

    for i=1:8
    
        net_pos_all = {}; R_pos = []; Predict_pos=[];
        net_neg_all = {}; R_neg = []; Predict_neg=[];
        clear model
        for j=1:6
        
           disp(['CogDimension',num2str(k),'---','Symptom',num2str(i),'--','Permute_component',num2str(j)])
           disp(models(number(i)+ (j-1)*8).name);
            model = load(fullfile(model_path,models(number(i)+ (j-1)*8).name));
           
            eval(['model_result = model.',model_id{k},';']);
            
            net_pos = model_result.pos_mask;
            net_neg = model_result.neg_mask;
            net_pos_all = [net_pos_all,net_pos];
            net_neg_all = [net_neg_all,net_neg];
            R_pos   = [R_pos;mean(model_result.permute_r_pos,2)];
            R_neg   = [R_neg;mean(model_result.permute_r_neg,2)];
            R_both  = [R_neg;mean(model_result.permute_r_both,2)];
            
            Predict_pos   = [Predict_pos,model_result.permute_predict_pos];
            Predict_neg   = [Predict_neg,model_result.permute_predict_neg];
            Predict_both  = [Predict_neg,model_result.permute_predict_both];
            
            model_name_id{j,i} = models(number(i)+ (j-1)*8).name;
        end 
        
        net_pos_permu{i} = net_pos_all;
        net_neg_permu{i} = net_neg_all;
        
        Perm_R_pos{i}  = R_pos;  Perm_Prdict_pos{i}  = Predict_pos; 
        Perm_R_neg{i}  = R_neg;  Perm_Prdict_neg{i}  = Predict_neg; 
        Perm_R_both{i} = R_both; Perm_Prdict_both{i} = Predict_both;
        
    end
    
    
   % overlap 
    p =0;
    for m = 1:8
         for n =1:8
%              if m~=n
             p = p+1;
            disp([m,n])
            
            net_1_pos = net_pos_permu{m};net_1_neg = net_neg_permu{m};
            net_2_pos = net_pos_permu{n};net_2_neg = net_neg_permu{n};
            
            pos_pos = cellfun(@(x,y) (intersect(x,y)),net_1_pos,net_2_pos,'un',0)';
            pos_neg = cellfun(@(x,y) (intersect(x,y)),net_1_pos,net_2_neg,'un',0)';
            neg_neg = cellfun(@(x,y) (intersect(x,y)),net_1_neg,net_2_neg,'un',0)';
            neg_pos = cellfun(@(x,y) (intersect(x,y)),net_1_neg,net_2_pos,'un',0)';
            
            pos_pos_all(:,p)= pos_pos';
            pos_neg_all(:,p)= pos_neg';
            neg_neg_all(:,p)= neg_neg';
            neg_pos_all(:,p)= neg_pos';
            
%             end
        end
    end
    
    % condition 
    
        Perm_R_pos_task{k}  = Perm_R_pos;  Perm_Prdict_pos_task{k}  = Perm_Prdict_pos; 
        Perm_R_neg_task{k}  = Perm_R_neg;  Perm_Prdict_neg_task{k}  = Perm_Prdict_neg; 
        Perm_R_both_task{k} = Perm_R_both; Perm_Prdict_both_task{k} = Perm_Prdict_both;
        Perm_Pos_edge{k}  = net_pos_permu;
        Perm_Neg_edge{k}  = net_neg_permu;
            
    for j = 1:1000
        disp([k,j])
        Exter_PosPos{j,k} = unique(cell2mat(pos_pos_all(j,exter_mask)));
        Exter_PosNeg{j,k} = unique(cell2mat(pos_neg_all(j,exter_mask)));  
        Exter_NegNeg{j,k} = unique(cell2mat(neg_neg_all(j,exter_mask)));
        Exter_NegPos{j,k} = unique(cell2mat(neg_pos_all(j,exter_mask)));
        
        Inter_PosPos{j,k} = unique(cell2mat(pos_pos_all(j,inter_mask)));
        Inter_PosNeg{j,k} = unique(cell2mat(pos_neg_all(j,inter_mask)));
        Inter_NegPos{j,k} = unique(cell2mat(neg_pos_all(j,inter_mask)));
        Inter_NegNeg{j,k} = unique(cell2mat(neg_neg_all(j,inter_mask)));
        
        Share_PosPos{j,k} = unique(cell2mat(pos_pos_all(j,across_mask)));
        Share_PosNeg{j,k} = unique(cell2mat(pos_neg_all(j,across_mask)));
        Share_NegPos{j,k} = unique(cell2mat(neg_pos_all(j,across_mask)));
        Share_NegNeg{j,k} = unique(cell2mat(neg_neg_all(j,across_mask)));
              
    end
    
end

% permute FC

TableS2_1_perm.Exter_PosPos = Exter_PosPos;
TableS2_1_perm.Exter_PosNeg = Exter_PosNeg;
TableS2_1_perm.Exter_NegNeg = Exter_NegNeg;
TableS2_1_perm.Exter_NegPos = Exter_NegPos;

TableS2_1_perm.Inter_PosPos = Inter_PosPos;
TableS2_1_perm.Inter_PosNeg = Inter_PosNeg;
TableS2_1_perm.Inter_NegNeg = Inter_NegNeg;
TableS2_1_perm.Inter_NegPos = Inter_NegPos;

TableS2_1_perm.Share_PosPos = Share_PosPos;
TableS2_1_perm.Share_PosNeg = Share_PosNeg;
TableS2_1_perm.Share_NegNeg = Share_NegNeg;
TableS2_1_perm.Share_NegPos = Share_NegPos;

TableS2_2_disorder.Perm_R_pos_task = Perm_R_pos_task;
TableS2_2_disorder.Perm_R_neg_task = Perm_R_neg_task;
TableS2_2_disorder.Perm_R_both_task = Perm_R_both_task;

TableS2_2_disorder.Perm_Prdict_pos_task = Perm_Prdict_pos_task;
TableS2_2_disorder.Perm_Prdict_neg_task = Perm_Prdict_neg_task;
TableS2_2_disorder.Perm_Prdict_both_task = Perm_Prdict_both_task;

save N2_1_2_make_Ex_In_factor_Permute.mat TableS2_1_perm TableS2_2_disorder -v7.3



































