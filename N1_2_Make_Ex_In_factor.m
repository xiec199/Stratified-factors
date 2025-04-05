% The analysis includes obtaining the results for parts d of Fig. 1:
% The goal is to construct externalizing and internalizing factors 
% based on shared and unique brain network connections (positive and negative). 

%% 
clear; clc; close all

% Define model paths
path_model = '/CPM_Models';

% Model names for the different tasks
model_name = {'_SST_stopsuc.mat', '_SST_stopfai.mat', 'SST_gowrong.mat', ...
    'feed_miss.mat', 'feed_hit.mat', 'antici_hit.mat', 'EFT_neutral.mat', 'EFT_angry.mat'};

% Task names for clarity in results
name = {'SST:stop success','SST:stop failure','SST:go wrong',...
       'MID:feedmiss','MID:feedhit','MID: antici hit', ...
       'EFT:neutral','EFT:angry'};

number = [4 2 5 7 3 1 6 8];
dis_names = {'asd','adhd','cd' ,'od','anxiety','dep','eat','speph'};

% Create masks for different regions and connections
mask_ind = tril(reshape(1:64,8,8), -1);
exter_mask = mask_ind(1:4, 1:4); exter_mask = exter_mask(exter_mask ~= 0);
inter_mask = mask_ind(5:8, 5:8); inter_mask = inter_mask(inter_mask ~= 0);
across_mask = mask_ind(5:8, 1:4); across_mask = across_mask(across_mask ~= 0);
mask_id_all = tril(reshape(1:268*268, 268, 268), -1);

% Loop through each task model
for k = 1:8

    models = dir(fullfile(path_model, ['*', model_name{k}]));

    % Iterate through different combinations of models
    for i = 1:8

        model1 = load(fullfile(path_model, models(number(i)).name));
        network1_neg = model1.CPM_Result.neg_mask;  
        network1_neg(network1_neg < 0.95) = 0; network1_neg(network1_neg > 0) = 1;
        network1_pos = model1.CPM_Result.pos_mask;  
        network1_pos(network1_pos < 0.95) = 0; network1_pos(network1_pos > 0) = 1;

        for j = 1:8
            disp([i, j])

            model2 = load(fullfile(path_model, models(number(j)).name));
            network2_neg = model2.CPM_Result.neg_mask; 
            network2_neg(network2_neg < 0.95) = 0; network2_neg(network2_neg > 0) = 1;
            network2_pos = model2.CPM_Result.pos_mask; 
            network2_pos(network2_pos < 0.95) = 0; network2_pos(network2_pos > 0) = 1;

            % Calculate overlapping and non-overlapping connections
            neg_pos = network1_neg .* network2_pos; 
            neg_pos_id = mask_id_all(neg_pos > 0); neg_pos_id(neg_pos_id == 0) = [];
            neg_neg = network1_neg .* network2_neg; 
            neg_neg_id = mask_id_all(neg_neg > 0); neg_neg_id(neg_neg_id == 0) = [];
            pos_pos = network1_pos .* network2_pos; 
            pos_pos_id = mask_id_all(pos_pos > 0); pos_pos_id(pos_pos_id == 0) = [];
            pos_neg = network1_pos .* network2_neg; 
            pos_neg_id = mask_id_all(pos_neg > 0); pos_neg_id(pos_neg_id == 0) = [];

            % Save results for later analysis
            neg_pos_all{i, j} = neg_pos;    
            neg_pos_all_id{i, j} = neg_pos_id;  
            neg_neg_all{i, j} = neg_neg;    
            neg_neg_all_id{i, j} = neg_neg_id; 
            pos_neg_all{i, j} = pos_neg;    
            pos_neg_all_id{i, j} = pos_neg_id;
            pos_pos_all{i, j} = pos_pos;    
            pos_pos_all_id{i, j} = pos_pos_id;

            % Count connections for each type
            neg_neg_num(i, j) = sum(neg_neg(:)) / 2;
            pos_pos_num(i, j) = sum(pos_pos(:)) / 2;
            pos_neg_num(i, j) = sum(pos_neg(:)) / 2;
            neg_pos_num(i, j) = sum(neg_pos(:)) / 2;

            % Calculate relationship between opposite and same-site connections
            opposite = pos_neg_num(i, j) + neg_pos_num(i, j);
            samesite = neg_neg_num(i, j) + pos_pos_num(i, j);
            all = (opposite + samesite);
        end
    end

    disp(['Overlapped edges'])

    % Save data for the different dimensions
    neg_pos_all_dimension_id{k} = neg_pos_all_id;     
    neg_neg_all_dimension_id{k} = neg_neg_all_id;    
    pos_neg_all_dimension_id{k} = pos_neg_all_id;     
    pos_pos_all_dimension_id{k} = pos_pos_all_id;    

    neg_pos_all_dimension{k} = neg_pos_all;    
    neg_neg_all_dimension{k} = neg_neg_all;    
    pos_neg_all_dimension{k} = pos_neg_all;    
    pos_pos_all_dimension{k} = pos_pos_all; 

    % Define externalizing and internalizing factors based on shared and unique connections
    share_neg_pos_all_dimension_id{k} = unique_cell(neg_pos_all_id(across_mask)); 
    share_neg_neg_all_dimension_id{k} = unique_cell(neg_neg_all_id(across_mask));  
    share_pos_neg_all_dimension_id{k} = unique_cell(pos_neg_all_id(across_mask));   
    share_pos_pos_all_dimension_id{k} = unique_cell(pos_pos_all_id(across_mask)); 

    exter_neg_pos_all_dimension_id{k} = unique_cell(neg_pos_all_id(exter_mask));  
    exter_neg_neg_all_dimension_id{k} = unique_cell(neg_neg_all_id(exter_mask)); 
    exter_pos_neg_all_dimension_id{k} = unique_cell(pos_neg_all_id(exter_mask));  
    exter_pos_pos_all_dimension_id{k} = unique_cell(pos_pos_all_id(exter_mask)); 

    inter_neg_pos_all_dimension_id{k} = unique_cell(neg_pos_all_id(inter_mask));  
    inter_neg_neg_all_dimension_id{k} = unique_cell(neg_neg_all_id(inter_mask));
    inter_pos_neg_all_dimension_id{k} = unique_cell(pos_neg_all_id(inter_mask));
    inter_pos_pos_all_dimension_id{k} = unique_cell(pos_pos_all_id(inter_mask));   
end

% Prepare results for externalizing and internalizing factor construction
TableS1_1 = table;
TableS1_1.name = name';

% Externalizing and Internalizing factor results
TableS1_1.share_pos_pos_all_dimension_id = share_pos_pos_all_dimension_id';
TableS1_1.share_pos_neg_all_dimension_id = share_pos_neg_all_dimension_id';
TableS1_1.share_neg_neg_all_dimension_id = share_neg_neg_all_dimension_id';
TableS1_1.share_neg_pos_all_dimension_id = share_neg_pos_all_dimension_id';

TableS1_1.exter_pos_pos_all_dimension_id = exter_pos_pos_all_dimension_id';
TableS1_1.exter_pos_neg_all_dimension_id = exter_pos_neg_all_dimension_id';
TableS1_1.exter_neg_neg_all_dimension_id = exter_neg_neg_all_dimension_id';
TableS1_1.exter_neg_pos_all_dimension_id = exter_neg_pos_all_dimension_id';

TableS1_1.inter_pos_pos_all_dimension_id = inter_pos_pos_all_dimension_id';
TableS1_1.inter_pos_neg_all_dimension_id = inter_pos_neg_all_dimension_id';
TableS1_1.inter_neg_neg_all_dimension_id = inter_neg_neg_all_dimension_id';
TableS1_1.inter_neg_pos_all_dimension_id = inter_neg_pos_all_dimension_id';

% Save the data tables
TableS1_2 = table;
TableS1_2.name = name';

TableS1_2.pos_pos_all_dimension = pos_pos_all_dimension';
TableS1_2.pos_neg_all_dimension = pos_neg_all_dimension';
TableS1_2.neg_neg_all_dimension = neg_neg_all_dimension';
TableS1_2.neg_pos_all_dimension = neg_pos_all_dimension';

TableS1_2.pos_pos_all_dimension_id = pos_pos_all_dimension_id';
TableS1_2.pos_neg_all_dimension_id = pos_neg_all_dimension_id';
TableS1_2.neg_neg_all_dimension_id = neg_neg_all_dimension_id';
TableS1_2.neg_pos_all_dimension_id = neg_pos_all_dimension_id';

% Save final output
save N1_1_1_make_Ex_In_factor_part1 TableS1_2 TableS1_1
