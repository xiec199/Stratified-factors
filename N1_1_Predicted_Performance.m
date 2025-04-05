
% We use the prediction results obtained from the NP-factor project 
% % (https://github.com/xiec199/NP-factor) for subsequent analysis

% The analysis includes obtaining the results for parts a-c of Fig. 1:
% a. Predictive performance for each task connectome and each symptom
% b. Correlation matrix between behavioral symptoms
% c. Correlation matrix between brain-predicted symptoms

% For more details of predictive process: https://github.com/xiec199/NP-factor

clear; clc;

%% Load symptom data
[BL_symptom, FU2_symptom, FU3_symptom] = xic_imagen_symptoms;

%% Brain Predicted Behaviours

% Define model parameters
path = '/CPM_Models/';
model_names = {'_SST_stopsuc.mat','_SST_stopfai.mat', ...
               'SST_gowrong.mat', 'feed_miss.mat', 'feed_hit.mat', ...
               'antici_hit.mat', 'EFT_neutral.mat', 'EFT_angry.mat'};

model_indices = [4 2 5 7 3 1 6 8];
dis_names = {'asd', 'adhd', 'cd', 'od', 'anxiety', 'dep', 'eat', 'speph'};

%% Subject selection
subject = BL_symptom.BL_subject;
for j = 1:8
    models = dir(fullfile(path, ['*', model_names{j}]));
    for i = 1:length(models)
        brain = load(models(model_indices(i)).name);
        subject = intersect(brain.CPM_Result.subject, subject);
    end
end

%% Predict Data
Pheno_predict = cell(1, 8);
Pheno_performance = zeros(8, length(models));

for j = 1:8
    models = dir(fullfile(path, ['*', model_names{j}]));
    
    for i = 1:length(models)
        fprintf('Loading subject data %02d %02d \n', j, i);
        
        brain = load(models(model_indices(i)).name);
        
        Pheno_subject{i,j} = brain.CPM_Result.subject;
        Pheno_performance(i,j) = brain.CPM_Result.both_r_mean;
        
        % Find the intersection of subjects
        [~, ind1] = intersect(brain.CPM_Result.subject, subject);
        predict = mean(brain.CPM_Result.both_predit, 2);
        
        Pheno_predict{j}(:,i) = predict(ind1);
        Pheno_predict_dat{i,j} = predict(ind1);
    end
end

%% Brain Predict Behavior
all_predict = NaN(size(Pheno_predict{1}, 1), 8);
for i = 1:8
    predict_ind = Pheno_performance(i,:) > 0.075;
    predict_beha = Pheno_predict_dat(i, predict_ind);
    all_predict(:,i) = mean(cell2mat(predict_beha), 2);
end

%% Correlation and Residual Analysis
[~, ind1, ind2] = intersect(CPM_subject, BL_symptom.BL_subject);
all_predict = all_predict(ind1,:);
bl_symptom = BL_symptom.BL_symptom_residual(ind2,:);

% Calculate correlation
brain_r = corr(all_predict);

% Extract and clean correlation matrices
brain_r = tril(brain_r, -1);
brain_ex_r = brain_r(1:4, 1:4);
brain_ex_r(brain_ex_r == 0) = [];
brain_in_r = brain_r(5:8, 5:8);
brain_in_r(brain_in_r == 0) = [];

% Display means
fprintf('Mean of external brain correlations: %.3f\n', mean(brain_ex_r));
fprintf('Mean of internal brain correlations: %.3f\n', mean(brain_in_r));
