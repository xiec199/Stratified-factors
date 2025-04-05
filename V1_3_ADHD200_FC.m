%% Analysis of Externalizing and Internalizing Factors in ADHD-200 Dataset

clear; clc;
load('N2_1_5_make_Ex_In_factor_eachFC.mat')
ADHD_dat = load('ADHD_FC_Map_268.mat');

%% Extract functional connectivity data
exter = exter_pospos_type;
inter = inter_negneg_type;
adhd_fc = reshape(ADHD_dat.fc_rest, [], size(ADHD_dat.fc_rest, 3));

% Calculate mean FC values for externalizing and internalizing factors
for i = 1:4
    exter_fc_rest(:, i) = mean(adhd_fc(exter{i}, :))';
    inter_fc_rest(:, i) = mean(adhd_fc(inter{i}, :))';
end

%% Load phenotypic data
pheno = importdata('adhd200_preprocessed_phenotypics.tsv.xlsx');
pheno_dat = pheno.data(:, 6);           % Diagnosis status
pheno_sec = pheno.textdata(2:end, 7);   % Secondary diagnosis
pheno_id = pheno.data(:, 1);            % Subject ID
pheno_age = pheno.data(:, 4);           % Age
pheno_site = pheno.data(:, 2);          % Acquisition site
pheno_gender = pheno.data(:, 3);        % Gender

%% Prepare data for analysis
fc_sub = str2double(string(ADHD_dat.fc_sub_name));
FC_mat = inter_fc_rest;                 % Using internalizing FC for analysis

% Match subjects between FC and phenotypic data
[all, ind1, ind2] = intersect(fc_sub, pheno_id);
ph_mat = pheno_dat(ind2, :);
ph_sec = pheno_sec(ind2, :);
ph_age = pheno_age(ind2, :);

% Create confound matrix for regression
confound = [dummyvar(pheno_site(ind2, :)), pheno_gender(ind2), pheno_age(ind2)];
confound = [ones(length(confound), 1), confound(:, 2:end)];

% Calculate mean FC and regress out confounds
FC_mats = mean(FC_mat(ind1, :), 2);
[~, ~, FC_mats] = regress(FC_mats, confound);

%% Define diagnostic groups
ph_sec_types = ph_sec(~cellfun(@isempty, ph_sec));

% Identify internalizing diagnosis subjects
inter = ph_mat ~= 0 & (contains(ph_sec, 'depres') + contains(ph_sec, 'anxiety') + ...
    contains(ph_sec, 'Depres') + contains(ph_sec, 'Anxie') + ...
    contains(ph_sec, 'phobia') + contains(ph_sec, 'Phobia') ~= 0);

% Identify externalizing diagnosis subjects
exter = ph_mat ~= 0 & contains(ph_sec, 'ODD');
exter(inter == 1) = 0;

% Identify control and ADHD subjects
control = ph_mat == 0;
ADHD = ph_mat ~= 0;

%% Conduct group comparisons
[~, p1, ~, t1] = ttest2(FC_mats(ADHD, :), FC_mats(control, :));
d1 = cohen_t2(t1.tstat, sum(ADHD), sum(control));

[~, p2, ~, t2] = ttest2(FC_mats(exter, :), FC_mats(control, :));
[~, p3, ~, t3] = ttest2(FC_mats(inter, :), FC_mats(control, :));

%% Compile results
Tvalues = table;
Tvalues.ADHD_Tvlaue = t1.tstat';
Tvalues.ADHD_Pvalue = p1/2';
Tvalues.ADHD_Dvalue = d1;
Tvalues.Exter_Tvlaue = t2.tstat';
Tvalues.Exter_Pvalue = p2/2';
Tvalues.Inter_Tvlaue = t3.tstat';
Tvalues.Inter_Pvalue = p3/2';

%% Save results
save('ADHD200_FC_Analysis.mat', 'Tvalues')