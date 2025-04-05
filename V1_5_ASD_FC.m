%% Meta-analysis of ABIDE II Dataset for ASD vs Control Comparison

clear; clc;

%% Load data
load('ABIDEII_FC_Map_268_DPABI.mat')
load('N2_1_5_make_Ex_In_factor_eachFC.mat')
load('asd2.mat')
fc_sub = split(fc_sub_name, '_');

%% Prepare FC data
fc_rest = reshape(fc_rest, [], size(fc_rest, 3));
exter = exter_pospos_type;
inter = inter_negneg_type;

%% Match subjects
[~, ind1, ind2] = intersect(asd2.dat_id, fc_sub(:, 3));
fc_rest = fc_rest(:, ind2);
sub_site = asd2.site(ind1);
sub_uniSite = unique(sub_site);
sub_dig = str2double(string(asd2.dig(ind1)));
sub_age = str2double(string(asd2.age(ind1)));

%% Analyze by site
for i = 1:length(sub_uniSite)
    % Identify subjects from current site
    con_site_ind = contains(sub_site, sub_uniSite{i}) & sub_dig == 2;
    asd_site_ind = contains(sub_site, sub_uniSite{i}) & sub_dig == 1;
    
    % Extract age information
    asd_age = sub_age(asd_site_ind);
    con_age = sub_age(con_site_ind);
    
    % Extract FC data
    con_site = fc_rest(:, con_site_ind);
    asd_site = fc_rest(:, asd_site_ind);
    
    % Calculate mean FC for externalizing factors
    clear asd_fc con_fc
    for j = 1:4
        asd_fc(:, i) = nanmean(asd_site(exter{j}, :))';
        con_fc(:, i) = nanmean(con_site(exter{j}, :))';
    end
    
    % Perform statistical comparison
    asd_fc_mean = mean(asd_fc, 2);
    con_fc_mean = mean(con_fc, 2);
    [h, p, c, t] = ttest2(asd_fc_mean, con_fc_mean);
    
    % Store results
    Tvalues(i, :) = t.tstat;
    Pvalues(i, :) = p;
    Age_con(i, 1) = nanmean(con_age);
    Age_asd(i, 1) = nanmean(asd_age);
    Nums(i, 1) = size(asd_site, 2);
    Nums(i, 2) = size(con_site, 2);
end

%% Exclude sites with insufficient subjects or older age
ex_ind = Nums(:, 2) < 20 | Nums(:, 1) < 20 | Age_con > 15 | Age_asd > 15;
sub_uniSite(ex_ind) = [];
Tvalues(ex_ind) = [];
Pvalues(ex_ind) = [];
Age_con(ex_ind) = [];
Age_asd(ex_ind) = [];
Nums(ex_ind, :) = [];

%% Perform meta-analysis
clear r z
for j = 1:size(Tvalues, 2)
    for i = 1:length(Nums)
        r = xic_t2r(Tvalues(i, j), (sum(Nums(i, :), 2) - 2));
        z(i, 1) = xic_fisherz(r);
    end
    
    % Calculate meta-analysis Z and p-values
    z_meta = sum(z .* (sum(Nums, 2) - 3)) / sqrt(sum(sum(Nums, 2) - 3));
    p_meta = (1 - cdf('Normal', abs(z_meta), 0, 1)) / 2;
    
    FC_meta.FC_metaZ(j, 1) = z_meta;
    FC_meta.FC_metaP(j, 1) = p_meta;
end

%% Calculate effect size
t_value = norminv(normcdf(FC_meta.FC_metaZ));
Cohen_d = cohen_t(t_value, 233 + 331);

%% Save results
save('ABIDEII_ASD_Meta_Analysis.mat', 'FC_meta', 'Cohen_d', 'sub_uniSite', 'Nums', 'Age_con', 'Age_asd');