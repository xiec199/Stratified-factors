# Hierarchical Neurocognitive Model of Externalizing and Internalizing Comorbidity

This repository contains the code for our paper: "Hierarchical Neurocognitive Model of Externalizing and Internalizing Comorbidity," which identified distinct neural factors underlying externalizing and internalizing symptoms, alongside a previously identified general neuropsychopathological factor.

## Overview

Our study leveraged a brain-behavior predictive framework with a 10-year longitudinal imaging-genetic cohort (IMAGEN, ages 14, 19, and 23, N = 1,750) to construct two neural factors underlying externalizing and internalizing symptoms, which were reproduced across six clinical and population-based datasets (ABCD, STRATIFY/ESTRA, ABIDE II, ADHD-200, and XiNan, from age 10 to age 36, N = 3,765).

These two neural factors exhibit distinct neural configurations:
- **Externalizing symptoms**: Hyperconnectivity in impulsivity-related circuits
- **Internalizing symptoms**: Hypoconnectivity in goal-directed circuits

Both factors also differ in their cognitive-behavior relevance, genetic substrates, and developmental profiles.

## Repository Structure

### Factor Construction and Analysis
- `N1_*.m`: Initial factor prediction performance and construction
- `N2_*.m`: Factor permutation testing and reliability analysis
- `N3_*.m`: Longitudinal analysis of factors
- `N4_*.m`: Factor specificity analysis and visualization

### Validation Across Datasets
- `V1_*.m`: Replication in external datasets (ABCD, ADHD-200, ASD, Depression, STRATIFY)

## Key Components

1. **Factor Construction**
   - `N1_1_Predicted_Performance.m`: Analyzes prediction performance for behavioral symptoms
   - `N1_2_Make_Ex_In_factor.m`: Constructs externalizing and internalizing factors

2. **Permutation Testing**
   - `N2_1_2_make_Ex_In_factor_Permute.m`: Permutation tests to validate factors
   - `N2_1_4_make_Ex_In_factor_EachType_sum.m`: Analysis of specific edge types

3. **Longitudinal Analysis**
   - `N3_1_1_make_Ex_In_factor_FC.m`: Factor stability analysis over time
   - `N3_2_2_Ex_In_factor_FC_Longitudinal_change.m`: Developmental trajectories

4. **External Validation**
   - `V1_1_ABCD_FC.m`: Replication in ABCD dataset
   - `V1_3_ADHD200_FC.m`: Analysis in ADHD-200 dataset
   - `V1_4_Xinan_FC.m`: Analysis in depression dataset
   - `V1_5_ASD_FC.m`: Analysis in ASD dataset
   - `V1_6_Striatify_FC_SST.m`: Analysis in STRATIFY dataset

