# Hierarchical Neurocognitive Model of Externalizing and Internalizing Comorbidity

This repository contains the code for our paper: "Hierarchical Neurocognitive Model of Externalizing and Internalizing Comorbidity," which identified distinct neural factors underlying externalizing and internalizing symptoms, alongside a previously identified general neuropsychopathological factor (A shared neural basis underlying psychiatric comorbidity, https://www.nature.com/articles/s41591-023-02317-4).

## Overview

Our study leveraged a brain-behavior predictive framework with a 10-year longitudinal imaging-genetic cohort (IMAGEN, ages 14, 19, and 23, N = 1,750) to construct two neural factors underlying externalizing and internalizing symptoms, which were reproduced across six clinical and population-based datasets (ABCD, STRATIFY/ESTRA, ABIDE II, ADHD-200, and XiNan, from age 10 to age 36, N = 3,765).

These two neural factors exhibit distinct neural configurations:
- **Externalizing symptoms**: Hyperconnectivity in impulsivity-related circuits
- **Internalizing symptoms**: Hypoconnectivity in goal-directed circuits

Both factors also differ in their cognitive-behavior relevance, genetic substrates, and developmental profiles.

## Software Requirements and Installation

### Prerequisites
- MATLAB R2018b or later
- Statistics and Machine Learning Toolbox

### Required Toolbox
- **CPM (Connectome-based Predictive Modeling) Toolbox**: https://github.com/YaleMRRC/CPM
  - Download or clone the CPM toolbox from the above repository
  - Add the CPM toolbox to your MATLAB path using `addpath('path/to/CPM')`

### Installation Steps
1. Clone this repository: `git clone https://github.com/xiec199/Stratified-factors`
2. Download and install the CPM toolbox (see link above)
3. Add both repositories to your MATLAB path
4. Ensure all required MATLAB toolboxes are installed

## Running Structure

The analysis scripts should be executed in the following order:

### Phase 1: Factor Construction and Analysis
Run scripts in numerical order within each category:

1. **Initial Factor Construction** (`N1_*.m`)
   - `N1_1_Predicted_Performance.m` → `N1_2_Make_Ex_In_factor.m`

2. **Permutation Testing** (`N2_*.m`)
   - `N2_1_2_make_Ex_In_factor_Permute.m` → `N2_1_4_make_Ex_In_factor_EachType_sum.m`

3. **Longitudinal Analysis** (`N3_*.m`)
   - `N3_1_1_make_Ex_In_factor_FC.m` → `N3_2_2_Ex_In_factor_FC_Longitudinal_change.m`

4. **Factor Specificity** (`N4_*.m`)
   - Run all `N4_*.m` scripts after completing N1-N3

### Phase 2: External Validation (`V1_*.m`)
These scripts can be run independently after Phase 1 completion:
- `V1_1_ABCD_FC.m`: ABCD dataset validation
- `V1_3_ADHD200_FC.m`: ADHD-200 dataset validation
- `V1_4_Xinan_FC.m`: Depression dataset validation
- `V1_5_ASD_FC.m`: ASD dataset validation
- `V1_6_Striatify_FC_SST.m`: STRATIFY dataset validation

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

## Data Requirements

Note: Due to data sharing agreements, the raw neuroimaging and behavioral data are not included in this repository. Users need to obtain access to the respective datasets (IMAGEN, ABCD, etc.) through their official channels.

## Contact

For questions about the code or methodology, please contact xiec199@gmail.com.
