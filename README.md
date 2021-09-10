# T2D_study
Diabetes Mellitus Type II sub-phenotyping and trajectory

## Files description

Code for the study on T2D clustering and trajectory modeling via MCMC organised as follow:
**1. File functions**
- functions.R: functions used across all files

**2. Folder patients_selection **
- patient_selection.R: main file processing of genes and hes data and for the pre-selection of patients for clustering and MCMC
- msm_selection.R: additional selection steps for multi-state models

**3. Folder clustering**
- diabetes_clustering.R: main file for T2D cluster analysis

**4. Folder MCMC**
- mcmc_v2.R, mcmc_v3.R, mcmc_v4.R, mcmc_v5.R, mcmc_v6.R: variables preparations for MCMC procedures studied and detailed below
- MCMC_recap.R: summary of results from MCMC procedures


## MCMC models description

**MCMC_V2**
- Model: PD -> PD/D/H
- Mathematical model: PD -> PD/(D+H) -> D/H
- Variables: All but LDL and DBP
- Scaling: Yes
- Parameters: 8 + 8 (all selected variables in both transitions)

**MCMC_V3**
- Model: PD -> PD/D/H
- Mathematical model: PD -> PD/(D+H) -> D/H
- Variables: All but LDL and DBP
- Scaling: Yes
- Parameters: 1 + 8 (only mu for first transition, all selected variables for second)

**MCMC_V4**
- Model: PD -> PD/D/H
- Mathematical model: PD -> H/(PD+D) -> PD/D
- Variables: All but LDL and DBP
- Scaling: Yes
- Parameters: 8 + 8 (all selected variables in both transitions)

**MCMC_V5**
- Model: D -> D/CVD
- Mathematical model: D/CVD
- Variables: Cluster membership based on T2D clusters from Clustering Model 2
- Scaling: Yes
- Parameters: 4

**MCMC_V6**
- Model: D -> D/CVD
- Mathematical model: D/CVD
- Variables: Cluster membership based on T2D clusters from Clustering Model 3
- Scaling: Yes
- Parameters: 5
