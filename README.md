# Context-of-Restoration

This repository conatisn script and data used in analyses related to the Context-dependent benefits of coral reef restoration study.

For this study, the fast version of ReefMod (Bozec et al. 2022), otherwise referred to as ReefMod Engine has been used to run various single reef customised intervention scenarios in a short time frame. The version of the ReefMod Engine used here is rme version 1.0.28 (folder of the model containing manual is provided)

Before running the scenarios using the Engine scenario, range of scenarios options need to be defined using parameters defined in Table 1 of the manuscript.

## Workflow
### Create input data

1) Script that set the intervention options= Inter_reefstate_environmentla_paper_script_030525.mat
Output files created using this script are:
- Spec_sub_1cyc_3init_3self_4larv_paper_test_150124.mat (folder 1cyc_1int_240124)
- Spec_sub_1dhw_3init_3self_4larv_paper_test_150124.mat (folder 1dhw_1int_240124)

Here, every row represents a unique scenario (2160 in total)
col 1= model start year
col 2= model end year
col 3= index initial coral cover (1-3)
col 4= index species composition (1-2)
col 5= initial cots setting (not used but important for model initiation)
col 6= initial rubble cover (%)
col 7= surface sedimentation (10,90 percentile of inshore reefs)
col 8= mid depth sedimentation (10,90 percentile of inshore reefs)
col 9= self retention rate
col 10= external larval supply
col 11= intervention year
col 12= disturbance indes (1= cyclone, 2= bleaching)
col 13-28= disturbance sequence from 2022-2037

2) Once option script has been produced, we can create the customised input data files to run the specific scenarios (more information in the "special_run.m" manual)= Preparation_syst_model_1dist_1inter_paper_030525.mat

3) in addition to the datafiles, we need to specify the "config" file using the file names of the input data files
Config files used corresponding to the specific intervention (cyclone or bleaching event) can be found in the 1cyc_1int_240124 or 1dhw_1int_240124 folder

### Run model
Once input data files have been created, we can run the intervention scenarios using script "HPC_Engine_intervention_model_1cyc_dhw_080824.mat"
This scripts describes activiting the model, creating an parameterising a scenario, set up intervention, model initialisation, model process and extraction of data
Model scenarios were run on the HPC using job array configuration (use matlab script and Batch script) but information is provided to run the code in matlab local computer

Output data is within reef and at reef level. Note in the output Cell_final, other data is present SITE/NOSITE but is not been used for this study.
Model output data can be provided upon request

### Extraction of desired model output data
For each individual model run, we extracted various ecological, connectivity, intervention and disturbance related data (see script HPC_Engine_intervention_model_1cyc_dhw_080824 for more information).

### Compile extracted output data and prep for downstream analyses
Every output data file is compiled for cyclone and bleaching parameterised scenarios seperately. Additionally, data is normalised and subset for further downstream analyses
(see script Compile_extracted_outputs_inter_scen_040525)

Following files can be requested upon request due to upload size limitations:
Compiled and normalised data for scenarios with no disturbance=Extract_final_OC_VR_1cyc_1int_3LS_nodist_norm_3SD_3PAR_170824.txt
Compiled and normalised data for scenarios with cyclone=Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.txt
Compiled and normalised data for scenarios with bleaching event=Extract_final_OC_VR_1dhw_1int_3LS_5_10_20_200824_cor_norm.txt

### Analysis of variable importance 
Variable importance is assessed for no disturbance/only cyclone or only bleaching scenarios in R. Additionally, script provides visualisation of variable impontance (Figure 3 and 5) analysis of pdp plots (Figure 4 and S3)
Script=Rscript_nodist_pdp_dist_randomforest_var_imp_040525.R

Following files can be requested upon request due to upload size limitations --> Randomforest models of the different disturbance scenario analyses.
Randomforest model for no disturbance scenarios=Nodist_cyc_nocorr_rf_3selfret_3larval_3stockdens_3PAR_norm_suppl_noHT_170824.R
Randomforest model for cyclone scenarios=Nodist_cyc_nocorr_rf_3selfret_3larval_3stockdens_3PAR_norm_suppl_noHT_170824.R
Randomforest model for bleaching event scenarios=Nodist_cyc_nocorr_rf_3selfret_3larval_3stockdens_3PAR_norm_suppl_noHT_170824.R

### Visualisation of results
Visualisation of remaining figures in main manuscript and supps
Script=Visualisation_results_figures_paper_040525.mat
