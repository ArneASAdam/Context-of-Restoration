% Script that compiles all extracted data csv to add columnnames and
% prepare data for analysis in R and figure visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for cyclone and bleaching scenarios go to the directory where the folder
% is that contains all extracted data files

%directory cyclones='VR_1dist_1int_cyc_080824'
%directory bleaching events='VR_1dist_1int_dhw_080824'
folder='Extract_080824_5_10_20'

extension='csv';
scenarios=dir([folder filesep '**' filesep '*extract.' extension]);
scenarios

cd Extract_080824_5_10_20

for i=1:36
    int_scen=scenarios(i).name;
    extract1=readmatrix(int_scen);
    if i==1
        extract_final=extract1;
    else
        extract_final=[extract_final;extract1];
    end
    i;
end

size(extract_final)
unique(extract_final(i,12))
extract_final2=extract_final;

% Edit no disturbance in file to zero as those are currently NA values

%for cyclone scenarios --> cyclone=1
%for bleaching scenarios --> cyclone=0
cyclone=1;
if cyclone==1
for i=1:size(extract_final,1)
    if isnan(extract_final(i,82))
        extract_final(i,82)=0;
    end
    if isnan(extract_final(i,83))
        extract_final(i,83)=1;
    end
    startyr=2023;
    censusyr=extract_final(i,102)+2023;

    if isnan(extract_final(i,84))
        extract_final(i,84)=(censusyr-startyr)+1;
    end
    if isnan(extract_final(i,85))
        last_int_before_censusyr=extract_final(i,79);
        extract_final(i,85)=last_int_before_censusyr+1; 
    end
    if isnan(extract_final(i,86))
        model_period=2022:1:2037;
        startyr_ind=find(model_period==startyr);
        censusyr_ind=find(model_period==censusyr);
        extract_final(i,86)=((censusyr_ind-last_int_before_censusyr)-startyr_ind)+1;
    end
    if isnan(extract_final(i,87))
        extract_final(i,87)=(censusyr-startyr)+1;
    end
end

else

for i=1:size(extract_final,1)
    if isnan(extract_final(i,91))

        extract_final(i,91)=0;
    end
    if isnan(extract_final(i,92))

        extract_final(i,92)=0;
    end
    if isnan(extract_final(i,93))

        extract_final(i,93)=2;
    end
    startyr=2023;
    censusyr=extract_final(i,102)+2023;


    if isnan(extract_final(i,94))
        extract_final(i,94)=(censusyr-startyr)+1;
    end

    if isnan(extract_final(i,95))
        last_int_before_censusyr=extract_final(i,79);
        extract_final(i,95)=last_int_before_censusyr+1;
    end
    if isnan(extract_final(i,96))
        model_period=2022:1:2037;
        startyr_ind=find(model_period==startyr);
        censusyr_ind=find(model_period==censusyr);
        extract_final(i,96)=((censusyr_ind-last_int_before_censusyr)-startyr_ind)+1;
    end
    if isnan(extract_final(i,97))
        extract_final(i,97)=(censusyr-startyr)+1;
    end
end
end

%add normalised data
% will need to be redone as data of 4 larval supply metrics are included.
% This is for the purpose of extending the matrix to 212 columns
extract_final(:,107:212)=normalize(extract_final(:,1:106));

[rows,cols]=size(extract_final);
tablesize=[rows,cols];

%with demographics including site data
varnames={'SSP','GCM','startyr','censusyr','reef index','self_ret','tot_larv_in_param_index',...
    'rubble_cover_bef_startyr','mean_larv_in_yr_bef_startyr','coral_cover_startyr','size_dist_index','species_comb_index','cover_prop_sp1_bef_styr','cover_prop_sp2_bef_styr',...
    'cover_prop_sp3_bef_styr','cover_prop_sp4_bef_styr','cover_prop_sp5_bef_styr','cover_prop_sp6_bef_styr','nb_juv_bef_startyr','nb_recr_bef_startyr','mean_nb_adol_bef_startyr_100m2','mean_nb_adult_bef_startyr_100m2'...
    'coral_cover_ref_censusyr_mean','coral_cover_ref_censusyr_stdev','coral_cover_int_censusyr_mean','coral_cover_int_censusyr_stdev','coral_cover_int_ref_censusyr_mean','coral_cover_int_ref_censusyr_stdev',...
    'coral_cover_nat_ref_censusyr_mean','coral_cover_nat_ref_censusyr_stdev','coral_cover_nat_int_censusyr_mean','coral_cover_nat_int_censusyr_stdev','coral_cover_outpl_int_censusyr_mean','coral_cover_outpl_int_censusyr_stdev','incr_rel_shelt_censusyr','incr_cover_prop_sp1_censusyr',...
    'incr_cover_prop_sp2_censusyr','incr_cover_prop_sp3_censusyr','incr_cover_prop_sp4_censusyr','incr_cover_prop_sp5_censusyr','incr_cover_prop_sp6_cenyr',...
    'mean_inc_larv_in_startyr_censusyr','mean_inc_larv_in_min5_censusyr_censusyr','mean_inc_rubble_min5_censusyr_censusyr',...
    'juv_reef_nat_ref_censusyr_100m2','adol_reef_nat_ref_censusyr_100m2','adult_reef_nat_ref_censusyr_100m2','juv_reef_nat_int_censusyr_100m2','adol_reef_nat_int_censusyr_100m2','adult_reef_nat_int_censusyr_100m2','juv_reef_outpl_int_censusyr_100m2','adol_reef_outpl_int_censusyr_100m2','adult_reef_outpl_int_censusyr_100m2',...
    'cover_site_ref_mean','cover_site_ref_stdev','cover_site_int_mean','cover_site_int_stdev','cover_site_int_ref_mean','cover_site_int_ref_stdev','cover_site_nat_ref_mean','cover_site_nat_ref_stdev','cover_site_nat_int_mean','cover_site_nat_int_stdev','cover_site_outpl_int_mean','cover_site_outpl_int_stdev',...
    'juv_site_nat_ref_censusyr_total_site','adol_site_nat_ref_censusyr_total_site','adult_site_nat_ref_censusyr_total_site','juv_site_nat_int_censusyr_total_site','adol_site_nat_int_censusyr_total_site','adult_site_nat_int_censusyr_total_site','juv_site_outpl_int_censusyr_total_site','adol_site_outpl_int_censusyr_total_site','adult_site_outpl_int_censusyr_total_site',...
    'prop_rest','stockdens','heattol',...
    'years_last_rest_rel_to_startyr','years_last_rest_rel_to_censusyr','freq_rest_events','tot_cycl_till_startyr','mean_magn_cycl_startyr_censusyr',...
    'last_cyc_bef_first_int','first_cyc_aft_first_int','first_cyc_aft_last_int','last_cyc_before_last_int','last_cyc_before_censusyr',...
    'cover_loss_cycl_ref','cover_loss_cycl_int','cover_loss_cycl_int_ref',...
    'tot_bleaching_till_startyr','mean_magn_bleaching_startyr_censusyr',...
    'last_bleaching_bef_first_int','first_bleaching_aft_first_int','first_bleaching_aft_last_int','last_bleaching_before_last_int','last_bleaching_before_censusyr',...
    'cover_loss_dhw_ref','cover_loss_dhw_int','cover_loss_dhw_int_ref',...
    'mean_numb_outplants_m2','years_postdeploy','int_year','species_comp_index','heattol_check','surface_WQ',...
    'zSSP','zGCM','zstartyr','zcensusyr','zreef index','zself_ret','ztot_larv_in_param_index',...
    'zrubble_cover_bef_startyr','zmean_larv_in_yr_bef_startyr','zcoral_cover_startyr','zsize_dist_index','zspecies_comb_index','zcover_prop_sp1_bef_styr','zcover_prop_sp2_bef_styr',...
    'zcover_prop_sp3_bef_styr','zcover_prop_sp4_bef_styr','zcover_prop_sp5_bef_styr','zcover_prop_sp6_bef_styr','znb_juv_bef_startyr','znb_recr_bef_startyr','zmean_nb_adol_bef_startyr_100m2','zmean_nb_adult_bef_startyr_100m2'...
    'zcoral_cover_ref_censusyr_mean','zcoral_cover_ref_censusyr_stdev','zcoral_cover_int_censusyr_mean','zcoral_cover_int_censusyr_stdev','zcoral_cover_int_ref_censusyr_mean','zcoral_cover_int_ref_censusyr_stdev',...
    'zcoral_cover_nat_ref_censusyr_mean','zcoral_cover_nat_ref_censusyr_stdev','zcoral_cover_nat_int_censusyr_mean','zcoral_cover_nat_int_censusyr_stdev','zcoral_cover_outpl_int_censusyr_mean','zcoral_cover_outpl_int_censusyr_stdev','zincr_rel_shelt_censusyr','zincr_cover_prop_sp1_censusyr',...
    'zincr_cover_prop_sp2_censusyr','zincr_cover_prop_sp3_censusyr','zincr_cover_prop_sp4_censusyr','zincr_cover_prop_sp5_censusyr','zincr_cover_prop_sp6_cenyr',...
    'zmean_inc_larv_in_startyr_censusyr','zmean_inc_larv_in_min5_censusyr_censusyr','zmean_inc_rubble_min5_censusyr_censusyr',...
    'zjuv_reef_nat_ref_censusyr_100m2','zadol_reef_nat_ref_censusyr_100m2','zadult_reef_nat_ref_censusyr_100m2','zjuv_reef_nat_int_censusyr_100m2','zadol_reef_nat_int_censusyr_100m2','zadult_reef_nat_int_censusyr_100m2','zjuv_reef_outpl_int_censusyr_100m2','zadol_reef_outpl_int_censusyr_100m2','zadult_reef_outpl_int_censusyr_100m2',...
    'zcover_site_ref_mean','zcover_site_ref_stdev','zcover_site_int_mean','zcover_site_int_stdev','zcover_site_int_ref_mean','zcover_site_int_ref_stdev','zcover_site_nat_ref_mean','zcover_site_nat_ref_stdev','zcover_site_nat_int_mean','zcover_site_nat_int_stdev','zcover_site_outpl_int_mean','zcover_site_outpl_int_stdev',...
    'zjuv_site_nat_ref_censusyr_total_site','zadol_site_nat_ref_censusyr_total_site','zadult_site_nat_ref_censusyr_total_site','zjuv_site_nat_int_censusyr_total_site','zadol_site_nat_int_censusyr_total_site','zadult_site_nat_int_censusyr_total_site','zjuv_site_outpl_int_censusyr_total_site','zadol_site_outpl_int_censusyr_total_site','zadult_site_outpl_int_censusyr_total_site',...
    'zprop_rest','zstockdens','zheattol',...
    'zyears_last_rest_rel_to_startyr','zyears_last_rest_rel_to_censusyr','zfreq_rest_events','ztot_cycl_till_startyr','zmean_magn_cycl_startyr_censusyr',...
    'zlast_cyc_bef_first_int','zfirst_cyc_aft_first_int','zfirst_cyc_aft_last_int','zlast_cyc_before_last_int','zlast_cyc_before_censusyr',...
    'zcover_loss_cycl_ref','zcover_loss_cycl_int','zcover_loss_cycl_int_ref',...
    'ztot_bleaching_till_startyr','zmean_magn_bleaching_startyr_censusyr',...
    'zlast_bleaching_bef_first_int','zfirst_bleaching_aft_first_int','zfirst_bleaching_aft_last_int','zlast_bleaching_before_last_int','zlast_bleaching_before_censusyr',...
    'zcover_loss_dhw_ref','zcover_loss_dhw_int','zcover_loss_dhw_int_ref',...
    'zmean_numb_outplants_m2','zyears_postdeploy','zint_year','zspecies_comp_index','zheattol_check','zsurface_WQ'}


vartypes={'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double'...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double'...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double','double','double','double','double','double','double'};

Extract_final_OC_VR_1cyc_1int_170824=table('Size',tablesize,'VariableTypes',vartypes,'VariableNames',varnames);
Extract_final_OC_VR_1cyc_1int_25_5_10_170824=table('Size',tablesize,'VariableTypes',vartypes,'VariableNames',varnames);
Extract_final_OC_VR_1cyc_1int_25_5_10_20_170824=table('Size',tablesize,'VariableTypes',vartypes,'VariableNames',varnames);
Extract_final_OC_VR_1dhw_1int_5_10_20_200824=table('Size',tablesize,'VariableTypes',vartypes,'VariableNames',varnames);
Extract_final_OC_VR_1dhw_1int_5_10_20_191124=table('Size',tablesize,'VariableTypes',vartypes,'VariableNames',varnames);


% PERCENTRANK = @(YourArray, TheProbes) reshape( mean( bsxfun(@le, YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) ); %This program determines the percentile of data popint X in vector Y

Extract_final_table=array2table(extract_final); % converts to a table

Extract_final_OC_VR_1cyc_1int_170824(:,:)=Extract_final_table(:,:);
Extract_final_OC_VR_1dhw_1int_5_10_20_200824(:,:)=Extract_final_table(:,:);
Extract_final_OC_VR_1dhw_1int_5_10_20_191124(:,:)=Extract_final_table(:,:);

% Below files will be used for further analyses
save('VR_1cyc_1int212col_36scen_prenorm_norm_all_5_10_20_170824.mat','Extract_final_OC_VR_1cyc_1int_170824','-v7.3') % for figures and etc.
save('VR_1dhw_1int212col_36scen_prenorm_norm_all_5_10_20_200824.mat','Extract_final_OC_VR_1dhw_1int_5_10_20_200824','-v7.3') % for figures and etc.

%first we normalise the data of scenarios
%to prepare the data for statistical analyses of data with no disturbances
%and measure variable importance for scenarios with disturbances

load('VR_1cyc_1int212col_36scen_prenorm_norm_all_5_10_20_170824.mat') % only if necessary
load('VR_1dhw_1int212col_36scen_prenorm_norm_all_5_10_20_200824.mat') % only if necessary

% cyclone
% Focus on 3 larval supplies only (1500/ 150,000/ 1.500.000)

%17/08/24 5_10_20 PAR
[i j]=find(Extract_final_OC_VR_1cyc_1int_170824.tot_larv_in_param_index == 1500);
[k l] =find(Extract_final_OC_VR_1cyc_1int_170824.tot_larv_in_param_index > 149000); 
m=[i;k];
Extract_final_OC_VR_1cyc_1int_3larvsup_170824=Extract_final_OC_VR_1cyc_1int_170824(m,:);

%normalise the data based on data with only 3 larval supply metrics
Extract_final_OC_VR_1cyc_1int_3larvsup_25_5_10_20_170824(:,107:end)=normalize(Extract_final_OC_VR_1cyc_1int_3larvsup_25_5_10_20_170824(:,1:106));
Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm=Extract_final_OC_VR_1cyc_1int_3larvsup_170824;

% this file will be used in R to measure variable importance under cyclone
% pressure
writetable(Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm);
% save('VR_1cyc_1int212col_36scen_prenorm_norm_fix_3LS_3SD_3PAR_170824.mat','Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm','-v7.3')

% now bleaching
% Focus on 3 larval supplies only (1500/ 150,000/ 1.500.000)

%20/08/24
[i j]=find(Extract_final_OC_VR_1dhw_1int_5_10_20_200824.tot_larv_in_param_index == 1500);
[k l] =find(Extract_final_OC_VR_1dhw_1int_5_10_20_200824.tot_larv_in_param_index > 149000); 
m=[i;k];
Extract_final_OC_VR_1dhw_1int_3LS_5_10_20_200824=Extract_final_OC_VR_1dhw_1int_5_10_20_200824(m,:);
Extract_final_OC_VR_1dhw_1int_3LS_5_10_20_200824(:,107:end)=normalize(Extract_final_OC_VR_1dhw_1int_3LS_5_10_20_200824(:,1:106));

Extract_final_OC_VR_1dhw_1int_3LS_5_10_20_200824_cor_norm=Extract_final_OC_VR_1dhw_1int_3LS_5_10_20_200824;
% this file will be used in R to measure variable importance under
% bleaching pressure
writetable(Extract_final_OC_VR_1dhw_1int_3LS_5_10_20_200824_cor_norm);

% Now prep the file containing only scenarios of no disturbance occurring
% post restoration. This can be either cyclone or bleaching compiled file

Nodist_scen=15:15:size(Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm,1);
Extract_final_OC_VR_1cyc_1int_170824_nodist=Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm(Nodist_scen,:);

% Focus on 3 larval supplies only (1500/ 150,000/ 1.500.000)
[i j]=find(Extract_final_OC_VR_1cyc_1int_170824_nodist.tot_larv_in_param_index == 1500);
[k l] =find(Extract_final_OC_VR_1cyc_1int_170824_nodist.tot_larv_in_param_index > 149000); 
m=[i;k];
Extract_final_OC_VR_1cyc_1int_nodist_3larvsup_170824=Extract_final_OC_VR_1cyc_1int_170824_nodist(m,:);

% normalise after subsetting the no disturbance scenarios only
Extract_final_OC_VR_1cyc_1int_nodist_3larvsup_170824(:,107:end)=normalize(Extract_final_OC_VR_1cyc_1int_nodist_3larvsup_170824(:,1:106));

% this file will be used in R to measure variable importance under
% no bleaching pressure
Extract_final_OC_VR_1cyc_1int_3LS_nodist_norm_3SD_3PAR_170824=Extract_final_OC_VR_1cyc_1int_nodist_3larvsup_170824;
writetable(Extract_final_OC_VR_1cyc_1int_3LS_nodist_norm_3SD_3PAR_170824);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
