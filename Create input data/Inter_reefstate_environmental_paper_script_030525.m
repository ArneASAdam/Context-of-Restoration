%%%%% Intervention scenario exploration- Virtual reef concept %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reef state, environmental parameters
shelf=1; % inshore
startyr=2022; % model initialisation year
endyr=2037; % model end year
rest_start=2023; % year when outplanting takes place


startyr_options=[startyr];
endyr_options=[endyr];

% geomorphology
% Reef Area, sand, rubble, shelf position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gbr_spatial=readtable("Reef_emulator_Jan_2024/rme_ml_2024_01_08/data_files/id/id_list_2023_03_30.csv"); % Can be found in folder data_files of the model

gbr_spatial.Properties.VariableNames={'Reef_ids','Reef_area_CH','prop_ungraz','Shelf_position'};
head(gbr_spatial);
gbr_spatial_shelf=(gbr_spatial(table2array(gbr_spatial(:,4))==1,:));

% Percentage sand on the grid which represents ungrazeable area
sandcover=0.1;  % as proportion: 10%

% Percentage rubble present on the grid
rubble_LB_M_UB=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ecological parameters (input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options for initial coral cover across all species
sp_init_cover_all=[1,2,3]; % 3 indexes, corresponding to 3 different initial coral cover states: 3, 6 and 9%

% size distribution
% Use 2 indexes to specify mean or sigma values for size distribution
% (normal and lognormal)

size_dist_all=[1,2]; % Option 1 corresponds to predominantly small individuals in the populations; options to corresponds to predominantly mid to larger individuals 

% Species composition indexes 
sp_comp=[1;2]; % option 1 corresponds to all species present on the reef versus only non-acroporidae

%Initial cots cover (still parameterised but will not be used for this study because the restored reef will be controlled for potential outbreaks, this is to avoid errors during model initialisation)
init_cover_cots_stats=[]
init_cots_gbr=readtable("Reef_emulator_Jan_2024/rme_ml_2024_01_08/data_files/initial/cots_2022",'Headerlines', 10);
init_cots_gbr.Properties.VariableNames={'Reef_ids',
    'cover_rep1',
    'cover_rep2',
    'cover_rep3',
    'cover_rep4',
    'cover_rep5',
    'cover_rep6',
    'cover_rep7',
    'cover_rep8',
    'cover_rep9',
    'cover_rep10',
    'cover_rep11',
    'cover_rep12',
    'cover_rep13',
    'cover_rep14',
    'cover_rep15',
    'cover_rep16',
    'cover_rep17',
    'cover_rep18',
    'cover_rep19',
    'cover_rep20'};

%isolate initial coral covers at inshore shelf position
init_cov_cots_shelf=init_cots_gbr(ismember(table2array(init_cots_gbr(:,1)),table2array(gbr_spatial_shelf(:,1))),:);
% Add shelf as a column
init_cov_cots_shelf.Shelf(1:size(init_cov_cots_shelf,1), 1) = shelf;

init_cov_cots_shelf_quantile = quantile(table2array(init_cov_cots_shelf(:,(2:21))),[0.05 0.25 0.5 0.75 0.95],"all");

cots_LB_M_UB = init_cov_cots_shelf_quantile(2:4); % interquartile range of cots found on the GBR in 2022 are not modelled because the restored reef will be controlled for potential outbreaks

%cots are not modelled in the simplistic intervention models

% cots connectivity
%use same as coral but deactivate cots because in normal case COTS will be
%controled for in the restored reef

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connectivity is replaced by larval supply in and self retention in thhe
% model itself but still need to set here to evaluate number of unique
% scenarios

%options of self retention
self_retention=[0.01; 0.07;0.28];

%External larval supply options
larval_supply_spec_m2_YM = [0.1e6/400; 1e6/400; 1e7/400 ;10e7/400]; % based on maximum recruitment threshold in the model (Bozec et al., 2022)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Environmental variables%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sedimentation
% For sedimentation we want to
sed_files=["SSC_surface_Baseline_q3b_WQ_mg_per_L.csv";"SSC_middepth_Baseline_q3b_WQ_mg_per_L.csv"];
sed_shelf_quantile=[]

for i=1:2
    shelf=1
% file_path_sed=append("Reef_emulator_virtual_June_2023/rme_ml_2023_06_27/data_files/water/",sed_files{i})
file_path_sed=append("Reef_emulator_Jan_2024/rme_ml_2024_01_08/data_files/water/",sed_files{i})

sed_gbr=readtable(file_path_sed,'Headerlines', 1);
sed_gbr.Properties.VariableNames={'Reef_ids',
    'sed_yr11',
    'sed_yr12',
    'sed_yr13',
    'sed_yr14',
    'sed_yr15',
    'sed_yr16'};

% isolate reefs based on shelf position (is important when modelling a
% virtual reef with different larval reduction factors linked by shelf
% position)

%isolate initial coral covers at certain shelf position
sed_shelf=sed_gbr(ismember(table2array(sed_gbr(:,1)),table2array(gbr_spatial_shelf(:,1))),:);
% Add shelf as a column
sed_shelf.Shelf(1:size(sed_shelf,1), 1) = shelf;

sed_shelf_quantile{i} = quantile(table2array(sed_shelf(:,(2:7))),[0.05 0.1 0.25 0.5 0.75 0.9 0.95],"all"); %GBR
end

% Time fluctuating variables --> disturbances



%intervention options
%option 1= only cyclones
%option 2= only bleaching
interv_option=[1;2]

%%%%%%%%%%%%%%%%%%%%%%%
% Building combinations
% start with startyr and endyr and remove incorrect combinations

% combinations (for 1 disturbance and 1 intervention) --> testing
% different options of self retention to capture the extreme values of
% recruitment rate (low/medium/high)
% First focus on scenarios with one disturbance as well as zero (case of no
% disturbance)
% We do this both for cyclone and bleaching event 

% Cyclone scenarios
E=array2table(allcomb(startyr_options(1:end), ...%only 1 start year
    endyr_options(1:end), ...% only 1 end years
    sp_init_cover_all([1 2 3]), ...%3 initial coral cover conditions
    sp_comp, ...%2 species composition options
    cots_LB_M_UB([2],1), ...%fixed cots setting
    rubble_LB_M_UB([1],1), ...%fixed rubble setting
    sed_shelf_quantile{1}([2 6],1), ...% surface sedimentation (depth sedimentation is linked with surface sedimentation afterwards to avoid introducing error scenarios) 
    self_retention, ... %3 self retention options
    larval_supply_spec_m2_YM([1 2 3 4],1), ...%4 options of larval supply 
    rest_start(1),...% start intervention in 2023
    interv_option([1])));

% Bleaching scenarios
E=array2table(allcomb(startyr_options(1:end), ...%only 1 start year
    endyr_options(1:end), ...% only 1 end years
    sp_init_cover_all([1 2 3]), ...%3 initial coral cover conditions
    sp_comp, ...%2 species composition options
    cots_LB_M_UB([2],1), ...%fixed cots setting
    rubble_LB_M_UB([1],1), ...%fixed rubble setting
    sed_shelf_quantile{1}([2 6],1), ...% surface sedimentation (depth sedimentation is linked with surface sedimentation afterwards to avoid introducing error scenarios) 
    self_retention, ... %3 self retention options
    larval_supply_spec_m2_YM([1 2 3 4],1), ...%4 options of larval supply
    rest_start(1),...% start intervention in 2023
    interv_option([2])));

% link surface sedimentation with depth sedimentation
size(E);
    WQ_mid=[];
    for i=1:size(E,1)
    if E.Var7(i)==sed_shelf_quantile{1}([2],1)
        WQ_mid(i,1)=sed_shelf_quantile{2}([2],1);
    else
WQ_mid(i,1)=sed_shelf_quantile{2}([6],1);
    end
    end

% Add middepth data to data
  E_final=[E(:,1:7) array2table(WQ_mid) E(:,8:end)];
  E_final.Properties.VariableNames=["Var1", "Var2", "Var3","Var4", "Var5", "Var6","Var7", "Var8", "Var9","Var10", "Var11", "Var12"];
E=E_final;

%%%% Disturbances activated after intervention started
F=[];
F1=[];

%option 1: only cyclones
for t=1:size(E,1)
    if E.Var12(t,:)==1 

        for i=1:endyr_options-startyr_options+1
            a{i}=[0 ;2]; % Post restoration years can only be no cyclone (cyclone value=0) or cyclone of magnitude 2 (Saffir-Simpson scale)
        end
        % if sum(v)<2
        v=allcomb(a{:});
        for b=1:size(v,1)
         if sum(v(b,:))>2 | v(b,1)>0 | v(b,2)>0 %this is for making sure
                % make sure that scenario is integrated where cyclone occurs when
                % restoration happens, which is in 2023
                % maximum number of cyclones is 1
                v(b,:)=nan;
            end
        end

        v(any(isnan(v),2),:)=[];
        v(1,:)=[];
        u=v;

        %remove NAs
        u(any(isnan(u),2),:)=[];
        w=repelem(0,endyr_options-startyr_options+1);
        u=[u;w];


        %replicate E row, n times size(u,1)
    E_copy={repmat(E(t,:),size(u,1),1)};

    Finalscen= [E_copy{1} array2table(u)];
                F1= cat(1,F1 ,Finalscen);

%option 2: only bleaching
%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif E.Var12(t,:)==2

        for i=1:endyr_options-startyr_options+1

            a{i}=[0 ;8];% bleaching options: no bleaching event or 1 bleaching event of 8 DHW (which is where we start to see impacts on coral cover)
        end
        v=allcomb(a{:});
        for b=1:size(v,1)            
        if sum(v(b,:))>16 | v(b,1)>0 | v(b,2)>0 | sum(v(b,:))<9 %remove option with more than 1 bleaching event of 8 DHW

                v(b,:)=nan;
            end
        end

        %remove NAs
        v(any(isnan(v),2),:)=[];
        v(1,:)=[];
        u=v;

        %remove NAs
        u(any(isnan(u),2),:)=[];
         w=repelem(0,endyr_options-startyr_options+1);
        u=[u;w];

        %replicate E row, n times size(u,1)
        E_copy={repmat(E(t,:),size(u,1),1)};
               Finalscen= [E_copy{1} array2table(u)];
        F1= cat(1,F1 ,Finalscen);
    end
end

special_run_subset2=F1;
scen_matrix=special_run_subset2;

% save option data set
save(['Spec_sub_1cyc_3init_3self_4larv_paper_test_150124.mat'],'special_run_subset2','-v7.3');
save(['Spec_sub_1dhw_3init_3self_4larv_paper_test_150124.mat'],'special_run_subset2','-v7.3');
