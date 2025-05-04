% Run intervention scenarios using the Engine for scenarios 
% For efficiency this was done on the UQ Bunya HPC using the matlab script
% and Batch script. To run the code on local computer make sure set
% different path to the model

function HPC_Engine_intervention_model_1cyc_dhw_080824(simul) % activate to run

addpath(genpath('rme_ml')) % activate to run on HPC

%%%%% Use below when you want to run scenarios with only cyclones
rme_init('','./rme_ml/1cyc_1int_240124/config_syst_1dist_1int_240124.xml') 
rme_init('','../rme_ml_2024_01_08/1int_1cyc_090124/config_syst_1dist_1int_150124.xml') % directory to where you saved the customised config file (needs to be where all input data files are) which need to be in the rme_ml_2024_01_08 folder

% Combination file need to be uploaded as wel

addpath('1cyc_1int_240124') % for cyclone scenarios
load('Spec_sub_1cyc_3init_3self_4larv_paper_test_150124.mat') % for cyclone scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Use below when you want to run scenarios with only bleaching events
rme_init('','./rme_ml/1dhw1int_240124/config_syst_1dist_1int_240124.xml') 
rme_init('','../rme_ml_2024_01_08/1int_1cyc_090124/config_syst_1dist_1int_150124.xml') % directory to where you saved the customised config file (needs to be where all input data files are) which need to be in the rme_ml_2024_01_08 folder

% Combination file need to be uploaded as wel
addpath('1dhw_1int_240124') % for bleaching scenarios

load('Spec_sub_1dhw_3init_3self_4larv_paper_test_150124.mat') % for bleaching scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

special_run_subset=special_run_subset2;

% For the HPC, remaining combinations of size distribution and intervention design still need to be defined outside the model data files 
size_dist_index=[1,2]; %2 indexes for size distribution
prop_rest_index=[5,10,20]; %3 options for percentage of reef restored
rest_freq_index=1; % outplanting only once
st_dens_index=[3.4/2,6.8/2,13.6/2]; % 3 options of stocking densities, need to be divided by 2 as the Engine outplants corals twice a year (see rme manual)
DHW_tolerance_index=[0,3]; %heat tolerance of outplant +0 DHW and +3 DHW
rest_end_index=2023;

% evaluate number of unique scenarios outside the model data files
Int_comb=allcomb(size_dist_index,prop_rest_index,rest_freq_index,st_dens_index,DHW_tolerance_index,rest_end_index); %36 unique combinations

% on the HPC everyone of the 36 scenarios are run in parallel unique job
% array id simul
size_dist_comb=Int_comb(simul,1);
prop_rest=Int_comb(simul,2);
rest_freq=Int_comb(simul,3);
st_dens=Int_comb(simul,4);
DHW_tolerance=Int_comb(simul,5);
rest_end=Int_comb(simul,6);

shelf=1

% before creating the model run, we need to set following options

rme_set_option('cyclone_year_one', 1); % these ensure the DHW and cyclone data are used in the order specified in the files
rme_set_option('reorder_dhw_years', 0); % don't reorder dhw years, use specific order
rme_set_option('reefs_have_same_seed', 1); % treat every row as unique reef with the same seed
rme_set_option('use_fixed_seed', 1); % fix seed for reproducability
rme_set_option('fixed_seed', 123); % set seed for reproducability
rme_set_option('thread_count', 12); %threat count
rme_set_option('cots_enabled', 0); %disable cots because cots are controlled for on the restored reefs

% size distribution

%normal distribution investigation
%use the species specific maximum size of adults with fixed ceiling of
%maximum cell size x axis = (100cm2)
% species_max_adult_size=[7853; 7853; 1963; 1256; 2827; 7853];
Lognorm_mean_size=[];

for i=1:6
sp_mean=[17.50380286;24.31965649;17.50380286;11.29180602;11.94162437;16.17448777]; %mean sizes (in diameter) of species 1-6 (see Bozec et al. 2022) based on Dietzel et al. 2020
sp_sigma=[16.86639;18.20781;16.86639;7.970298;14.3261;24.56004]; %sd sizes of species 1-6 (in diameter) (see Bozec et al. 2022) based on Dietzel et al. 2020

Lognorm_mean_size{i}= [log(pi*(sp_mean(i)/2)^2);log(pi*((sp_mean(i)*5)/2)^2)]; % first value is represtative of population predominantly occupied by small individuals, second value is predominantly mid to larger individuals
Lognorm_sigma_size{i}= [log(sp_sigma(i));log(sp_sigma(i)/5)]; % first value is represtative of population predominantly occupied by small individuals, second value is predominantly mid to larger individuals

end

size_mean_dist=Lognorm_mean_size;
size_sigma_dist=Lognorm_sigma_size;

% specify size distribution based on specified setting by simul --> size_dist_comb
    for i=1:6
       rme_set_option('initial_coral_size_dbn', 'LOGNORMAL');
        rme_set_option(['initial_coral_size_mean_sp',num2str(i),''],size_mean_dist{i}(size_dist_comb))
        rme_set_option(['initial_coral_size_sigma_sp',num2str(i),''],size_sigma_dist{i}(size_dist_comb))
    end

    % Run Id name
if size_dist_comb==1
run_id=['VR_OC_SD_LN1_EYR_',num2str(rest_end),'_PRA',num2str(prop_rest),'_RF',num2str(rest_freq),'_SD',num2str(st_dens*2),'_HT',num2str(DHW_tolerance)]
else
run_id=['VR_OC_SD_LN2_EYR_',num2str(rest_end),'_PRA',num2str(prop_rest),'_RF',num2str(rest_freq),'_SD',num2str(st_dens*2),'_HT',num2str(DHW_tolerance)]
end

startyr=2022;
endyr=2037;
% Create model run, setting arbitrary climate scenarios as we use our
% customised sequence and magnitude of disturbances. Run 10 replicates

rme_run_create(run_id,startyr ,endyr, 'ssp_test', 'gsm_test', 10);


% this would be the place to add interventions
%integrate intervention related parameters
rme_set_option('restoration_dhw_tolerance_outplants',DHW_tolerance)
rme_iv_add('VR_outplant', 'outplant', 'All reefs', rest_end, rest_end,rest_freq, prop_rest,st_dens.*(prop_rest/100)) %for more information see manual

% track years of restoration
iv_years=[];
iv_list=rme_iv_list;
for iv_numb=1:rme_iv_count;
    iv_first= rme_iv_first_year(iv_list{iv_numb});
    iv_last=rme_iv_last_year(iv_list{iv_numb});
    if iv_first<iv_last;
        freq=rme_iv_year_step(iv_list{iv_numb});
        iv_multiyear=iv_first:freq:iv_last;
        iv_years=[iv_years iv_multiyear];
    else
        iv_multiyear=[];
        iv_years=[iv_years iv_first];
    end
end

Int_meta=[];
Int_meta.iv_years=iv_years;

y=rme_run_first_year:1:rme_run_last_year;
years_rest_begins=zeros(1,size(y,2));
cumulative_rest_events=zeros(1,size(y,2));
years_last_rest_event=zeros(1,size(y,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise the run at year 2022 
rme_run_init
year=1;

years_rest_begins(1,year)==0;
cumulative_rest_events(1,year)==0;
years_last_rest_event(1,year)==0;

processed_year = rme_run_current_year() % so we can track at what year to model run is at


CELL_final=[];
reef_index_restored=[];
restored_cell_indices=[];

storing_indices=[];
unrestored_cell_indices=[];
unrestored_cell_indices_test=[];


% To track what is happening on the reef at timepoint 1= 2022
reef_index_restored=rme_reefset_get_as_id_list('All reefs');


                        % the number of corals in each cell of the test reef
                        % with and without interventions applied
                        CELL=[];                       
                        REEF=[];
                        SITE=[];
                        NOSITE=[];

                        % Record coral sizes in reference and intervention
                        % model run
                        for repeat=1:10
                            for reef_id=1:size(reef_index_restored,2)

                                %
                                REEF.cc_ref(reef_id,:, repeat) = rme_run_get_cell_data_ref('coral_cm2', reef_id, repeat);
                                REEF.cc_ref(reef_id,:, repeat)=single(REEF.cc_ref(reef_id,:, repeat));
                                REEF.cc_int(reef_id,:,repeat) = rme_run_get_cell_data_int('coral_cm2', reef_id, repeat);
                                REEF.cc_int(reef_id,:,repeat) =single(REEF.cc_int(reef_id,:,repeat));
                                REEF.cc_int_ref(reef_id,:, repeat)= REEF.cc_int(reef_id,:,repeat)- REEF.cc_ref(reef_id,:,repeat);
                                REEF.cc_int_ref(reef_id,:, repeat)=single(REEF.cc_int_ref(reef_id,:, repeat));

                                species_cell=['species_1_cm2'; 'species_2_cm2'; 'species_3_cm2';'species_4_cm2';'species_5_cm2';'species_6_cm2'];
                                for sp=1:6
                                    REEF.cc_sp_ref(reef_id,:,sp, repeat)=rme_run_get_cell_data_ref(species_cell(sp,:), reef_id, repeat);
                                    REEF.cc_sp_ref(reef_id,:,sp, repeat)=single( REEF.cc_sp_ref(reef_id,:,sp, repeat));
                                    REEF.cc_sp_int(reef_id,:,sp, repeat)=rme_run_get_cell_data_int(species_cell(sp,:),reef_id, repeat);
                                    REEF.cc_sp_int(reef_id,:,sp, repeat)=single(REEF.cc_sp_int(reef_id,:,sp, repeat));
                                    REEF.cc_sp_int_ref(reef_id,:,sp, repeat)=REEF.cc_sp_int(reef_id,:,sp, repeat)-REEF.cc_sp_ref(reef_id,:,sp, repeat);
                                    REEF.cc_sp_int_ref(reef_id,:,sp, repeat)=single(REEF.cc_sp_int_ref(reef_id,:,sp, repeat));
                                end

                                REEF.coral_size_ref{repeat}{reef_id} = rme_run_get_coral_sizes_cm2_ref(reef_id, repeat);
                                REEF.coral_size_ref{repeat}{reef_id} =single(REEF.coral_size_ref{repeat}{reef_id});
                                REEF.coral_size_int{repeat}{reef_id} = rme_run_get_coral_sizes_cm2_int(reef_id, repeat);
                                REEF.coral_size_int{repeat}{reef_id}=single(REEF.coral_size_int{repeat}{reef_id});
                                REEF.coral_size_int_outplant{repeat}{reef_id} = rme_run_get_coral_sizes_cm2_outplant(reef_id, repeat);
                                REEF.coral_size_int_outplant{repeat}{reef_id}=single(REEF.coral_size_int_outplant{repeat}{reef_id});

                                % Track demographic composition (more
                                % details see Bozec et al., 2022)


                                     CORAL.size_threshold_pm = 250; % Size threshold for partial mortality (natural and due to cyclones) (Edwards et al. 2011)

                                     CORAL.juv_max_size = 13 ; % at least 2 years old, escaped the post-settlement bottlenecks (Doropoulos et al. 2016)
                                     % Used for cyclone mortality, natural mortality and growth

                                     % Only for tracking size distribution: threshold diameters for transitions juv -> adol -> adults
                                     CORAL.adol_max_diam = floor(sqrt(CORAL.size_threshold_pm/pi)*2); % maximum diameter of adolescents = lower bound of mature corals in the Caribbean
                                     CORAL.adult_max_diam = sqrt(100*100); % defined by cell size

                                     % Define bin size for the size distributions for (1) juvenile, (2) adolescent and (3) adult colonies
                                     CORAL.size_bins = [3 ; 20 ; 500]; % in cm2
                                     CORAL.juv_max_diam = 5; % every diameter <5cm is a juvenile
                                     CORAL.diam_bins = [2 ; 4 ; 14]; % in cm

% Juveniles, sub adults and adults are derived from sizes using
% f_count_sizefreq function
                                     [ref_size_count,class,ref_size_sum] = f_count_sizefreq(REEF.coral_size_ref{repeat}{reef_id}, CORAL);
                                     REEF.ref_size_count{repeat}{reef_id}=ref_size_count;
                                     REEF.ref_class_total{repeat}{reef_id}=ref_size_sum;
                                     [int_size_count,class,int_size_sum] = f_count_sizefreq(REEF.coral_size_int{repeat}{reef_id}, CORAL);
                                     REEF.int_size_count{repeat}{reef_id}=int_size_count;
                                     REEF.int_class_total{repeat}{reef_id}=int_size_sum;
                                     [int_outplant_size_count,class,int_outplant_size_sum] = f_count_sizefreq(REEF.coral_size_int_outplant{repeat}{reef_id}, CORAL);
                                     REEF.int_outplant_size_count{repeat}{reef_id}=int_outplant_size_count;
                                     REEF.int_outplant_class_total{repeat}{reef_id}=int_outplant_size_sum;

                                     %difference between intervention and
                                     %counterfactual model run
                                     fn = fieldnames(REEF.ref_class_total{repeat}{reef_id});
                                     for i = 1 : numel(fn)
                                         REEF.int_ref_size_count{repeat}{reef_id}.(fn{i}) = REEF.int_size_count{repeat}{reef_id}.(fn{i})- REEF.ref_size_count{repeat}{reef_id}.(fn{i});
                                         REEF.int_ref_class_total{repeat}{reef_id}.(fn{i}) = REEF.int_class_total{repeat}{reef_id}.(fn{i})- REEF.ref_class_total{repeat}{reef_id}.(fn{i});
                                       REEF.int_natural_size_count{repeat}{reef_id}.(fn{i})=REEF.int_size_count{repeat}{reef_id}.(fn{i})-REEF.int_outplant_size_count{repeat}{reef_id}.(fn{i});
                                     REEF.int_natural_class_total{repeat}{reef_id}.(fn{i})=REEF.int_class_total{repeat}{reef_id}.(fn{i})-REEF.int_outplant_class_total{repeat}{reef_id}.(fn{i});

                                     end
                                        
                                     SITE.indices=restored_cell_indices;
                                     NOSITE.indices=unrestored_cell_indices;
                                 clear REEF.coral_size_ref{repeat}{reef_id};
                                clear REEF.coral_size_ref{repeat}{reef_id};
                                clear REEF.coral_size_int{repeat}{reef_id};
                                clear REEF.coral_size_int{repeat}{reef_id};
                                clear REEF.coral_size_int_outplant{repeat};
                                clear REEF.coral_size_int_outplant{repeat}{reef_id}
                     
                            end
                        end

                        CELL.REEF=REEF;
                        CELL.SITE=SITE;
                        CELL.NOSITE=NOSITE;
                        CELL_final{year}=CELL;

% Now we can parameterise the specific larval supply, here we use that the
% supply is the same for all 6 species in this study
larvae_in = (special_run_subset.Var10);

rme_run_set_add_coral_larvae_in(1, larvae_in);
rme_run_set_add_coral_larvae_in(2, larvae_in);
rme_run_set_add_coral_larvae_in(3, larvae_in);
rme_run_set_add_coral_larvae_in(4, larvae_in);
rme_run_set_add_coral_larvae_in(5, larvae_in);
rme_run_set_add_coral_larvae_in(6, larvae_in);

% Now let's process the remaining 15 years, modelling year by year to
% extract within reef processes
while (rme_run_current_year() < endyr)

    rme_run_process_years(1);

    processed_year = rme_run_current_year();
    year=year+1

    REEF=[];

                            for repeat=1:10 % for every repeat
                                  for reef_id=1:size(reef_index_restored,2)
                                     
                                    REEF.cc_ref(reef_id,:, repeat) = rme_run_get_cell_data_ref('coral_cm2', reef_id, repeat);
                                     REEF.cc_ref(reef_id,:, repeat)=single(REEF.cc_ref(reef_id,:, repeat));
                                     REEF.cc_int(reef_id,:,repeat) = rme_run_get_cell_data_int('coral_cm2', reef_id, repeat);
                                     REEF.cc_int(reef_id,:,repeat) =single(REEF.cc_int(reef_id,:,repeat));
                                     REEF.cc_int_ref(reef_id,:,repeat)= REEF.cc_int(reef_id,:,repeat)- REEF.cc_ref(reef_id,:,repeat);
                                     REEF.cc_int_ref(reef_id,:,repeat)=single(REEF.cc_int_ref(reef_id,:, repeat));

                                    species_cell=['species_1_cm2'; 'species_2_cm2'; 'species_3_cm2';'species_4_cm2';'species_5_cm2';'species_6_cm2'];
                                    for sp=1:6
                                        REEF.cc_sp_ref(reef_id,:,sp, repeat)=rme_run_get_cell_data_ref(species_cell(sp,:),reef_id, repeat);
                                        REEF.cc_sp_ref(reef_id,:,sp, repeat)=single(REEF.cc_sp_ref(reef_id,:,sp, repeat));
                                        REEF.cc_sp_int(reef_id,:,sp, repeat)=rme_run_get_cell_data_int(species_cell(sp,:),reef_id, repeat);
                                        REEF.cc_sp_int(reef_id,:,sp, repeat)=single(REEF.cc_sp_int(reef_id,:,sp, repeat));
                                        REEF.cc_sp_int_ref(reef_id,:,sp, repeat)=REEF.cc_sp_int(reef_id,:,sp, repeat)-REEF.cc_sp_ref(reef_id,:,sp, repeat);
                                        REEF.cc_sp_int_ref(reef_id,:,sp, repeat)=single(REEF.cc_sp_int_ref(reef_id,:,sp, repeat));
                                    end

                                     REEF.coral_size_ref{repeat}{reef_id} = rme_run_get_coral_sizes_cm2_ref(reef_id, repeat);
                                     REEF.coral_size_ref{repeat}{reef_id} =single(REEF.coral_size_ref{repeat}{reef_id});
                                     REEF.coral_size_int{repeat}{reef_id} = rme_run_get_coral_sizes_cm2_int(reef_id, repeat);
                                     REEF.coral_size_int{repeat}{reef_id}=single(REEF.coral_size_int{repeat}{reef_id});
                                      REEF.coral_size_int_outplant{repeat}{reef_id} = rme_run_get_coral_sizes_cm2_outplant(reef_id, repeat);
                                REEF.coral_size_int_outplant{repeat}{reef_id}=single(REEF.coral_size_int_outplant{repeat}{reef_id});

                                [ref_size_count,class,ref_size_sum] = f_count_sizefreq(REEF.coral_size_ref{repeat}{reef_id}, CORAL);
                                     REEF.ref_size_count{repeat}{reef_id}=ref_size_count;
                                     REEF.ref_class_total{repeat}{reef_id}=ref_size_sum;
                                     [int_size_count,class,int_size_sum] = f_count_sizefreq(REEF.coral_size_int{repeat}{reef_id}, CORAL);
                                     REEF.int_size_count{repeat}{reef_id}=int_size_count;
                                     REEF.int_class_total{repeat}{reef_id}=int_size_sum;
                                      [int_outplant_size_count,class,int_outplant_size_sum] = f_count_sizefreq(REEF.coral_size_int_outplant{repeat}{reef_id}, CORAL);
                                     REEF.int_outplant_size_count{repeat}{reef_id}=int_outplant_size_count;
                                     REEF.int_outplant_class_total{repeat}{reef_id}=int_outplant_size_sum;

                                     % Difference between intervention and
                                     % counterfactual model run
                                      fn = fieldnames(REEF.ref_class_total{repeat}{reef_id});
                                     for i = 1 : numel(fn)
                                         REEF.int_ref_size_count{repeat}{reef_id}.(fn{i}) = REEF.int_size_count{repeat}{reef_id}.(fn{i})- REEF.ref_size_count{repeat}{reef_id}.(fn{i});
                                         REEF.int_ref_class_total{repeat}{reef_id}.(fn{i}) = REEF.int_class_total{repeat}{reef_id}.(fn{i})- REEF.ref_class_total{repeat}{reef_id}.(fn{i});
                                         REEF.int_natural_size_count{repeat}{reef_id}.(fn{i})=REEF.int_size_count{repeat}{reef_id}.(fn{i})-REEF.int_outplant_size_count{repeat}{reef_id}.(fn{i});
                                     REEF.int_natural_class_total{repeat}{reef_id}.(fn{i})=REEF.int_class_total{repeat}{reef_id}.(fn{i})-REEF.int_outplant_class_total{repeat}{reef_id}.(fn{i});

                                     end
                                      clear REEF.coral_size_ref{repeat}{reef_id};
                                clear REEF.coral_size_ref{repeat}{reef_id};
                                clear REEF.coral_size_int{repeat}{reef_id};
                                clear REEF.coral_size_int{repeat}{reef_id};
                                clear REEF.coral_size_int_outplant{repeat};
                                clear REEF.coral_size_int_outplant{repeat}{reef_id}
                                  end


                            CELL.REEF=REEF;
                            CELL_final{year}=CELL; 
                            end

% Also for the next year continue the supply of external larval supply
rme_run_set_add_coral_larvae_in(1, larvae_in);
rme_run_set_add_coral_larvae_in(2, larvae_in);
rme_run_set_add_coral_larvae_in(3, larvae_in);
rme_run_set_add_coral_larvae_in(4, larvae_in);
rme_run_set_add_coral_larvae_in(5, larvae_in);
rme_run_set_add_coral_larvae_in(6, larvae_in);

% Track when the intervention took place compared to current model year.
% This is important to track the timepoint a disturbance occurs relative to
% the intervention
years_before_intervention=rme_run_first_year:iv_years(1)-1;

if ismember(processed_year,iv_years)==1 & ismember(processed_year-1,years_before_intervention)==1;
    years_rest_begins(1,year)=1;
    cumulative_rest_events(1,year)=1;
    years_last_rest_event(1,year)=0;
elseif ismember(processed_year,iv_years)==0 & ismember(processed_year,years_before_intervention)==1;
    years_rest_begins(1,year)=0;
    cumulative_rest_events(1,year)=0;
else
    years_rest_begins(1,year)=years_rest_begins(1,year-1)+1;
    if ismember(processed_year-1,iv_years)==1 & ismember(processed_year,iv_years)==0;
        cumulative_rest_events(1,year)=cumulative_rest_events(1,year-1);
    elseif ismember(processed_year-1,iv_years)==0 & ismember(processed_year,iv_years)==0 & processed_year>iv_years(1);
        cumulative_rest_events(1,year)=cumulative_rest_events(1,year-1);
    else
        cumulative_rest_events(1,year)=cumulative_rest_events(1,year-1)+1;
    end

    if cumulative_rest_events(1,year)==cumulative_rest_events(1,year-1) & years_last_rest_event(1,year-1)==0;
        years_last_rest_event(1,year)=1;
    elseif cumulative_rest_events(1,year)==cumulative_rest_events(1,year-1) & years_last_rest_event(1,year-1)>0;
        years_last_rest_event(1,year)=years_last_rest_event(1,year-1)+1;
    end
end
end



 %% Extract results at individual REEF level at the end of the model run (reefs that have been restored)
 % for both counterfactual as well as intervention model run. Also the
 % difference between the two is calculated (Restoration benefit)

                            y = 2022:1:endyr

                            INT_GBR=[]
                            COVER_GBR=[]
                            COTS_GBR=[]
                            DIST_GBR=[]
                            GEO_GBR=[]
                            CONN_GBR=[]
                            DEMO_GBR=[]
                            SHELTER_GBR=[]

                            for repeat=1:10
                                for position=1:16

                                    %%%%%% GBR-scale ecology
                                    %coral cover
                                    COVER_GBR.coral_pct_ref(:,position,repeat)= rme_run_get_data_ref('coral_pct', y(position),'All reefs',repeat);
                                    COVER_GBR.coral_pct_ref(:,position,repeat)=single(COVER_GBR.coral_pct_ref(:,position,repeat));
                                    COVER_GBR.coral_pct_int(:,position,repeat)= rme_run_get_data_int('coral_pct', y(position),'All reefs',repeat);
                                    COVER_GBR.coral_pct_int(:,position,repeat)=single(COVER_GBR.coral_pct_int(:,position,repeat));
                                    COVER_GBR.coral_pct_int_ref(:,position,repeat)= COVER_GBR.coral_pct_int(:,position,repeat)- COVER_GBR.coral_pct_ref(:,position,repeat);
                                    COVER_GBR.coral_pct_int_ref(:,position,repeat)=single(COVER_GBR.coral_pct_int_ref(:,position,repeat));
                                    species_cover=['species_1_pct'; 'species_2_pct'; 'species_3_pct';'species_4_pct';'species_5_pct';'species_6_pct'];
                                    for sp=1:6
                                        COVER_GBR.coral_pct_sp_ref(:,position,sp,repeat)=rme_run_get_data_ref(species_cover(sp,:), y(position),'All reefs',repeat);
                                        COVER_GBR.coral_pct_sp_ref(:,position,sp,repeat)=single(COVER_GBR.coral_pct_sp_ref(:,position,sp,repeat));
                                        COVER_GBR.coral_pct_sp_int(:,position,sp,repeat)=rme_run_get_data_int(species_cover(sp,:), y(position),'All reefs',repeat);
                                        COVER_GBR.coral_pct_sp_int(:,position,sp,repeat)=single(COVER_GBR.coral_pct_sp_int(:,position,sp,repeat));
                                        COVER_GBR.coral_pct_sp_int_ref(:,position,sp,repeat)=COVER_GBR.coral_pct_sp_int(:,position,sp,repeat)-COVER_GBR.coral_pct_sp_ref(:,position,sp,repeat);
                                        COVER_GBR.coral_pct_sp_int_ref(:,position,sp,repeat)=single(COVER_GBR.coral_pct_sp_int_ref(:,position,sp,repeat));
                                    end

                                    % reef-scale Disturbance
                                    DIST_GBR.cycl_cat_ref(:,position,repeat)= rme_run_get_data_ref('cyclone_cat', y(position),'All reefs',repeat);
                                    DIST_GBR.cycl_cat_ref(:,position,repeat)=single(DIST_GBR.cycl_cat_ref(:,position,repeat));
                                    DIST_GBR.cycl_cat_int(:,position,repeat)= rme_run_get_data_int('cyclone_cat', y(position),'All reefs',repeat);
                                    DIST_GBR.cycl_cat_int(:,position,repeat)=single(DIST_GBR.cycl_cat_int(:,position,repeat));
                                    DIST_GBR.cycl_cat_int_ref(:,position,repeat)= DIST_GBR.cycl_cat_int(:,position,repeat)- DIST_GBR.cycl_cat_ref(:,position,repeat);
                                    DIST_GBR.cycl_cat_int_ref(:,position,repeat)=single(DIST_GBR.cycl_cat_int_ref(:,position,repeat));
                                    DIST_GBR.cor_loss_cyc_ref(:,position,repeat)= rme_run_get_data_ref('cyclone_loss_pct', y(position),'All reefs',repeat);
                                    DIST_GBR.cor_loss_cyc_ref(:,position,repeat)=single(DIST_GBR.cor_loss_cyc_ref(:,position,repeat));
                                    DIST_GBR.cor_loss_cyc_int(:,position,repeat)= rme_run_get_data_int('cyclone_loss_pct', y(position),'All reefs',repeat);
                                    DIST_GBR.cor_loss_cyc_int(:,position,repeat)=single(DIST_GBR.cor_loss_cyc_int(:,position,repeat));
                                    DIST_GBR.cor_loss_cyc_int_ref(:,position,repeat)= DIST_GBR.cor_loss_cyc_int(:,position,repeat)-DIST_GBR.cor_loss_cyc_ref(:,position,repeat);
                                    DIST_GBR.cor_loss_cyc_int_ref(:,position,repeat)=single(DIST_GBR.cor_loss_cyc_int_ref(:,position,repeat));
                                    DIST_GBR.dhw_mag_ref(:,position,repeat)= rme_run_get_data_ref('max_dhw', y(position),'All reefs',repeat);
                                    DIST_GBR.dhw_mag_ref(:,position,repeat)=single(DIST_GBR.dhw_mag_ref(:,position,repeat));
                                    DIST_GBR.dhw_mag_int(:,position,repeat)= rme_run_get_data_int('max_dhw', y(position),'All reefs',repeat);
                                    DIST_GBR.dhw_mag_int(:,position,repeat)=single(DIST_GBR.dhw_mag_int(:,position,repeat));
                                    DIST_GBR.dhw_mag_int_ref(:,position,repeat)=DIST_GBR.dhw_mag_int(:,position,repeat)-DIST_GBR.dhw_mag_ref(:,position,repeat);
                                    DIST_GBR.dhw_mag_int_ref(:,position,repeat)=single(DIST_GBR.dhw_mag_int_ref(:,position,repeat));
                                    DIST_GBR.cor_loss_dhw_ref(:,position,repeat)= rme_run_get_data_ref('dhw_loss_pct', y(position),'All reefs',repeat);
                                    DIST_GBR.cor_loss_dhw_ref(:,position,repeat)=single(DIST_GBR.cor_loss_dhw_ref(:,position,repeat));
                                    DIST_GBR.cor_loss_dhw_int(:,position,repeat)= rme_run_get_data_int('dhw_loss_pct', y(position),'All reefs',repeat);
                                    DIST_GBR.cor_loss_dhw_int(:,position,repeat)=single(DIST_GBR.cor_loss_dhw_int(:,position,repeat));
                                    DIST_GBR.cor_loss_dhw_int_ref(:,position,repeat)= DIST_GBR.cor_loss_dhw_int(:,position,repeat)-DIST_GBR.cor_loss_dhw_ref(:,position,repeat);
                                    DIST_GBR.cor_loss_dhw_int_ref(:,position,repeat)=single(DIST_GBR.cor_loss_dhw_int_ref(:,position,repeat));
                                    DIST_GBR.WQ_mid_ref(:,position,repeat)= rme_run_get_data_ref('ssc_mg_per_l', y(position),'All reefs',repeat);
                                    DIST_GBR.WQ_mid_ref(:,position,repeat)=single(DIST_GBR.WQ_mid_ref(:,position,repeat));
                                    DIST_GBR.WQ_mid_int(:,position,repeat)= rme_run_get_data_int('ssc_mg_per_l', y(position),'All reefs',repeat);
                                    DIST_GBR.WQ_mid_int(:,position,repeat)=single(DIST_GBR.WQ_mid_int(:,position,repeat));
                                    DIST_GBR.WQ_surf_ref(:,position,repeat)= rme_run_get_data_ref('ssc_surface_mg_per_l', y(position),'All reefs',repeat);
                                    DIST_GBR.WQ_surf_ref(:,position,repeat)=single(DIST_GBR.WQ_surf_ref(:,position,repeat));
                                    DIST_GBR.WQ_surf_int(:,position,repeat)= rme_run_get_data_int('ssc_surface_mg_per_l', y(position),'All reefs',repeat);
                                    DIST_GBR.WQ_surf_int(:,position,repeat)=single(DIST_GBR.WQ_surf_int(:,position,repeat));

                                    % Connectivity and geomorphology
                                    GEO_GBR.sand_pct_ref(:,position,repeat)= rme_run_get_data_ref('sand_pct', y(position),'All reefs',repeat);
                                    GEO_GBR.sand_pct_ref(:,position,repeat)=single(GEO_GBR.sand_pct_ref(:,position,repeat));
                                    GEO_GBR.sand_pct_int(:,position,repeat)= rme_run_get_data_int('sand_pct', y(position),'All reefs',repeat);
                                    GEO_GBR.sand_pct_int(:,position,repeat)=single(GEO_GBR.sand_pct_int(:,position,repeat));
                                    GEO_GBR.sand_pct_int_ref(:,position,repeat)=GEO_GBR.sand_pct_int(:,position,repeat)-GEO_GBR.sand_pct_ref(:,position,repeat);
                                    GEO_GBR.sand_pct_int_ref(:,position,repeat)=single(GEO_GBR.sand_pct_int_ref(:,position,repeat));

                                    GEO_GBR.rubble_pct_ref(:,position,repeat)= rme_run_get_data_ref('rubble_pct', y(position),'All reefs',repeat);
                                    GEO_GBR.rubble_pct_ref(:,position,repeat)=single(GEO_GBR.rubble_pct_ref(:,position,repeat));
                                    GEO_GBR.rubble_pct_int(:,position,repeat)= rme_run_get_data_int('rubble_pct', y(position),'All reefs',repeat);
                                    GEO_GBR.rubble_pct_int(:,position,repeat)=single(GEO_GBR.rubble_pct_int(:,position,repeat));
                                    GEO_GBR.rubble_pct_int_ref(:,position,repeat)=GEO_GBR.rubble_pct_int(:,position,repeat)-GEO_GBR.rubble_pct_ref(:,position,repeat);
                                    GEO_GBR.rubble_pct_int_ref(:,position,repeat)=single(GEO_GBR.rubble_pct_int_ref(:,position,repeat));
                                    CONN_GBR.conn_matrix_year(repeat,:)= rme_run_data_years(repeat);
                                    CONN_GBR.conn_matrix_year(repeat,:)=single(CONN_GBR.conn_matrix_year(repeat,:));

                                    CONN_GBR.cor_larv_in_total_ref(:,position,repeat)= rme_run_get_data_ref('coral_larvae_in_per_m2', y(position),'All reefs',repeat);
                                    CONN_GBR.cor_larv_in_total_ref(:,position,repeat)=single(CONN_GBR.cor_larv_in_total_ref(:,position,repeat));
                                    CONN_GBR.cor_larv_in_total_int(:,position,repeat)= rme_run_get_data_int('coral_larvae_in_per_m2', y(position),'All reefs',repeat);
                                    CONN_GBR.cor_larv_in_total_int(:,position,repeat)=single(CONN_GBR.cor_larv_in_total_int(:,position,repeat));
                                    CONN_GBR.cor_larv_in_total_int_ref(:,position,repeat)=CONN_GBR.cor_larv_in_total_int(:,position,repeat)-CONN_GBR.cor_larv_in_total_ref(:,position,repeat);
                                    CONN_GBR.cor_larv_in_total_int_ref(:,position,repeat)=single(CONN_GBR.cor_larv_in_total_int_ref(:,position,repeat));
                                    CONN_GBR.cor_larv_out_total_ref(:,position,repeat)= rme_run_get_data_ref('coral_larvae_out_per_m2', y(position),'All reefs',repeat);
                                    CONN_GBR.cor_larv_out_total_ref(:,position,repeat)=single(CONN_GBR.cor_larv_out_total_ref(:,position,repeat));
                                    CONN_GBR.cor_larv_out_total_int(:,position,repeat)= rme_run_get_data_int('coral_larvae_out_per_m2', y(position),'All reefs',repeat);
                                    CONN_GBR.cor_larv_out_total_int(:,position,repeat)=single(CONN_GBR.cor_larv_out_total_int(:,position,repeat));
                                    CONN_GBR.cor_larv_out_total_int_ref(:,position,repeat)=CONN_GBR.cor_larv_out_total_int(:,position,repeat)-CONN_GBR.cor_larv_out_total_ref(:,position,repeat);
                                    CONN_GBR.cor_larv_out_total_int_ref(:,position,repeat)=single(CONN_GBR.cor_larv_out_total_int_ref(:,position,repeat));
                                    species_larv_in=['species_1_larvae_in_per_m2'; 'species_2_larvae_in_per_m2'; 'species_3_larvae_in_per_m2';'species_4_larvae_in_per_m2';'species_5_larvae_in_per_m2';'species_6_larvae_in_per_m2'];
                                    for sp=1:6
                                        CONN_GBR.cor_larv_in_sp_ref(:,position,sp,repeat)=rme_run_get_data_ref(species_larv_in(sp,:), y(position),'All reefs',repeat);
                                        CONN_GBR.cor_larv_in_sp_ref(:,position,sp,repeat)=single(CONN_GBR.cor_larv_in_sp_ref(:,position,sp,repeat));
                                        CONN_GBR.cor_larv_in_sp_int(:,position,sp,repeat)=rme_run_get_data_int(species_larv_in(sp,:), y(position),'All reefs',repeat);
                                        CONN_GBR.cor_larv_in_sp_int(:,position,sp,repeat)=single(CONN_GBR.cor_larv_in_sp_int(:,position,sp,repeat));
                                        CONN_GBR.cor_larv_in_sp_int_ref(:,position,sp,repeat)=CONN_GBR.cor_larv_in_sp_int(:,position,sp,repeat)-CONN_GBR.cor_larv_in_sp_ref(:,position,sp,repeat);
                                        CONN_GBR.cor_larv_in_sp_int_ref(:,position,sp,repeat)=single(CONN_GBR.cor_larv_in_sp_int_ref(:,position,sp,repeat));
                                    end
                                    species_larv_out=['species_1_larvae_out_per_m2'; 'species_2_larvae_out_per_m2'; 'species_3_larvae_out_per_m2';'species_4_larvae_out_per_m2';'species_5_larvae_out_per_m2';'species_6_larvae_out_per_m2'];
                                    for sp=1:6
                                        CONN_GBR.cor_larv_out_sp_ref(:,position,sp,repeat)=rme_run_get_data_ref(species_larv_out(sp,:), y(position),'All reefs',repeat);
                                        CONN_GBR.cor_larv_out_sp_ref(:,position,sp,repeat)=single(CONN_GBR.cor_larv_out_sp_ref(:,position,sp,repeat));
                                        CONN_GBR.cor_larv_out_sp_int(:,position,sp,repeat)=rme_run_get_data_int(species_larv_out(sp,:), y(position),'All reefs',repeat);
                                        CONN_GBR.cor_larv_out_sp_int(:,position,sp,repeat)=single(CONN_GBR.cor_larv_out_sp_int(:,position,sp,repeat));
                                        CONN_GBR.cor_larv_out_sp_int_ref(:,position,sp,repeat)=CONN_GBR.cor_larv_out_sp_int(:,position,sp,repeat)-CONN_GBR.cor_larv_out_sp_ref(:,position,sp,repeat);
                                        CONN_GBR.cor_larv_out_sp_int_ref(:,position,sp,repeat)=single(CONN_GBR.cor_larv_out_sp_int_ref(:,position,sp,repeat));
                                    end
                                    %

                                    %%Demographics
                                    % juvenile
                                    DEMO_GBR.juv_total_ref(:,position,repeat)= rme_run_get_data_ref('coral_juvenile_count_per_m2', y(position),'All reefs',repeat);
                                    DEMO_GBR.juv_total_ref(:,position,repeat)=single(DEMO_GBR.juv_total_ref(:,position,repeat));
                                    DEMO_GBR.juv_total_int(:,position,repeat)= rme_run_get_data_int('coral_juvenile_count_per_m2', y(position),'All reefs',repeat);
                                    DEMO_GBR.juv_total_int(:,position,repeat)=single(DEMO_GBR.juv_total_int(:,position,repeat));
                                    DEMO_GBR.juv_total_int_ref(:,position,repeat)=DEMO_GBR.juv_total_int(:,position,repeat)-DEMO_GBR.juv_total_ref(:,position,repeat);
                                    DEMO_GBR.juv_total_int_ref(:,position,repeat)=single(DEMO_GBR.juv_total_int_ref(:,position,repeat));

                                    species_juv=['sp1_juvenile_count_per_m2'; 'sp2_juvenile_count_per_m2'; 'sp3_juvenile_count_per_m2';'sp4_juvenile_count_per_m2';'sp5_juvenile_count_per_m2';'sp6_juvenile_count_per_m2'];
                                    for sp=1:6
                                        DEMO_GBR.juv_total_sp_ref(:,position,sp,repeat)=rme_run_get_data_ref(species_juv(sp,:), y(position),'All reefs',repeat);
                                        DEMO_GBR.juv_total_sp_ref(:,position,sp,repeat)=single(DEMO_GBR.juv_total_sp_ref(:,position,sp,repeat));
                                        DEMO_GBR.juv_total_sp_int(:,position,sp,repeat)=rme_run_get_data_int(species_juv(sp,:), y(position),'All reefs',repeat);
                                        DEMO_GBR.juv_total_sp_int(:,position,sp,repeat)=single(DEMO_GBR.juv_total_sp_int(:,position,sp,repeat));
                                        DEMO_GBR.juv_total_sp_int_ref(:,position,sp,repeat)=DEMO_GBR.juv_total_sp_int(:,position,sp,repeat)-DEMO_GBR.juv_total_sp_ref(:,position,sp,repeat);
                                        DEMO_GBR.juv_total_sp_int_ref(:,position,sp,repeat)=single(DEMO_GBR.juv_total_sp_int_ref(:,position,sp,repeat));
                                    end

                                    % recruits
                                    DEMO_GBR.recr_total_ref(:,position,repeat)= rme_run_get_data_ref('coral_recruit_count_per_m2', y(position),'All reefs',repeat);
                                    DEMO_GBR.recr_total_ref(:,position,repeat)=single(DEMO_GBR.recr_total_ref(:,position,repeat));
                                    DEMO_GBR.recr_total_int(:,position,repeat)= rme_run_get_data_int('coral_recruit_count_per_m2', y(position),'All reefs',repeat);
                                    DEMO_GBR.recr_total_int(:,position,repeat)=single(DEMO_GBR.recr_total_int(:,position,repeat));
                                    DEMO_GBR.recr_total_int_ref(:,position,repeat)=DEMO_GBR.recr_total_int(:,position,repeat)-DEMO_GBR.recr_total_ref(:,position,repeat);
                                    DEMO_GBR.recr_total_int_ref(:,position,repeat)=single(DEMO_GBR.recr_total_int_ref(:,position,repeat));

                                    species_recr=['sp1_recruit_count_per_m2'; 'sp2_recruit_count_per_m2'; 'sp3_recruit_count_per_m2';'sp4_recruit_count_per_m2';'sp5_recruit_count_per_m2';'sp6_recruit_count_per_m2'];
                                    for sp=1:6
                                        DEMO_GBR.recr_total_sp_ref(:,position,sp,repeat)=rme_run_get_data_ref(species_recr(sp,:), y(position),'All reefs',repeat);
                                        DEMO_GBR.recr_total_sp_ref(:,position,sp,repeat)=single(DEMO_GBR.recr_total_sp_ref(:,position,sp,repeat));
                                        DEMO_GBR.recr_total_sp_int(:,position,sp,repeat)=rme_run_get_data_int(species_recr(sp,:), y(position),'All reefs',repeat);
                                        DEMO_GBR.recr_total_sp_int(:,position,sp,repeat)=single(DEMO_GBR.recr_total_sp_int(:,position,sp,repeat));
                                        DEMO_GBR.recr_total_sp_int_ref(:,position,sp,repeat)=DEMO_GBR.recr_total_sp_int(:,position,sp,repeat)-DEMO_GBR.recr_total_sp_ref(:,position,sp,repeat);
                                        DEMO_GBR.recr_total_sp_int_ref(:,position,sp,repeat)=single(DEMO_GBR.recr_total_sp_int_ref(:,position,sp,repeat));
                                    end



                                    % restoration specific parameters
                                    INT_GBR.outplant_m2_ref(:,position,repeat)= rme_run_get_data_ref('outplant_count_per_m2', y(position),'All reefs',repeat);
                                    INT_GBR.outplant_m2_int(:,position,repeat)= rme_run_get_data_int('outplant_count_per_m2', y(position), 'All reefs',repeat);
                                    INT_GBR.outplant_m2_int_ref(:,position,repeat)=INT_GBR.outplant_m2_int(:,position,repeat)-INT_GBR.outplant_m2_ref(:,position,repeat);
                                end
                            end

             %  Intervention time specific metadata
             INT_GBR.years_rest_begins=years_rest_begins;
             INT_GBR.cumulative_rest_events=cumulative_rest_events;
             INT_GBR.years_last_rest_event=years_last_rest_event;
             INT_GBR.heattol=rme_get_option('restoration_dhw_tolerance_outplants');


             GBR.COVER_GBR=COVER_GBR;
             GBR.COTS_GBR=COTS_GBR;
             GBR.DIST_GBR=DIST_GBR;
             GBR.GEO_GBR=GEO_GBR;
             GBR.CONN_GBR=CONN_GBR;
             GBR.DEMO_GBR=DEMO_GBR;
             GBR.SHELTER_GBR=SHELTER_GBR;
             GBR.INT_GBR=INT_GBR;


             DENSITY.GBR=GBR;
             DENSITY.CELL_final=CELL_final;

%Save all output files             
save([run_id,'_demo_080824.mat'],'DENSITY','-v7.3');

end
    