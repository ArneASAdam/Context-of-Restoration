% Script below is used to extract various data parameters from the model outcomes necessary for
% downstream analyses
function HPC_extractdata_1cyc_dhw_paper_080824(scen)

addpath(genpath('rme_ml'))
format shortG

%%%%% Use below when extracting data related to scenarios with only cyclones

addpath('1cyc_1int_240124') % for cyclone scenarios
load('Spec_sub_1cyc_3init_3self_4larv_paper_test_150124.mat') % for cyclone scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Use below when extracting data related to scenarios with only cyclones
addpath('1dhw_1int_240124') % for bleaching scenarios

load('Spec_sub_1dhw_3init_3self_4larv_paper_test_150124.mat') % for bleaching scenarios

special_run_subset=special_run_subset2;
startyr=2023;

%Directory where output files are
cd '/QRISdata/Q5785/Virtual_reef/version_080124/'

folder='VR_1dist_1int_cyc_080824' %for cyclone scenarios
folder='VR_1dist_1int_dhw_080824' %for bleaching scenarios

extension='mat';
scenarios=dir([folder filesep '**' filesep '*080824.' extension]);

scenarios
cd(folder);

int_scen=scenarios(scen).name;

%Load file
load(int_scen);
out2=zeros(1,106);% to collate data from multiple simulations

    for t=1:15
censusyr_options=[2023:2037]; %extract data from every timepoint on or after the restoration took place
censusyr=censusyr_options(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out=zeros(1,106); % for simulation specific normalisation

model_period=2023:1:2037;
startyr_ind=find(model_period==startyr);
censusyr_ind=find(model_period==censusyr);
[numreefs,years,sims]=size(DENSITY.GBR.COVER_GBR.coral_pct_int); %should be 2160 unique scenarios,16 years,10 replicates

s=1;
count=1;
for r=1:numreefs
r;
out(count,1)=1; %non specific, represent SSPs
out(count,2)=1; %non specific, represent GCMs
out(count,3)=startyr; % start year
out(count,4)=censusyr; % census year
out(count,5)=r; % scenario row

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ecology and demographics-at initiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out(count,6)=special_run_subset.Var9(r); % self retention rate unique scenario
out(count,7)=(mean(DENSITY.GBR.CONN_GBR.cor_larv_in_total_ref(r,2,:),3)); %mean larval supply (ref) when restoration took place, used as check
out(count,8)=mean(DENSITY.GBR.GEO_GBR.rubble_pct_ref(r,1,:),3); % rubble cover at model initiation
larv_BSY=mean(DENSITY.GBR.CONN_GBR.cor_larv_in_total_int(r,1:startyr_ind,:),2); %mean larval supply (int) when restoration took place, used as check
out(count,9)= mean(larv_BSY,3);
out(count,10)=mean(DENSITY.GBR.COVER_GBR.coral_pct_ref(r,startyr_ind-1,:),3); % coral cover state at model initialisation

% size distribution index used is in the name
idx = strfind(int_scen,'SD_LN');
if str2num(int_scen(idx(1)+5))==2
out(count,11)=2;
else
 out(count,11)=1;
end

%species composition index
if DENSITY.GBR.COVER_GBR.coral_pct_sp_int(r,1,1, 1)>0
    out(count,12)=1;
else
     out(count,12)=2;   
end

% Following is proportion of species represented in total coral cover at
% model initialisation
col_index_startyr_cover=[13 14 15 16 17 18];
for species=1:6
    out(count,col_index_startyr_cover(species))=mean(DENSITY.GBR.COVER_GBR.coral_pct_sp_int(r,startyr_ind-1,species, :),4)/mean(DENSITY.GBR.COVER_GBR.coral_pct_int(r,startyr_ind-1, :),3);
end

%Number of juveniles at model initialisation
out(count,19)=mean(DENSITY.GBR.DEMO_GBR.juv_total_ref(r,startyr_ind-1,:),3);
%Number of recruits at model initialisation (2022)
out(count,20)=mean(DENSITY.GBR.DEMO_GBR.recr_total_ref(r,startyr_ind-1,:),3);

for rep=1:10
    adol_reef(rep,1)=DENSITY.CELL_final{startyr_ind-1}.REEF.ref_class_total{rep}{r}.adol;
    adult_reef(rep,1)=DENSITY.CELL_final{startyr_ind-1}.REEF.ref_class_total{rep}{r}.adult;
end

out(count,21)=mean(adol_reef); %number of subadults and adults at start year
out(count,22)=mean(adult_reef);%number of adults and adults at start year

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reef level model outcomes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out(count,23)=mean(DENSITY.GBR.COVER_GBR.coral_pct_ref(r,censusyr_ind,:),3); % mean total coral cover CF at census timepoint (mean across 10 replicates)
out(count,24)=std(DENSITY.GBR.COVER_GBR.coral_pct_ref(r,censusyr_ind,:),0,3);% sd total coral cover at CF census timepoint (mean across 10 replicates)

out(count,25)=mean(DENSITY.GBR.COVER_GBR.coral_pct_int(r,censusyr_ind,:),3); % mean total coral cover intervention at census timepoint (mean across 10 replicates)
out(count,26)=std(DENSITY.GBR.COVER_GBR.coral_pct_int(r,censusyr_ind,:),0,3); % sd total coral cover at CF census timepoint (mean across 10 replicates)

out(count,27)=mean(DENSITY.GBR.COVER_GBR.coral_pct_int_ref(r,censusyr_ind,:),3); % mean total coral cover (INT minus CF) at census timepoint (mean across 10 replicates)
out(count,28)=std(DENSITY.GBR.COVER_GBR.coral_pct_int_ref(r,censusyr_ind,:),0,3); % sd total coral cover at (INT minus CF) at census timepoint (mean across 10 replicates)

% Not used for this study but kept in the data files
% sum_int_outplanted_area_m2_reef=[];
% sum_int_all_area_m2_reef=[];
% sum_int_natural_area_m2_reef=[];
% int_natural_cover_reef=[];
% int_outplant_cover_reef=[];
% 
% first_int_year_index=find(DENSITY.GBR.INT_GBR.years_rest_begins==1);
% 
% for site_year=1:16
% if site_year>=first_int_year_index
% 
%     for repeat=1:10
% 
% %         for reef=1:48
%             %cyclone
%             sum_int_outplanted_area_m2_reef=sum(DENSITY.CELL_final{site_year}.REEF.coral_size_int_outplant{repeat}{r}./10000,2);
%             sum_int_all_area_m2_reef=sum(DENSITY.CELL_final{site_year}.REEF.coral_size_int{repeat}{r}./10000,2);
%             sum_int_natural_area_m2_reef=sum_int_all_area_m2_reef-sum_int_outplanted_area_m2_reef;
%             int_natural_cover_reef(r,site_year,repeat)=(sum_int_natural_area_m2_reef/100)*100;
%             int_outplant_cover_reef(r,site_year,repeat)=(sum_int_outplanted_area_m2_reef/100)*100;
%     end
% %         end
% 
% else
%      for repeat=1:10
%     sum_int_outplanted_area_m2_reef=str2double('NA');
%     sum_int_all_area_m2_reef=str2double('NA');
%     sum_int_natural_area_m2_reef=str2double('NA');
%     int_natural_cover_reef(r,site_year,repeat)=str2double('NA');
%     int_outplant_cover_reef(r,site_year,repeat)=str2double('NA');
%      end
%         end
%             end
% % end
% 
% % mean and stdev of the natural population in reference run
%         out(count,29)=mean(DENSITY.GBR.COVER_GBR.coral_pct_ref(r,censusyr_ind,:),3); %natural
%         out(count,30)=std(DENSITY.GBR.COVER_GBR.coral_pct_ref(r,censusyr_ind,:),0,3);%stdev
% 
% % mean and stdev of the natural population in intervention run at census
% % year
% if  isnan(int_natural_cover_reef)
%         out(count,31)=int_natural_cover_reef;
%         out(count,32)=int_natural_cover_reef;%stdev    
% else
%         out(count,31)=mean(int_natural_cover_reef(r,censusyr_ind,:),3);
%         out(count,32)=std(int_natural_cover_reef(r,censusyr_ind,:),0,3);%stdev
% end 
%         % mean and stdev of the natural population in intervention run
%         if  isnan(int_outplant_cover_reef)
%         out(count,33)=int_outplant_cover_reef;
%         out(count,34)=int_outplant_cover_reef;%stdev
%         else
%         out(count,33)=mean(int_outplant_cover_reef(r,censusyr_ind,:),3);
%         out(count,34)=std(int_outplant_cover_reef(r,censusyr_ind,:),0,3);%stdev
%         end

%         out(count,35)=mean(DENSITY.GBR.SHELTER_GBR.shelter_total_int_ref(r,censusyr_ind,:),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coral cover benefit at census timepoint for individual species
        col_index_censusyr_cover=[36 37 38 39 40 41];

        for species=1:6
                out(count,col_index_censusyr_cover(species))=mean(DENSITY.GBR.COVER_GBR.coral_pct_sp_int_ref(r,censusyr_ind,species, :),4)/mean(DENSITY.GBR.COVER_GBR.coral_pct_int_ref(r,censusyr_ind, :),3);
        end

        % Mean larval supply ben between restoration event and census
        % year
        Larve_all_styr_cyr=mean(DENSITY.GBR.CONN_GBR.cor_larv_in_total_int_ref(r,startyr_ind:censusyr_ind,:),2);
        out(count,42)=mean(Larve_all_styr_cyr,3);

        % Mean larval supply ben 5 years prior to census
        % year
if censusyr_ind-5>0
    Larve_all_cyr5_cyr=mean(DENSITY.GBR.CONN_GBR.cor_larv_in_total_int_ref(r,censusyr_ind-5:censusyr_ind,:),2);
    out(count,43)=mean(Larve_all_cyr5_cyr,3);

    Rubble_all_cyr5_cyr=mean(DENSITY.GBR.GEO_GBR.rubble_pct_int_ref(r,censusyr_ind-5:censusyr_ind,:),2);
    out(count,44)=mean(Rubble_all_cyr5_cyr,3);

else
  out(count,43)=str2double('NA');
  out(count,44)=str2double('NA');
end


%Demographic info at reef scale, natural vs outplanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Natural and outplanted corals are added for supplementary analysis
for rep=1:10
        juv_reef_nat_ref(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.ref_class_total{rep}{r}.juv_CF;
        adol_reef_nat_ref(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.ref_class_total{rep}{r}.adol;
        adult_reef_nat_ref(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.ref_class_total{rep}{r}.adult;
        if censusyr_ind>=first_int_year_index
        juv_reef_nat_int(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.int_natural_class_total{rep}{r}.juv_CF;
        adol_reef_nat_int(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.int_natural_class_total{rep}{r}.adol;
        adult_reef_nat_int(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.int_natural_class_total{rep}{r}.adult;
        juv_reef_out_int(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.int_outplant_class_total{rep}{r}.juv_CF;
        adol_reef_out_int(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.int_outplant_class_total{rep}{r}.adol;
        adult_reef_out_int(rep,1)=DENSITY.CELL_final{censusyr_ind}.REEF.int_outplant_class_total{rep}{r}.adult;
        else
        juv_reef_nat_int(rep,1)=str2double('NA');
        adol_reef_nat_int(rep,1)=str2double('NA');
        adult_reef_nat_int(rep,1)=str2double('NA');
        juv_reef_out_int(rep,1)=str2double('NA');
        adol_reef_out_int(rep,1)=str2double('NA');
        adult_reef_out_int(rep,1)=str2double('NA');
        end
           
end

        out(count,45)=mean(juv_reef_nat_ref);
        out(count,46)=mean(adol_reef_nat_ref);
        out(count,47)=mean(adult_reef_nat_ref);
        out(count,48)=mean(juv_reef_nat_int);
        out(count,49)=mean(adol_reef_nat_int);
        out(count,50)=mean(adult_reef_nat_int);
        out(count,51)=mean(juv_reef_out_int);
        out(count,52)=mean(adol_reef_out_int);
        out(count,53)=mean(adult_reef_out_int);

% restoration site level
sum_int_outplanted_area_m2_site=[];
sum_int_all_area_m2_site=[];
sum_int_natural_area_m2_site=[];
int_natural_cover_site=[];
int_outplant_cover_site=[];
sum_ref_all_area_m2_site=[];
ref_all_cover_site=[];
sum_ref_all_area_m2_site=[];

% check restoration area
idx = strfind(int_scen,'PRA');
allcell_test=str2num(int_scen(idx(1)+4))
if isempty(allcell_test)
allcell_total=str2num(int_scen(idx(1)+3))
else
allcell_total=str2num(int_scen([idx(1)+3 idx(1)+4]))
end

%Not used for this study but still recorded
% for site_year=1:16
%           for repeat=1:10
%             ref_all_cover_site=str2double('NA');
%              sum_int_ref_all_area_m2_site=str2double('NA');
%             sum_int_natural_area_m2_site=str2double('NA');
%             int_all_cover_site(r,site_year,repeat)=str2double('NA');
%             int_ref_all_cover_site(r,site_year,repeat)=str2double('NA');
%             int_natural_cover_site(r,site_year,repeat)=str2double('NA');
%             int_outplant_cover_site(r,site_year,repeat)=str2double('NA');
%           end
% %         end
%     end
% % end
% 
%     out(count,54)=str2double('NA'); %natural
%         out(count,55)=str2double('NA');%stdev
% 
% % mean and stdev of site level in intervention run
%         out(count,56)=str2double('NA');
%         out(count,57)=str2double('NA');%stdev
% 
% % mean and stdev of site level in difference
%         out(count,58)=str2double('NA');
%         out(count,59)=str2double('NA');%stdev
% 
% % mean and stdev of natural population at site level reference run
%         out(count,60)=str2double('NA'); %natural
%         out(count,61)=str2double('NA');%stdev
% 
% % mean and stdev of natural population at site level intervention run
%         out(count,62)=str2double('NA'); %natural
%         out(count,63)=str2double('NA');%stdev
% 
% % mean and stdev of outplanting population at site level intervention run
%         out(count,64)=str2double('NA'); %natural
%         out(count,65)=str2double('NA');%stdev

%Demographic info at site scale, natural vs outplanted
first_int_year_index=find(DENSITY.GBR.INT_GBR.years_rest_begins==1);

% for startyr_ind=1:16
%                 if censusyr_ind>=first_int_year_index %& Disturbance==1

% for repeat=1:10
% %             if startyr_ind>=first_int_year_index
%         juv_site_nat_ref(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_ref_total{repeat}{r}(3);
%         adol_site_nat_ref(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_ref_total{repeat}{r}(4);
%         adult_site_nat_ref(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_ref_total{repeat}{r}(5);
% %         if startyr_ind>=first_int_year_index
%         juv_site_nat_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_natural_total{repeat}{r}(3);
%         adol_site_nat_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_natural_total{repeat}{r}(4);
%         adult_site_nat_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_natural_total{repeat}{r}(5);
%         juv_site_out_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_outplant_total{repeat}{r}(3);
%         adol_site_out_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_outplant_total{repeat}{r}(4);
%         adult_site_out_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_outplant_total{repeat}{r}(5);
%         juv_site_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_total{repeat}{r}(3);
%         adol_site_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_total{repeat}{r}(4);
%         adult_site_int(repeat,1)=DENSITY.CELL_final{censusyr_ind}.SITE.coral_size_int_total{repeat}{r}(5);
% end
% %                 else
%                     for repeat=1:10
%         juv_site_nat_ref(repeat,1)=str2double('NA');
%         adol_site_nat_ref(repeat,1)=str2double('NA');
%         adult_site_nat_ref(repeat,1)=str2double('NA');
% 
%         juv_site_nat_int(repeat,1)=str2double('NA');
%         adol_site_nat_int(repeat,1)=str2double('NA');
%         adult_site_nat_int(repeat,1)=str2double('NA');
%         juv_site_out_int(repeat,1)=str2double('NA');
%         adol_site_out_int(repeat,1)=str2double('NA');
%         adult_site_out_int(repeat,1)=str2double('NA');
%         juv_site_int(repeat,1)=str2double('NA');
%         adol_site_int(repeat,1)=str2double('NA');
%         adult_site_int(repeat,1)=str2double('NA');
%         end
%            
% % end
% % end
%         out(count,66)=mean(juv_site_nat_ref);
%         out(count,67)=mean(adol_site_nat_ref);
%         out(count,68)=mean(adult_site_nat_ref);
%         out(count,69)=mean(juv_site_nat_int);
%         out(count,70)=mean(adol_site_nat_int);
%         out(count,71)=mean(adult_site_nat_int);
%         out(count,72)=mean(juv_site_out_int);
%         out(count,73)=mean(adol_site_out_int);
%         out(count,74)=mean(adult_site_out_int);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intervention related variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


idx = strfind(int_scen,'PRA'); % percentage reef restored
area_deploy=str2double(int_scen(idx(1)+3:idx(1)+4))
if isnan(area_deploy)
out(count,75)=str2double(int_scen(idx(1)+3));
else
out(count,75)=str2double(int_scen(idx(1)+3:idx(1)+4));
end


idx = strfind(int_scen,'SD'); % stocking density used (recorded in the name)
out(count,76)=str2double(int_scen(idx(2)+2:idx(2)+4));
idx = strfind(int_scen,'HT'); % heat tolerance used (recorded in the name)
out(count,77)=str2double(int_scen(idx(1)+2));

% Record of when intervention took place in the model (check)
    first_int_year_index=find(DENSITY.GBR.INT_GBR.years_rest_begins==1);
    last_int_startyr=DENSITY.GBR.INT_GBR.years_last_rest_event(1,startyr_ind);

    if first_int_year_index<=startyr_ind;

        out(count,78)=first_int_year_index-startyr_ind; % there hasn't been an intervention
else
    last_int_startyr=DENSITY.GBR.INT_GBR.years_last_rest_event(1,startyr_ind);
last_int_startyr=last_int_startyr+1
        out(count,78)=last_int_startyr;
end

%Time between census year and intervention
last_int_census=DENSITY.GBR.INT_GBR.years_last_rest_event(1,censusyr_ind);
out(count,79)=last_int_census;

idx = strfind(int_scen,'RF'); % restoration frequency
out(count,80)=str2double(int_scen(idx(1)+2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disturbances --> record time position between restoration event, census
% year and time since disturbance

cyclones=DENSITY.GBR.DIST_GBR.cycl_cat_int(r,startyr_ind:censusyr_ind,s);
out(count,81)=nnz(cyclones); %number of cyclones
if nnz(cyclones)>0
    cyclonesizes=cyclones(1,cyclones(1,:)>0);
    out(count,82)=mean(cyclonesizes);
    [empty,cycloneyrs]=find(DENSITY.GBR.DIST_GBR.cycl_cat_int(r,:,s)>0);
    cycloneyrs_last_int =[];
    if cycloneyrs(1)>=first_int_year_index;
        out(count,83)=(first_int_year_index); % there hasn't been a cyclone before the intervention at or before startyear
    else
        for i=1:length(cycloneyrs);
            if cycloneyrs(1,i)<first_int_year_index;
                cycloneyrs_last_int=(first_int_year_index)-cycloneyrs(1,:);
                out(count,83)=min(cycloneyrs_last_int(cycloneyrs_last_int>=0));

            end

        end
    end

 % First cyclone after intervention relative to censusy year
     first_int_year_index=find(DENSITY.GBR.INT_GBR.years_rest_begins==1);

post_first_rest_cycl=[];
        if cycloneyrs(end)<=first_int_year_index;
            out(count,84)=(censusyr_ind-first_int_year_index); % there hasn't been a cyclone event since the census timepoint

        elseif cycloneyrs(1,1)>first_int_year_index && cycloneyrs(1,1)>censusyr_ind; % there was a cyclone event relative to the census year
             out(count,84) = (censusyr_ind-first_int_year_index);
        else
            for i=1:length(cycloneyrs);
                if cycloneyrs(1,i)>first_int_year_index && cycloneyrs(1,i)<=censusyr_ind;
                    if cycloneyrs(1,i)==censusyr_ind
                        post_first_rest_cycl(1,i)=cycloneyrs(1,i)-first_int_year_index;
                    else
                        post_first_rest_cycl(1,i)=cycloneyrs(1,i)-first_int_year_index;
                
                end
                end
            end
            out(count,84) = min(post_first_rest_cycl(post_first_rest_cycl>0));
        end
    
[empty,cycloneyrs]=find(DENSITY.GBR.DIST_GBR.cycl_cat_int(r,:,s)>0);

% next cyclone event after restoration event, relative to census year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
post_rest_cycl=[];

    if cycloneyrs(1,1)<censusyr_ind-last_int_census && cycloneyrs(1,end)>=censusyr_ind;

        out(count,85)=last_int_census; % there hasn't been a cyclone event since census year
    elseif cycloneyrs(end)<censusyr_ind-last_int_census
        out(count,85)=last_int_census;


                else
                    for i=1:length(cycloneyrs);
   if cycloneyrs(1,i)>=censusyr_ind-last_int_census && cycloneyrs(1,i)<=censusyr_ind;
        if cycloneyrs(1,i)==censusyr_ind
            post_rest_cycl(1,i)=last_int_census;
        elseif cycloneyrs(1,i)==censusyr_ind-last_int_census
                post_rest_cycl(1,i)=0;

        else
            post_rest_cycl(1,i)=cycloneyrs(1,i)-(censusyr_ind-last_int_census);

        end
                    if any(post_rest_cycl)>0
        out(count,85) = min(post_rest_cycl(post_rest_cycl>0));
                    else
                    out(count,45) = 0;
                    end
    end
end
        end
    
        last_cycl_bef_last_rest=[];

 % last cyclone event before last intervention --> not relevant for this
 % study but still recorded
%  if cycloneyrs(1,1)>censusyr_ind-last_int_census;
% %             out(count,46)=(censusyr_ind-last_int_census)-(startyr_ind); % there hasn't been a cyclone
%             out(count,86)=(censusyr_ind-last_int_census); % there hasn't been a cyclone
% 
% %             out(count,75)=sqrt(out(count,37));
%         else
%             for i=1:length(cycloneyrs);
%                 if cycloneyrs(1,i)>=startyr_ind && cycloneyrs(1,i)<=censusyr_ind-last_int_census;
%                     last_cycl_bef_last_rest(1,i)=(censusyr_ind-last_int_census)-cycloneyrs(1,i);
% %                     out(count,36)
%                 end
%             end
% %             out(count,39) = min(last_cyclone_event(last_cyclone_event>0));
%             out(count,86) = min(last_cycl_bef_last_rest);
% 
% %             out(count,75)=sqrt(out(count,37));
%  end

 % last cyclone event before census year
 last_cyclone_event=[];
 if cycloneyrs(1,1)>=censusyr_ind;
     out(count,87)=(censusyr_ind-startyr_ind)+1; % there hasn't been a cyclone since the census timepoint
 else
     for i=1:length(cycloneyrs);
                      if cycloneyrs(1,i)>=startyr_ind && cycloneyrs(1,i)<censusyr_ind;

             last_cyclone_event(1,i)=censusyr_ind-cycloneyrs(1,i);
         end
     end
 out(count,87) = min(last_cyclone_event);

        end
        else

            % incase there are no cyclone events at all (so specific for
            % bleaching scenarios
            out(count,82)=str2double('NA');
            out(count,83)=str2double('NA'); 
            out(count,84)=str2double('NA'); 
            out(count,85)=str2double('NA'); 
            out(count,86)=str2double('NA'); 
            out(count,87)=str2double('NA'); 
        end

        % mean coral mortality cyclone census year reference, intervention and
        % difference
out(count,88)=mean(DENSITY.GBR.DIST_GBR.cor_loss_cyc_ref(r,censusyr_ind,:),3);
out(count,89)=mean(DENSITY.GBR.DIST_GBR.cor_loss_cyc_int(r,censusyr_ind,:),3);
out(count,90)=mean(DENSITY.GBR.DIST_GBR.cor_loss_cyc_int_ref(r,censusyr_ind,:),3);


            %   DHW - Heatwave output results --> similar as for cyclones
        dhw=DENSITY.GBR.DIST_GBR.dhw_mag_int(r,startyr_ind:censusyr_ind,s);
        [empt,allbleachingevents]=find(DENSITY.GBR.DIST_GBR.dhw_mag_int(r,:,s)>0);
        dhws_period=dhw(dhw(1,:)>0);
        if nnz(dhws_period)>0
        out(count,91)=length(dhws_period);
        out(count,92)=mean(dhws_period);

        % time since last heat wave before intervention started--> is the
        % same for all scenarios
        allbleachingevents_last_int=[];
        if allbleachingevents(1)>=first_int_year_index;
            out(count,93)=(first_int_year_index); 
        else
            for i=1:length(allbleachingevents)
                if allbleachingevents(1,i)<first_int_year_index;
                     allbleachingevents_last_int=(first_int_year_index)-allbleachingevents(1,:);
                     out(count,93)=min(allbleachingevents_last_int(allbleachingevents_last_int>=0));
 end

            end
        end


 % First bleaching after intervention took place relative to censusy year
     first_int_year_index=find(DENSITY.GBR.INT_GBR.years_rest_begins==1);

post_first_rest_bleaching=[];
        if allbleachingevents(end)<=first_int_year_index; 
            out(count,94)=((censusyr_ind-first_int_year_index)); % there hasn't been a bleaching event since census year
        elseif allbleachingevents(1,1)>first_int_year_index && allbleachingevents(1,1)>censusyr_ind;
            out(count,94) = censusyr_ind-first_int_year_index;
        else
            for i=1:length(allbleachingevents);
                if allbleachingevents(1,i)>first_int_year_index && allbleachingevents(1,i)<=censusyr_ind;
                    if allbleachingevents(1,i)==censusyr_ind
                        post_first_rest_bleaching(1,i)=allbleachingevents(1,i)-(first_int_year_index);
                    else
                        post_first_rest_bleaching(1,i)=allbleachingevents(1,i)-(first_int_year_index);
                
                end
                end
            end
            out(count,94) = min(post_first_rest_bleaching(post_first_rest_bleaching>0));
        end
    

% next bleaching event after restoration event relative to census year
post_rest_bleaching=[];
        if allbleachingevents(1,1)<censusyr_ind-last_int_census && allbleachingevents(1,end)>=censusyr_ind
            out(count,95)=(last_int_census); % there hasn't been a bleaching event since census year

        elseif allbleachingevents(end)<censusyr_ind-last_int_census;
            out(count,95)=(last_int_census); 
        else
            for i=1:length(allbleachingevents);
                if allbleachingevents(1,i)>=censusyr_ind-last_int_census && allbleachingevents(1,i)<=censusyr_ind;
                    if allbleachingevents(1,i)==censusyr_ind
                        post_rest_bleaching(1,i)=last_int_census;
                     elseif allbleachingevents(1,i)==censusyr_ind-last_int_census
                post_rest_bleaching(1,i)=0;
   
                    else
                        post_rest_bleaching(1,i)=allbleachingevents(1,i)-(censusyr_ind-last_int_census);
                
                end
            if any(post_rest_bleaching)>0
            out(count,95) = min(post_rest_bleaching(post_rest_bleaching>0));
            else
             out(count,95) = 0;
            end

                end
            end
        end
    

 % last bleaching event before last intervention --> double data in this
 % case
 if allbleachingevents(1,1)>censusyr_ind-last_int_census;
            out(count,96)=(censusyr_ind-last_int_census); 
        else
            for i=1:length(allbleachingevents);
                if allbleachingevents(1,i)>=startyr_ind && allbleachingevents(1,i)<=censusyr_ind-last_int_census
                    last_bleaching_bef_last_rest(1,i)=(censusyr_ind-last_int_census)-allbleachingevents(1,i);
                end
            end
            out(count,96) = min(last_bleaching_bef_last_rest);
 end

 % last heatwave event before census year
 last_allbleachingevents_event=[];
      if allbleachingevents(1,1)>=censusyr_ind;

            out(count,97)=(censusyr_ind-startyr_ind)+1; % there hasn't been a bleaching event before census year
        else
            for i=1:length(allbleachingevents);
                if allbleachingevents(1,i)>=startyr_ind && allbleachingevents(1,i)<=censusyr_ind;
                    last_allbleachingevents_event(1,i)=censusyr_ind-allbleachingevents(1,i);
                end
            end
            out(count,97) = min(last_allbleachingevents_event);
 end
        else % in case the scenarios are only cyclones
            out(count,91)=str2double('NA');
           out(count,92)=str2double('NA');
            out(count,93)=str2double('NA'); 

            out(count,94)=str2double('NA'); 
            out(count,95)=str2double('NA'); 
            out(count,96)=str2double('NA'); 
            out(count,97)=str2double('NA'); 
        end
out(count,98)=mean(DENSITY.GBR.DIST_GBR.cor_loss_dhw_ref(r,censusyr_ind,:),3);
out(count,99)=mean(DENSITY.GBR.DIST_GBR.cor_loss_dhw_int(r,censusyr_ind,:),3);
out(count,100)=mean(DENSITY.GBR.DIST_GBR.cor_loss_dhw_int_ref(r,censusyr_ind,:),3);

% Mean deployment of coral outplants across modelling period
for i=1:10
outpl_dens_index{i}=find(DENSITY.GBR.INT_GBR.outplant_m2_int(r,:,i)>0);
outpl_dens{i}=mean(DENSITY.GBR.INT_GBR.outplant_m2_int(r,outpl_dens_index{i},:,:),3);
end
outpl_dens_comb=cat(2,outpl_dens{:});

out(count,101)=mean(outpl_dens_comb);
out(count,102)=censusyr-2023; %index census year relative to restoration year

out(count,103)=model_period(first_int_year_index); % double check when restoration took place

    out(count,104)=special_run_subset.Var4(r); % double check when restoration took place
    out(count,105)=DENSITY.GBR.INT_GBR.heattol;% heat tolerance appliede
    out(count,106)=special_run_subset.Var7(r); % larval supply double check

        count=count+1;
  
end

        if t==1
            out2=out;
        else
            out2=[out2;out];
        end


    end

% where to store output data
cd '/QRISdata/Q5785/Virtual_reef/version_080124/VR_1dist_1int_dhw_080824/Extract_080824'

writematrix(out2,[int_scen(1:end-11),'_censusyr_',num2str(censusyr_options(t)),'_extract.csv'])
save([int_scen(1:end-11),'_censusyr_',num2str(censusyr_options(t)),'_extract.mat'],'out2')
end
