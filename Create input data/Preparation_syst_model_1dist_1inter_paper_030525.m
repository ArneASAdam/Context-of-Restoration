%%%%%%% Preparation of file

% load file with all options, for this you need to be in the folder where this folder has been saved 
addpath('1cyc_1int_240124') % for cyclone scenarios
addpath('1dhw_1int_240124') % for bleaching scenarios

load('Spec_sub_1cyc_3init_3self_4larv_paper_test_150124.mat') % for cyclone scenarios
load('Spec_sub_1dhw_3init_3self_4larv_paper_test_150124.mat') % for bleaching scenarios

special_run_subset=special_run_subset2;

%To set unique scenarios, we need to create specific input data files using
%the "special_run.m" approach (see manual)

%Reef_ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
special_reef_ID=table();
shelf=1;

% geomorphology
% Reef Area, sand, shelf position
gbr_spatial=readtable("Reef_emulator_virtual_June_2023/rme_ml_2023_06_27/data_files/id/id_list_2023_03_30.csv");

gbr_spatial.Properties.VariableNames={'Reef_ids','Reef_area_CH','prop_ungraz','Shelf_position'};
head(gbr_spatial);
gbr_spatial_shelf=(gbr_spatial(table2array(gbr_spatial(:,4))==1,:));

% reef_area_stats
reef_area_quantile = quantile(table2array(gbr_spatial_shelf(:,2)),[0.05 0.25 0.5 0.75 0.95]); %median value will be used for the model (fixed value)

% Preparation of the ID input file
for a=1:size(special_run_subset,1)
special_reef_ID(a,:)={['REEF',num2str(a)], 0,0,0};
special_reef_ID(a,2)={reef_area_quantile(1,3)}; 
special_reef_ID(a,3)={0.10}; % sand cover at 10%
special_reef_ID(a,4)={1}; % shelf number 
end

writetable(special_reef_ID,['special_reef_ID_240124.csv'],'Delimiter',',','WriteVariableNames',false);  % is not actually produced but is here used as a script design for the paper

%%%%% coral connectivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%coral connectivity
special_con_self_ret=diag(special_run_subset.Var9); %self retention in diagonal, remaining values are 0 as connectivity is defined by external larval supply
writematrix(special_con_self_ret,['special_con_self_ret_240124.csv'],'Delimiter',','); 

% larval supply will be supply outside the data files produced and will be
% based on the combination file "spec...."

%Initial coral cover settings will depend on species composition
init_shelf_quantile=[];
for sp=1:6
init_shelf_quantile{sp}=[0.5;1;1.5];

end

special_coral_init=[];
for sp=1:6
    special_coral_init{sp}=table;
    for a=1:size(special_run_subset,1)

        if special_run_subset.Var4(a)==2 & sp==1 % if special_run_subset.Var4(a)==2, only species 4, 5 and 6 are present (no acroporiidae)

            special_coral_init{sp}(a,:)={['REEF',num2str(a)],0};
            special_coral_init{sp}(a,2)={0};

        elseif special_run_subset.Var4(a)==2 & sp==2 % if special_run_subset.Var4(a)==2, only species 4, 5 and 6 are present (no acroporiidae)

            special_coral_init{sp}(a,:)={['REEF',num2str(a)],0};
            special_coral_init{sp}(a,2)={0};

        elseif special_run_subset.Var4(a)==2 & sp==3 % if special_run_subset.Var4(a)==2, only species 4, 5 and 6 are present (no acroporiidae)

            special_coral_init{sp}(a,:)={['REEF',num2str(a)],0};
            special_coral_init{sp}(a,2)={0};
            
        elseif special_run_subset.Var4(a)==2 & sp==4 % if special_run_subset.Var4(a)==2, only species 4, 5 and 6 are present (no acroporiidae) in double amount to get the parameterised total coral cover to either 3, 6 or 9

            special_coral_init{sp}(a,:)={['REEF',num2str(a)],0};
            special_coral_init{sp}(a,2)={init_shelf_quantile{sp}(special_run_subset.Var3(a)).*2};

        elseif special_run_subset.Var4(a)==2 & sp==5 % if special_run_subset.Var4(a)==2, only species 4, 5 and 6 are present (no acroporiidae) in double amount to get the parameterised total coral cover to either 3, 6 or 9

            special_coral_init{sp}(a,:)={['REEF',num2str(a)],0};
            special_coral_init{sp}(a,2)={init_shelf_quantile{sp}(special_run_subset.Var3(a)).*2};

        elseif special_run_subset.Var4(a)==2 & sp==6 % if special_run_subset.Var4(a)==2, only species 4, 5 and 6 are present (no acroporiidae) in double amount to get the parameterised total coral cover to either 3, 6 or 9

            special_coral_init{sp}(a,:)={['REEF',num2str(a)],0};
            special_coral_init{sp}(a,2)={init_shelf_quantile{sp}(special_run_subset.Var3(a)).*2};

        else

            special_coral_init{sp}(a,:)={['REEF',num2str(a)],0};
            special_coral_init{sp}(a,2)={init_shelf_quantile{sp}(special_run_subset.Var3(a))};
        end
    end

    writetable(special_coral_init{sp},['special_sp',num2str(sp),'_init_cover_240124.csv'],'Delimiter',',','WriteVariableNames',false);

end

% initial cots cover (not used as assumption is made to remove cots if in
% outbreak number on the restored reefs)
%%%%%%%%%%%%%%%%%%%%%%%%%
special_cots_init=table;
for a=1:size(special_run_subset,1)

    special_cots_init(a,:)={['REEF',num2str(a)],0};
    special_cots_init(a,2)={(special_run_subset.Var5(a)/100)*400}; % because in reefmod Engine, proportion needs to be converted from gridsize source code
end
writetable(special_cots_init,['special_cots_init_cover_240124.csv'],'Delimiter',',','WriteVariableNames',false);

% initial rubble cover
%%%%%%%%%%%%%%%%%%%%%%%
special_rubble_init=table;
for a=1:size(special_run_subset,1)

    special_rubble_init(a,:)={['REEF',num2str(a)],0};
    special_rubble_init(a,2)={special_run_subset.Var6(a)};
end
writetable(special_rubble_init,['special_rubble_init_cover_22024.csv'],'Delimiter',',','WriteVariableNames',false);

%suspended sedimentation file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
special_surface_sed_init=table();
for a=1:size(special_run_subset,1)

    special_surface_sed_init(a,:)={['REEF',num2str(a)],0};
    special_surface_sed_init(a,2)={special_run_subset.Var7(a)};


end
special_surface_sed_init.Properties.VariableNames={'id','2011'}; % necessary for parameterisation, year 2011 is arbitrary in this case

writetable(special_surface_sed_init,['special_sed_surface_240124.csv'],'Delimiter',',');

%middepth
special_middepth_sed_init=table();
for a=1:size(special_run_subset,1)

    special_middepth_sed_init(a,:)={['REEF',num2str(a)],0};
        special_middepth_sed_init(a,2)={special_run_subset.Var8(a)};

%     end
end
special_middepth_sed_init.Properties.VariableNames={'id','2011'};  % necessary for parameterisation, year 2011 is arbitrary in this case

writetable(special_middepth_sed_init,['special_sed_middepth_240124.csv'],'Delimiter',',');

%COTS larval reduction file (not used in this study but still need to be
%parameterised)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shelf==1
    special_COTS_larval_red=table();
    for a=1:size(special_run_subset,1)

        special_COTS_larval_red(a,:)={['REEF',num2str(a)],0};
        special_COTS_larval_red(a,2)={1};
    end
end
special_COTS_larval_red.Properties.VariableNames={'id','2011'};
writetable(special_COTS_larval_red,['special_COTS_red_factor_240124.csv'],'Delimiter',','); % is not used in this case because COTS are turned off

%Disturbances sequence file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Only cyclones
special_cyc=table();
if special_run_subset.Var12==1

    for a=1:size(special_run_subset,1)
        special_cyc(a,:)={['REEF',num2str(a)],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        for i=1:16
            special_cyc(a,i+1)={table2array(special_run_subset(a,i+12))};
        end

    end

else

    for a=1:size(special_run_subset,1)
        special_cyc(a,:)={['REEF',num2str(a)],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        for i=1:16
            special_cyc(a,i+1)={0};
        end
    end
end

writetable(special_cyc,['special_cyc_240124.csv'],'Delimiter',',','WriteVariableNames',false);

% only bleaching events (DHW trajectory)
special_DHW=table();
if special_run_subset.Var12==2

    for a=1:size(special_run_subset,1)

        special_DHW(a,:)={['REEF',num2str(a)],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        for i=1:16
            special_DHW(a,i+1)={table2array(special_run_subset(a,i+12))};

        end

    end

    special_DHW.Properties.VariableNames(1)={'id'}; %years are necessary to know what DHW to implement (0 or 8 DHW)
    for i=1:size(special_DHW,2)-1
        special_DHW.Properties.VariableNames{i+1}=model_yrs_char(:,:,i);
    end
else
    for a=1:size(special_run_subset,1)
        special_DHW(a,:)={['REEF',num2str(a)],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        for i=1:16
            special_DHW(a,i+1)={0};
        end

    end
    special_DHW.Properties.VariableNames(1)={'id'};
    for i=1:size(special_DHW,2)-1
        special_DHW.Properties.VariableNames{i+1}=model_yrs_char(:,:,i);
    end
end
writetable(special_DHW,['special_DHW_240124.csv'],'Delimiter',',');
