% Script is used to visualise figure 1, 2 and 6 in the manuscript
% As well as Figure S4, S5 an S6 in supps

% Load scenarios with cyclone (can also be used for figures related to no
% disturbances
load('VR_1cyc_1int212col_36scen_prenorm_norm_fix_3LS_3SD_3PAR_170824.mat')
% Load scenarios with bleaching
load('VR_1dhw_1int212col_36scen_prenorm_norm_all_5_10_20_200824.mat')


% For figure 1 and 2, we need to isolate results related to scenarios in
% the absence of disturbances

% Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%
Dist_1int_data=[];
Demo_matrix_ref=[];
Demo_table_ref=[];
Demo_matrix_int=[];
Demo_table_int=[];

figure
tiledlayout(4,2, 'Padding','Compact')
unique_larv_sup=unique(Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.tot_larv_in_param_index);

for l=1:3 %plot trends for every larval supply setting
%     Pick 1 specific scenario of env, reef, int settings (see Table 1,
%     bold values)
Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm_1scen=Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm(Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.self_ret==0.28 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.tot_larv_in_param_index==unique_larv_sup(l) & ...
     Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.surface_WQ<0.5 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.size_dist_index==1 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.coral_cover_startyr==3 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.species_comb_index==1 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.prop_rest==10 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.heattol==3 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.stockdens==6.8,:);

% Only show results for scenarios without disturbance
Nodist_scen=15:15:size(Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm_1scen,1);

Dist_1int_data{l}=Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm_1scen(Nodist_scen,:);

% from reef to m2 for both reference and intervention (100 cells on 1 grid
% in ReefMod)
ben_juv_ref=(Dist_1int_data{l}.juv_reef_nat_ref_censusyr_100m2)./100;
ben_adol_ref=(Dist_1int_data{l}.adol_reef_nat_ref_censusyr_100m2)./100;
ben_adult_ref=(Dist_1int_data{l}.adult_reef_nat_ref_censusyr_100m2)./100;

Demo_table_ref{l}=table(Dist_1int_data{l}.years_postdeploy,ben_juv_ref,ben_adol_ref,ben_adult_ref);
Demo_matrix_ref{l}=[Dist_1int_data{l}.years_postdeploy,ben_juv_ref,ben_adol_ref,ben_adult_ref];

ben_juv_int=(Dist_1int_data{l}.juv_reef_nat_int_censusyr_100m2+Dist_1int_data{l}.juv_reef_outpl_int_censusyr_100m2)./100;
ben_adol_int=(Dist_1int_data{l}.adol_reef_nat_int_censusyr_100m2+Dist_1int_data{l}.adol_reef_outpl_int_censusyr_100m2)./100;
ben_adult_int=(Dist_1int_data{l}.adult_reef_nat_int_censusyr_100m2+Dist_1int_data{l}.adult_reef_outpl_int_censusyr_100m2)./100;

Demo_table_int{l}=table(Dist_1int_data{l}.years_postdeploy,ben_juv_int,ben_adol_int,ben_adult_int);
Demo_matrix_int{l}=[Dist_1int_data{l}.years_postdeploy,ben_juv_int,ben_adol_int,ben_adult_int];

end

%Figure 1 A
nexttile
errorbar(Dist_1int_data{1}.years_postdeploy,Dist_1int_data{1}.coral_cover_ref_censusyr_mean,Dist_1int_data{1}.coral_cover_ref_censusyr_stdev,'Color','r','LineWidth',1)
hold on
errorbar(Dist_1int_data{2}.years_postdeploy,Dist_1int_data{2}.coral_cover_ref_censusyr_mean,Dist_1int_data{2}.coral_cover_ref_censusyr_stdev,'--','Color','r','LineWidth',1)
hold on
errorbar(Dist_1int_data{3}.years_postdeploy,Dist_1int_data{3}.coral_cover_ref_censusyr_mean,Dist_1int_data{3}.coral_cover_ref_censusyr_stdev,'-.','Color','r','LineWidth',1)

hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
hold on
ylabel('Mean total coral cover (%)')
hold on
xlabel('Years post deployment')
axis([0 14 0 100])

%Figure 1 B
    nexttile
errorbar(Dist_1int_data{1}.years_postdeploy,Dist_1int_data{1}.coral_cover_int_censusyr_mean,Dist_1int_data{1}.coral_cover_int_censusyr_stdev,'Color','b','LineWidth',1)
hold on
errorbar(Dist_1int_data{2}.years_postdeploy,Dist_1int_data{2}.coral_cover_int_censusyr_mean,Dist_1int_data{2}.coral_cover_int_censusyr_stdev,'--','Color','b','LineWidth',1)
hold on
errorbar(Dist_1int_data{3}.years_postdeploy,Dist_1int_data{3}.coral_cover_int_censusyr_mean,Dist_1int_data{3}.coral_cover_int_censusyr_stdev,'-.','Color','b','LineWidth',1)

hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
hold on
ylabel('Mean total coral cover (%)')
hold on
xlabel('Years post deployment')
axis([0 14 0 100]) 
yticks(0:20:100)

%Figure 1 C
% Counterfactual juveniles
    nexttile
plot(Demo_table_ref{1}.Var1,Demo_table_ref{1}.Var2,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
plot(Demo_table_ref{2}.Var1,Demo_table_ref{2}.Var2,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
plot(Demo_table_ref{3}.Var1,Demo_table_ref{3}.Var2,'-.','Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
hold on
ylabel('Mean # individual/m2')
hold on
xlabel('Years post deployment')
axis([0 14 0 25]) 

%Figure 1 D
% intervention juveniles
    nexttile
plot(Demo_table_int{1}.Var1,Demo_table_int{1}.Var2,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
plot(Demo_table_int{2}.Var1,Demo_table_int{2}.Var2,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
plot(Demo_table_int{3}.Var1,Demo_table_int{3}.Var2,'-.','Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
hold on
ylabel('Mean # individual/m2')
hold on
xlabel('Years post deployment')
axis([0 14 0 25])

%Figure 1 E
% Counterfactual subadults
    nexttile
plot(Demo_table_ref{1}.Var1,Demo_table_ref{1}.Var3,'Color',[0.8 0.1 0.8],'LineWidth',1)
hold on
plot(Demo_table_ref{2}.Var1,Demo_table_ref{2}.Var3,'--','Color',[0.8 0.1 0.8],'LineWidth',1)
hold on
plot(Demo_table_ref{3}.Var1,Demo_table_ref{3}.Var3,'-.','Color',[0.8 0.1 0.8],'LineWidth',1)
hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
hold on
ylabel('Mean # individual/m2')
hold on
xlabel('Years post deployment')
axis([0 14 0 25]) %species1


%Figure 1 F
% intervention subadults
    nexttile
plot(Demo_table_int{1}.Var1,Demo_table_int{1}.Var3,'Color',[0.8 0.1 0.8],'LineWidth',1)
hold on
plot(Demo_table_int{2}.Var1,Demo_table_int{2}.Var3,'--','Color',[0.8 0.1 0.8],'LineWidth',1)
hold on
plot(Demo_table_int{3}.Var1,Demo_table_int{3}.Var3,'-.','Color',[0.8 0.1 0.8],'LineWidth',1)
hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
hold on
ylabel('Mean # individual/m2')
hold on
xlabel('Years post deployment')
axis([0 14 0 25])

%Figure 1 G
% counterfactual adults
    nexttile
plot(Demo_table_ref{1}.Var1,Demo_table_ref{1}.Var4,'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
hold on
plot(Demo_table_ref{2}.Var1,Demo_table_ref{2}.Var4,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',1)
hold on
plot(Demo_table_ref{3}.Var1,Demo_table_ref{3}.Var4,'-.','Color',[0.4660 0.6740 0.1880],'LineWidth',1)
hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
hold on
ylabel('Mean # individual/m2')
hold on
xlabel('Years post deployment')
axis([0 14 0 25])

%Figure 1 H
% intervention adults
    nexttile
plot(Demo_table_int{1}.Var1,Demo_table_int{1}.Var4,'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
hold on
plot(Demo_table_int{2}.Var1,Demo_table_int{2}.Var4,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',1)
hold on
plot(Demo_table_int{3}.Var1,Demo_table_int{3}.Var4,'-.','Color',[0.4660 0.6740 0.1880],'LineWidth',1)
hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
hold on
ylabel('Mean # individual/m2')
hold on
xlabel('Years post deployment')
axis([0 14 0 25]) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 2
Nodist_scen=15:15:size(Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm,1);
nodist_scen_all=Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm(Nodist_scen,:);


ben_nodist=[];
year=2023:1:2037;
for yr=1:15
censusyr=find(nodist_scen_all.censusyr == year(yr));
ben_nodist(:,yr)=table2array(nodist_scen_all(censusyr,"coral_cover_int_ref_censusyr_mean"));
end

y=0:14;
figure
boxplot(ben_nodist,y, 'Colors',color("ForestGreen"),'PlotStyle','compact','symbol', '');

ylim([-2.7 7])
xticks(1:1:15)

test=[{"0", "1", "2", "3", "4","5", "6", "7", "8", "9","10", "11", "12", "13", "14"}]
xticklabels(test)
set(gca,'FontSize',10,'XTickLabelRotation',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6

% load specific scenarios to go the raw data as extracted data based on
% average data over replicates

%directory of raw data of cyclone scenarios= VR_1dist_1int_cyc_080824
load('VR_OC_SD_LN2_EYR_2023_PRA10_RF1_SD6.8_HT3_demo_080824.mat') % figure 6 top 4 panels (+3DHW outplants)
load('VR_OC_SD_LN2_EYR_2023_PRA10_RF1_SD6.8_HT0_demo_080824.mat') % figure S4 top 4 panels (+0DHW outplants)

count=1;
total_benefit=[];
sequence=[1:15; 121:135; 166:180; 46:60]; % represent scenarios with low LS (larval supply)- low SR (self retention).
% low LS- high SR, high LS - low SR, high LS- high SR, respectively
for panel=1:4
    for reef=sequence(panel,:)

        for rep=1:10
            total_benefit(count,1)=reef;
            total_benefit(count,2)=rep;
            total_benefit(count,3:18)=DENSITY.GBR.COVER_GBR.coral_pct_int_ref(reef,:,rep);
            count=count+1;
        end
    end
end
writematrix(total_benefit,'Ben_ind_rep_1cyc_default_setting_220425.csv')   
writematrix(total_benefit,'Ben_ind_rep_1cyc_default_setting_0DHW_final_290425.csv') 


% same for bleaching
%directory of raw data of cyclone scenarios= VR_1dist_1int_dhw_080824
load('VR_OC_SD_LN2_EYR_2023_PRA10_RF1_SD6.8_HT3_demo_080824.mat') % figure 6 top 4 panels (+3DHW outplants)
load('VR_OC_SD_LN2_EYR_2023_PRA10_RF1_SD6.8_HT0_demo_080824.mat') % figure S4 top 4 panels (+0DHW outplants)


count=1;
total_benefit=[];
sequence=[1:15; 121:135; 166:180; 46:60]; % represent scenarios with low LS (larval supply)- low SR (self retention).
% low LS- high SR, high LS - low SR, high LS- high SR, respectively
for panel=1:4
    for reef=sequence(panel,:)

        for rep=1:10
            total_benefit(count,1)=reef;
            total_benefit(count,2)=rep;
            total_benefit(count,3:18)=DENSITY.GBR.COVER_GBR.coral_pct_int_ref(reef,:,rep);
            count=count+1;
        end
    end
end
writematrix(total_benefit,'Ben_ind_rep_1dhw_default_setting_220425.csv')  
writematrix(total_benefit,'Ben_ind_rep_1dhw_default_setting_0DHW_final_290425.csv')    

% Afterwards plot data
%%%%%%%%
% corresponding to the scenarios under certain larval supply and self
% retention.
cyc=1;
if cyc==1 % first plot cyclones, than change "cyc=0", and run bleaching results to have them on the same figure
    total_ben=load("Ben_ind_rep_1cyc_default_setting_220425.csv")
    total_ben=load("Ben_ind_rep_1cyc_default_setting_0DHW_final_290425.csv")

else
    total_ben=load("Ben_ind_rep_1dhw_default_setting_220425.csv")
    total_ben=load("Ben_ind_rep_1dhw_default_setting_0DHW_final_290425.csv")

end

% Select rows (scenarios that have disturbance event occurs 1-13 years post
% deployment, Y values correspond to benefits after disturbance (2-14).
% Other vectos are down the list of 2160 scenarios but correspond to the
% same disturbance sequence only different LS and SR setting
LS_L_SR_L=2:14;
LS_L_SR_H=122:134;
LS_H_SR_H=167:179;
LS_H_SR_L=47:59;

plots=[LS_L_SR_L;LS_L_SR_H;LS_H_SR_L;LS_H_SR_H];

for m=1:4
Dist_1int_data=table();
Dist_1int_data.Position=(plots(m,:))';
Dist_1int_data.YearDisturbancePostDeployment=flip(1:13)'
Dist_1int_data.censusyr2=flip(2025:2037)';
Dist_1int_data.censusyr=flip(6:18)';

for i=1:13 %calculate stats
        
Dist_1int_data_sub=total_ben(total_ben(:,1)==Dist_1int_data.Position(i),:);
Dist_1int_data.coral_cover_int_ref_censusyr_mean(i)=mean(Dist_1int_data_sub(:,Dist_1int_data.censusyr(i)));
Dist_1int_data.coral_cover_int_ref_censusyr_median(i)=median(Dist_1int_data_sub(:,Dist_1int_data.censusyr(i)));
Dist_1int_data.coral_cover_int_ref_censusyr_25perc(i)=prctile(Dist_1int_data_sub(:,Dist_1int_data.censusyr(i)),25);
Dist_1int_data.coral_cover_int_ref_censusyr_75perc(i)=prctile(Dist_1int_data_sub(:,Dist_1int_data.censusyr(i)),75);

end

nexttile
%plot only median and interquartile interval
if cyc==1
h1=plot(Dist_1int_data.YearDisturbancePostDeployment,Dist_1int_data.coral_cover_int_ref_censusyr_median,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
plotColor = get(h1, 'Color');
hold on
curve1 = Dist_1int_data.coral_cover_int_ref_censusyr_25perc';
curve2 = Dist_1int_data.coral_cover_int_ref_censusyr_75perc';

h2=plot(Dist_1int_data.YearDisturbancePostDeployment,curve1,"--");
set(h2, 'Color', plotColor);
hold on
x= (13:-1:1);
x1 = 1:13;

x2 = [x, x1];
inBetween1 = [curve1(1,:), fliplr(curve2(1,:))];
fill(x2, inBetween1, [plotColor], 'FaceAlpha',.15,'EdgeAlpha',.15);
hold on
h3=plot(Dist_1int_data.YearDisturbancePostDeployment,curve2,"--");
set(h3, 'Color', plotColor);

hold on
plot(Dist_1int_data.YearDisturbancePostDeployment,Dist_1int_data.coral_cover_int_ref_censusyr_median,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5)

set(h3, 'Color', plotColor);
hold on
ylabel('Mean coral cover benefit after cyclone event (%)');
xlabel('Year of cyclone event occuring post deployment (years)');

else %bleaching in different color
h1=plot(Dist_1int_data.YearDisturbancePostDeployment,Dist_1int_data.coral_cover_int_ref_censusyr_median,'Color','b','LineWidth',1.5);
plotColor = get(h1, 'Color');
hold on
curve1 = Dist_1int_data.coral_cover_int_ref_censusyr_25perc';
curve2 = Dist_1int_data.coral_cover_int_ref_censusyr_75perc';

h2=plot(Dist_1int_data.YearDisturbancePostDeployment,curve1,"--");
set(h2, 'Color', plotColor);
hold on
x= (13:-1:1);
x1 = 1:13;

x2 = [x, x1];
inBetween1 = [curve1(1,:), fliplr(curve2(1,:))];
fill(x2, inBetween1, [plotColor], 'FaceAlpha',.15,'EdgeAlpha',.15);
hold on
h3=plot(Dist_1int_data.YearDisturbancePostDeployment,curve2,"--");
set(h3, 'Color', plotColor);

hold on
plot(Dist_1int_data.YearDisturbancePostDeployment,Dist_1int_data.coral_cover_int_ref_censusyr_median,'Color','b','LineWidth',1.5)

set(h3, 'Color', plotColor);
hold on

ylabel('Mean coral cover benefit after bleaching event (%)');
xlabel('Year of bleaching event occuring post deployment (years)');
end
hold on
plot(13:-1:0, repelem(0,14), 'Color', 'black', 'LineWidth',0.5,'HandleVisibility','off');
hold on
plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10); % marks timepoint of outplanting
axis([0 13 -3.5 15]) 
end
text(1,1,'Random variables'); % used to label x axis vertically, remove all other y axis labels for editing (keep only 1 x axis label)
H=findobj(gca,'Type','text');
 set(H,'Rotation',90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure S1

figure
tiledlayout(1,2, 'Padding','Compact')
small=1; % predominantly small individuals

    for i=2:2 % plot species 2


sp_mean=[17.50380286;24.31965649;17.50380286;11.29180602;11.94162437;16.17448777]; % mean diameter size based on Dietzel et al., 2020
sp_sigma=[16.86639;18.20781;16.86639;7.970298;14.3261;24.56004];

% AA: coral group selection for transect measurements (see excel file Intersept_data), only focussing on recent measurements of years 2016 and 2017 (Dietzel et al., 2021):
%     1.  in Dietzel: other Acropora
%     2.  in Dietzel: tabular Acropora
%     3. in Dietzel: other Acropora
%     4. in Dietzel: P. damicornis, other Pocillopora, Stylophora
%     5. in Dietzel: Faviidae
%     6. in Dietzel: Porites
% 
species_max_adult_size=[7853; 7853; 1963; 1256; 2827; 7853]; % see Bozec et al., 2022

nexttile
%option 1
if small==1
m=sp_mean(i) % skewed to right
v= sp_sigma(i)
else % predominantly mid to large individuals
m=sp_mean(i)*5 % skewed to left
v= sp_sigma(i)/5
end

rng(123)
s=lognrnd(log(pi*(m/2)^2), log(v),1,2000);
s(s > species_max_adult_size(i)) = []; % bound by max size of species
s(s < 2) = 2;
d=sqrt(s/pi)*2;
hist(d,0:10:100);;
xlabel('Diameter size (cm2)');
ylabel('Count');
end
exportgraphics(gcf,['Lognormal_dist_mean_log10000_std_log2_species_',num2str(i),'.png']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure S2

% Test to visualise recruitment: data include=larval supply, SR,
% #larval in, # larvae out, # recruits, # juv, # coral cover (model scripts
% and extraction script can be provided)
larv_supp_recr2=readtable('VR_OC_SD_LN1_EYR_2023_PRA5_RF1_SD3.4_HT0_demo_150124_extract_300124.csv'); % need to add larval supply 2023 and coral cover state 2022 from Reference model run in order visualise the data below

% Visualisation

y=0:14

self_ret=[0.01,0.07, 0.28]
figure
tiledlayout(3,1, 'Padding','Compact')

for i=1:3
    % Counterfactual
Y2_1500_allspecies=larv_supp_recr2(larv_supp_recr2.selfRetention==self_ret(i) & ...
    larv_supp_recr2.sizeDistributionIndex==1 & ...
    larv_supp_recr2.larvalSupply2023==1500 & ...
    larv_supp_recr2.WQ_surface<0.05 & ...
     larv_supp_recr2.Species_comp==1 & ...
    larv_supp_recr2.coral_cover_Ref_2022==3,:);


Y2_150000_allspecies=larv_supp_recr2(larv_supp_recr2.selfRetention==self_ret(i) & ...
    larv_supp_recr2.sizeDistributionIndex==1 & ...
    larv_supp_recr2.larvalSupply2023==150000 & ...
    larv_supp_recr2.WQ_surface<0.05 & ...
     larv_supp_recr2.Species_comp==1 & ...
    larv_supp_recr2.coral_cover_Ref_2022==3,:);

Y2_1500000_allspecies=larv_supp_recr2(larv_supp_recr2.selfRetention==self_ret(i) & ...
    larv_supp_recr2.sizeDistributionIndex==1 & ...
    larv_supp_recr2.larvalSupply2023==1500000 & ...
    larv_supp_recr2.WQ_surface<0.05 & ...
     larv_supp_recr2.Species_comp==1 & ...
    larv_supp_recr2.coral_cover_Ref_2022==3,:);

recr1=Y2_1500_allspecies(:,[37:52]); %recruitment under low LS

recr2=Y2_150000_allspecies(:,[37:52]); %recruitment under mid LS

recr3=Y2_1500000_allspecies(:,[37:52]); %recruitment under highLS

nexttile

plot(y, table2array(recr1(15,2:end)),'LineWidth',1,'Color',[0 1 1]); % scenario without disturbance
hold on
plot(y, table2array(recr5(15,2:end)),'LineWidth',1,'Color',[0 0.5 0.5]); %juveniles
hold on
plot(y, table2array(recr6(15,2:end)),'LineWidth',1,'Color',[0 0 1]); %juveniles

hold on
plot(y(1),0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
axis([0 14 0 12])

xlabel('Year post deployment')
ylabel('# recruits (ind/m2)')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S3 in R
% Figure S4 (see Figure 6)

% Figure S5

%14/04/24
Extract_final_OC_VR_1dhw_1int_5_10_20_200824_1scen=Extract_final_OC_VR_1dhw_1int_5_10_20_200824(Extract_final_OC_VR_1dhw_1int_5_10_20_200824.self_ret==0.28 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.tot_larv_in_param_index==1500000 & ...
     Extract_final_OC_VR_1dhw_1int_5_10_20_200824.surface_WQ<0.5 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.size_dist_index==1 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.coral_cover_startyr==3 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.species_comb_index==1 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.prop_rest==10 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.heattol==3 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.stockdens==6.8,:);


Dist_1int_data=Extract_final_OC_VR_1dhw_1int_5_10_20_200824_1scen;

YearDisturbancePostDeployment=repmat(flip(0:14),1,15);
Dist_1int_data.YearDisturbancePostDeployment=YearDisturbancePostDeployment';

% Figure S5 A
figure
h=heatmap(Dist_1int_data,'years_postdeploy','YearDisturbancePostDeployment','ColorVariable','coral_cover_int_ref_censusyr_mean','ColorLimits',[0 6])
colormap(flipud(summer(7)))
h.CellLabelFormat = '%.1f';
h.YDisplayData = flipud(h.YDisplayData) %change orientation of yaxis
ylabel('Year of bleaching event occuring post deployment')
xlabel ('Years post deployment')

% Figure S5 B
figure
h=heatmap(Dist_1int_data,'years_postdeploy','YearDisturbancePostDeployment','ColorVariable','coral_cover_ref_censusyr_mean','ColorLimits',[0 90])
colormap(flipud(summer(7)))
h.CellLabelFormat = '%.1f';
h.YDisplayData = flipud(h.YDisplayData) %change orientation of yaxis
ylabel('Year of bleaching event occuring post deployment')
xlabel ('Years post deployment')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S6
% Demographic plots
% Cyclone
Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm

Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm_1scen=Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm(Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.self_ret==0.28 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.tot_larv_in_param_index==1500000 & ...
     Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.surface_WQ<0.5 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.size_dist_index==1 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.coral_cover_startyr==3 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.species_comb_index==1 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.prop_rest==10 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.heattol==3 & ...
    Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.stockdens==6.8,:);

Dist_1int_data=Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm_1scen;

% Bleaching
Extract_final_OC_VR_1dhw_1int_5_10_20_200824_1scen=Extract_final_OC_VR_1dhw_1int_5_10_20_200824(Extract_final_OC_VR_1dhw_1int_5_10_20_200824.self_ret==0.28 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.tot_larv_in_param_index==1500000 & ...
     Extract_final_OC_VR_1dhw_1int_5_10_20_200824.surface_WQ<0.5 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.size_dist_index==1 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.coral_cover_startyr==3 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.species_comb_index==1 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.prop_rest==10 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.heattol==3 & ...
    Extract_final_OC_VR_1dhw_1int_5_10_20_200824.stockdens==6.8,:);


Dist_1int_data=Extract_final_OC_VR_1dhw_1int_5_10_20_200824_1scen;

time_dist=[15 10 6]; % 15 corresponding to scenario without bleaching
% 10 corresponds to scenario with bleaching on year 9
% 6 corresponds to scenario with bleaching on year 5
figure
tiledlayout(3,1, 'Padding','Compact')

for dist=1:3
dist_scen10=time_dist(dist):15:225 %no bleaching event

Dist_1int_HLS_HSR_dhw_sub=Dist_1int_data(dist_scen10,:);

y=unique(Dist_1int_HLS_HSR_dhw_sub.years_postdeploy)
ben_juv=[];
ben_adol=[];
ben_adult=[];

% Combine demographic classes for natural and outplanted corals for the
% intervention model run (INT). Counterfactual only consists of natural
% corals
ben_juv=((Dist_1int_HLS_HSR_dhw_sub.juv_reef_nat_int_censusyr_100m2+Dist_1int_HLS_HSR_dhw_sub.juv_reef_outpl_int_censusyr_100m2)-Dist_1int_HLS_HSR_dhw_sub.juv_reef_nat_ref_censusyr_100m2)./100;
ben_adol=((Dist_1int_HLS_HSR_dhw_sub.adol_reef_nat_int_censusyr_100m2+Dist_1int_HLS_HSR_dhw_sub.adol_reef_outpl_int_censusyr_100m2)-Dist_1int_HLS_HSR_dhw_sub.adol_reef_nat_ref_censusyr_100m2)./100;
ben_adult=((Dist_1int_HLS_HSR_dhw_sub.adult_reef_nat_int_censusyr_100m2+Dist_1int_HLS_HSR_dhw_sub.adult_reef_outpl_int_censusyr_100m2)-Dist_1int_HLS_HSR_dhw_sub.adult_reef_nat_ref_censusyr_100m2)./100;

Demo_matrix=[];
Demo_table=[];
Demo_table=table(ben_juv,ben_adol,ben_adult);
Demo_matrix=[ben_juv,ben_adol,ben_adult];
x = 0:14;
z = 1:3;
nexttile
plot(y, repelem(0,15), 'Color', 'black', 'LineWidth',0.5,'HandleVisibility','off'); % plot horizontal line

for i=1:3
plot(x,table2array(Demo_table(:,i)),'LineWidth',2) % plot in 3D space
hold on
end
axis([0 14 -0.5 2.5])

xlabel('Years post deployment')
zlabel('Benefit in # individual/m2')
hold on

plot(0,0,"-p",'MarkerEdgeColor','#008744','MarkerFaceColor','#008744','MarkerSize',10)
plot(y, repelem(0,15), 'Color', 'black', 'LineWidth',0.5,'HandleVisibility','off');
hold on
if dist==2
plot(5,0,"-p",'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
elseif dist==3
plot(9,0,"-p",'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%