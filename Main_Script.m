%% Main Script
% This script functions as main script. The script loads the model
% parameters and runs the simulink models. Then it makes plots
clear all 
close all
clc
% delete(findall(0));

%% Explanation
% This script runs the models, calibrates the model, and can be used to
% make plots. 
% In the Simulink files one can also see the plots, however those are not
% callibrated to the real-world values.
% The simulink scripts serve as a tool with wich you can play around
% First the three networks are simulated and plottes seperately
% Then I've included some plots that compare the different models. 
% At the end I've included the comparison of the controlled and
% uncontrolled core network

% IMPORTANT: Before running this scipt check the following:
% 1) Make sure that the path to find the
% Excel file is correct in the Import_KNMI_Wind_Data_2019() function
% 2) Make sure that the path to find the
% Excel file is correct in the mport_TNO_Price_Calibration_Data() function
% 3) Make sure the Simulink files model the scenario you want to model
% (standard it is set to model the one year wind distribution. You can add
% shocks, use other inputs, etc.

%% Stop Time

years=1; % years of simulation
Z=years*365; % in days

%% Get model parameters
run("Model_Parameters.m"); % Load model parameters

%% Run Base Model SS
Core_Cal=sim('Core_Model_Cal',Z); 
X=Core_Cal.time; 
%% Fitting Fourier Series on average wind distribution
%See Fourier_Series_Wind.m
%This is integrated as a function in Simulink
%% Import one year data wind distribution KNMI 2019 

Wind_data_2019 = Import_KNMI_Wind_Data_2019(); %Function that extracts data from excel
Wind_Distr_2019 = Wind_data_2019{:,3};

SF_wind=1.1135; %Scaling factor to get the wind speeds at hub hight (116m)

Wind_Distr_2019_Scaled=Wind_Distr_2019.*SF_wind;
Wind_Distr_2019_input=[X,Wind_Distr_2019_Scaled()]; %Input Wind distribution Simulink

%% Run Models for Scenario analysis

% Base Model WD
Core=sim('Core_Model',Z); %Wind Distribution

% Other Market Designs
%% Model Calibration

% Flow calibration
timesteps=numel(X);
Coulomb_365=Core_Cal.total_supply(8761); 

% Scaling
D_year=280; %yearly demand in PJ
HHV=141.9; %HHV H2 MJ/kg
PJMWH=0.27778e6;  %1 PJ equals 277.78 GWh

MWh_year=D_year*PJMWH; %in MWh
GWh_year=MWh_year*1e-3; %in GWh


Coulomb_sf=(GWh_year/8760)/(Coulomb_365/Z); %1 Coulomb = 0.676 GWh

SF_flow_test=Coulomb_sf;

capacitycheck=GWh_year/8760;

hpy=numel(X);%hours


% Price calibration
% Import data from spreadsheet
IELGAS_cal = Import_TNO_Price_Calibration_Data(); 
Time_cal= IELGAS_cal{:,1};
H2_35_cal= IELGAS_cal{:,2}; %X
elec_35_cal= IELGAS_cal{:,3};
NG_35_cal= IELGAS_cal{:,4};

H2_price_MWh_cal=Core_Cal.H2_price(:);
H2_price_MWh_sort_cal_b=sort(H2_price_MWh_cal,'descend');%F


avg_H2_price_IELGAS=mean(H2_35_cal);
avg_H2_price_model=mean(H2_price_MWh_sort_cal_b);%F
avg_NG_price_IELGAS=mean(NG_35_cal);

SF_price_test=avg_H2_price_IELGAS/avg_H2_price_model;
SF_GWh=round(SF_price_test*1e-3);
SF_MWh=SF_GWh*1e3;
% SF_MWh=SF_GWh*1e3;

SF_kg= 0.039; %To get price per kg hydrogen

% SF_GWh=SF_MWh*1e-3;
%Calibration of storage:
SL_sc=1000*0.6; %Desired storage level of Salt Caverns: 60% (in GWh)
SL_ic=45*0.6; %Desired storage level of local storages: 60% (in MWh)

% %Check cummulative supply:
% cum=base.total_supply*SF;
% 
% % Supply comparison
% figure()
% plot(X,cum)
% hold off
%% 

% 1) Seperate Markets for renewable and low-carbon hydrogen
Green_Blue=sim('Green_Blue_Model',Z); 

% 2) Including ammonia market to include ammonia import
Ammonia=sim('Ammonia_Model',Z); 

% 3) Controlled core network
Core_Control=sim('Core_Model_Controlled',Z); 

delete(findall(0));

%% Calculations Core Calibration Model (Seasonal wind distrution)

%Supply
I_ship_BSS=Core_Cal.ship*SF_GWh; % Regassified LH2 (GWh/h) 
I_pipe_BSS=Core_Cal.pipe*SF_GWh; % Import by pipeline (GWh/h) 
P_OWE_BSS=Core_Cal.green*SF_GWh; % Green hydrogen production (GWh/h) 
P_SMR_BSS=Core_Cal.blue*SF_GWh; % Blue hydrogen production

S_total_BSS=I_ship_BSS+I_pipe_BSS+P_OWE_BSS+P_SMR_BSS; % Total H2 supply (GWh/h)
S_total_BSS_y=Core_Cal.total_supply(8761)*SF_GWh; % Total H2 supply (GWh/year)

%Demand
D_GR_BSS=Core_Cal.D_GR(:)*SF_GWh; % Groningen (GWh/h)
D_AMS_BSS=Core_Cal.D_AMS(:)*SF_GWh; % Amsterdam (GWh/h
D_RDAM_BSS=Core_Cal.D_RDAM(:)*SF_GWh; % Rotterdam (GWh/h)
D_ZL_BSS=Core_Cal.D_ZL(:)*SF_GWh; % Zeeland (GWh/h)
D_LIM_BSS=Core_Cal.D_LIM(:)*SF_GWh; % Limburg (GWh/h)

%Prices
H2_price_MWh_BSS=Core_Cal.H2_price(:)*SF_MWh; % H2 Price (per MWh)
avg_H2_price_MWh_BSS=mean(H2_price_MWh_BSS); % Average H2 Price (per MWh)

H2_price_kg_BSS=H2_price_MWh_BSS*SF_kg; % H2 Price (per kg)

NG_price_MWh_BSS=Core_Cal.NG_price(:)*SF_MWh; % NG price (per MWh)
avg_NG_price_MWh_BSS=mean(NG_price_MWh_BSS); % Average NG price (per MWh)

LH2_price_MWh_BSS=Core_Cal.LH2_price(:)*SF_MWh; % LH2 price (per MWh)

% Stocks
CS_stock_BSS=Core_Cal.CS_stock(:)*SF_GWh; % Salt Caverns stock (GWh)
NG_stock_BSS=Core_Cal.NG_stock(:)*SF_GWh; % NG stock (GWh)
LH2_stock_BSS=Core_Cal.LH2_stock(:)*SF_GWh; % LH2 stock (GWh)

%% Calculaions Core Model (standard=Wind Distribution)

%Supply
I_ship_BWD=Core.ship*SF_GWh; % RegaWDified LH2 (GWh/h) 
I_pipe_BWD=Core.pipe*SF_GWh; % Import by pipeline (GWh/h) 
P_OWE_BWD=Core.green*SF_GWh; % Green hydrogen production (GWh/h) 
P_SMR_BWD=Core.blue*SF_GWh; % Blue hydrogen production

S_total_BWD=I_ship_BWD+I_pipe_BWD+P_OWE_BWD+P_SMR_BWD; % Total H2 supply (GWh/h)
S_total_BWD_y=Core.total_supply(8761)*SF_GWh; % Total H2 supply (GWh/year)

%Demand
D_GR_BWD=Core.D_GR(:)*SF_GWh; % Groningen (GWh/h)
D_AMS_BWD=Core.D_AMS(:)*SF_GWh; % Amsterdam (GWh/h
D_RDAM_BWD=Core.D_RDAM(:)*SF_GWh; % Rotterdam (GWh/h)
D_ZL_BWD=Core.D_ZL(:)*SF_GWh; % Zeeland (GWh/h)
D_LIM_BWD=Core.D_LIM(:)*SF_GWh; % Limburg (GWh/h)

%Prices
H2_price_MWh_BWD=Core.H2_price(:)*SF_MWh; % H2 Price (per MWh)
avg_H2_price_MWh_BWD=mean(H2_price_MWh_BWD); % Average H2 Price (per MWh)

H2_price_kg_BWD=H2_price_MWh_BWD*SF_kg; % H2 Price (per kg)

NG_price_MWh_BWD=Core.NG_price(:)*SF_MWh; % NG price (per MWh)
avg_NG_price_MWh_BWD=mean(NG_price_MWh_BWD); % Average NG price (per MWh)

LH2_price_MWh_BWD=Core.LH2_price(:)*SF_MWh; % LH2 price (per MWh)

% Stocks
CS_stock_BWD=Core.CS_stock(:)*SF_GWh; % Salt Caverns stock (GWh)
NG_stock_BWD=Core.NG_stock(:)*SF_GWh; % NG stock (GWh)
LH2_stock_BWD=Core.LH2_stock(:)*SF_GWh; % LH2 stock (GWh)


%% Calculations Ammonia Model (Standard=Wind Distribution)
%Supply
I_ship_AWD=Ammonia.ship*SF_GWh; % RegaWDified LH2 (GWh/h) 
I_pipe_AWD=Ammonia.pipe*SF_GWh; % Import by pipeline (GWh/h) 
P_OWE_AWD=Ammonia.green*SF_GWh; % Green hydrogen production (GWh/h) 
P_SMR_AWD=Ammonia.blue*SF_GWh; % Blue hydrogen production

I_ammonia_AWD=Ammonia.ammonia*SF_GWh; % Ammonia import (GWh/h)

S_total_AWD=I_ship_AWD+I_pipe_AWD+P_OWE_AWD+P_SMR_AWD; % Total H2 supply (GWh/h)
SA_total_AWD=S_total_AWD+I_ammonia_AWD;% Total H2 + ammonia supply (GWh/h)
S_total_AWD_y=Ammonia.total_supply(8761)*SF_GWh; % Total H2 supply (GWh/year)

%Demand Hydrogen
D_GR_AWD=Ammonia.D_GR(:)*SF_GWh; % Groningen (GWh/h)
D_AMS_AWD=Ammonia.D_AMS(:)*SF_GWh; % Amsterdam (GWh/h
D_RDAM_AWD=Ammonia.D_RDAM(:)*SF_GWh; % Rotterdam (GWh/h)
D_ZL_AWD=Ammonia.D_ZL(:)*SF_GWh; % Zeeland (GWh/h)
D_LIM_AWD=Ammonia.D_LIM(:)*SF_GWh; % Limburg (GWh/h)

%Demand Ammonia
Da_RDAM_AWD=Ammonia.D_RDAM1(:)*SF_GWh; % Rotterdam (GWh/h)
Da_ZL_AWD=Ammonia.D_ZL1(:)*SF_GWh; % Zeeland (GWh/h)
Da_LIM_AWD=Ammonia.D_LIM1(:)*SF_GWh; % Limburg (GWh/h)

%Prices
H2_price_MWh_AWD=Ammonia.H2_price(:)*SF_MWh; % H2 Price (per MWh)
avg_H2_price_MWh_AWD=mean(H2_price_MWh_AWD); % Average H2 Price (per MWh)

H2_price_kg_AWD=H2_price_MWh_AWD*SF_kg; % H2 Price (per kg)

NG_price_MWh_AWD=Ammonia.NG_price(:)*SF_MWh; % NG price (per MWh)
avg_NG_price_MWh_AWD=mean(NG_price_MWh_AWD); % Average NG price (per MWh)

LH2_price_MWh_AWD=Ammonia.LH2_price(:)*SF_MWh; % LH2 price (per MWh)

Ammonia_price_MWh_AWD=Ammonia.ammonia_price(:)*SF_MWh; % Ammonia price (per MWh)
Ammonia_price_kg_AWD=Ammonia_price_MWh_AWD*SF_kg; % Ammonia price (per kg)

% Stocks
CS_stock_AWD=Ammonia.CS_stock(:)*SF_GWh; % Salt Caverns stock (GWh)
NG_stock_AWD=Ammonia.NG_stock(:)*SF_GWh; % NG stock (GWh)
LH2_stock_AWD=Ammonia.LH2_stock(:)*SF_GWh; % LH2 stock (GWh)
Ammonia_stock_AWD=Ammonia.ammonia_stock(:)*SF_GWh; % Ammonia stock (GWh)




%% Calculations Renewable vs Low-Carbon Model (standard=Wind Distribution)
%Supply
I_ship_GBWD=Green_Blue.ship*SF_GWh; % RegGBWDified LH2 (GWh/h) 
I_pipe_GBWD=Green_Blue.pipe*SF_GWh; % Import by pipeline (GWh/h) 
P_OWE_GBWD=Green_Blue.green*SF_GWh; % Green hydrogen production (GWh/h) 
P_SMR_GBWD=Green_Blue.blue*SF_GWh; % Blue hydrogen production

S_total_GBWD=I_ship_GBWD+I_pipe_GBWD+P_OWE_GBWD+P_SMR_GBWD; % Total H2 supply (GWh/h)
S_total_B_GBWD=I_ship_GBWD+I_pipe_GBWD+P_SMR_GBWD; %Total Blue H2 supply (GWh/h)
S_total_GBWD_y=Green_Blue.total_supply(8761)*SF_GWh; % Total H2 supply (GWh/year)

%Demand Blue
Db_GR_GBWD=Green_Blue.D_GR(:)*SF_GWh; % Groningen (GWh/h)
Db_AMS_GBWD=Green_Blue.D_AMS(:)*SF_GWh; % Amsterdam (GWh/h
Db_RDAM_GBWD=Green_Blue.D_RDAM(:)*SF_GWh; % Rotterdam (GWh/h)
Db_ZL_GBWD=Green_Blue.D_ZL(:)*SF_GWh; % Zeeland (GWh/h)
Db_LIM_GBWD=Green_Blue.D_LIM(:)*SF_GWh; % Limburg (GWh/h)

%Demand Green
Dg_GR_GBWD=Green_Blue.D_GR1(:)*SF_GWh; % Groningen (GWh/h)
Dg_AMS_GBWD=Green_Blue.D_AMS1(:)*SF_GWh; % Amsterdam (GWh/h
Dg_RDAM_GBWD=Green_Blue.D_RDAM1(:)*SF_GWh; % Rotterdam (GWh/h)
Dg_ZL_GBWD=Green_Blue.D_ZL1(:)*SF_GWh; % Zeeland (GWh/h)
Dg_LIM_GBWD=Green_Blue.D_LIM1(:)*SF_GWh; % Limburg (GWh/h)

%Prices
B_H2_price_MWh_GBWD=Green_Blue.LC_H2_price(:)*SF_MWh; % Blue H2 Price (per MWh)
avg_B_H2_price_MWh_GBWD=mean(B_H2_price_MWh_GBWD); % Average H2 Price (per MWh)

B_H2_price_kg_GBWD=B_H2_price_MWh_GBWD*SF_kg; % Blue H2 Price (per kg)

G_H2_price_MWh_GBWD=Green_Blue.R_H2_price(:)*SF_MWh; % Green H2 Price (per MWh)
avg_G_H2_price_MWh_GBWD=mean(G_H2_price_MWh_GBWD); % Average H2 Price (per MWh)
G_H2_price_kg_GBWD=B_H2_price_MWh_GBWD*SF_kg; %Green H2 Price (per kg)

NG_price_MWh_GBWD=Green_Blue.NG_price(:)*SF_MWh; % NG price (per MWh)
avg_NG_price_MWh_GBWD=mean(NG_price_MWh_GBWD); % Average NG price (per MWh)

LH2_price_MWh_GBWD=Green_Blue.LH2_price(:)*SF_MWh; % LH2 price (per MWh)

% Stocks
CS_stock_G_GBWD=Green_Blue.CS_stock_G(:)*SF_GWh; % Green
CS_stock_B_GBWD=Green_Blue.CS_stock_B(:)*SF_GWh; % Blue
CS_stock_GBWD=CS_stock_G_GBWD+CS_stock_B_GBWD; % Salt Caverns stock (GWh)
NG_stock_GBWD=Green_Blue.NG_stock(:)*SF_GWh; % NG stock (GWh)
LH2_stock_GBWD=Green_Blue.LH2_stock(:)*SF_GWh; % LH2 stock (GWh)

%% Sorting price curves (If one is interested in price duration curves)

H2_price_MWh_BSS_sort=sort(H2_price_MWh_BSS,'descend');
NG_price_MWh_BSS_sort=sort(NG_price_MWh_BSS,'descend');
LH2_price_MWh_BSS_sort=sort(LH2_price_MWh_BSS,'descend');

H2_price_MWh_BWD_sort=sort(H2_price_MWh_BWD,'descend');
NG_price_MWh_BWD_sort=sort(NG_price_MWh_BWD,'descend');
LH2_price_MWh_BWD_sort=sort(LH2_price_MWh_BWD,'descend');

H2_price_MWh_AWD_sort=sort(H2_price_MWh_AWD,'descend');
NG_price_MWh_AWD_sort=sort(NG_price_MWh_AWD,'descend');
LH2_price_MWh_AWD_sort=sort(LH2_price_MWh_AWD,'descend');
Ammonia_price_MWh_AWD_sort=sort(Ammonia_price_MWh_AWD,'descend');

G_H2_price_MWh_GBWD_sort=sort(G_H2_price_MWh_GBWD,'descend');
B_H2_price_MWh_GBWD_sort=sort(B_H2_price_MWh_GBWD,'descend');
NG_price_MWh_GBWDS_sort=sort(NG_price_MWh_GBWD,'descend');
LH2_price_MWh_GBWD_sort=sort(LH2_price_MWh_GBWD,'descend');
%% Government price control by subsidy 

%Supply
I_ship_BWD_C=Core_Control.ship*SF_GWh; %GWh/h = GW (capaciteit)
I_pipe_BWD_C=Core_Control.pipe*SF_GWh;
P_OWE_BWD_C=Core_Control.green*SF_GWh;
P_SMR_BWD_C=Core_Control.blue*SF_GWh;
H2_Sub_BWD_C=Core_Control.H2_subsidy*SF_GWh;
P_SMRT_BWD_C=P_SMR_BWD_C+H2_Sub_BWD_C;

S_total_BWD_C=I_ship_BWD_C+I_pipe_BWD_C+P_OWE_BWD_C+P_SMR_BWD_C+H2_Sub_BWD_C;
S_total_BWD_C_y=Core_Control.total_supply(8761)*SF_GWh;

%Demand
D_total_BWD_C=Core_Control.Total_demand(:)*SF_GWh;


%H2 Price (MWh)
H2_price_MWh_BWD_C=Core_Control.H2_price(:)*SF_MWh;
avg_H2_price_MWh_BWD_C=mean(H2_price_MWh_BWD_C);

%H2 Price (kg)
H2_price_kg_BWD_C=H2_price_MWh_BWD_C*SF_kg;

% % Price dynamics
% H2_price_dyn_BWD_C=base.H2_price_dyn(:)*SF_GWh*SF_MWh;

% NG price (MWh)
NG_price_MWh_BWD_C=Core_Control.NG_price(:)*SF_MWh;
avg_NG_price_MWh_BWD_C=mean(NG_price_MWh_BWD_C);

% LH2 price (MWh)
LH2_price_MWh_BWD_C=Core_Control.LH2_price(:)*SF_MWh;

%Salt Caverns stock
CS_stock_BWD_C=Core_Control.CS_stock(:)*SF_GWh;
LH2_stock_BWD_C=Core_Control.LH2_stock(:)*SF_GWh;
NG_stock_BWD_C=Core_Control.NG_stock(:)*SF_GWh;

H2_price_MWh_BWD_C_sort=sort(H2_price_MWh_BWD_C,'descend');
NG_price_MWh_BWD_C_sort=sort(NG_price_MWh_BWD_C,'descend');
LH2_price_MWh_BWD_C_sort=sort(LH2_price_MWh_BWD_C,'descend');

% Calculate subsidy costs
SC_h_BWD_C=Core_Control.subsidy_cost_h(:);
CSC_y_BWD_C=Core_Control.subsidy_costs_y(:);
TSC_y_BWD_C=Core_Control.subsidy_costs_y(8761);
display(TSC_y_BWD_C)

% Effect on prices
H2_price_decrease_BWD_C=-(avg_H2_price_MWh_BWD_C-avg_H2_price_MWh_BWD)/avg_H2_price_MWh_BWD_C*100;
NG_price_increase_BWD_C=(avg_NG_price_MWh_BWD_C-avg_NG_price_MWh_BWD)/avg_NG_price_MWh_BWD_C*100;
display(H2_price_decrease_BWD_C)
display(NG_price_increase_BWD_C)

% Effect on supply:
Supply_increase_BWD_C=(S_total_BWD_C_y-S_total_BWD_y)/S_total_BWD_y*100;
display(Supply_increase_BWD_C)

%% Define colors for plotting:

%Hydrogen
color_cyan=1/255*[123,214,214]; % General H2
color_green= 1/255*[147, 192 ,145]; % Green H2
color_blue= 1/255*[108,142,191]; % Blue H2
color_livid= 1/255*[146, 174, 195]; % Regassified LH2
color_dark_gray= 1/255*[170, 170, 170]; % H2 import by pipeline
color_orange=1/255*[244, 166, 87]; % H2 production by government subsidy

%Other
color_purple=1/255*[150, 115, 166]; % LH2
color_brown=1/255*[184, 157, 127]; % NG
color_yellow=1/255*[255, 230,150]; % OWE
color_pink=1/255*[255, 150, 150]; % Ammonia

color_light_gray= 1/255*[95, 95, 95]; % H2+ ammonia



%% Base Scenario plots (standard=wind distribution)
% Supply comparison
SF_time=365/8760;

figure()
plot(X,P_OWE_BWD, 'Color', color_green)
hold on
plot(X,P_SMR_BWD, 'Color', color_blue)
hold on
plot(X,I_ship_BWD,'Color', color_livid)
hold on
plot(X, I_pipe_BWD,'Color', color_dark_gray)
legend('Renewable H2 Production', 'Low-carbon H2 production','LH2 import by ship ', 'H2 import by pipe')
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh/h');
title('Supply Comparison');
 xlim([0 365*years])
 ylim([0 6])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

% Total Supply
figure()
plot(X,S_total_BWD,'Color', color_cyan)
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh/h');
title('Total Supply');
% axis([0 365*years 8 10.5]);
xlim([0 365*years])
ylim([6 14])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})


% Demand comparison
figure()
plot(X,D_GR_BWD,'Color',color_blue)
hold on
plot(X,D_AMS_BWD,'Color',color_purple )
hold on
plot(X,D_RDAM_BWD,'Color',color_orange)
hold on
plot(X,D_ZL_BWD,'Color',color_livid)
hold on
plot(X,D_LIM_BWD,'Color',color_green)
legend on
legend('Groningen', 'Amsterdam', 'Rotterdam', 'Zeeland', 'Limburg')
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh/h');
title('Demand Comparison');
xlim([0 365*years])
ylim([0 4.5])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

% Price H2 per MWh
figure()
plot(X,H2_price_MWh_BWD,'Color', color_light_gray)
grid on
grid minor
xlabel('Time (hours)');
ylabel('€/MWh');
title('Hydrogen price (MWh)');
xlim([0 365*years])
ylim([0 100])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

% Price comparison per MWh
figure()
plot(X,H2_price_MWh_BWD,'Color',color_cyan)
hold on
plot(X,NG_price_MWh_BWD,'Color',color_brown)
hold on
plot(X,LH2_price_MWh_BWD,'Color',color_purple)
legend('Hydrogen','Natural gas','Liquid hydrogen')
grid on
grid minor
xlabel('Time (hours)');
ylabel('€/MWh');
xlabel('hours');
title('Price comparison (MWh)');
xlim([0 365*years])
ylim([0 100])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

% Price H2 kg
figure()
plot(X,H2_price_kg_BWD,'Color', color_light_gray)
grid on
grid minor
xlabel('Time (hours)');
ylabel('€/kg');
title('Price (kg)');
xlim([0 365*years])
ylim([0 9])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

% Salt Caverns Hydrogen stock
figure()
plot(X,CS_stock_BWD,'Color', color_light_gray)
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh');
title('Hydrogen stock');
xlim([0 365*years])
ylim([-1.6 0.6])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off


% Stock Comparison
figure()
plot(X,CS_stock_BWD,'Color', color_cyan)
hold on
plot(X,NG_stock_BWD,'Color',color_brown)
hold on
plot(X,LH2_stock_BWD,'Color',color_purple)
legend on
legend ('Hydrogen','Natural Gas','Liquid Hydrogen')
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh');
title('Stock comparison');
xlim([0 365*years])
ylim([-3 3])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off


%Price Duration Curves MWh
figure()
plot(X*24,H2_price_MWh_BWD_sort,'Color',color_cyan, 'LineWidth',2)
hold on
plot(X*24,NG_price_MWh_BWD_sort,'Color',color_brown, 'LineWidth',2)
hold on
plot(X*24,LH2_price_MWh_BWD_sort,'Color',color_purple, 'LineWidth',2)
legend('Hydrogen','Natural Gas','Liquid hydrogen')
grid on
grid minor
xlabel('Time (hours)');
ylabel('€/MWh');
axis([0 365*24*years 0 100]);
title('Price Duration Curves (MWh)');
hold off




%% Plots Ammonia Model (standard=wind distribution)
% Supply comparison
figure()
plot(X,P_OWE_AWD, 'Color', color_green)
hold on
plot(X,P_SMR_AWD, 'Color', color_blue)
hold on
plot(X,I_ship_AWD,'Color', color_livid)
hold on
plot(X, I_pipe_AWD,'Color', color_light_gray)
hold on
plot(X, I_ammonia_AWD,'Color', color_pink)
legend on
legend('Renewable H2 Production', 'Low-carbon H2 production','LH2 import by ship ', 'H2 import by pipe','Ammonia')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh/h');
title('Supply Comparison');
xlim([0 365*years])
ylim([0 6])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

% Total Supply
figure()
plot(X,S_total_AWD,'Color', color_dark_gray)
hold on
plot(X,SA_total_AWD,'Color', color_dark_gray)
legend on
legend('Hydrogen only', 'Hydrogen and ammonia together')
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh/h');
title('Total Supply');
axis([0 365*years 8 10.5]);
xlim([0 365*years])
ylim([6 14])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})


% Total Supply
figure()
plot(X,S_total_AWD,'Color', color_cyan)
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh/h');
title('Total Supply');
axis([0 365*years 8 10.5]);
xlim([0 365*years])
ylim([6 14])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})


% Demand comparison
figure()
plot(X,D_GR_AWD,'Color',color_blue)
hold on
plot(X,D_AMS_AWD,'Color',color_purple )
hold on
plot(X,D_RDAM_AWD,'Color',color_orange)
hold on
plot(X,D_ZL_AWD,'Color',color_livid)
hold on
plot(X,D_LIM_AWD,'Color',color_green)
hold on
plot(X,Da_RDAM_AWD,':','Color',color_orange)
hold on
plot(X,Da_ZL_AWD,':','Color',color_livid)
hold on
plot(X,Da_LIM_AWD,':','Color',color_green)
legend on
legend('Groningen', 'Amsterdam', 'Rotterdam', 'Zeeland', 'Limburg')
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh/h');
title('Demand Comparison');
xlim([0 365*years])
ylim([0 4])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

% Price H2 per MWh
figure()
plot(X,H2_price_MWh_AWD,'Color', color_light_gray)
grid on
grid minor
xlabel('Time (hours)');
ylabel('€/MWh');
title('Hydrogen price (MWh)');
xlim([0 365*years])
ylim([0 75])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

%Price comparison MWh
figure()
plot(X,H2_price_MWh_AWD,'Color',color_cyan)
hold on
plot(X,Ammonia_price_MWh_AWD, 'Color', color_pink)
hold on
plot(X,NG_price_MWh_AWD,'Color',color_brown)
hold on
plot(X,LH2_price_MWh_AWD,'Color',color_livid)
legend('H2','Ammonia','NG','LH2')
grid on
grid minor
xlabel('hours');
ylabel('€/MWh');
title('Price comparison H2/Ammonia (MWh)');
xlim([0 365*years])
ylim([0 100])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off


% Stock Comparison
figure(9)
plot(X,CS_stock_AWD,'Color', color_cyan)
hold on
plot(X,NG_stock_AWD,'Color',color_brown)
hold on
plot(X,LH2_stock_AWD,'Color',color_purple)
hold on
plot(X,Ammonia_stock_AWD,'Color',color_pink)
legend('H2','Ammonia','NG','LH2')
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh');
title('Stock comparison');
xlim([0 365*years])
ylim([-1.2 0.6])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off


% Price duration curves
figure()
plot(X,H2_price_MWh_AWD_sort,'Color',color_cyan, 'LineWidth',2)
hold on
plot(X,Ammonia_price_MWh_AWD_sort, 'Color', color_pink, 'LineWidth',2)
hold on
plot(X,NG_price_MWh_AWD_sort,'Color',color_brown, 'LineWidth',2)
hold on
plot(X,LH2_price_MWh_AWD_sort,'Color',color_livid, 'LineWidth',2)
legend('H2','Ammonia','NG','LH2')
grid on
grid minor
xlabel('Time (day)');
ylabel('€/MWh');
title('Price Duration comparison H2/Ammonia (MWh)');
xlim([0 365*years])
ylim([0 100])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off




%% Renewable vs Low-Carbon model (standard=wind_distribution)

% Supply comparison
figure()
plot(X,P_OWE_GBWD, 'Color', color_green)
hold on
plot(X,P_SMR_GBWD, 'Color', color_blue)
hold on
plot(X,I_ship_GBWD,'Color', color_livid)
hold on
plot(X, I_pipe_GBWD,'Color', color_light_gray)
legend on
legend('Renewable H2 Production', 'Low-carbon H2 production','LH2 import by ship ', 'H2 import by pipe')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh/h');
title('Supply Comparison');
xlim([0 365*years])
ylim([0 6])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
% xlim([72 79])
% ylim([0 6])
% xlim([82 89])
% ylim([0 6])
% xlim([172 179])
% ylim([0 6])
% xlim([182 189])
% ylim([0 6])
hold off


% Total Supply
figure()
plot(X,S_total_GBWD,'Color', color_cyan)
hold on
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh/h');
title('Total Supply');
xlim([0 365*years])
ylim([6 14])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})



% Demand comparison
figure()
plot(X,Db_GR_GBWD,'Color',color_blue)
hold on
plot(X,Db_AMS_GBWD,'Color',color_purple )
hold on
plot(X,Db_RDAM_GBWD,'Color',color_orange)
hold on
plot(X,Db_ZL_GBWD,'Color',color_livid)
hold on
plot(X,Db_LIM_GBWD,'Color',color_green)
hold on
plot(X,Dg_GR_GBWD,':','Color',color_blue)
hold on
plot(X,Dg_AMS_GBWD,':','Color',color_purple )
hold on
plot(X,Dg_RDAM_GBWD,':','Color',color_orange)
hold on
plot(X,Dg_ZL_GBWD,':','Color',color_livid)
hold on
plot(X,Dg_LIM_GBWD,':','Color',color_green)
legend on
legend('Groningen', 'Amsterdam', 'Rotterdam', 'Zeeland', 'Limburg')
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh/h');
title('Demand Comparison');
xlim([0 365*years])
ylim([0 3.5])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
hold off

%Price comparison MWh
figure()
plot(X,B_H2_price_MWh_GBWD,'Color',color_blue)
hold on
plot(X,G_H2_price_MWh_GBWD, 'Color', color_green)
hold on
plot(X,NG_price_MWh_GBWD, 'Color', color_brown)
hold on
plot(X,LH2_price_MWh_GBWD, 'Color', color_purple)
legend('Blue H2','Green H2','Natural Gas','Liquid Hydrogen')
grid on
grid minor
xlabel('hours');
ylabel('€/MWh');
title('Price comparison H2 (MWh)');
xlim([0 365*years])
ylim([0 100])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
% xlim([72 79])
% ylim([0 90])
% ylim([0 6])
% xlim([82 89])
% ylim([0 6])
% xlim([172 179])
% ylim([0 6])
% xlim([182 189])
% ylim([0 6])
hold off
% 
%Price comparison kg
figure()
plot(X,B_H2_price_kg_GBWD,'Color',color_blue)
hold on
plot(X,G_H2_price_kg_GBWD,'Color', color_green)
legend('Blue H2','Green H2')
grid on
grid minor
xlabel('Time (day)');
ylabel('€/MWh');
title('Price comparison (kg)');
hold off


% Stock Comparison
figure()
plot(X,CS_stock_B_GBWD,'Color', color_blue)
hold on
plot(X,CS_stock_G_GBWD,'Color', color_green)
hold on
plot(X,NG_stock_GBWD,'Color',color_brown)
hold on
plot(X,LH2_stock_GBWD,'Color',color_purple)
grid on
grid minor
xlabel('Time (hours)');
ylabel('GWh');
title('Stock comparison');
xlim([0 365*years])
ylim([-1.2 0.6])
xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
% xlim([72 79])
% ylim([-0.2 0.6])
% ylim([0 6])
% xlim([82 89])
% ylim([0 6])
% xlim([172 179])
% ylim([0 6])
% xlim([182 189])
% ylim([0 6])

%Price Duration Curves MWh
figure()
% plot(X*24,H2_price_MWh_sort,'Color',color_cyan)
% hold on
plot(X*24,B_H2_price_MWh_GBWD_sort,'Color',color_blue,'Linewidth',2)
hold on
plot(X*24,G_H2_price_MWh_GBWD_sort, 'Color', color_green,'Linewidth',2)
legend('Blue H2','Green H2')
grid on
grid minor
xlabel('Time (hours)');
ylabel('€/MWh');
axis([0 365*24*years 0 100]);
title('Price Duration Curves (MWh)');
hold off


%% Three network comparisons
% One can choose to simulate the

% The stormy spring period: xlim([82 87])
% The low-wind summer period: xlim([192 199])
% The volatile winter period: xlim([342 349])



%H2 price comparison
figure()
plot(X,H2_price_MWh_BWD,'Color',color_cyan)
hold on
plot(X,G_H2_price_MWh_GBWD,'Color',color_green)
hold on
plot(X,B_H2_price_MWh_GBWD,'Color',color_blue)
hold on
plot(X,H2_price_MWh_AWD, 'Color', color_pink)
legend('Base','Renewable','Low-carbon','Ammonia')
grid on
grid minor
xlabel('Time (day)');
ylabel('€/MWh');
title('Price comparison H2 (MWh) (Wind distribution)');
axis([0 365*years 0 100]);
hold off



%Comparison Price Duration Curves MWh
figure()
plot(X*24,H2_price_MWh_BWD_sort,'Color',color_cyan, 'LineWidth',2)
hold on
plot(X*24,G_H2_price_MWh_GBWD_sort,'Color',color_green, 'LineWidth',2)
hold on
plot(X*24,B_H2_price_MWh_GBWD_sort,'Color',color_blue, 'LineWidth',2)
hold on
plot(X*24,H2_price_MWh_AWD_sort,'Color',color_pink, 'LineWidth',2)
legend('Base','Renewable','Low-carbon','Ammonia')
grid on
grid minor
xlabel('Time (hours)');
ylabel('€/MWh');
axis([0 365*24*years 0 100]);
title('Price Duration Curves (MWh) (Wind distribution)');
hold off


% Total hydrogen Supply
figure()
plot(X,S_total_BWD,'Color', color_cyan)
hold on
plot(X,S_total_GBWD,'Color', color_livid)
hold on
plot(X,S_total_AWD,'Color', color_pink)
legend('Base','Renewable/Low-carbon','Ammonia')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh/h');
title('Total hydrogen Supply');
axis([0 365*years 5 13]);

% Supply comparison
figure()
plot(X,P_OWE_BWD, 'Color', color_green, 'LineWidth',1.2)
hold on
plot(X,P_SMR_BWD, 'Color', color_blue, 'LineWidth',1.2)
hold on
plot(X,I_ship_BWD,'Color', color_livid, 'LineWidth',1.2)
hold on
plot(X, I_pipe_BWD,'Color', color_light_gray, 'LineWidth',1.2)
hold on
plot(X,P_SMR_GBWD,':', 'Color', color_blue, 'LineWidth',1.2)
hold on
plot(X,I_ship_GBWD,':','Color', color_livid, 'LineWidth',1.2)
hold on
plot(X, I_pipe_GBWD,':','Color', color_light_gray, 'LineWidth',1.2)
hold on
plot(X,P_SMR_AWD,'--', 'Color', color_blue, 'LineWidth',1.2)
hold on
plot(X,I_ship_AWD,'--','Color', color_livid, 'LineWidth',1.2)
hold on
plot(X, I_pipe_AWD,'--','Color', color_light_gray, 'LineWidth',1.2)
hold on
plot(X, I_ammonia_AWD,'--','Color', color_pink, 'LineWidth',1.2)
legend on
legend('Renewable H2 Production', 'Low-carbon H2 production','LH2 import by ship ', 'H2 import by pipe')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh/h');
title('Supply Comparison');
% xlim([0 365*years])
% ylim([0 6])
% xticks([0 1000*SF_time 2000*SF_time 3000*SF_time 4000*SF_time 5000*SF_time 6000*SF_time 7000*SF_time 8000*SF_time])
% xticklabels({'0','1000','2000','3000','4000','5000', '6000','7000', '8000'})
% xlim([82 87])
ylim([-0.5 6])
% xlim([192 199])
xlim([342 349])
hold off

% Total hydrogen Supply
figure()
plot(X,S_total_BWD,'Color', color_cyan, 'LineWidth',1.2)
hold on
plot(X,S_total_GBWD,':','Color', color_cyan, 'LineWidth',1.2)
hold on
plot(X,SA_total_AWD,'--','Color', color_dark_gray, 'LineWidth',1.2)
legend('Base','Renewable/Low-carbon','Ammonia')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh/h');
title('Total hydrogen Supply');
xlim([82 87])
ylim([5 13])
% xlim([192 199])
xlim([342 349])

%Stock comparison
figure()
plot(X,CS_stock_BWD,'Color', color_cyan, 'LineWidth',1.2)
hold on
plot(X,NG_stock_BWD,'Color',color_brown, 'LineWidth',1.2)
hold on
plot(X,LH2_stock_BWD,'Color',color_purple, 'LineWidth',1.2)
hold on
plot(X,CS_stock_B_GBWD,':','Color', color_blue, 'LineWidth',1.2)
hold on
plot(X,CS_stock_G_GBWD,':','Color', color_green, 'LineWidth',1.2)
hold on
plot(X,NG_stock_GBWD,':','Color',color_brown, 'LineWidth',1.2)
hold on
plot(X,LH2_stock_GBWD,':','Color',color_purple, 'LineWidth',1.2)
hold on
plot(X,CS_stock_AWD,'--','Color', color_cyan, 'LineWidth',1.2)
hold on
plot(X,NG_stock_AWD,'--','Color',color_brown, 'LineWidth',1.2)
hold on
plot(X,LH2_stock_AWD,'--','Color',color_purple, 'LineWidth',1.2)
hold on
plot(X,Ammonia_stock_AWD,'--','Color',color_pink, 'LineWidth',1.2)
legend('Base','Renewable','Low-carbon','Ammonia')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh');
title('Price comparison H2 (MWh) (Wind distribution)');
% xlim([82 87])
ylim([-0.6 0.6])
% xlim([192 199])
xlim([342 349])
hold off

%Price comparison
figure()
plot(X,H2_price_MWh_BWD,'Color',color_cyan, 'LineWidth',1.2)
hold on
plot(X,NG_price_MWh_BWD, 'Color', color_brown, 'LineWidth',1.2)
hold on
plot(X,LH2_price_MWh_BWD, 'Color', color_purple, 'LineWidth',1.2)
hold on
plot(X,G_H2_price_MWh_GBWD,':','Color',color_green, 'LineWidth',1.2)
hold on
plot(X,B_H2_price_MWh_GBWD,':','Color',color_blue, 'LineWidth',1.2)
hold on
plot(X,NG_price_MWh_GBWD,':', 'Color', color_brown, 'LineWidth',1.2)
hold on
plot(X,LH2_price_MWh_GBWD,':', 'Color', color_purple, 'LineWidth',1.2)
hold on
plot(X,H2_price_MWh_AWD,'--', 'Color', color_cyan, 'LineWidth',1.2)
hold on
plot(X,NG_price_MWh_AWD,'--', 'Color', color_brown, 'LineWidth',1.2)
hold on
plot(X,LH2_price_MWh_AWD,'--', 'Color', color_purple, 'LineWidth',1.2)
hold on
plot(X,Ammonia_price_MWh_AWD,'--', 'Color', color_pink, 'LineWidth',1.2)
legend('Base','Renewable','Low-carbon','Ammonia')
grid on
grid minor
xlabel('hours');
ylabel('€/MWh');
title('Price comparison H2 (MWh) (Wind distribution)');
xlim([82 87])
ylim([-20 100])
% xlim([192 199])
% xlim([342 349])
hold off


%% Government price control by subsidy 
% Plots comparison between uncontrolled and controlled core network

% Supply comparison
figure()
plot(X,P_OWE_BWD, 'Color', color_green, 'LineWidth',1.2)
hold on
plot(X,P_SMR_BWD, 'Color', color_blue, 'LineWidth',1.2)
hold on
plot(X,I_ship_BWD,'Color', color_livid, 'LineWidth',1.2)
hold on
plot(X, I_pipe_BWD,'Color', color_light_gray, 'LineWidth',1.2)
hold on
legend on
plot(X,P_OWE_BWD_C, ':','Color', color_green, 'LineWidth',1.2)
hold on
plot(X,P_SMRT_BWD_C, ':','Color', color_blue, 'LineWidth',1.2)
hold on
plot(X,I_ship_BWD_C,':','Color', color_livid, 'LineWidth',1.2)
hold on
plot(X, I_pipe_BWD_C,':','Color', color_light_gray, 'LineWidth',1.2)
legend('Renewable H2 Production', 'Low-carbon H2 production','LH2 import by ship ', 'H2 import by pipe')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh/h');
title('Supply Comparison');
% axis([0 365*years 0 4]);
xlim([82 87])
ylim([-0.5 6])

hold off
% Total Supply
figure()
plot(X,S_total_BWD,'Color', color_light_gray, 'LineWidth',1.2)
hold on
plot(X,S_total_BWD_C,':','Color', color_dark_gray, 'LineWidth',1.2)
legend ('Total Supply', 'Total Supply (C)')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh/h');
title('Total Supply');
% axis([0 365*years 8.7 10.2]);
xlim([82 87])
ylim([6.5 12.5])



% Price NG per MWh
figure()
plot(X,H2_price_MWh_BWD,'Color',color_cyan, 'LineWidth',1.2)
hold on
plot(X,NG_price_MWh_BWD,'Color',color_brown, 'LineWidth',1.2)
hold on
plot(X,LH2_price_MWh_BWD,'Color',color_purple, 'LineWidth',1.2)
hold on
plot(X,H2_price_MWh_BWD_C,":",'Color',color_cyan, 'LineWidth',1.2)
hold on
plot(X,NG_price_MWh_BWD_C,":",'Color',color_brown, 'LineWidth',1.2)
hold on
plot(X,LH2_price_MWh_BWD_C,":",'Color',color_purple, 'LineWidth',1.2)
legend('Hydrogen','Natural gas','Liquid hydrogen','Hydrogen (C)','Natural gas (C)','Liquid hydrogen (C)')
grid on
grid minor
xlabel('Time (day)');
ylabel('€/MWh');
title('Price comparison (MWh)');
xlim([82 87])
ylim([0 100])
hold off


% Salt Caverns Hydrogen stock
figure()
plot(X,CS_stock_BWD,'Color', color_light_gray, 'LineWidth',1.2)
hold on
plot(X,CS_stock_BWD_C,":",'Color',  color_light_gray, 'LineWidth',1.2)
legend ('Hydrogen stock', 'Hydrogen stock (C)')
grid on
grid minor
xlabel('Time (day)');
ylabel('GWh');
title('Hydrogen stock');
% axis([0 365*years -5 5]);
xlim([82 87])
ylim([-0.4 0.4])
hold off


%Price Duration Curves MWh
figure()
plot(X*24,H2_price_MWh_BWD_sort,'Color',color_cyan, 'LineWidth',2)
hold on
plot(X*24,NG_price_MWh_BWD_sort,'Color',color_brown, 'LineWidth',2)
hold on
plot(X*24,LH2_price_MWh_BWD_sort,'Color',color_livid, 'LineWidth',2)
hold on
plot(X*24,H2_price_MWh_BWD_C_sort,":",'Color',color_cyan, 'LineWidth',2)
hold on
plot(X*24,NG_price_MWh_BWD_C_sort,":",'Color',color_brown, 'LineWidth',2)
hold on
plot(X*24,LH2_price_MWh_BWD_C_sort,":",'Color',color_livid, 'LineWidth',2)
legend('Hydrogen','Natural Gas','Liquid hydrogen', 'Hydrogen (C)','Natural Gas (C)','Liquid hydrogen (C)')
grid on
grid minor
xlabel('Time (hours)');
ylabel('€/MWh');
axis([0 365*24*years 0 100]);
title('Price Duration Curves (MWh)');
hold off


