%% Model Parameters
% delete(findall(0));
%% Time vector
% years=1;   % Specify how many years you want to simulate
% n=1/24;
% t=0:years*365-n:n; 
% 

%% Base Model
% First the base model paraters are listed, from which many are also used
% in the other two market designs. 
%% LH2 Market
% LH2 import by ships
I_LH2=0.275;

% Demand LH2
id_lh2=0.275;
L_lh2= 0.3827*0.5;

% Demand LH2 for H2 productionI j
id_lh2h=0.009;
L_lh2h=0.8862*0.5;

% Supply in H2 market
L_sih=0.9293*0.5;

% Price effects through mutual inductance
c_lh212=0.965;
c_lh213= -0.3062;
c_lh223=-0.2228;

% LH2 storage
C_lh=300*0.5;

% Brokerage LH2
R_lh=0.025; 
%% Gasueous hydrogen import by Pipeline
L_i=0.05*0.5; % elasticity

%% Offshore wind Modeling
% See function in simulink

%% Offshore wind energy market

% Demand electricity (L2)
L_e=0.5*0.5;
id_e=0.1*0.5;

% Demand electricity for H2 production (L1)SF_flow
L_eh=0.3*0.5;
id_eh=0.16;

% Price effects through mutual inductance
cc_e=-0.85;

% Electrolyzer efficiency 
c_eh=0.85;

%% Natural Gas market

% Supply NG
L_sng=0.009*0.5; %variable
v_sng=0.4; %inelastic supply

% Storage NG
C_ng=500*0.5;

% Brokerage NG
R_ng=0.025;
% R_ng=0.025*0.3; %underdemped

% Demand NG (L3)
L_ng=L_lh2;
id_ng_IO=0.5; %sinusoidal offset IO
id_ng_IA=0.08; %amplitude
id_ng_FREQ= 1/365; %yearly distribution, period = 1 year
id_ng_TD= -0.25*365; %to make sure seasonality holds

% Demand NG for H2 production (L2)heb j
L_ngh=0.5862*0.5;
id_ngh=0.02;

% Supply low-carbon H2 in H2 market (L1)
L_sbh=0.56*0.5;

% Price effects through mutual inductance
c_ng12=c_lh212;
c_ng13=c_lh213;
c_ng23=c_lh223;
%% H2 Brokerage
R_H2=0.05;
% R_H2=0.05*0.1; %ovverdemped
%% Salt Caverns storage
C_sc=750*0.5;
%% Industry Clusters as demanders

% Groningen
L_g=0.1;
% 2e-1*0.5;
id_g=0.8*0.072*1.2;
C_g=1;

% Amsterdam
L_a=0.25;
% 9e-1*0.5;
id_a=0.8*0.028*1.2;
C_a=1;

% Rotterdam
L_r=0.05;
% 0.85e-2*0.5;
id_r=0.4*0.428*1.2;
C_r=1;

% Zeeland
L_z=0.08;
% 0.95e-1*0.5;
id_z=0.4*0.33*1.2;
C_z=1;

% Limburg
L_l=0.2;
% 8e-1*0.5;
id_l=0.4*0.139*1.2;
C_l=1;

%% Design 1: Seperate markets for renewable and low-carbon hydrogen

% 
S_g_g= 0.29; %share green hydrogen
S_g_a=0.3;
S_g_r=0.19;
S_g_z=0.19;
S_g_l=0.3;


% id_frac=1;
id_frac=0.69;

id_g_BG=id_g*id_frac;
L_LCH2_g= 0.9*0.5;
L_RH2_g= 0.4*0.5;
cc_H2_g=0.85;

id_a_BG=id_a*id_frac;
L_LCH2_a= 0.8*0.5;
L_RH2_a= 0.6*0.5;
cc_H2_a=0.6;

id_r_BG=id_r*id_frac;
L_LCH2_r= 0.5*0.5;
L_RH2_r= 0.2*0.5;
cc_H2_r=0.9;

id_z_BG=id_z*id_frac;
L_LCH2_z= 0.5*0.5;
L_RH2_z= 0.2*0.5;
cc_H2_z=0.85;

id_l_BG=id_l*id_frac;
L_LCH2_l= 0.95*0.5;
L_RH2_l= 0.4*0.5;
cc_H2_l=0.75;



%% Design 2: Ammonia
% H2+ammonia in line with the rest
id_frac_A=0.72;
id_r_A=id_r*id_frac_A;
id_z_A=id_z*id_frac_A;
id_l_A=id_l*id_frac_A;
I_ammonia=0.063; %Inelastic supply
L_ammonia=0.2; %Elastic supply
R_ammonia=0.025;
% R_ammonia=0.025*0.3;%underdemped
C_ammonia=500*0.5;
% S_a=0.4; %share ammonia

S_a_r= 0.2; %0.3;
L_H2_r= 0.45;
L_A_r= 0.35;
cc_AH2_r=0.75;

S_a_z=0.25; %0.35;
L_H2_z= 0.45;
L_A_z= 0.35;
cc_AH2_z=0.65;

S_a_l=0.3; %0.45;
L_H2_l= 0.45;
L_A_l= 0.3;
cc_AH2_l=0.6;


