%% Price Dynamics Analytically

%% Get model parameters
%% Stop Time

years=1; % years of simulation
Z=years*365; % in days
Time=[0:1/24:Z]';
run("Model_Parameters.m");


%% Hydrogen



%% LH2 Mutual Inductor
% Calculate the determinant term in the denominator
denominator_LH2 = 1 / (1 - c_lh212^2 - c_lh213^2 - c_lh223^2 + 2 * c_lh212 * c_lh213 * c_lh223);

% Calculate individual elements of the inverse matrix
M_LH2_inv = denominator_LH2 * [
    (1 - c_lh223^2) / L_sih,            (c_lh223*c_lh213 - c_lh212) / sqrt(L_sih*L_lh2h),   (c_lh212*c_lh223 - c_lh213) / sqrt(L_sih*L_lh2);
    (c_lh223*c_lh213 - c_lh212) / sqrt(L_sih*L_lh2h),   (1 - c_lh213^2) / L_lh2h,             (c_lh212*c_lh213 - c_lh223) / sqrt(L_lh2h*L_lh2);
    (c_lh212*c_lh223 - c_lh213) / sqrt(L_sih*L_lh2),    (c_lh212*c_lh213 - c_lh223) / sqrt(L_lh2h*L_lh2), (1 - c_lh212^2) / L_lh2
];

% Display the symbolic inverse matrix
disp(M_LH2_inv);

%% OWE Mutual Inductor
denominator_OWE=1/(L_eh * L_e*(1-cc_e^2));
M_OWE_inv = denominator_OWE * [L_e, -cc_e * sqrt(L_eh * L_e);
                                -cc_e * sqrt(L_eh * L_e), L_eh];
% M_OWE_inv = inv(M_OWE);

%% NG Mutual Inductor
% Calculate the determinant term in the denominator
denominator_NG = 1 / (1 - c_ng12^2 - c_ng13^2 - c_ng23^2 + 2 * c_ng12 * c_ng13 * c_ng23);

% Calculate individual elements of the inverse matrix
M_NG_inv = denominator_NG * [
    (1 - c_ng23^2) / L_sbh,            (c_ng23*c_ng13 - c_ng12) / sqrt(L_sbh*L_ngh),   (c_ng12*c_ng23 - c_ng13) / sqrt(L_sbh*L_ng);
    (c_ng23*c_ng13 - c_ng12) / sqrt(L_sbh*L_ngh), (1 - c_ng13^2) / L_ngh,             (c_ng12*c_ng13 - c_ng23) / sqrt(L_ngh*L_ng);
    (c_ng12*c_ng23 - c_ng13) / sqrt(L_sbh*L_ng), (c_ng12*c_ng13 - c_ng23) / sqrt(L_ngh*L_ng), (1 - c_ng12^2) / L_ng
];

% Display the symbolic inverse matrix
disp(M_NG_inv);
% % 
% % %% Excess flows
% % 
% % f_H2_e=
% % f_OWE_e=
% % f_LH2_e=
% % f_NG_e=

%% Variables

%LH2
eps_11_LH2= M_LH2_inv(1,1) ;
eps_22_LH2 =M_LH2_inv(2,2) ;
eps_23_LH2 = M_LH2_inv(2,3) ;
eps_33_LH2 = M_LH2_inv(3,3) ;
eps_12_LH2 = M_LH2_inv(1,2) ; 
eps_13_LH2 = M_LH2_inv(1,3) ;
% eps_13_LH2 = 0 ;

w_LH2=I_LH2;
u2_LH2=id_lh2h;
u3_LH2=id_lh2;
k_LH2 =1/C_lh ;
b_LH2 = R_lh;

%OWE
eps_11_OWE = M_OWE_inv(1,1);
eps_12_OWE = M_OWE_inv(1,2);
eps_22_OWE = M_OWE_inv(2,2);

%Choose wind input
u1_OWE=id_eh;
u2_OWE=id_e;
c_OWE = c_eh;
% SW=Seasonal_wind_distribution(); 
% w_OWE=OWE_Generation(); %Seasonal wind distribution
% wd_OWE=diff(w_OWE)/(1/24); %Seasonal wind distribution
w_OWE=OWE_Generation2(); %Measured wind distribution
wd_OWE=diff(w_OWE)/(1/24); %Measured wind distribution


%NG
eps_11_NG = M_LH2_inv(1,1);
eps_22_NG = M_LH2_inv(2,2);
eps_23_NG = M_LH2_inv(2,3);
eps_33_NG = M_LH2_inv(3,3);
eps_12_NG = M_LH2_inv(1,2);
eps_13_NG = M_LH2_inv(1,3);
eps_s_NG = 1/L_sng;

w_NG=v_sng;
u2_NG=id_ngh;
u3_NG= id_ng_IO+id_ng_IA*sin(2*pi*id_ng_FREQ*(Time-id_ng_TD));
k_NG = 1/C_ng ;
b_NG =R_ng;
g_NG=0; %controlled insatiable want



%H2
eps_s_2_H2 = 1/L_i;
eps_d_H2 =1/L_g+1/L_a+1/L_r+1/L_z+1/L_l;

u_H2= id_g+id_a+id_r+id_z+id_l;
b_H2 =R_H2;
k_H2 = 1/C_sc ;
g_H2=0; %controlled insatiable want

% syms eps_22_LH2 eps_23_LH2 eps_33_LH2 k_LH2 b_LH2 eps_12_LH2 eps_13_LH2
% syms eps_22_NG eps_23_NG eps_33_NG eps_s_NG k_NG b_NG eps_12_NG eps_13_NG
% syms eps_11_OWE eps_12_OWE eps_22_OWE c_OWE
% syms eps_d_H2 eps_11_LH2 eps_s_2_H2 eps_11_NG
% syms b_H2 k_H2

%% States and inputs
% 
% x=[q_LH2; p_LH2; p_OWE ; q_NG ; p_NG ; qH2; p_H2 ]; %states

u= [w_LH2; u2_LH2; u3_LH2; 0; 0 ; u1_OWE; u2_OWE; w_NG; u2_NG; 0; u_H2]; %inputs


% Assuming the constants are scalar values
% constant_values = [w_LH2; u2_LH2; u3_LH2; w_OWE(1); diff(w_OWE); u1_OWE; u2_OWE; w_NG; u2_NG; u3_NG; u_H2];

% Create u_matrix with constant values for all inputs
u_matrix = repmat(u, 1, length(w_OWE));

% Replace the column corresponding to w_OWE with its time-varying values
u_matrix(4, :) = w_OWE';
u_matrix(5, :) = [0 wd_OWE'];
u_matrix(10, :)= u3_NG';






%% Matrices

A = [0, -(eps_22_LH2 + 2*eps_23_LH2 + eps_33_LH2), 0, 0, 0, 0, -(eps_12_LH2 + eps_13_LH2);
     k_LH2, -b_LH2*(eps_22_LH2 + 2*eps_23_LH2 + eps_33_LH2), 0, 0, 0, 0, -b_LH2*(eps_12_LH2 + eps_13_LH2);
     0,0,0,0,0,0,0;
      0, 0, 0, 0, -(eps_22_NG + 2*eps_23_NG + eps_33_NG + eps_s_NG), 0, -(eps_12_NG + eps_13_NG);
      0, 0, 0, k_NG, -b_NG*(eps_22_NG + 2*eps_23_NG + eps_33_NG + eps_s_NG), 0, -b_NG*(eps_12_NG + eps_13_NG);
     0, -(eps_12_LH2 + eps_13_LH2), c_OWE*(eps_11_OWE + eps_12_OWE), 0, -(eps_12_NG + eps_13_NG), 0, -(eps_d_H2 + eps_11_LH2 + eps_s_2_H2 + eps_11_NG);
     0, -b_H2*(eps_12_LH2 + eps_13_LH2), b_H2*c_OWE*(eps_11_OWE + eps_12_OWE), 0, -b_H2*(eps_12_NG + eps_13_NG), k_H2, -b_H2*(eps_d_H2 + eps_11_LH2 + eps_s_2_H2 + eps_11_NG)];

B = [
    -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    -b_LH2, b_LH2, b_LH2, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, -1 / (eps_11_OWE + 2 * eps_12_OWE + eps_22_OWE), 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 0;
    0, 0, 0, 0, 0, 0, 0, -b_NG, b_NG, b_NG, 0;
    0, 0, 0, 0, 0, -c_OWE, 0, 0, 0, 0, 1,;
    0, 0, 0, 0, 0, -b_H2 * c_OWE, 0, 0, 0, 0, b_H2
];

C = [
    1, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 1;
];

D = zeros(7, 11);

uncontrollable_states=length(A)-rank(ctrb(A,B));

states = {'q_LH2' ,'p_LH2' , 'p_OWE', 'q_NG', 'p_NG' ,'qH2', 'p_H2' };
inputs = {'w_LH2' ,'u2_LH2' ,'u3_LH2', 'w_OWE' ,'wd_OWE' ,'u1_OWE' ,'u2_OWE', 'w_NG', 'u2_NG', 'u3_NG' ,'u_H2 '}; %without controller
% inputs = {'w_LH2' ,'u2_LH2' ,'u3_LH2', 'w_OWE' ,'wd_OWE' ,'u1_OWE' ,'u2_OWE', 'w_NG', 'u2_NG', 'u3_NG' ,'g_NG' ,'u_H2 ','g_H2 '}; %with controller
outputs = {'q_LH2' ,'p_LH2' , 'p_OWE', 'q_NG', 'p_NG', 'qH2' ,'p_H2'};

% sys_mimo = ss(A,B,C,D,'statename',states,...
% 'inputname',inputs,...
% 'outputname',outputs);

sys_MIMO = ss(A,B,C,D,'StateName',states,'InputName',inputs,'OutputName',outputs);

tf_MIMO=tf(sys_MIMO);
pole(tf_MIMO)



figure()
pzmap(tf_MIMO)
grid on
ylim([-0.4 0.4])
xlim([-10 0])
hold off


%% Plotting

% Define initial conditions
initial_conditions = zeros(7, 1); 



% Simulate the system
[y, t, x] = lsim(sys_MIMO, u_matrix, Time, initial_conditions);

% Plot the results
figure();
subplot(2, 1, 1);
plot(t, SF_MWh*y(:, 2)); % Plot LH2 
hold on
plot(t, SF_MWh*y(:, 5)); % Plot NG 
hold on
plot(t, SF_MWh*y(:, 7)); % Plot H2 
% hold on
% plot(t, y(:, 3)); % Plot OWE 
grid on
grid minor
title('Prices');
ylabel('€/MWh');
xlim([0 365])
legend('p_{LH2}', 'p_{NG}', 'p_{H2}', 'p_{OWE}');

subplot(2, 1, 2);
plot(t, -SF_GWh*y(:, 1)); % Plot LH2 
hold on
% plot(t, y(:, 4)); % Plot OWE 
% hold on
plot(t, -SF_GWh*y(:, 4)); % Plot NG 
hold on
plot(t, -SF_GWh*y(:, 6)); % Plot H2 
title('Commodity Stock');
legend('q_{LH2}', 'q_{NG}', 'q_{H2}');
grid on
grid minor
xlim([0 365])
xlabel('Time');
ylabel('GWh');
hold off

