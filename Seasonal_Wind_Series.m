% clear all
% clc
% close all

%% wind data: 
%Data TNO
X=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]; % per month *365/11
Y=[11 12 11.1 8.7 8.5 8.7 8.2 8 8.6 10.3 10.1 11 12 11.1 8.7 8.5 8.7 8.2 8 8.6 10.3 10.1]; % in m/s

[f1,g1]=fit(X',Y','fourier5');


color_yellow=[1, 0.9,0.6];
color_dark_blue= [0.73, 0.78, 0.83];
color_light_blue= [0.8, 0.85, 0.93];
color_dark_gray= [0.74, 0.74, 0.74];
color_light_gray= [0.86, 0.86, 0.86];
color_green=[0.77, 0.84, 0.77];
figure

figure
% plot(out.yearly,'Color',color_yellow)
% hold on

hold on
% plot(f)
% hold on
% plot(f2)
% hold on
% figure


% color_blue= [186, 200, 211];
% hexColor_blue='#BAC8D3',  'Color', color_blue
plot(Wind_Distr_2019_input(:,1)/365*12,Wind_Distr_2019_input(:,2), 'Color',color_dark_blue)
hold on
scatter(X,Y,"filled",'MarkerFaceColor', color_yellow )

hold on
plot(f1,'black');
legend off
grid on
grid minor
gridColor='gray';
xlabel('Month');
ylabel('Wind speed (m/s)');
title('Wind distribution ');
axis([0 12 0 20]);
hold off
%%

Wind_plot=[Wind_Distr_2019_Scaled(:,1)/365*12,Wind_Distr_2019_Scaled(:,2)];

%%
     % Coefficients
       a0 =       9.564  ;
       a1 =       1.636  ;
       b1 =      0.4314  ;
       a2 =      0.2612  ;
       b2 =      0.2799  ;
       a3 =     -0.4869  ;
       b3 =      0.4023  ;
       a4 =      -0.141  ;
       b4 =     0.06985  ;
       a5 =      0.1669  ;
       b5 =     0.03972  ;
       w =     0.01721  ;

yearly_w =  a0 + a1*cos(X*w) + b1*sin(X*w) + a2*cos(2*X*w) + b2*sin(2*X*w) + a3*cos(3*X*w) + b3*sin(3*X*w) + a4*cos(4*X*w) + b4*sin(4*X*w) + a5*cos(5*X*w) + b5*sin(5*X*w);
figure
plot(X,yearly_w)

cut_in=3;
cut_out=4;




%% Daily new:
X3=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]/25; % per hour
Y3=[10.3 10.2 10.15 10.1 10.05 10 10 10 10 10 10 10 10 10 10.1 10.2 10.3 10.4 10.5 10.6 10.65 10.6 10.5 10.4 10.3]; % in m/s

figure
plot(X3,Y3);

figure
scatter(X3,Y3,"filled",'MarkerFaceColor', color_yellow )
hold on
[f3,g3]=fit(X3',Y3','fourier3');
hold on
plot(f3,'blue')

% Create the second subplot
figure
scatter(X3,Y3,"filled",'MarkerFaceColor', color_yellow )
hold on
[f3,g3]=fit(X3',Y3','fourier3');
hold on
plot(f3,'black')
legend off
grid on
grid minor
gridColor='gray';
xticks([0 4 8 12 16 20 24])
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','00:00'})
xlabel('Hour');
ylabel('Wind speed (m/s)');
title('Daily wind distribution');
axis([0 24 6 14]);

% % Save the second subplot
% saveas(gcf, 'winddistr_daily.eps'); 



       aa0 =       10.21  ;
       aa1 =      0.1631 ;
       bb1 =     -0.2515 ;
       aa2 =    -0.05622  ;
       bb2 =    -0.08132  ;
       aa3 =    -0.01446  ;
       bb3 =     0.01009  ;
       ww =      0.2618  ;
X=out.X(:); 
daily_w=aa0 + aa1*cos(X*ww) + bb1*sin(X*ww) + aa2*cos(2*X*ww) + bb2*sin(2*X*ww) + aa3*cos(3*X*ww) + bb3*sin(3*X*ww);


% figure
% pot(X,daily_w)







%%
time=1:0.00001:365;
wind=yearly_w+daily_w;
figure
plot(time,yearly_w)
% hold on
% plot(X,wind)

%% Prodution pattern:

total=integral(yearly_w,0,365)


%% power curve

color_yellow=1/255*[255, 230,150];

% Power curve
    wind_speeds = [0, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 25, 26]; % Wind speed values
    power_output = [0, 0, 0, 0.5, 1.8, 3.3, 5.1, 7.1, 9.5, 9.5, 9.5, 9.5, 9.5, 9.5, 9.5, 0]; % Corresponding power output values 
    %powercurve: https://www.researchgate.net/figure/Power-curve-of-the-Vestas-V164-95-turbine-including-some-technical-data-left-side-and_fig4_351600993

    % Interpolate the power curve
    % power_curve = interp1(wind_speeds, power_output, wind, 'linear', 'extrap');  % 'extrap' for extrapolationRR
    % 
    % % Limit power output to positive values
    % power_curve(power_curve < 0) = 0;

figure
plot(wind_speeds,power_output,'black','LineWidth',1.2)
ylabel('Power (MW)');
xlabel('Wind speed (m/s)');
grid on
grid minor