
% AE419 Spring 2017 HW1
% Siyun He 
%date: 02/03/2017
% netid: siyunhe2

%HW1-Problem2
%Use the Average method to predict SG+SR(With headwind = 15knots / No slope condition)
%Result SG+SR = 2.3086e+03


%The purpose of this Matlab code is to use a)average acceleration method
%and b)Numerical method(time-stepping method) to calculate the 
%ground run distance (SG + SR) of  the A-4C/L Skyhawk and compare the  
%results with Section XI in the flight manual. Use a
%gross weight of 18,000 lb for all cases.

%run time stepping method over a range of altitudes
% h varies from sea level to 10,000 ft in increments of 1,000ft 
%and compare with figure 11-8 in flight manual

%case1
%Assumptions
    %with headwind, Vw = 15kn
    %No slope,alpha = 0
    %Cd0 ~0.015-0025
    %eo ~0.80-0.90
    %standard atmosphere
    %tire coefficienty of friciton mu = 0.025

%case1
%Assumptions
    %headwind of 15kn(knots), Vw = 15kn
    %No slope,alpha = 0
    %Cd0 ~0.015-0025
    %eo ~0.80-0.90
    % W = 18,000 ib
    
%Type: A4 Skyhawk    
%Parameters
    %Span = 27.5 ft 
    %Wing Area = 260 ft^2  
    %AR=2.909
    %e_0 = 0.8
    %Cd_o = 0.018
    
    %Empty weight = 10,602 lb
    %gross weight = 15,783 lb
    %max weight = 24,3300 lb
    %max uninstalled thust at sea level = 11,187 lb
    
%Average Accleleration method
%S_G = (1/a_avg)*(V_G)^2
%V_G = V_LOE +/- V_wind
%a_avg = (g/W)*(T_avg - D_avg - W_alpha - mu*(W - L_avg))

%V_avg = V_G/(sqrt(2)) = 0.707 V_LOF

%Since lift and drag are function of Cl, L,D(Cl)
%Define to maximize a_avg with respect to Cl
%d(a_avg)/d(Cl) = d((g/W)*(T_avg - D_avg - W_alpha - mu*(W -
%L_avg)))/d(Cl)= 0

%declare a array of velocity
clear all
clc

W = 18000;
e_0 = 0.90;
AR = 2.909;
k = 1/(pi * e_0 * AR);
mu = 0.025;
b_Ft = 27.5;
Sref_Ft2 = 260;
Cd_o = 0.025;
g_fts2 = 32.1174;
Rho_Slugsft3 = 0.002377; %1 kg/m3  =  0.062428 lb/ft3
alpha = 0;
%Velocity
kn_to_fts = 1.68791; 
%1knot = 1.68781 ft/s
%V_stall = 116;

%Unit Conversion
m_to_foot = 3.28084;
foot_to_m = 1 / 3.28084;


%from take off distance chart  
%V_stall = sqrt(2*W/(Rho*S_ref*Cl_max));
V_wind = 15 *kn_to_fts;
V_LOF = 136 *kn_to_fts; %Knots
V_TAS = V_LOF - V_wind;
V_G = 0.707 * V_TAS;
V_avg = V_G;

%Coefficient
%Cl_opt = (2*W)./(Rho_Slugsft3 * Sref_Ft2 * (V_TAS.^2));
Cl_opt = mu/(2*k);

C_d = Cd_o + k * Cl_opt.^2;


Temp_R = 518.67;
gam_R = 1.4; %gamma 
R_gas = 1716.59; % ft2/(s2·°R)
a_SL= sqrt(gam_R*R_gas*Temp_R);
M = V_avg/a_SL;

q_avg = 0.5*(Rho_Slugsft3*(V_G).^2);


%All avg forces    
W_alpha = W*alpha;
L_avg = q_avg*Sref_Ft2*Cl_opt;
D_avg = q_avg*Sref_Ft2*C_d;
N_avg = W - L_avg;

%sigma = rho_h / rho_SL;
%Advanced turbojet(dry) Thrust
T_o = 11187;
Tx= (0.76 * (0.907 + 0.262 *((abs(M-0.5)).^1.5)))*(1.^0.7)* T_o;
T_avg = Tx;

a_avg = (g_fts2/W)*(T_avg - D_avg - mu*(W - L_avg));
%Acceleration 
%a = (g_fts2/W)*(T_avg - D - W_alpha - mu*(W - L));
S_G = (1/a_avg)*(V_avg).^2*0.5;