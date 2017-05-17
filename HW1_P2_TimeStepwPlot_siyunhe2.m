
% AE419 Spring 2017 HW1
% Siyun He 
%date: 02/03/2017
% netid: siyunhe2

%Problem2 time stepping method by using Heun's method

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
    %No wind, Vw = 0
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
    %max weight = 24,500 lb
    %max uninstalled thust at sea level = 11,187 lb
    
%Average Accleleration method
%S_G = (1/a_avg)*(V_G)^2
%V_G = V_LOE +/- V_wind
%a_avg = (g/W)*(T_avg - D_avg - W_alpha - mu*(W - L_avg))

%in case 1, with no wind conditon,
%V_avg = V_G/(sqrt(2)) = 0.707 V_LOF

%Since lift and drag are function of Cl, L,D(Cl)
%Define to maximize a_avg with respect to Cl
%d(a_avg)/d(Cl) = d((g/W)*(T_avg - D_avg - W_alpha - mu*(W -
%L_avg)))/d(Cl)= 0

clear all
clc
%declare a array of velocity

W = 18000;
e_0 = 0.90;
AR = 2.909;
k = 1/(pi*e_0*AR);
mu = 0.025;
b_Ft = 27.5;
Sref_Ft2 = 260;
Cd_o = 0.025;
g_fts2 = 32.1174;
Rho_Slugsft3 = 0.002377;
alpha = 0;
%Velocity
%Unit Conversion
kn_to_fts = 1.68791;

V_LOF = 136*kn_to_fts; %Knots
V_G = 0.707*V_LOF;
V_avg = V_G;

%Coefficient
Cl_opt = mu/(2*k);
%Cl_opt = 2*W/(Rho_Slugsft3*Sref_Ft2*(V_LOF).^2);
C_d = Cd_o + k*Cl_opt.^2;

Temp_R = 518.67;
gam_R = 1.4;  
R_gas = 1716.59; % ft2/(s2·°R)
a_SL= sqrt(gam_R*R_gas*Temp_R);
M = V_avg/a_SL;

%Advanced turbojet(dry) Thrust
T_o = 11187;
Tx= (0.76*(0.907+0.262*((abs(M-0.5)).^1.5)))*(1.^0.7)*T_o;


%time stepping method
%Calculate initial condition when t = 0, v= 0

S_o = 0;
V_o = 0;
M_o = V_o/a_SL;
T_vo = (0.76*(0.907+0.262*((abs(M_o-0.5)).^1.5)))*(1.^0.7)*T_o;

T_t = zeros(330,1);
M_t = zeros(330,1);
t = zeros(330,1);
a = zeros(330,1);
V = zeros(330,1);
S = zeros(330,1);
delta_t = 0.05;
a_o = (g_fts2/W)*(T_vo - mu*W);
V_1 = V_o+ a_o*delta_t;
S_1 = S_o + (V_o+V_1)*0.5*delta_t;
V(1) = V_1;
t(1) = 0;
a(1) = a_o;
S(1) = S_1;
M_t(1) = V(1)/a_SL;
T_t(1) = (0.76*(0.907+0.262*((abs(M_t(1)-0.5)).^1.5)))*(1.^0.7)*T_o;
i = 1;
for i = i: 330
   if (V(i,:) <= V_LOF)
        t(i+1) = t(i)+delta_t;
        a(i+1,:) =(g_fts2/W)*(T_t(i) - mu* W -(C_d-mu*Cl_opt)*0.5*Rho_Slugsft3*V(i,:).^2);
        V(i+1,:) = V(i,:)+ a(i+1,:)*delta_t;
        M_t(i+1,:) = V(i+1,:)/a_SL;
        T_t(i+1,:) = (0.76*(0.907+0.262*((abs(M_t(i+1,:)-0.5)).^1.5)))*(1.^0.7)*T_o;
        S(i+1,:) = S(i,:) + (V(i,:) + V(i+1,:))*0.5*delta_t;   
    else 
        break;
    end;
end;


q_t = 1/2*(Rho_Slugsft3*(V).^2);

%All forces    
W_alpha = W*alpha;
W_t = zeros(330,1);
N_t = zeros(330,1);
W_t(:) = 18000;
L_t = q_t*Sref_Ft2*Cl_opt;
D_t = q_t*Sref_Ft2*C_d;
N_t = W_t-L_t;

%The TA provided a much eaiser funciton for plotting

%Plot1:Co-plot all forces vs. Time, Distance, and V^2 for the timestepping method
x1 = V.^2;                                           
  

y1 = D_t;                                                                               
y2 = L_t;                                                  
y3 = T_t;   
y4 = W_t;
y5 = N_t;

Plotmeinput.x              = {x1, x1, x1, x1, x1};   
Plotmeinput.y              = {y1, y2, y3, y4, y5};  
Plotmeinput.xlabelname     = 'V.^2';                
Plotmeinput.ylabelname     = 'F';                    
Plotmeinput.title          = 'No wind Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'No wind  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2 = t;                                           

z1 = D_t;                                                                               
z2 = L_t;                                                  
z3 = T_t;   
z4 = W_t;
z5 = N_t;

Plotmeinput.x              = {x2, x2, x2, x2, x2};   
Plotmeinput.y              = {z1, z2, z3, z4, z5};  
Plotmeinput.xlabelname     = 't';                
Plotmeinput.ylabelname     = 'F';                    
Plotmeinput.title          = 'No wind Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'No wind  Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance
x3 = S;                                           

f1 = D_t;                                                                               
f2 = L_t;                                                  
f3 = T_t;   
f4 = W_t;
f5 = N_t;

Plotmeinput.x              = {x3, x3, x3, x3, x3};   
Plotmeinput.y              = {f1, f2, f3, f4, f5};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F';                    
Plotmeinput.title          = 'No wind Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'No wind Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4 = S;                                           

k1 = S;                                                                               
k2 = V;                                                  
k3 = a;   


Plotmeinput.x              = {x4, x4, x4};   
Plotmeinput.y              = {k1, k2, k3};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'No wind Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'No wind Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5 = t;                                           

n1 = S;                                                                               
n2 = V;                                                  
n3 = a;   


Plotmeinput.x              = {x5, x5, x5};   
Plotmeinput.y              = {n1, n2, n3};  
Plotmeinput.xlabelname     = 't';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'No wind Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'No wind Time stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 



