
% AE419 Spring 2017 HW1
% Siyun He 
%date: 02/03/2017
%resubmission date: 02/20/2017
% netid: siyunhe2

%The purpose of this Matlab code is to use 
%a)average acceleration method
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
    
%Type: A4 Skyhawk    
%Parameters
    %Span = 27.5 ft 
    %Wing Area = 260 ft^2  
    %AR=2.909
    %e_0 = 0.8
    %Cd_o = 0.018
    
    %Empty weight = 10,602 lb
    %gross weight = 15,783 lb
    %max weight = 24,3310 lb
    %max uninstalled thust at sea level = 11,187 lb

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
Rho_Slugsft3 = 0.002377;
alpha = 0;

kn = 1.68791; 
%1knot = 1.68781 ft/s
%V_stall = 116;

%Unit Conversion
m_to_foot = 3.28084;
foot_to_m = 1 / 3.28084;

f = 0;
b = 1000* foot_to_m;
c = 10000* foot_to_m;

ms_to_fts = (1./(3.048*10.^-1));
%1ft/s = 3.048×10?1 m/s 
kgm3_to_kbft3 = 1/16.01846337;
%1	lb/ft3	? 16.01846337 kg/m3

[T_x, a_x, P_x, rho_x] = atmosisa(f : b: c);
% Return of Temp in Kelvin
% Speed of Sound m/s
% Pressure in Pascal
% Density in kilograms per meter cubed

%Convert from International units to SI Units
temp_K = T_x;
temp_R_h = temp_K*(9/5);
a_ms = a_x;
%a_knot_h = a_ms*1.94384*kn;  
a_knot_h = a_ms * ms_to_fts;
p_pa = P_x; 
p_psi_h = p_pa*(0.0208854);
rho_kgm3 = rho_x;
%1 kg/m3 = 0.000036127298147753 lb/in3
rho_Ibft3_h = rho_kgm3*0.002377/1.225;  %unit in slug/ft3
%rho_Ibft3_h = rho_kgm3 * kgm3_to_kbft3 ;


%from take off distance chart  
%V_stall = sqrt(2*W/(Rho*S_ref*Cl_max));
V_wind = 0 * kn;
V_LOF = 136 * kn; %Knots
V_G = 0.707 * ( V_LOF - V_wind);
V_avg = V_G;

%ssl property
Temp_R = 518.67;
gam_R = 1.4; %gamma 
R_gas = 1716.59; % ft2/(s2·°R)
a_SL= sqrt(gam_R*R_gas*Temp_R);
M = V_avg/a_SL;
T_o = 11187;


%Initilization of Standard Atmosphere properties from built in function
h = 0:1000:10000;
rho_ssl = Rho_Slugsft3;
g = 32.174;
p_ssl = 2116.2;  %lbf/ft2




temp_h_k =  temp_R_h;
p_h =  p_psi_h;
rho_h = rho_Ibft3_h;
a_h =  a_knot_h;
 
z = 1;
for z = z: 11
 %Properties change with altitude from 0 - 10,000ft, ie h(j) =
 %0:10000:10000 
 %temp_h_k(:,j) = (59 + ((-3.57) * (h(j)/1000) )) + 459.67;
 %p_h(:,j) = (p_ssl * (g*h(j)) )/(R_gas * temp_h_k(:,j) );
 %rho_h(:,j) = p_h(:,j)/(R_gas * temp_h_k(:,j) );
    sigma(:,z) = rho_h(:,z) / rho_ssl;        
    sigma_reverse(:,z) = sqrt(rho_ssl/ rho_h(:,z));
    TAS(:,z) = V_LOF * sigma_reverse(:,z);
    square_TAS(:,z) = TAS(:,z).^2;
%Corrected EAS read directly from flight manual
%Transform to TAS by EAS = TAS * sqrt(rho/rho_ssl)
    
    
    
%    rho_h_transpose(:,z) = (rho_h(:,z)).';
%    Cl_opt(:,z) = rho_h_transpose(:,z).' * (square_TAS(:,z))* Sref_Ft2;
%    Cl_opt_t_h(:,z) = 2*W/Cl_opt(:,z);
%    C_d_t_h(:,z) = Cd_o + k * Cl_opt_t_h(:,z).^2;
 end;
 
rho_h_transpose = rho_h.';
k = 1;
for k = k:11
    
    Cl_opt(:,k) = rho_h_transpose(k,:) * (square_TAS(:,k))* Sref_Ft2;
    Cl_opt_t_h(:,k) = (2*W)./Cl_opt(:,k);
    C_d_t_h(:,k) = Cd_o + k * Cl_opt_t_h(:,k).^2;
end;


%time stepping method
%Calculate initial condition when t = 0, v= 0, 
T_t_h = zeros(331,11);
M_t_h = zeros(331,11);
t = zeros(331,11);
a = zeros(331,11);
V = zeros(331,11);
S = zeros(331,11);
S(1,:) = 0;
V(1,:) = 0;
t(1,:) = 0;
a(1,:) = 0;


%calcualte for initial condition for at t = 0 and t = 0.05s; then update
% there using time stepping method 
delta_t = 0.05;

%Set altitude range from 0ft to 10,000ft with 1000ft increment
%Thus we have all the data mapping in X x 11 matrix.
%each column represent the parameter change in according altitude
%First Calculate all initial conditions for alittude at
% 0ft and 1000ft
%Then use time stepping method 

x = 1;
for x = x:11
    M_t_h(1,x)= V(1,x)./a_knot_h(1,x); 
    T_t_h(1,x) = (0.76*(0.907+0.262*((abs(M_t_h(1,x))-0.5)).^1.5))*(sigma(1,x).^0.7)*T_o;
    t(2,:) = t(1,:) + delta_t;
    a(2,x) = (g_fts2./W)*(T_t_h(1,x) - mu*W);
    V(2,x) = V(1,x)+ a(2,x)*delta_t;
    S(2,x) = S(1,x) + (V(1,x)+V(2,x))*0.5*delta_t;
    M_t_h(2,x) = V(2,x)./a_knot_h(1,x);
    T_t_h(2,x) = (0.76*(0.907+0.262*((abs(M_t_h(2,x)-0.5)).^1.5)))*(sigma(1,x).^0.7)*T_o;
end;

%Use time stepping method

for j = 1 : 11
    for i = 2: 510
        if (V(i,j) >= TAS(:,j)) %break when the True lift off speed is reached
            break;
        elseif (V(i,j) == 0)
            break;    
        else
             t(i+1,:) = t(i,:)+delta_t;
             a(i+1,j) =(g_fts2/W)*(T_t_h(i,j) - mu* W -(C_d_t_h(:,j) - mu*Cl_opt_t_h(:,j) )* 0.5 * rho_h(:,j) * V(i,j).^2);
             V(i+1,j) = V(i,j)+ a(i+1,j)*delta_t;
             M_t_h(i+1,j) = V(i+1,j)/a_h(:,j);
             T_t_h(i+1,j) = (0.76 *( 0.907 +0.262*((abs(M_t_h(i+1,j)-0.5)).^1.5)))*(sigma(:,j).^0.7)*T_o;
             S(i+1,j) = S(i,j) + (V(i,j) + V(i+1,j))*0.5*delta_t; 
         end;
     end;
end;

W_t_h =  8000;

%Calculate the dynamic pressure q, 
 % Lift, Drag
 % Weight, Normal
 % forces gererated by the time stepping method
 
for j = 1: 11
    for i = 1: 509
        q_t_h(i,j) = 0.5*(rho_h(:,j)*V(i,j).^2); 
        L_t_h(i,j) = q_t_h(i,j)*Sref_Ft2*Cl_opt_t_h(:,j);
        D_t_h(i,j) = q_t_h(i,j)*Sref_Ft2*C_d_t_h(:,j);
        N_t_h(i,j) = W_t_h - L_t_h(i,j);
    end;
end;
 

 %Since now V,a,S are 509x11 matrix, which means at h=0, t=0, you can map
 %the according V
V_0 = nonzeros(V(:,1));
V_1 = nonzeros(V(:,2));
V_2 = nonzeros(V(:,3));
V_3 = nonzeros(V(:,4));
V_4 = nonzeros(V(:,5));
V_5 = nonzeros(V(:,6));
V_6 = nonzeros(V(:,7));
V_7 = nonzeros(V(:,8));
V_8 = nonzeros(V(:,9));
V_9 = nonzeros(V(:,10));
V_10 = nonzeros(V(:,11));

a_0 = nonzeros(a(:,1));
a_1 = nonzeros(a(:,2));
a_2 = nonzeros(a(:,3));
a_3 = nonzeros(a(:,4));
a_4 = nonzeros(a(:,5));
a_5 = nonzeros(a(:,6));
a_6 = nonzeros(a(:,7));
a_7 = nonzeros(a(:,8));
a_8 = nonzeros(a(:,9));
a_9 = nonzeros(a(:,10));
a_10 = nonzeros(a(:,11));

S_0 = nonzeros(S(:,1));
S_1 = nonzeros(S(:,2));
S_2 = nonzeros(S(:,3));
S_3 = nonzeros(S(:,4));
S_4 = nonzeros(S(:,5));
S_5 = nonzeros(S(:,6));
S_6 = nonzeros(S(:,7));
S_7 = nonzeros(S(:,8));
S_8 = nonzeros(S(:,9));
S_9 = nonzeros(S(:,10));
S_10 = nonzeros(S(:,11));

Ti = t(:,1);

D_t_h_n1(:,1) = nonzeros(D_t_h(:,1));  
L_t_h_n1(:,1) = nonzeros(L_t_h(:,1));  
T_t_h_n1(:,1) = nonzeros(T_t_h(:,1));  
N_t_h_n1(:,1) = nonzeros(N_t_h(:,1));  

%Plot all forces(Thust, Lift,Drag,weight and Normal forces)
%as function of time t, distance s, and V^2 for the take off

%The TA provided a much eaiser funciton for plotting
% h = 0 ft
%Plot1:Co-plot all forces vs. Time, Distance, and V^2 for the timestepping method

x1 = V_0.^2;                                           
W_t_h_n1 (1:331, 1:1) = 18000;  
N_t_h_n1 = W_t_h_n1 - L_t_h_n1;
out1 =1;
out1 = T_t_h_n1;
out1(1,:) = [];
T_t_h_n1 = out1;
 
y1 = D_t_h_n1(:,1);                                                                               
y2 = L_t_h_n1(:,1);                                                  
y3 = T_t_h_n1(:,1);   
y4 = W_t_h_n1;
y5 = N_t_h_n1;

Plotmeinput.x              = {x1, x1, x1, x1, x1};   
Plotmeinput.y              = {y1, y2, y3, y4, y5};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/m^2)';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H0ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H0ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2 = Ti(1:331,1:1);                                           

z1 = D_t_h_n1(:,1);                                                                                
z2 = L_t_h_n1(:,1);                                              
z3 = T_t_h_n1(:,1);   
z4 = W_t_h_n1;
z5 = N_t_h_n1;

Plotmeinput.x              = {x2, x2, x2, x2, x2};   
Plotmeinput.y              = {z1, z2, z3, z4, z5};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H0ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H0ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance
x3 = S_0;                                           

f1 = D_t_h_n1(:,1);                                                                               
f2 = L_t_h_n1(:,1);                                                  
f3 = T_t_h_n1(:,1);   
f4 = W_t_h_n1;
f5 = N_t_h_n1;

Plotmeinput.x              = {x3, x3, x3, x3, x3};   
Plotmeinput.y              = {f1, f2, f3, f4, f5};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H0ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H0ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4 = S_0;                                           

k1 = S_0;                                                                               
k2 = V_0;                                                  
k3 = a_0;   


Plotmeinput.x              = {x4, x4, x4};   
Plotmeinput.y              = {k1, k2, k3};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H0ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H0ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5 = Ti(1:331,1:1);                                           

n1 = S_0;                                                                               
n2 = V_0;                                                  
n3 = a_0;   


Plotmeinput.x              = {x5, x5, x5};   
Plotmeinput.y              = {n1, n2, n3};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H0ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H0ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%h=1000ft

x1_1000 = V_1.^2;  
W_t_h_n2 (1:344, 1:1) = 18000;  

D_t_h_n2(:,2) = nonzeros(D_t_h(:,2));  
L_t_h_n2(:,2) = nonzeros(L_t_h(:,2));  
T_t_h_n2(:,2) = nonzeros(T_t_h(:,2));  

out2 = 1;
out2 = T_t_h_n2;
out2(1,:) = [];
T_t_h_n2 = out2;

N_t_h_n2 = W_t_h_n2(:,1) - L_t_h_n2(:,2);

y1_1000 = D_t_h_n2(:,2);                                                                               
y2_1000 = L_t_h_n2(:,2);                                                  
y3_1000 = T_t_h_n2(:,2);   
y4_1000 = W_t_h_n2;
y5_1000 = N_t_h_n2;

Plotmeinput.x              = {x1_1000, x1_1000, x1_1000, x1_1000, x1_1000};   
Plotmeinput.y              = {y1_1000, y2_1000, y3_1000, y4_1000, y5_1000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H1000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H1000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_1000 = Ti(1:344, 1:1);                                           

z1_1000 = D_t_h_n2(:,2);                                                                               
z2_1000 = L_t_h_n2(:,2);                                                  
z3_1000 = T_t_h_n2(:,2);   
z4_1000 = W_t_h_n2;
z5_1000 = N_t_h_n2;

Plotmeinput.x              = {x2_1000, x2_1000, x2_1000, x2_1000, x2_1000};   
Plotmeinput.y              = {z1_1000, z2_1000, z3_1000, z4_1000, z5_1000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H1000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H1000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_1000 = S_1;                                           

f1_1000 = D_t_h_n2(:,2);                                                                               
f2_1000 = L_t_h_n2(:,2);                                                  
f3_1000 = T_t_h_n2(:,2);   
f4_1000 = W_t_h_n2;
f5_1000 = N_t_h_n2;

Plotmeinput.x              = {x3_1000, x3_1000, x3_1000, x3_1000, x3_1000};   
Plotmeinput.y              = {f1_1000, f2_1000, f3_1000, f4_1000, f5_1000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H1000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H1000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_1000 = S_1;                                           

k1_1000 = S_1;                                                                               
k2_1000 = V_1;                                                  
k3_1000 = a_1;   


Plotmeinput.x              = {x4_1000, x4_1000, x4_1000};   
Plotmeinput.y              = {k1_1000, k2_1000, k3_1000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H1000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H1000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_1000 = Ti(1:344, 1:1);                                           

n1_1000 = S_1;                                                                               
n2_1000 = V_1;                                                  
n3_1000 = a_1;   


Plotmeinput.x              = {x5_1000, x5_1000, x5_1000};   
Plotmeinput.y              = {n1_1000, n2_1000, n3_1000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H1000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H1000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%h=2000ft


x1_2000 = V_2.^2;  
W_t_h_n3 (1:359, 1:1) = 18000;  

D_t_h_n3 = nonzeros(D_t_h(:,3));  
L_t_h_n3 = nonzeros(L_t_h(:,3));  
T_t_h_n3 = nonzeros(T_t_h(:,3));  

out3 = 1;
out3 = T_t_h_n3;
out3(1,:) = [];
T_t_h_n3 = out3;

N_t_h_n3 = W_t_h_n3 - L_t_h_n3;
                                    
  
y1_2000 = D_t_h_n3;                                                                               
y2_2000 = L_t_h_n3;                                                  
y3_2000 = T_t_h_n3;   
y4_2000 = W_t_h_n3;
y5_2000 = N_t_h_n3;

Plotmeinput.x              = {x1_2000, x1_2000, x1_2000, x1_2000, x1_2000};   
Plotmeinput.y              = {y1_2000, y2_2000, y3_2000, y4_2000, y5_2000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H2000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H2000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_2000 = Ti(1:359, 1:1);                                          

z1_2000 = D_t_h_n3;                                                                               
z2_2000 = L_t_h_n3;                                                  
z3_2000 = T_t_h_n3;      
z4_2000 = W_t_h_n3;
z5_2000 = N_t_h_n3;

Plotmeinput.x              = {x2_2000, x2_2000, x2_2000, x2_2000, x2_2000};   
Plotmeinput.y              = {z1_2000, z2_2000, z3_2000, z4_2000, z5_2000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H2000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H2000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_2000 = S_2;                                           

f1_2000 = D_t_h_n3;                                                                          
f2_2000 = L_t_h_n3;                                      
f3_2000 = T_t_h_n3;      
f4_2000 = W_t_h_n3;   
f5_2000 = N_t_h_n3;   

Plotmeinput.x              = {x3_2000, x3_2000, x3_2000, x3_2000, x3_2000};   
Plotmeinput.y              = {f1_2000, f2_2000, f3_2000, f4_2000, f5_2000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H2000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H2000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_2000 = S_2;                                           

k1_2000 = S_2;                                                                               
k2_2000 = V_2;                                                  
k3_2000 = a_2;   


Plotmeinput.x              = {x4_2000, x4_2000, x4_2000};   
Plotmeinput.y              = {k1_2000, k2_2000, k3_2000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H2000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H2000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_2000 = Ti(1:359, 1:1);                                            

n1_2000 = S_2;                                                                               
n2_2000 = V_2;                                                  
n3_2000 = a_2;   


Plotmeinput.x              = {x5_2000, x5_2000, x5_2000};   
Plotmeinput.y              = {n1_2000, n2_2000, n3_2000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H2000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H2000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%h = 3000ft
x1_3000 = V_3.^2;    

W_t_h_n4 (1:374, 1:1) = 18000;  

D_t_h_n4 = nonzeros(D_t_h(:,4));  
L_t_h_n4 = nonzeros(L_t_h(:,4));  
T_t_h_n4 = nonzeros(T_t_h(:,4));  

out4 = 1;
out4 = T_t_h_n4;
out4(1,:) = [];
T_t_h_n4 = out4;

N_t_h_n4 = W_t_h_n4 - L_t_h_n4;
                                    
y1_3000 = D_t_h_n4;                                                                               
y2_3000 = L_t_h_n4;                                                  
y3_3000 = T_t_h_n4;   
y4_3000 = W_t_h_n4;
y5_3000 = N_t_h_n4;

Plotmeinput.x              = {x1_3000, x1_3000, x1_3000, x1_3000, x1_3000};   
Plotmeinput.y              = {y1_3000, y2_3000, y3_3000, y4_3000, y5_3000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H3000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H3000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_3000 = Ti(1:374, 1:1);                                           

z1_3000 = D_t_h_n4;                                                                               
z2_3000 = L_t_h_n4;                                                  
z3_3000 = T_t_h_n4;   
z4_3000 = W_t_h_n4;
z5_3000 = N_t_h_n4;

Plotmeinput.x              = {x2_3000, x2_3000, x2_3000, x2_3000, x2_3000};   
Plotmeinput.y              = {z1_3000, z2_3000, z3_3000, z4_3000, z5_3000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H3000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H3000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_3000 = S_3;                                           

f1_3000 = D_t_h_n4;                                                                               
f2_3000 = L_t_h_n4;                                                  
f3_3000 = T_t_h_n4;   
f4_3000 = W_t_h_n4;
f5_3000 = N_t_h_n4;

Plotmeinput.x              = {x3_3000, x3_3000, x3_3000, x3_3000, x3_3000};   
Plotmeinput.y              = {f1_3000, f2_3000, f3_3000, f4_3000, f5_3000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H3000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H3000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_3000 = S_3;                                           

k1_3000 = S_3;                                                                               
k2_3000 = V_3;                                                  
k3_3000 = a_3;   


Plotmeinput.x              = {x4_3000, x4_3000, x4_3000};   
Plotmeinput.y              = {k1_3000, k2_3000, k3_3000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H3000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H3000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_3000 = Ti(1:374, 1:1);                                           

n1_3000 = S_3;                                                                               
n2_3000 = V_3;                                                  
n3_3000 = a_3;   


Plotmeinput.x              = {x5_3000, x5_3000, x5_3000};   
Plotmeinput.y              = {n1_3000, n2_3000, n3_3000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H3000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H3000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%h = 4000ft

x1_4000 = V_4.^2;                                           
 
W_t_h_n5 (1:390, 1:1) = 18000;  

D_t_h_n5 = nonzeros(D_t_h(:,5));  
L_t_h_n5 = nonzeros(L_t_h(:,5));  
T_t_h_n5 = nonzeros(T_t_h(:,5));  

out4 = 1;
out5 = T_t_h_n5;
out5(1,:) = [];
T_t_h_n5 = out5;

N_t_h_n5 = W_t_h_n5 - L_t_h_n5;
                                    
y1_4000 = D_t_h_n5;                                                                               
y2_4000 = L_t_h_n5;                                                  
y3_4000 = T_t_h_n5;   
y4_4000 = W_t_h_n5;
y5_4000 = N_t_h_n5;

Plotmeinput.x              = {x1_4000, x1_4000, x1_4000, x1_4000, x1_4000};   
Plotmeinput.y              = {y1_4000, y2_4000, y3_4000, y4_4000, y5_4000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H4000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H4000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_4000 = Ti(1:390, 1:1);                                           

z1_4000 = D_t_h_n5;                                                                               
z2_4000 = L_t_h_n5;                                                  
z3_4000 = T_t_h_n5;   
z4_4000 = W_t_h_n5;
z5_4000 = N_t_h_n5;

Plotmeinput.x              = {x2_4000, x2_4000, x2_4000, x2_4000, x2_4000};   
Plotmeinput.y              = {z1_4000, z2_4000, z3_4000, z4_4000, z5_4000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H4000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H4000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_4000 = S_4;                                           

f1_4000 = D_t_h_n5;                                                                               
f2_4000 = L_t_h_n5;                                                  
f3_4000 = T_t_h_n5;
f4_4000 = W_t_h_n5;
f5_4000 = N_t_h_n5;

Plotmeinput.x              = {x3_4000, x3_4000, x3_4000, x3_4000, x3_4000};   
Plotmeinput.y              = {f1_4000, f2_4000, f3_4000, f4_4000, f5_4000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H4000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H4000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_4000 = S_4;                                           

k1_4000 = S_4;                                                                               
k2_4000 = V_4;                                                  
k3_4000 = a_4;   


Plotmeinput.x              = {x4_4000, x4_4000, x4_4000};   
Plotmeinput.y              = {k1_4000, k2_4000, k3_4000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H4000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H4000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_4000 = Ti(1:390,1:1);                                           

n1_4000 = S_4;                                                                               
n2_4000 = V_4;                                                  
n3_4000 = a_4;   


Plotmeinput.x              = {x5_4000, x5_4000, x5_4000};   
Plotmeinput.y              = {n1_4000, n2_4000, n3_4000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H4000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H4000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%h = 5000ft

x1_5000 = V_5.^2; 
                                        
W_t_h_n6 (1:407, 1:1) = 18000;  

D_t_h_n6 = nonzeros(D_t_h(:,6));  
L_t_h_n6 = nonzeros(L_t_h(:,6));  
T_t_h_n6 = nonzeros(T_t_h(:,6));  

out6 = 1;
out6 = T_t_h_n6;
out6(1,:) = [];
T_t_h_n6 = out6;

N_t_h_n6 = W_t_h_n6 - L_t_h_n6;
                                   
y1_5000 = D_t_h_n6;                                                                               
y2_5000 = L_t_h_n6;                                                  
y3_5000 = T_t_h_n6;   
y4_5000 = W_t_h_n6;
y5_5000 = N_t_h_n6;

Plotmeinput.x              = {x1_5000, x1_5000, x1_5000, x1_5000, x1_5000};   
Plotmeinput.y              = {y1_5000, y2_5000, y3_5000, y4_5000, y5_5000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H5000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H5000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_5000 = Ti(1:407, 1:1);                                           

z1_5000 = D_t_h_n6;                                                                               
z2_5000 = L_t_h_n6;                                                  
z3_5000 = T_t_h_n6;   
z4_5000 = W_t_h_n6;
z5_5000 = N_t_h_n6;

Plotmeinput.x              = {x2_5000, x2_5000, x2_5000, x2_5000, x2_5000};   
Plotmeinput.y              = {z1_5000, z2_5000, z3_5000, z4_5000, z5_5000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H5000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H5000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_5000 = S_5;                                           

f1_5000 = D_t_h_n6;                                                                               
f2_5000 = L_t_h_n6;                                                  
f3_5000 = T_t_h_n6;   
f4_5000 = W_t_h_n6;
f5_5000 = N_t_h_n6;

Plotmeinput.x              = {x3_5000, x3_5000, x3_5000, x3_5000, x3_5000};   
Plotmeinput.y              = {f1_5000, f2_5000, f3_5000, f4_5000, f5_5000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H5000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H5000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_5000 = S_5;                                           

k1_5000 = S_5;                                                                               
k2_5000 = V_5;                                                  
k3_5000 = a_5;   


Plotmeinput.x              = {x4_5000, x4_5000, x4_5000};   
Plotmeinput.y              = {k1_5000, k2_5000, k3_5000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H5000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H5000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_5000 = Ti(1:407, 1:1);                                                  

n1_5000 = S_5;                                                                               
n2_5000 = V_5;                                                  
n3_5000 = a_5;   


Plotmeinput.x              = {x5_5000, x5_5000, x5_5000};   
Plotmeinput.y              = {n1_5000, n2_5000, n3_5000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H5000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H5000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%h = 6000ft

x1_6000 = V_6.^2;    
                                        
W_t_h_n7 (1:425, 1:1) = 18000;  

D_t_h_n7 = nonzeros(D_t_h(:,7));  
L_t_h_n7 = nonzeros(L_t_h(:,7));  
T_t_h_n7 = nonzeros(T_t_h(:,7));  

out7 = 1;
out7 = T_t_h_n7;
out7(1,:) = [];
T_t_h_n7 = out7;

N_t_h_n7 = W_t_h_n7 - L_t_h_n7;
                                     
y1_6000 = D_t_h_n7;                                                                               
y2_6000 = L_t_h_n7;                                                  
y3_6000 = T_t_h_n7;   
y4_6000 = W_t_h_n7;
y5_6000 = N_t_h_n7;

Plotmeinput.x              = {x1_6000, x1_6000, x1_6000, x1_6000, x1_6000};   
Plotmeinput.y              = {y1_6000, y2_6000, y3_6000, y4_6000, y5_6000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H6000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H6000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_6000 = Ti(1:425, 1:1);                                           

z1_6000 = D_t_h_n7;                                                                               
z2_6000 = L_t_h_n7;                                                  
z3_6000 = T_t_h_n7;   
z4_6000 = W_t_h_n7;
z5_6000 = N_t_h_n7;

Plotmeinput.x              = {x2_6000, x2_6000, x2_6000, x2_6000, x2_6000};   
Plotmeinput.y              = {z1_6000, z2_6000, z3_6000, z4_6000, z5_6000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H6000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H6000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_6000 = S_6;                                           

f1_6000 = D_t_h_n7;                                                                               
f2_6000 = L_t_h_n7;                                                  
f3_6000 = T_t_h_n7;   
f4_6000 = W_t_h_n7;
f5_6000 = N_t_h_n7;

Plotmeinput.x              = {x3_6000, x3_6000, x3_6000, x3_6000, x3_6000};   
Plotmeinput.y              = {f1_6000, f2_6000, f3_6000, f4_6000, f5_6000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H6000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H6000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_6000 = S_6;                                           

k1_6000 = S_6;                                                                               
k2_6000 = V_6;                                                  
k3_6000 = a_6;   


Plotmeinput.x              = {x4_6000, x4_6000, x4_6000};   
Plotmeinput.y              = {k1_6000, k2_6000, k3_6000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H6000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H6000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_6000 = Ti(1:425, 1:1);                                           

n1_6000 = S_6;                                                                               
n2_6000 = V_6;                                                  
n3_6000 = a_6;   


Plotmeinput.x              = {x5_6000, x5_6000, x5_6000};   
Plotmeinput.y              = {n1_6000, n2_6000, n3_6000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H6000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H6000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

% h =7000ft

x1_7000 = V_7.^2;       
                                        
W_t_h_n8 (1:444, 1:1) = 18000;  

D_t_h_n8 = nonzeros(D_t_h(:,8));  
L_t_h_n8 = nonzeros(L_t_h(:,8));  
T_t_h_n8 = nonzeros(T_t_h(:,8));  

out8 = 1;
out8 = T_t_h_n8;
out8(1,:) = [];
T_t_h_n8 = out8;

N_t_h_n8 = W_t_h_n8 - L_t_h_n8;
                                      
y1_7000 = D_t_h_n8;                                                                               
y2_7000 = L_t_h_n8;                                                  
y3_7000 = T_t_h_n8;   
y4_7000 = W_t_h_n8;
y5_7000 = N_t_h_n8;

Plotmeinput.x              = {x1_7000, x1_7000, x1_7000, x1_7000, x1_7000};   
Plotmeinput.y              = {y1_7000, y2_7000, y3_7000, y4_7000, y5_7000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H7000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H7000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_7000 = Ti(1:444, 1:1);                                           

z1_7000 = D_t_h_n8;                                                                               
z2_7000 = L_t_h_n8;                                                  
z3_7000 = T_t_h_n8;   
z4_7000 = W_t_h_n8;
z5_7000 = N_t_h_n8;

Plotmeinput.x              = {x2_7000, x2_7000, x2_7000, x2_7000, x2_7000};   
Plotmeinput.y              = {z1_7000, z2_7000, z3_7000, z4_7000, z5_7000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H7000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H7000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_7000 = S_7;                                           

f1_7000 = D_t_h_n8;                                                                               
f2_7000 = L_t_h_n8;                                                  
f3_7000 = T_t_h_n8;   
f4_7000 = W_t_h_n8;
f5_7000 = N_t_h_n8;

Plotmeinput.x              = {x3_7000, x3_7000, x3_7000, x3_7000, x3_7000};   
Plotmeinput.y              = {f1_7000, f2_7000, f3_7000, f4_7000, f5_7000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H7000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H7000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_7000 = S_7;                                           

k1_7000 = S_7;                                                                               
k2_7000 = V_7;                                                  
k3_7000 = a_7;   


Plotmeinput.x              = {x4_7000, x4_7000, x4_7000};   
Plotmeinput.y              = {k1_7000, k2_7000, k3_7000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H7000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H7000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_7000 = Ti(1:444, 1:1);                                           

n1_7000 = S_7;                                                                               
n2_7000 = V_7;                                                  
n3_7000 = a_7;   


Plotmeinput.x              = {x5_7000, x5_7000, x5_7000};   
Plotmeinput.y              = {n1_7000, n2_7000, n3_7000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H7000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H7000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%h = 8000ft

x1_8000 = V_8.^2;      
                                            
W_t_h_n9 (1:464, 1:1) = 18000;  

D_t_h_n9 = nonzeros(D_t_h(:,9));  
L_t_h_n9 = nonzeros(L_t_h(:,9));  
T_t_h_n9 = nonzeros(T_t_h(:,9));  

out9 = 1;
out9 = T_t_h_n9;
out9(1,:) = [];
T_t_h_n9 = out9;

N_t_h_n9 = W_t_h_n9 - L_t_h_n9;
                                      
y1_8000 = D_t_h_n9;                                                                               
y2_8000 = L_t_h_n9;                                                  
y3_8000 = T_t_h_n9;   
y4_8000 = W_t_h_n9;
y5_8000 = N_t_h_n9;

Plotmeinput.x              = {x1_8000, x1_8000, x1_8000, x1_8000, x1_8000};   
Plotmeinput.y              = {y1_8000, y2_8000, y3_8000, y4_8000, y5_8000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H8000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H8000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_8000 = Ti(1:464, 1:1);                                           

z1_8000 = D_t_h_n9;                                                                               
z2_8000 = L_t_h_n9;                                                  
z3_8000 = T_t_h_n9;   
z4_8000 = W_t_h_n9;
z5_8000 = N_t_h_n9;

Plotmeinput.x              = {x2_8000, x2_8000, x2_8000, x2_8000, x2_8000};   
Plotmeinput.y              = {z1_8000, z2_8000, z3_8000, z4_8000, z5_8000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H8000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H8000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_8000 = S_8;                                           

f1_8000 = D_t_h_n9;                                                                               
f2_8000 = L_t_h_n9;                                                  
f3_8000 = T_t_h_n9;   
f4_8000 = W_t_h_n9;
f5_8000 = N_t_h_n9;

Plotmeinput.x              = {x3_8000, x3_8000, x3_8000, x3_8000, x3_8000};   
Plotmeinput.y              = {f1_8000, f2_8000, f3_8000, f4_8000, f5_8000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H8000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H8000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_8000 = S_8;                                           

k1_8000 = S_8;                                                                               
k2_8000 = V_8;                                                  
k3_8000 = a_8;   


Plotmeinput.x              = {x4_8000, x4_8000, x4_8000};   
Plotmeinput.y              = {k1_8000, k2_8000, k3_8000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H8000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H8000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_8000 = Ti(1:464, 1:1);                                           

n1_8000 = S_8;                                                                               
n2_8000 = V_8;                                                  
n3_8000 = a_8;   


Plotmeinput.x              = {x5_8000, x5_8000, x5_8000};   
Plotmeinput.y              = {n1_8000, n2_8000, n3_8000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H8000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H8000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%h = 9000ft
x1_9000 = V_9.^2;       
                                            
W_t_h_n10 (1:485, 1:1) = 18000;  

D_t_h_n10 = nonzeros(D_t_h(:,10));  
L_t_h_n10 = nonzeros(L_t_h(:,10));  
T_t_h_n10 = nonzeros(T_t_h(:,10));  

out10 = 1;
out10 = T_t_h_n10;
out10(1,:) = [];
T_t_h_n10 = out10;

N_t_h_n10 = W_t_h_n10 - L_t_h_n10;

y1_9000 = D_t_h_n10;                                                                               
y2_9000 = L_t_h_n10;                                                  
y3_9000 = T_t_h_n10;   
y4_9000 = W_t_h_n10;
y5_9000 = N_t_h_n10;

Plotmeinput.x              = {x1_9000, x1_9000, x1_9000, x1_9000, x1_9000};   
Plotmeinput.y              = {y1_9000, y2_9000, y3_9000, y4_9000, y5_9000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H9000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H9000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_9000 = Ti(1:485, 1:1);                                           

z1_9000 = D_t_h_n10;                                                                               
z2_9000 = L_t_h_n10;                                                  
z3_9000 = T_t_h_n10;   
z4_9000 = W_t_h_n10;
z5_9000 = N_t_h_n10;

Plotmeinput.x              = {x2_9000, x2_9000, x2_9000, x2_9000, x2_9000};   
Plotmeinput.y              = {z1_9000, z2_9000, z3_9000, z4_9000, z5_9000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H9000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H9000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_9000 = S_9;                                           

f1_9000 = D_t_h_n10;                                                                               
f2_9000 = L_t_h_n10;                                                  
f3_9000 = T_t_h_n10;   
f4_9000 = W_t_h_n10;
f5_9000 = N_t_h_n10;

Plotmeinput.x              = {x3_9000, x3_9000, x3_9000, x3_9000, x3_9000};   
Plotmeinput.y              = {f1_9000, f2_9000, f3_9000, f4_9000, f5_9000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H9000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H9000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_9000 = S_9;                                           

k1_9000 = S_9;                                                                               
k2_9000 = V_9;                                                  
k3_9000 = a_9;   


Plotmeinput.x              = {x4_9000, x4_9000, x4_9000};   
Plotmeinput.y              = {k1_9000, k2_9000, k3_9000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H9000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H9000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_9000 = Ti(1:485, 1:1);                                           

n1_9000 = S_9;                                                                               
n2_9000 = V_9;                                                  
n3_9000 = a_9;   


Plotmeinput.x              = {x5_9000, x5_9000, x5_9000};   
Plotmeinput.y              = {n1_9000, n2_9000, n3_9000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H9000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H9000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%h = 10,000ft
x1_10000 = V_10.^2;                                                
                                            
W_t_h_n11 (1:508, 1:1) = 18000;  

D_t_h_n11 = nonzeros(D_t_h(:,11));  
L_t_h_n11 = nonzeros(L_t_h(:,11));  
T_t_h_n11 = nonzeros(T_t_h(:,11));  

out11 = 1;
out11 = T_t_h_n11;
out11(1,:) = [];
T_t_h_n11 = out11;

N_t_h_n11 = W_t_h_n11 - L_t_h_n11;

y1_10000 = D_t_h_n11;                                                                               
y2_10000 = L_t_h_n11;                                                  
y3_10000 = T_t_h_n11;   
y4_10000 = W_t_h_n11;
y5_10000 = N_t_h_n11;

Plotmeinput.x              = {x1_10000, x1_10000, x1_10000, x1_10000, x1_10000};   
Plotmeinput.y              = {y1_10000, y2_10000, y3_10000, y4_10000, y5_10000};  
Plotmeinput.xlabelname     = 'V.^2(ft^2/s^2)^';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H10000ft Time stepping Co-plotting forces as a function of V^2 '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(V^2)', 'Lift(V^2)', 'Thrust(V^2)','Weight(V^2)', 'Normal(V^2)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H10000ft  Time stepping force(V^2) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 


%Plot2,Co-plot all forces vs. Time
x2_10000 = Ti(1:508, 1:1);                                           

z1_10000 = D_t_h_n11;                                                                               
z2_10000 = L_t_h_n11;                                                  
z3_10000 = T_t_h_n11;   
z4_10000 = W_t_h_n11;
z5_10000 = N_t_h_n11;
Plotmeinput.x              = {x2_10000, x2_10000, x2_10000, x2_10000, x2_10000};   
Plotmeinput.y              = {z1_10000, z2_10000, z3_10000, z4_10000, z5_10000};  
Plotmeinput.xlabelname     = 't/sec';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H10000ft Time stepping Co-plotting forces as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(t)', 'Lift(t)', 'Thrust(t)','Weight(t)', 'Normal(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H10000ft Time steppingforce(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot3.Co-plot all forces vs.  Distance

x3_10000 = S_10;                                           

f1_10000 = D_t_h_n11;                                                                               
f2_10000 = L_t_h_n11;                                                  
f3_10000 = T_t_h_n11;   
f4_10000 = W_t_h_n11;
f5_10000 = N_t_h_n11;

Plotmeinput.x              = {x3_10000, x3_10000, x3_10000, x3_10000, x3_10000};   
Plotmeinput.y              = {f1_10000, f2_10000, f3_10000, f4_10000, f5_10000};  
Plotmeinput.xlabelname     = 'S';                
Plotmeinput.ylabelname     = 'F/Ib';                    
Plotmeinput.title          = 'H10000ftn Time stepping Co-plotting forces as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'Drag(S)', 'Lift(S)', 'Thrust(S)','Weight(S)', 'Normal(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H10000ft Time stepping Method force(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 5;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot4. Coplot d(s), V(s), a(s)
x4_10000 = S_10;                                           

k1_10000 = S_10;                                                                               
k2_10000 = V_10;                                                  
k3_10000 = a_10;   


Plotmeinput.x              = {x4_10000, x4_10000, x4_10000};   
Plotmeinput.y              = {k1_10000, k2_10000, k3_10000};  
Plotmeinput.xlabelname     = 'S/ft';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H10000ft Time stepping Co-plotting S,V,a as a function of S '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(S)', 'V(S)', 'a(S)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H10000ft Time stepping Method for d(S), V(S), a(S) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 

%Plot5. Coplot d(s), V(s), a(s)
x5_10000 = Ti(1:508, 1:1);                                           

n1_10000 = S_10;                                                                               
n2_10000 = V_10;                                                  
n3_10000 = a_10;   


Plotmeinput.x              = {x5_10000, x5_10000, x5_10000};   
Plotmeinput.y              = {n1_10000, n2_10000, n3_10000};  
Plotmeinput.xlabelname     = 't/s';                
Plotmeinput.ylabelname     = 'S,V,a';                    
Plotmeinput.title          = 'H10000ft Time stepping Co-plotting S,V,a as a function of t '; % Specify the plot title
Plotmeinput.legendoption   = 'On';                % Use 'Off' when plotting a single curve. Use 'On' when plotting two or more curves on the same graph
Plotmeinput.legendnames    = {'d(t)', 'V(t)', 'a(t)'};   % Provide the Legend entries here
Plotmeinput.filename       = 'H10000ftTime stepping Method for d(t), V(t), a(t) compare';   % Provide a filename to the figure being saved
Plotmeinput.nVar           = 3;     % Number of curves being plotted on the same graph
Plotme(Plotmeinput) 















