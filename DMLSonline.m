% Program DMLSonline
%
% Purpose: The function simulates the isotherm migration method (IMM)
% equations to find the response when laser power is decreased
% instantaneously to 50%.
%
% Inputs:
% ======
%
% (a) Power (P, in W)
% (b) Speed (v, in m/s)
% (c) Absorption coefficient (A)
% (d) Number of isotherms between T_0 and T_m (n1, integer)
% (e) plot_flag (boolean indicating whether to plot convergence or not)
%
% Outputs:
% =======
%
% (a) Time (t_array_, in s)
% (b) Isotherms (y_array_, in m)
% (c) Grid granularity (d_T, in ºC)
%
% Example:
% [t_array_,y_array_,d_T] = DMLSonline(195,0.8,0.6,5,true)
%
% Developed by: Felipe Lopez, based on Devesse's IMM model
%
% Date: 10/15/2015
%
% Ref: W. Devesse et al., Int. Journal of Heat and Mass Transfer 75(2015),
% pp. 726-735.
%
% Revision history
% ================
%
% 09/03/2015    Convergence was verified with order p = 3.
% 10/15/2015    First version deployed to GitHub.
% 10/16/2015    Adjusted maximum width instead of laser-height width
%

function [t_array_,y_array_,d_T] = DMLSonline(P,v,A,n1,plot_flag)
%% Definition of global variables
global T_m m alpha_0
global T dT T_0 T_max hl t_sim
%% Read thermophysical properties
T_m = 1320.0; % Melting temperature (ºC)
T_0 = 80.0; % Room temperature (ºC)
max_T = 2560.0; % Maximum (user-defined) temperature in the grid (ºC) (Li Ma said 2263.0)
hl = 2.97e+5; % Latent heat of fusion (J/kg)
%% Definition of the temperature grid
dT = -(T_m-T_0)/n1; % Temperature increment in isotherms (dT < 0, ºC)
T_max = T_0-floor((T_0-max_T)/dT)*dT; % Maximum temperature in grid (ºC)
n = round((T_0-T_max)/dT); % Number of gridpoints (integer)
T = T_max:dT:(T_0-dT); % Temperature grid (n-long array of temperatures)
m = length(T)-n1+1; % Location of melting isotherm (integer)
%% Definition of initial state
alpha_0 = k(T_max)/rho(T_max)/Cp(T_max);
y_nom = zeros(n,1);
% Rosenthal's solution
for i=1:n
   c1 = (T(i)-T_0)*2*pi*rho(T_max)*Cp(T_max)*alpha_0/A/P;
   c2 = 2*alpha_0/v;
   y_nom(i,1) = c2*lambertw(1/c1/c2); % Solution by Lambert's product-log
end
%% Simulate to find true steady-state
S = 4.0*alpha_0/v/v; % Characteristic time (s)
t_sim = 1.0e+5*S; % Simulation time: Looong time (s)
[t_array,y_array] = ode23s(@(t,x)SLM_Rate(t,x,[P,v],A,'steady'),[0 t_sim],y_nom);
%% Simulate transient
[t_array2,y_array2] = ode23s(@(t,x)SLM_Rate(t,x,[0.5*P,v],A,'transient'),[t_sim 2.0*t_sim],y_array(end,:)');
%% Concatenate
t_array_ = [t_array;t_array2];
y_array_ = [y_array;y_array2];
P_array_ = [P*ones(size(t_array)); 0.5*P*ones(size(t_array2))];
% Compute maximum width
L = 2.0*ThermalDiffusivity(T_m)/v;
C = y_array_(:,m).*exp(y_array_(:,m)/L);
a = C/2 + L/4*lambertw(2.0*C/L);
b = sqrt(L^2*lambertw(C./L.*exp((C-a)/L)).^2 - (C-a).^2);
width_vector = 2.0e+6*b;
%% Plot
if plot_flag
    figure(1)
    [ax,p1,p2] = plotyy(1.0e+3*t_array_,P_array_,1.0e+3*t_array_,width_vector,'plot','plot');
    xlim(ax(1),1.0e+3*[(t_sim-2*S) (t_sim+4*S)]);
    xlim(ax(2),1.0e+3*[(t_sim-2*S) (t_sim+4*S)]);
    grid(ax(1),'on'); p1.LineWidth = 2.0; p2.LineWidth = 2.0;
    xlabel(ax(1), 'Time (ms)'); ylabel(ax(1), 'Laser power (W)');
    ylabel(ax(2), 'Melt pool width (\mum)'); 
end
d_T = dT;
end

function alpha_transition = SmoothThermalDiffusivity(T,tau)
global alpha_0
alpha_transition = alpha_0+(ThermalDiffusivity(T)-alpha_0)*tau;
end

function alpha = ThermalDiffusivity(T)
alpha = k(T)./Cp(T)./rho(T); % Thermal diffusivity
end

function k=k(T)
global T_max T_0
% Linear interpolation for thermal conductivity given temperature-dependent
% data
T_k_array = [T_0 21.0 38.0 93.0 204.0 316.0 427.0 538.0 649.0 760.0 ...
    871.0 982.0 1314.0 1369.0 2000.0 T_max]; % % Array of temperatures (ºC)
%k_array = [9.8 9.8 10.1 10.8 12.5 14.1 15.7 17.5 19.0 20.8 22.8 25.2 ...
%    31.0 25.0 27.0 27.0]; % Array of thermal conductivities (W/m-C)
k_array = 27.0*ones(size(T_k_array));
k = interp1(T_k_array,k_array,T,'linear')'; % 1-D data interpolation
%k = 27.0*ones(size(T));
end

function Cp=Cp(T)
global T_max T_0
% Linear interpolation for specific heat given temperature-dependent
% data
T_Cp_array = [T_0 21.0 93.0 204.0 316.0 427.0 538.0 649.0 760.0 871.0 ...
    982.0 1093.0 T_max]; % % Array of temperatures (ºC)
Cp_array = [410 410 427 456 481 511 536 565 590 620 645 ...
    670 670]; % Array of specific heat
Cp = interp1(T_Cp_array,Cp_array,T,'linear')'; % 1-D data interpolation
%Cp = 670*ones(size(T));
end

function rho=rho(T)
global T_max T_0
% Linear interpolation for density given temperature-dependent
% data
T_rho_array = [T_0 27.0 227.0 727.0 1314.0 1317.0 1327.0 1337.0 1347.0 ...
    1537.0 1367.0 1369.0 1727.0 2227.0 T_max]; % % Array of temperatures (ºC)
rho_array = [8620 8620 8536 8316 8022 8014 7977 7928 7863 7777 7660 ...
    7636 7313 6860 6860]; % Array of densities (kg/m^3)
rho = interp1(T_rho_array,rho_array,T,'linear')'; % 1-D data interpolation
%rho = 6860*ones(size(T));
end

function y_dot=SLM_Rate(t,x,u,A,mode)
global alpha hl m
global T dT T_0 T_max T_m t_sim
%% Read inputs
y = x; % Read vertical isotherms (states) (m)
P = u(1); % Laser power (W)
v = u(2); % Scan speed (m/s)
%% Interpolate thermal diffusivity
% Make smooth transition to 'real' values during the first fifth of
% simulation
if strcmp(mode,'steady')
    if (t < t_sim/5.0)
        alpha = SmoothThermalDiffusivity(T,t*5.0/t_sim);
    else
        alpha = ThermalDiffusivity(T);
    end
else
    alpha = ThermalDiffusivity(T);
end
%% Construct rate equations
n = length(x); y_dot = zeros(n,1);
% For i=1
Cp_ = Cp(T(1)); rho_ = rho(T(1)); alpha_ = alpha(1);
L1 = 2*alpha_/v;
A_l = T_0-hl/(Cp_+2*alpha_*Cp_/v/y(m)+2*alpha_*alpha_*Cp_/y(m)/y(m)/v/v);
y_dot(1) = -2*alpha_*y(2)/y(1)/(y(2)-y(1))-(A*P/pi/rho_/Cp_/dT/y(1)/y(1))...
    *(1+y(1)/L1)*exp(-y(1)/L1)+(v*v/4/alpha_)*(T_max-A_l)*(y(2)-y(1))/dT;
% For 1 < i < m-1
for i=2:m-1
    Cp_ = Cp(T(i)); alpha_ = alpha(i);
    A_l = T_0-hl/(Cp_+2*alpha_*Cp_/v/y(m)+2*alpha_*alpha_*Cp_/y(m)/y(m)/v/v);
    y_dot(i) = -(alpha_/y(i))*(y(i+1)/(y(i+1)-y(i))-y(i-1)/(y(i)...
        -y(i-1)))+(v*v/8/alpha_)*(T(i)-A_l)*(y(i+1)-y(i-1))/dT;
end
% For i = m
Cp_ = Cp(T(m)); alpha_ = alpha(m);
A_l = T_0-hl/(Cp_+2*alpha_*Cp_/v/y(m)+2*alpha_*alpha_*Cp_/y(m)/y(m)/v/v);
Ste = -Cp_*dT/hl; % Stefan number
C1 = -(A_l-2*T_m+T_0)/dT;
C2 = (1/(y(m+1)-y(m))+1/(y(m)-y(m-1)))^-1;
C3 = 1+Ste^-1;
y_dot(m) = -(alpha_/y(m)/C3)*(y(m+1)/(y(m+1)-y(m)) - y(m-1)/(y(m)-y(m-1))) ...
    + (C1*C2*v*v/4/alpha_/C3)*(1+2*(Ste^-1)*(1+(C1*C2*v/2/alpha_)^2)^-1);
% For m < i < n-1
for i=m+1:n-1
    alpha_ = alpha(i);
    y_dot(i) = -(alpha_/y(i))*(y(i+1)/(y(i+1)-y(i))-y(i-1)/(y(i)...
        -y(i-1)))+(v*v/8/alpha_)*(T(i)-T_0)*(y(i+1)-y(i-1))/dT;
end
% For i = n
alpha_ = alpha(n);
y_dot(n) = -alpha_/y(n)*((y(n)-2*y(n-1))/(y(n)-y(n-1)))-(v*v/4/alpha_)*...
    (y(n)-y(n-1));
end