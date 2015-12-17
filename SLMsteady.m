% PROGRAM SLMsteady
%
% PURPOSE: The function simulates the isotherm migration method (IMM)
% equations to find the  steady-state location of isotherms.
%
% INPUTS:
% (a) Power (P, in W)
% (b) Speed (v, in m/s)
% (c) Absorption coefficient (A)
% (d) Diffusion efficiency (mu)
% (e) Number of isotherms between ambient T_0 and melting T_m (n1, integer)
% (f) plot_flag (boolean, plot convergence or not?)
%
% OUTPUTS:
% (a) Temperature grid (T_ss, in ºC)
% (b) Steady-state isotherm (half-width) location (y_ss, in m)
% (c) Grid granularity (d_T, in ºC)
% (d) Widths at laser location and maximum melt pool width
%    (width.laser and width.max, both in microns)
%
% Example:
% [T_ss,y_ss,d_T,width] = SLMsteady(195,0.800,0.6,1.0,5,false)
%
% Developed by: Felipe Lopez, based on Devesse's IMM model
% Date: 12/09/2015
% References:
% (a) W. Devesse et al., IJHMT 75(2015), pp. 726-735.
% (b) F. Lopez et al., ASME MSEC (2016).
%
% Revision history
% ================
%
% 09/03/2015    Convergence was verified with order p = 3.
% 10/15/2015    First deployable version uploaded to GitHub.
% 10/16/2015    Adjusted maximum width instead of laser-height width
% 12/10/2015    Included mu as an extra parameter.
%

function [T_ss,y_ss,d_T,width] = SLMsteady(P,v,A,mu,n1,plot_flag)


%% Definition of global variables
global T_m m alpha_0 alpha_T
global T dT T_0 T_max hl t_sim


%% Read thermophysical properties
T_m = 1320.0; % Melting temp. (ºC)
T_0 = 80.0; % Room temp. (ºC)
max_T = 2560.0; % Suggested max. temp. in grid (ºC), to be redefined
hl = 2.97e+5; % Latent heat of fusion (J/kg)


%% Definition of the temperature grid
dT = -(T_m-T_0)/n1; % Temp. increment in isotherms (dT < 0, ºC)
T_max = T_0-floor((T_0-max_T)/dT)*dT; % Adjusted max. temp. in grid (ºC)
T = T_max:dT:(T_0-dT); % Temperature grid (n-long array of temp.)
m = length(T)-n1+1; % Location of melting isotherm (integer): T(m) = T_m
alpha_T = ThermalDiffusivity(T); % Temp. dependent thermal diffusivity


%% Definition of initial state (initial guess)
% Determined assumed constant thermo-physical properties evaluated at T_max
alpha_0 = ThermalDiffusivity(T_max); % Thermal diffusivity at T_max
rho_ = rho(T_max); Cp_ = Cp(T_max); % Density and Cp at T_max

% Rosenthal's solution for initial guess of steady-state
c1 = (T-T_0)*2*pi*rho_*Cp_*mu*alpha_0/A/P;
c2 = 2*mu*alpha_0/v;
y_nom = c2*lambertw(1./c1/c2); % Guessed half widths:Lambert's product-log


%% Simulate to find true steady-state
S = 4.0*mu*alpha_0/v/v; % Characteristic time (s)
t_sim = 1.0e+5*S; % Simulation time: Looong time (s)
[t_array,y_array] = ode23s(@(t,x)SLM_Rate(t,x,[P,v],A,mu),[0 t_sim],y_nom);

% Compute maximum melt pool width: eq. (36) of paper (a)
L = 2.0*ThermalDiffusivity(T_m)/v; % Characteristic length (m)
C = y_array(:,m).*exp(y_array(:,m)/L);
a = C/2 + L/4*lambertw(2.0*C/L);
b = sqrt(L^2*lambertw(C./L.*exp((C-a)/L)).^2 - (C-a).^2);
width_vector = 2.0e+6*b; % Time history of max. melt pool width


%% Plotting routine
if plot_flag
    figure (1)
    plot(1.0e+3*t_array,width_vector,'r','LineWidth',2.0); grid on;
    xlabel('Time (ms)'); ylabel('Melt pool width (\mum)');
    xlim([0 1.0+3*t_sim]); ylim(width_vector(end)*[0.8 1.2]);
end


%% Output results
T_ss = T; % Array of temperature grid
y_ss = y_array(end,:); % Array of steady-state half-widths
d_T = dT; % Temp. grid granularity
width.laser = 2.0e+6*y_array(end,m); % Melt pool widths under laser
width.max = width_vector(end); % Max. melt pool width (behind laser)


end


function alpha_transition = SmoothThermalDiffusivity(T,tau)

global alpha_0
% Assume a smooth transition to actual thermal diffusivity
% 0<tau<1 : alpha(tau=0) = alpha_0, alpha(tau=1) = ThermalDiffusivity(T)
alpha_transition = alpha_0 + ...
    (ThermalDiffusivity(T)-alpha_0)*tau;

end


function alpha = ThermalDiffusivity(T)

alpha = k(T)./Cp(T)./rho(T); % Thermal diffusivity

end


function k=k(T)

global T_max

% Linear interpolation for thermal conductivity given temperature-dependent
% data
T_k_array = [21.0 38.0 93.0 204.0 316.0 427.0 538.0 649.0 760.0 ...
    871.0 982.0 1314.0 1369.0 2000.0 T_max]; % % Array of temperatures (ºC)
k_array = [9.8 10.1 10.8 12.5 14.1 15.7 17.5 19.0 20.8 22.8 25.2 ...
    31.0 25.0 27.0 27.0]; % Array of thermal conductivities (W/m-C)
k = interp1(T_k_array,k_array,T,'linear')'; % 1-D data interpolation

end


function Cp=Cp(T)

global T_max

% Linear interpolation for specific heat given temperature-dependent
% data.
T_Cp_array = [21.0 93.0 204.0 316.0 427.0 538.0 649.0 760.0 871.0 ...
    982.0 1093.0 T_max]; % % Array of temperatures (ºC)
Cp_array = [410 427 456 481 511 536 565 590 620 645 ...
    670 670]; % Array of specific heat
Cp = interp1(T_Cp_array,Cp_array,T,'linear')'; % 1-D data interpolation

end


function rho=rho(T)

global T_max

% Linear interpolation for density given temperature-dependent
% data
T_rho_array = [27.0 227.0 727.0 1314.0 1317.0 1327.0 1337.0 1347.0 ...
    1537.0 1367.0 1369.0 1727.0 2227.0 T_max]; % % Array of temperatures (ºC)
rho_array = [8620 8536 8316 8022 8014 7977 7928 7863 7777 7660 ...
    7636 7313 6860 6860]; % Array of densities (kg/m^3)
rho = interp1(T_rho_array,rho_array,T,'linear')'; % 1-D data interpolation

end


function y_dot=SLM_Rate(t,x,u,A,mu)

global hl m alpha_T
global T dT T_0 T_max T_m t_sim

%% Read inputs
y = x; % Read half-widths (states) (m)
P = u(1); % Laser power (W)
v = u(2); % Scan speed (m/s)


%% Interpolate thermal diffusivity
% Make smooth transition to 'real' values during the first fifth of
% simulation. Implemented for numerical stability.
if (t < t_sim/5.0)
    alpha = SmoothThermalDiffusivity(T, t*5.0/t_sim );
else
    alpha = alpha_T;
end


%% Construct rate equations
n = length(x); % Size of state vector
y_dot = zeros(n,1); % Vector of rate equations

% First case
% y1_dot: Inner-most point of the grid
Cp_ = Cp(T(1)); rho_ = rho(T(1)); alpha_ = mu*alpha(1);
L1 = 2*alpha_/v;
A_l = T_0-hl/(Cp_+2*alpha_*Cp_/v/y(m)+2*alpha_*alpha_*Cp_/y(m)/y(m)/v/v);
y_dot(1) = -2*alpha_*y(2)/y(1)/(y(2)-y(1))-(A*P/pi/rho_/Cp_/dT/y(1)/y(1))...
    *(1+y(1)/L1)*exp(-y(1)/L1)+(v*v/4/alpha_)*(T_max-A_l)*(y(2)-y(1))/dT;

% Second case
% yi_dot: Grid points in the liquid phase
for i=2:m-1
    Cp_ = Cp(T(i)); alpha_ = mu*alpha(i);
    A_l = T_0-hl/(Cp_+2*alpha_*Cp_/v/y(m)+2*alpha_*alpha_*Cp_/y(m)/y(m)/v/v);
    y_dot(i) = -(alpha_/y(i))*(y(i+1)/(y(i+1)-y(i))-y(i-1)/(y(i)...
        -y(i-1)))+(v*v/8/alpha_)*(T(i)-A_l)*(y(i+1)-y(i-1))/dT;
end

% Third case
% ym_dot: Melting front
Cp_ = Cp(T(m)); alpha_ = mu*alpha(m);
A_l = T_0-hl/(Cp_+2*alpha_*Cp_/v/y(m)+2*alpha_*alpha_*Cp_/y(m)/y(m)/v/v);
Ste = -Cp_*dT/hl; % Stefan number
C1 = -(A_l-2*T_m+T_0)/dT;
C2 = (1/(y(m+1)-y(m))+1/(y(m)-y(m-1)))^-1;
C3 = 1+Ste^-1;
y_dot(m) = -(alpha_/y(m)/C3)*(y(m+1)/(y(m+1)-y(m)) - y(m-1)/(y(m)-y(m-1))) ...
    + (C1*C2*v*v/4/alpha_/C3)*(1+2*(Ste^-1)*(1+(C1*C2*v/2/alpha_)^2)^-1);

% Fourth case
% yi_dot: Grid points in the solid phase
for i=m+1:n-1
    alpha_ = mu*alpha(i);
    y_dot(i) = -(alpha_/y(i))*(y(i+1)/(y(i+1)-y(i))-y(i-1)/(y(i)...
        -y(i-1)))+(v*v/8/alpha_)*(T(i)-T_0)*(y(i+1)-y(i-1))/dT;
end

% Fifth case
% yN_dot: Outer-most point of the grid
alpha_ = mu*alpha(n);
y_dot(n) = -alpha_/y(n)*((y(n)-2*y(n-1))/(y(n)-y(n-1)))-(v*v/4/alpha_)*...
    (y(n)-y(n-1));

end