% Program SLMstochastic
%
% Purpose: The function simulates the isotherm migration method (IMM)
% equations to find the  steady-state location of isotherms.
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
% (a) Temperature grid (T_ss, in ºC)
% (b) Steady-state isotherm (half-width) location (y_ss, in m)
% (c) Grid granularity (d_T, in ºC)
% (d) Widths at laser location and maximum melt pool width
%    (width.laser and width.max, both in m)
%
% Example:
% SLMstochastic(195,0.800,0.6,5,true)
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
% 10/15/2015    First deployable version uploaded to GitHub.
% 10/16/2015    Adjusted maximum width instead of laser-height width
%

function SLMstochastic(P,v,A,n1,plot_flag)
%% Definition of global variables
global T_m m alpha Ste
global T dT T_0 T_max hl
%% Read thermophysical properties
T_m = 1320.0; % Melting temperature (ºC)
T_0 = 80.0; % Room temperature (ºC)
max_T = 2560.0; % Maximum (user-defined) temperature in the grid (ºC)
hl = 2.97e+5; % Latent heat of fusion (J/kg)
%% Definition of the temperature grid
dT = -(T_m-T_0)/n1; % Temperature increment in isotherms (dT < 0, ºC)
T_max = T_0-floor((T_0-max_T)/dT)*dT; % Maximum temperature in grid (ºC)
N = round((T_0-T_max)/dT); % Number of gridpoints (integer)
n = N+1; % Number of state variables
T = T_max:dT:(T_0-dT); % Temperature grid (n-long array of temperatures)
m = length(T)-n1+1; % Location of melting isotherm (integer)
y_nom = 1.0e-6*[29.853626291529350; 30.698226486579170; 31.626762614199848;  ...
    32.657551883897540; 33.819868158807290; 35.167370363015770; ...
    37.062352155064346; 39.457359095032785; 42.857585450430754; ...
    50.162997800647120]; % Obtained with DMLSoffline
alpha = ThermalDiffusivity(T); % Thermal diffusivity array
Ste = -Cp(T_m)*dT/hl; % Stefan number
%% Simulate to find true steady-state
S = 4.0*alpha(1)/v/v; % Characteristic time (s)
dt = S/10; % Time stepping for forward Euler routine
N_sim = 399; % Number of simulations
x_array = zeros(n,N_sim+1);
t_array = 0:dt:dt*N_sim;
% Assume bulk material firt, then overhang, then bulk material
mu_array = [1.0*ones(100,1); 2.2*ones(200,1); 1.0*ones(100,1)];
x_array(:,1) = [y_nom; mu_array(1)];
for i=1:N_sim
    [~,int_array] = ode23s(@(t,x)SLM_Stochastic(t,x,[P,v],A),...
        [t_array(i) t_array(i+1)], x_array(:,i));
    x_array(:,i+1) = int_array(end,:);
    x_array(end,i+1) = mu_array(i+1);
end
% Compute maximum width
L = 2.0*ThermalDiffusivity(T_m)/v;
C = x_array(m,:).*exp(x_array(m,:)/L);
a = C/2 + L/4*lambertw(2.0*C/L);
b = sqrt(L^2*lambertw(C./L.*exp((C-a)/L)).^2 - (C-a).^2);
width_vector = 2.0e+6*b;
% Plot
if plot_flag
    figure (1)
    plot(1.0e+3*t_array,width_vector,'r','LineWidth',2.0); grid on;
    xlabel('Time (ms)'); ylabel('Melt pool width (\mum)');
    ylim([0 300]); xlim([0 1.45]);
end
%% Simulate measurements for temperatures between 500ºC and 1100ºC
% Assuming measurements are taken in line with heat source
v_meas = 26.0; % 26 microns for std. dev.
meas = 2.0e+6*x_array(m+1:m+3,:)+normrnd(0,v_meas,3,N_sim+1); % laser width = 2 * isotherm width
if plot_flag
   figure(2)
   plot(1.0e+3*t_array,meas(1,:),'r','LineWidth',1.0); hold on;
   plot(1.0e+3*t_array,meas(2,:),'b','LineWidth',1.0);
   plot(1.0e+3*t_array,meas(3,:),'k','LineWidth',1.0); grid on;
   xlabel('Time (ms)'); ylabel('Width (\mum)');
   legend([num2str(T(m+1)),'ºC'],[num2str(T(m+2)),'ºC'],[num2str(T(m+3)),'ºC']);
end
% Compute matrices
x_nom = x_array(:,end); u_nom = [P; v];
[A_m,B_m,G_m,C_m] = ComputeEvolutionMatrices(x_nom,u_nom,A);
conSys = ss([A_m G_m; zeros(1,n)], [B_m; 0 0], [C_m zeros(3,1)], 0);
disSys = c2d(conSys,dt,'zoh');
Phi = disSys.a; Lambda = disSys.b; H = disSys.c;
y_nom = H*x_nom;
Q = zeros(n,n); Q(n,n) = dt*(1.0)^2; % 100% variation in a time step
R = v_meas*eye(3);
dxest_array = zeros(n,N_sim+1); Pest_array = zeros(n,n,N_sim+1);
Pest_array(:,:,1) = diag([R(1,1)*ones(1,N)  Q(n,n)]);
for i =1:N_sim
    % Prediction
    dx_pred = Phi*dxest_array(:,i) + Lambda*([P; v]-u_nom);
    P_pred = Phi*Pest_array(:,:,i)*Phi' + Q;
    % Update
    innov = meas(:,i) - (y_nom+H*dx_pred);
    S = H*P_pred*H' + R;
    K = P_pred*H'/S;
    dxest_array(:,i+1) = dx_pred + K*innov;
    Pest_array(:,:,i+1) = P_pred - K*H*P_pred;
end
xest_array = dxest_array + x_nom*ones(1,N_sim+1);
%% Plot melt pool width estimates
% Estimated mu
mu_est = xest_array(end,:)';
% Lower and upper bound, remove singleton dimensions
mu_min = xest_array(end,:)'-2.0*sqrt(squeeze(sqrt(Pest_array(end,end,:))));
mu_max = xest_array(end,:)'+2.0*sqrt(squeeze(sqrt(Pest_array(end,end,:))));
% Estimated melt pool
ym_est = xest_array(m,:)';
sigma_ym = squeeze(sqrt(Pest_array(m,m,:)));
sigma_ym(1:5) = 0; % Ignore the first five elements
% Lower and upper bound, remove singleton dimensions
ym_min = xest_array(m,:)'-2.0*sigma_ym;
ym_max = xest_array(m,:)'+2.0*sigma_ym;
% Compute maximum width
C_est = ym_est.*exp(ym_est/L);
a_est = C_est/2 + L/4*lambertw(2.0*C_est/L);
b_est = sqrt(L^2*lambertw(C_est./L.*exp((C_est-a_est)/L)).^2 - (C_est-a_est).^2);
width_est = 2.0e+6*b_est;
C_min = ym_min.*exp(ym_min/L);
a_min = C_min/2 + L/4*lambertw(2.0*C_min/L);
b_min = sqrt(L^2*lambertw(C_min./L.*exp((C_min-a_min)/L)).^2 - (C_min-a_min).^2);
width_min = 2.0e+6*b_min;
C_max = ym_max.*exp(ym_max/L);
a_max = C_max/2 + L/4*lambertw(2.0*C_max/L);
b_max = sqrt(L^2*lambertw(C_max./L.*exp((C_max-a_max)/L)).^2 - (C_max-a_max).^2);
width_max = 2.0e+6*b_max;
figure(3)
subplot(2,1,1)
plot(1.0e+3*t_array,mu_est,'b','LineWidth',2.0); hold on;
plot(1.0e+3*t_array,mu_min,'r--','LineWidth',1.0);
plot(1.0e+3*t_array,mu_max,'r--','LineWidth',1.0); grid on;
xlabel('Time (ms)'); ylabel('Diffusion efficiency (\mu)'); xlim([0 1.45]);
subplot(2,1,2)
plot(1.0e+3*t_array,width_est,'b','LineWidth',2.0); hold on;
plot(1.0e+3*t_array,width_min,'r--','LineWidth',1.0);
plot(1.0e+3*t_array,width_max,'r--','LineWidth',1.0); grid on;
xlabel('Time (ms)'); ylabel('Melt pool width (\mum)'); xlim([0 1.45]);
ylim([0 300]);
end

function y_dot = SLM_Stochastic(t,x,u,A)
global alpha hl m Ste T dT T_0 T_max T_m
%% Read inputs
n  = length(x); % Number of states
N = n-1; % Number of isotherms
y = x(1:N); % Isotherms in temperature grid (m)
mu = x(n); % Diffusion efficiency
P = u(1); % Laser power (W)
v = u(2); % Scan speed (m/s)
Ste_ = 1/Ste;
%% Construct rate equations
y_dot = zeros(n,1);
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
Tm_bar = ( A_l-2*T_m+T_0 )/dT;
g1 = -( alpha_/y(m)/(1+Ste_) ) * ( y(m+1)/(y(m+1)-y(m)) ...
    - y(m-1)/(y(m)-y(m-1)) );
g2 = 1/( y(m+1)-y(m) ) + 1/( y(m)-y(m-1) );
g3 = 2 * Ste_ *( 1 + ( Tm_bar*v / ( 2*alpha(m)*g2 ) )^2 )^-1;
y_dot(m) = g1 - v^2/(4*alpha_*(1+Ste_))*Tm_bar*g2^-1*(1+g3);
% Fourth case
% yi_dot: Grid points in the solid phase
for i=m+1:N-1
    alpha_ = mu*alpha(i);
    y_dot(i) = -(alpha_/y(i))*(y(i+1)/(y(i+1)-y(i))-y(i-1)/(y(i)...
        -y(i-1)))+(v*v/8/alpha_)*(T(i)-T_0)*(y(i+1)-y(i-1))/dT;
end
% Fifth case
% yN_dot: Outer-most point of the grid
alpha_ = mu*alpha(N);
y_dot(N) = -alpha_/y(N)*((y(N)-2*y(N-1))/(y(N)-y(N-1)))-(v*v/4/alpha_)*...
    (y(N)-y(N-1));
% Diffusion efficiency
% Assume no variation
y_dot(n) = 0;
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