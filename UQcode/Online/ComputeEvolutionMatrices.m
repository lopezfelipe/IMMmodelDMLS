% Prediction (before measurement is gathered) of maximum likelihood and
% covariance using Extended Kalman Filter (EKF)
%function [x_bar,P_bar] = Prediction(x_old,P_old,u)
function [A_m,B_m,G_m,C_m] = ComputeEvolutionMatrices(x,u,A)
global T T_0 T_m hl Ste alpha m dT
% Read inputs
n = length(x); % Number of states
N = n-1; % Number of isotherms in temperaturre grid
y = x(1:N); % Isotherms come from state vector (m)
mu = x(n); % Diffusion efficiency
P = u(1); % Laser power (W)
v = u(2); % Scan speed (m/s)
Ste_ = 1/Ste; % Inverse of Stefan number
%% Population of Jacobian matrices
%% Jacobian matrix
%% State matrix: A_m = d(x_dot)/dx
A_m = zeros(N); % Initialize matrix with zeros
Al = zeros(N,1); % Apparent ambient temperature in the liquid phase
dAl_dym = zeros(N,1); dAl_dv = zeros(N,1); dAl_dmu = zeros(N,1);
for i=1:N
    Cp_ = Cp(T(i));
    Al(i) = T_0 - hl/( Cp_ + 2*mu*alpha(i)*Cp_/v/y(m) + 2*mu^2*alpha(i)^2*Cp_/(v*y(m))^2 );
    dAl_dv(i) = - ( 2*hl*v*y(m)^2/Cp_ ) * mu*alpha(i) * ( 2*mu*alpha(i)+v*y(m) ) ...
        /( 2*mu^2*alpha(i)^2+2*mu*alpha(i)*v*y(m)+v^2*y(m)^2 )^2;
    dAl_dym(i) = - ( 2*hl*v^2*y(m)/Cp_ ) * mu*alpha(i) * ( 2*mu*alpha(i)+v*y(m) ) ...
        /( 2*mu^2*alpha(i)^2+2*mu*alpha(i)*v*y(m)+v^2*y(m)^2 )^2;
    dAl_dmu(i) = ( 2*hl*v^2*y(m)^2/Cp_ ) * alpha(i) * ( 2*mu*alpha(i)+v*y(m) ) ...
        /( 2*mu^2*alpha(i)^2+2*mu*alpha(i)*v*y(m)+v^2*y(m)^2 )^2;
end
g1 = -( mu*alpha(m)/y(m)/(1+Ste_) ) * ( y(m+1)/(y(m+1)-y(m)) ...
    - y(m-1)/(y(m)-y(m-1)) );
Tm_bar = ( Al(m)-2*T_m+T_0 )/dT;
g2 = 1/( y(m+1)-y(m) ) + 1/( y(m)-y(m-1) );
g3 = 2 * Ste_ *( 1 + ( Tm_bar*v / ( 2*mu*alpha(m)*g2 ) )^2 )^-1;
dg1_dymm = mu*alpha(m) / ( 1+Ste_ ) / ( y(m)-y(m-1) )^2;
dg1_dym = ( mu*alpha(m)/( 1+Ste_ ) ) * ( (1/y(m)^2) * ( y(m+1)/(y(m+1)-y(m)) ...
    - y(m-1)/(y(m)-y(m-1)) ) - (1/y(m)) * ( y(m+1)/(y(m+1)-y(m))^2 ...
    + y(m-1)/(y(m)-y(m-1))^2 ) );
dg1_dymp = mu*alpha(m)/( 1+Ste_ )/( y(m+1)-y(m) )^2;
dg1_dmu = -( alpha(m)/y(m)/(1+Ste_) ) * ( y(m+1)/(y(m+1)-y(m)) ...
    - y(m-1)/(y(m)-y(m-1)) );
dTm_bar_dv = dAl_dv(m)/dT;
dTm_bar_dym = dAl_dym(m)/dT;
dTm_bar_dmu = dAl_dmu(m)/dT;
dg2_dymm = 1/(y(m)-y(m-1))^2;
dg2_dym = 1/(y(m+1)-y(m))^2 - 1/(y(m)-y(m-1))^2;
dg2_dymp = -1/(y(m+1)-y(m))^2; % Checked until here
den_g3 = ( 1.0 + (Tm_bar*v/(2*mu*alpha(m)*g2))^2 )^-2;
dg3_dv = -den_g3 * (Tm_bar*v/mu^2/alpha(m)^2) * (Ste_/g2^2) * (v*dTm_bar_dv+Tm_bar);
dg3_dmu = den_g3 * (Tm_bar*v^2/mu^2/alpha(m)^2) * (Ste_/g2^2) * (Tm_bar/mu-dTm_bar_dmu);
dg3_dymm = den_g3 * (Tm_bar*v/mu/alpha(m))^2 * (Ste_/g2^3) * dg2_dymm;
dg3_dym = -den_g3 * (Tm_bar*v^2/mu^2/alpha(m)^2) * (Ste_/g2^2) * (v*dTm_bar_dym ...
    -(Tm_bar/g2)*dg2_dym);
dg3_dymp = den_g3 * (Tm_bar*v/mu/alpha(m))^2 * (Ste_/g2^3) * dg2_dymp;
% First case
% y1_dot: Inner-most point of the grid
Cp_ = Cp(T(1)); rho_ = rho(T(1));
A_m(1,1) = -2.0*mu*alpha(1)*y(2)*( 2*y(1)-y(2) ) / ( y(1) * (y(2)-y(1)) )^2 ...
    + A*P*exp( -v*y(1)/2/mu/alpha(1) ) / ( 4*pi*rho_*Cp_*dT*y(1)^3*mu^2*alpha(1)^2 ) ...
    *( 8*mu^2*alpha(1)^2+4*mu*alpha(1)*v*y(1)+v^2*y(1)^2 ) ...
    - v^2*( T(1)-Al(1) )/( 4*mu*alpha(1)*dT );
A_m(1,2) = 2*mu*alpha(1) / ( y(2)-y(1) )^2 + v^2*( T(1)-Al(1) )/( 4*mu*alpha(1)*dT );
A_m(1,m) = -v^2*( y(2)-y(1) )*dAl_dym(1) / ( 4*mu*alpha(1)*dT );
% Second case
% yi_dot: Grid points in the liquid phase
for i=2:m-1
    A_m(i,i-1) = mu*alpha(i) / ( y(i)-y(i-1) )^2 - v^2*( T(i)-Al(i) )/( 8*mu*alpha(i)*dT );
    A_m(i,i) = ( mu*alpha(i)/y(i)^2 )*( y(i+1)/(y(i+1)-y(i)) - y(i-1)/(y(i)-y(i-1)) ) ...
        -(  mu*alpha(i)/y(i) )*( y(i+1)/(y(i+1)-y(i))^2 + y(i-1)/(y(i)-y(i-1))^2 );
    A_m(i,i+1) = mu*alpha(i) / ( y(i+1)-y(i) )^2 + v^2*( T(i)-Al(i) )/( 8*mu*alpha(i)*dT );
    A_m(i,m) = A_m(i,m) - v^2*( y(i+1)-y(i-1) ) * dAl_dym(i) / ( 8*mu*alpha(i)*dT );        
end
% Third case
% ym_dot: Melting front
A_m(m,m-1) = dg1_dymm - ( (v^2*Tm_bar) / (4*mu*alpha(m)*g2*(1+Ste_)) ) ...
    *( dg3_dymm - (1+g3)*dg2_dymm/g2 );
A_m(m,m) = dg1_dym - ( (v^2) / (4*mu*alpha(m)*(1+Ste_)*g2) ) * ( (1+g3) ...
    * dTm_bar_dym - (1+g3)*Tm_bar/g2 * dg2_dym + Tm_bar*dg3_dym );
A_m(m,m+1) = dg1_dymp - ( (v^2*Tm_bar) / (4*mu*alpha(m)*g2*(1+Ste_)) ) ...
    *( dg3_dymp - (1+g3)*dg2_dymp/g2 );
% Fourth case
% yi_dot: Grid points in the solid phase
for i=m+1:N-1
    A_m(i,i-1) = mu*alpha(i)/(y(i)-y(i-1))^2 - v^2*( T(i)-T_0 )/dT/( 8*mu^2*alpha(i) );
    A_m(i,i) = (mu*alpha(i)/y(i)^2) * ( y(i+1)/(y(i+1)-y(i)) ...
        - y(i-1)/(y(i)-y(i-1)) ) - (mu*alpha(i)/y(i)) ...
        * ( y(i+1)/(y(i+1)-y(i))^2 + y(i-1)/(y(i)-y(i-1))^2 );
    A_m(i,i+1) = mu*alpha(i)/(y(i+1)-y(i))^2 + v^2*( T(i)-T_0 )/dT/( 8*mu^2*alpha(i) );
end
% Fifth case
% yN_dot: Outer-most point of the grid
A_m(N,N-1) = mu*alpha(N)/(y(N)-y(N-1))^2 + v^2/(4*mu*alpha(N));
A_m(N,N) = -mu*alpha(N)/(y(N)-y(N-1))/y(N) * ( y(N-1)/(y(N)-y(N-1)) ...
    - (y(N)-2*y(N-1))/y(N) ) - v^2/(4*mu*alpha(N));
%% Input matrix
%% Input matrix: B_m = d(x_dot)/du
B_m = zeros(N,2); % Initialize matrix with zeros
% First case
% y1_dot: Inner-most point of the grid
Cp_ = Cp(T(1)); rho_ = rho(T(1));
B_m(1,1) = -A*exp( -v*y(1)/(2*mu*alpha(1)) )*( 1 ...
    + v*y(1)/(2*mu*alpha(1)) )/(pi*rho_*Cp_*dT*y(1)^2);
B_m(1,2) = A*P*v*exp( -v*y(1)/(2*mu*alpha(1)) )/(4*pi*rho_*Cp_*dT*mu^2*alpha(1)^2) ...
    + v*(y(2)-y(1))*( 2*(T(1)-Al(1)) - v*dAl_dv(1) )/(4*mu*alpha(1)*dT); % Checked up to here
% Second case
% yi_dot: Grid points in the liquid phase
for i=2:m-1
    B_m(i,2) = (v/(4*mu*alpha(i)*dT)) * ( y(i+1)-y(i-1) ) ...
        * ( T(i)-Al(i)-(v/2)*dAl_dv(i) );    
end
% Third case
% ym_dot: Melting front
B_m(m,2) = -v/( 4*mu*alpha(m)*(1+Ste_)*g2 ) * ( (1+g3)*(2*Tm_bar+v*dTm_bar_dv) ...
    + v*Tm_bar*dg3_dv );
% Fourth case
% yi_dot: Grid points in the solid phase
for i=m+1:N-1
    B_m(i,2) = v*(T(i)-T_0)*(y(i+1)-y(i-1))/(4*mu*alpha(i)*dT);
end
% Fifth case
% yN_dot: Outer-most point of the grid
B_m(N,2) = -v*(y(N)-y(N-1))/(2*mu*alpha(N));
%% Perturbation matrix
%% Perturbation matrix: G_m = d(x_dot)/dmu
G_m = zeros(N,1); % Initialize matrix with zeros
% First case
% y1_dot: Inner-most point of the grid
Cp_ = Cp(T(1)); rho_ = rho(T(1));
G_m(1,1) = -2*alpha(1)*y(2)/y(1)/(y(2)-y(1)) ...
    - A*P*v^2*exp( -v*y(1)/(2*mu*alpha(1)) )/( 4*pi*rho_*Cp_*dT*mu^3*alpha(1)^2 ) ...
    - ( v^2*(y(2)-y(1)) ) * ( mu*dAl_dmu(1)+T(1)-Al(1) ) / ( 4*mu^2*dT*alpha(1) );
% Second case
% yi_dot: Grid points in the liquid phase
for i=2:m-1
    G_m(i,1) = -(alpha(i)/y(i))*( y(i+1)/(y(i+1)-y(i)) - y(i-1)/(y(i)-y(i-1)) ) ...
        -v^2*(y(i+1)-y(i-1)) * ( mu*dAl_dmu(i)+T(i)-Al(i) ) / ( 8*mu^2*dT*alpha(i) );
end
% Third case
% ym_dot: Melting front
G_m(m,1) = dg1_dmu - v^2/(4*mu*alpha(m)*(1+Ste_)*g2)*( (dTm_bar_dmu ...
    -Tm_bar/mu)*(1+g3) + Tm_bar*dg3_dmu );
% Fourth case
% yi_dot: Grid points in the solid phase
for i=m+1:N-1
    G_m(i,1) = -(alpha(i)/y(i))*( y(i+1)/(y(i+1)-y(i)) - y(i-1)/(y(i)-y(i-1)) ) ...
        - v^2*( T(i)-T_0 )*(y(i+1)-y(i-1))/(8*mu^2*alpha(i)*dT);
end
% Fifth case
% yN_dot: Outer-most point of the grid
G_m(N,1) = -(alpha(N)/y(N))*( (y(N)-2*y(N-1))/(y(N)-y(N-1)) )...
    + v^2*(y(N)-y(N-1))/(4*mu^2*alpha(N));
%% Observation matrix
%% Observation matrix: C_m = d(z)/dx
C_m = zeros(3,N); % Initialize matrix with zeros
C_m(1,m+1) = 2.0e+6;
C_m(2,m+2) = 2.0e+6;
C_m(3,m+3) = 2.0e+6;
end

%% Material properties
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
