% Prediction (before measurement is gathered) of maximum likelihood and
% covariance using Extended Kalman Filter (EKF)
function [x_bar,P_bar] = Prediction(x,P,u)
global Cp rho alpha T A T0 Tm hl Ste_
% Read inputs
n = length(x); % Number of states
N = n-1; % Number of isotherms in temperaturre grid
y = x(1:N); % Isotherms come from state vector (m)
mu = x(n); % Diffusion efficiency
P = u(1); % Laser power (W)
v = u(2); % Scan speed (m/s)
% Initialize prediction
x_bar = zeros(length(x),1);
P_bar = zeros(length(x));
%% Population of Jacobian matrices
%% Jacobian matrix
%% State matrix: A_m = d(x_dot)/dx
A_m = zeros(n); % Initialize matrix with zeros
Al = zeros(n,1); % Equivalent ambient temperature in the liquid phase
dAl_dym = zeros(n,1); dAl_dv = zeros(n,1);
for i=1:n
    Al(i) = T0 - hl/( Cp+2*alpha(i)*Cp/v/y(m) + 2*alpha(i)^2*Cp/(v*y(m))^2 );
    dAl_dv(i) = - ( 2*hl*v*y(m)^2/Cp ) * alpha(i) * ( 2*alpha(i)+v*y(m) ) ...
        /( 2*alpha(i)^2+2*alpha(i)*v*y(m)+v^2*y(m)^2 )^2;
    dAl_dym(i) = - ( 2*hl*v^2*y(m)/Cp ) * alpha(i) * ( 2*alpha(i)+v*y(m) ) ...
        /( 2*alpha(i)^2+2*alpha(i)*v*y(m)+v^2*y(m)^2 )^2;
end
g1 = -( alpha(m)/y(m)/(1+Ste_) ) * ( y(m+1)/(y(m+1)-y(m)) ...
    - y(m-1)/(y(m)-y(m-1)) );
Tm_bar = ( Al(m)-2*Tm+T0 )/dT;
g2 = 1/( y(m+1)-y(m) ) + 1/( y(m)-y(m-1) );
g3 = 2 * Ste_ *( 1 + ( Tm_bar*v / ( 2*dT*alpha(m)*g2 ) )^2 )^-1;
dg1_dymm = alpha(m) / ( 1+Ste_ ) / ( y(m)-y(m-1) )^2;
dg1_dym = ( 1/( 1+Ste_ ) ) * ( (alpha(m)/y(m)^2) * ( y(m+1)/(y(m+1)-y(m)) ...
    - y(m-1)/(y(m)-y(m-1)) ) - (alpha(m)/y(m)) * ( y(m+1)/(y(m+1)-y(m))^2 ...
    + y(m-1)/(y(m)-y(m-1))^2 ) );
dg1_dyp = alpha(m)/( 1+Ste_ )/( y(m+1)-y(m) )^2;
dTm_bar_dv = dAl_dv(m)/dT;
dTm_bar_dym = dAl_dv(m)/dT;
dg2_dymm = 1/(y(m)-y(m-1))^2;
dg2_dym = 1/(y(m+1)-y(m))^2 - 1/(y(m)-y(m-1))^2;
dg2_dymp = -1/(y(m+1)-y(m))^2;
dg3_dv = -Ste_*Tm_bar*v / (dT*alpha(m)*g2)^2 * ( dTm_bar_dv*v + Tm_bar ) ...
    *( 1 + ( Tm_bar*v / (2*dT*alpha(m)*g2) )^2 )^-2;
dg3_dymm = 



% First case
% y1_dot: Inner-most point of the grid
A_m(1,1) = -2.0*alpha(1)*y(2)*( 2*y(1)-y(2) ) / ( y(1) * (y(2)-y(1)) )^2 ...
    + A*P*exp( -v*y(1)/2/alpha(1) ) / ( 4*pi*rho*Cp*dT*y(1)^3*alpha(1)^2 ) ...
    *( 8*alpha(1)^2+4*alpha(1)*v*y(1)+v^2*y(1)^2 ) ...
    - v^2*( T1-Al(1) )/( 4*alpha(1)*dT );
A_m(1,2) = 2*alpha(1) / ( y(2)-y(1) )^2 + v^2*( T1-Al(1) )/( 4*alpha(1)*dT );
A_m(1,m) = -v^2*( y(2)-y(1) )*dAl_dym(1) / ( 4*alpha(1)*dT );
% Second case
% yi_dot: Grid points in the liquid phase
for i=2:m-1
    A_m(i,i-1) = alpha(i) / ( y(i)-y(i-1) )^2 - v^2*( T(i)-Al(i) )/( 8*alpha(i)*dT );
    A_m(i,i) = ( alpha(i)/y(i)^2 )*( y(i+1)/(y(i+1)-y(i)) - y(i-1)/(y(i)-y(i-1)) ) ...
        -(  alpha(i)/y(i) )*( y(i+1)/(y(i+1)-y(i))^2 + y(i-1)/(y(i)-y(i-1))^2 );
    A_m(i,i+1) = alpha(i) / ( y(i+1)-y(i) )^2 + v^2*( T(i)-Al(i) )/( 8*alpha(i)*dT );
    A_m(i,m) = A_m(i,m) - v^2*( y(i+1)-y(i-1) ) * dAl_dym(i) / ( 8*alpha(i)*dT );        
end
% Third case
% ym_dot: Melting front
A_m(m,m-1) = dg1_dymm + 




end
