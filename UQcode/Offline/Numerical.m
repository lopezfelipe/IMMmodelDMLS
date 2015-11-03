function GCI = Numerical(plot_convergence)
if nargin < 1
   plot_convergence = false; 
end
% Select three grids
% h3 = 620ºC (n1=2), h2 = 310ºC (n1=4), and 248ºC (n1=5)
n1 = [5, 4, 2];
for i=1:length(n1)
   [~,~,d_T(i),width(i)] = SLMsteady(195,0.800,0.6,1.0,n1(i),false) ;
end
phi = [width.max]; % Predictions
epsilon_32 = phi(3)-phi(2); r_32 = d_T(3)/d_T(2);
epsilon_21 = phi(2)-phi(1); r_21 = d_T(2)/d_T(1);
s = sign(epsilon_32/epsilon_21);
p = 1.0:0.1:3.0;
for j=1:length(p)
   q = log((r_21^p(j)-s)/(r_32^p(j)-s));
   p_c(j) = (1/log(r_21))*(log(epsilon_32/epsilon_21)+q);
end
% Find minimum
[~,i]=min(abs(p-p_c));
% Order of accuracy is
p = p(i);
% Now, calculate extrapolated values
%phi_ext = (r_21^p*phi(1))/abs(r_21^p-1);
e21_a = abs(phi(1)-phi(2));
GCI = 1.25*e21_a/(r_21^p-1); % Absolute value of 95% CI
% Computation resulted in +/- 2.74 microns
if (plot_convergence)
    % Plot convergence plot
    n1_c = [8, 5, 4, 2];
    for i=1:length(n1_c)
       [~,~,d_T(i),width(i)] = DMLSoffline(195,0.800,0.6,n1_c(i),false) ;
    end
    phi = [width.max]; % Predictions
    n1_max = 40;
    [~,~,~,width_max] = DMLSoffline(195,0.800,0.6,n1_max,false) ;
    error = [width.max]-width_max.max;
    loglog(abs(d_T),error,'rs','LineWidth',2.0); grid on;
    xlabel('Mesh size, \DeltaT (ºC)'); ylabel('Error, Solution-Exact (\mum)');
    xlim([100 1000]); ylim([1.0e-1 2.0e+1]); hold on;
end
end