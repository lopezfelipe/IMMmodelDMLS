% Code snippet CompareCMU_v
%
% Purpose: Compare melt pool widhts predicted for DMLS of Alloy 625 in an
% EOS M270 machine with experimental data provided by Colt Montgomery
% (cmontogom@andrew.cmu.edu). Comparisons are shown for different laser
% powers at 800 mm/s.

% Simulation parameters
abs = 4.0;
P = [50 100 150 195]; % Laser power (W)
v = [0.2 0.4 0.6 0.8 1.0 1.2]; % Scan speed (m/s)

% Experimental data
e = 10.4*ones(size(P)); % 95% confidence interval for measurement error
widths = zeros(length(P),length(v)); % Array of widths
widths(1,1) = 115; widths(1,2) = 95; widths(1,3) = 85;
widths(1,4) = NaN; widths(1,5) = 70; widths(1,6) = 60;
widths(2,1) = 160; widths(2,2) = NaN; widths(2,3) = 110;
widths(2,4) = 100; widths(2,5) = 100; widths(2,6) = 90;
widths(3,1) = 250; widths(3,2) = 140; widths(3,3) = 120;
widths(3,4) = 105; widths(3,5) = 110; widths(3,6) = 100;
widths(4,1) = 285; widths(4,2) = 215; widths(4,3) = 150;
widths(4,4) = 125; widths(4,5) = NaN; widths(4,6) = 115;

% Model predictions
pred_widths = zeros(length(P)-1,length(v)-2); % Array of predicted widths
for i=2:length(P)
    for j=2:length(v)-1
        [T_ss,y_ss,y_m] = DMLSoffline(P(i),v(j),abs,4,false);
        pred_widths(i-1,j-1) = 2.0e+6*y_m;
    end
end
UQ_error = 21.4; % 95% confidence interval in model predictions
figure(1)
errorbar(P,widths(:,4),e,'bd','LineWidth',2.0); hold on;
plot(P(2:end),pred_widths(:,3),'b');
plot(P(2:end),pred_widths(:,3)-(pred_widths(3,4)-widths(3,4)+UQ_error)*ones(3,1),'b');
plot(P(2:end),pred_widths(:,3)-(pred_widths(3,4)-widths(3,4)-UQ_error)*ones(3,1),'b');
xlabel('Power (W)'); ylabel('Width (\mum)');
xlim([80 200]); ylim([0 300]); grid on;