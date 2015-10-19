% Code snippet CompareCMU_Pv
%
% Purpose: Compare melt pool widhts predicted for DMLS of Alloy 625 in an
% EOS M270 machine with experimental data provided by Colt Montgomery
% (cmontogom@andrew.cmu.edu). Comparisons are shown for different laser
% powers at 800 mm/s and for different speeds at 195 W.

% Simulation parameters
abs = 0.60;
P = [50 100 150 195]; % Laser power (W)
v = [0.2 0.4 0.6 0.8 1.0 1.2]; % Scan speed (m/s)

% Experimental error, reported by Montgomery (10.4 microns, 95% C.I.)
e_v = 10.4*ones(2,1); % Only two powers (constant speed)
e_P = 10.4*ones(3,1); % Only three speeds (constant power)
% Array of widths
widths = zeros(length(P),length(v));
widths(1,1) = 115; widths(1,2) = 95; widths(1,3) = 85;
widths(1,4) = NaN; widths(1,5) = 70; widths(1,6) = 60;
widths(2,1) = 160; widths(2,2) = NaN; widths(2,3) = 110;
widths(2,4) = 100; widths(2,5) = 100; widths(2,6) = 90;
widths(3,1) = 250; widths(3,2) = 140; widths(3,3) = 120;
widths(3,4) = 105; widths(3,5) = 110; widths(3,6) = 100;
widths(4,1) = 285; widths(4,2) = 215; widths(4,3) = 150;
widths(4,4) = 125; widths(4,5) = NaN; widths(4,6) = 115;

% Model predictions
pred_widths = zeros(length(P),length(v)); % Array of predicted widths
for i=3:length(P) % Do only 150W and 195W
    for j=3:length(v)-1 % Do only 600mm/s, 800mm/s, and 1000mm/s
        [T_ss,y_ss,d_T,width] = DMLSoffline(P(i),v(j),abs,5,false);
        pred_widths(i,j) = width.max;
    end
end

UQ_error = 38.8; % 95% confidence interval in model predictions
% Constant speed
figure(1)
subplot(2,1,1)
errorbar(P(3:4),widths(3:4,4),e_v,'bd','LineWidth',2.0); hold on;
plot(P(3:4),pred_widths(3:4,4),'b');
plot(P(3:4),pred_widths(3:4,4)-(pred_widths(4,4)-widths(4,4)+UQ_error)*ones(2,1),'b--');
plot(P(3:4),pred_widths(3:4,4)-(pred_widths(4,4)-widths(4,4)-UQ_error)*ones(2,1),'b--');
xlabel('Power (W)'); ylabel('Width (\mum)');
title('Melt pool widths at 800mm/s');
xlim([80 200]); ylim([0 300]); grid on;
% Constant power
subplot(2,1,2)
errorbar(1.0e+3*v(3:5),widths(4,3:5),e_P,'bd','LineWidth',2.0); hold on;
plot(1.0e+3*v(3:5),pred_widths(4,3:5),'b');
plot(1.0e+3*v(3:5),pred_widths(4,3:5)-(pred_widths(4,4)-widths(4,4)+UQ_error)*ones(1,3),'b--');
plot(1.0e+3*v(3:5),pred_widths(4,3:5)-(pred_widths(4,4)-widths(4,4)-UQ_error)*ones(1,3),'b--');
xlabel('Scan speed (mm/s)'); ylabel('Width (\mum)');
title('Melt pool widths at 195W');
xlim([0 1200]); ylim([0 300]); grid on;