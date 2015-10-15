% Code snippet CompareCMU
%
% Purpose: Compare melt pool widhts predicted for DMLS of Alloy 625 in an
% EOS M270 machine with experimental data provided by Colt Montgomery
% (cmontogom@andrew.cmu.edu).

% Compare steady state to CMU's width data
% Array of power used
abs = 4.0;
P = [50 100 150 195]; % Laser power (W)
v = [0.2 0.4 0.6 0.8 1.0 1.2]; % Scan speed (m/s)
e = 10.4*ones(size(v)); % Error
widths = zeros(length(P),length(v)); % Array of widths
widths(1,1) = 115; widths(1,2) = 95; widths(1,3) = 85;
widths(1,4) = NaN; widths(1,5) = 70; widths(1,6) = 60;
widths(2,1) = 160; widths(2,2) = NaN; widths(2,3) = 110;
widths(2,4) = 100; widths(2,5) = 100; widths(2,6) = 90;
widths(3,1) = 250; widths(3,2) = 140; widths(3,3) = 120;
widths(3,4) = 105; widths(3,5) = 110; widths(3,6) = 100;
widths(4,1) = 285; widths(4,2) = 215; widths(4,3) = 150;
widths(4,4) = 125; widths(4,5) = NaN; widths(4,6) = 115;
pred_widths = zeros(length(P)-1,length(v)-2); % Array of predicted widths
tic;
for i=2:length(P)
    for j=2:length(v)-1
        [T_ss,y_ss,y_m] = DMLSoffline(P(i),v(j),abs,4,false);
        pred_widths(i-1,j-1) = 2.0e+6*y_m;
    end
end
UQ_error = 21.4; % From UQ
toc;
figure(1)
errorbar(1.0e+3*v,widths(2,:),e,'g^','LineWidth',2.0); hold on;
errorbar(1.0e+3*v,widths(3,:),e,'rs','LineWidth',2.0);
errorbar(1.0e+3*v,widths(4,:),e,'bd','LineWidth',2.0);
plot(1.0e+3*v(2:end-1),pred_widths(1,:),'g');
plot(1.0e+3*v(2:end-1),pred_widths(2,:),'r');
plot(1.0e+3*v(2:end-1),pred_widths(3,:),'b');
xlabel('Power (W)'); ylabel('Width (\mum)');
legend('100 W', '150 W', '195 W'); xlim([0 1400]); ylim([0 300]); grid on;