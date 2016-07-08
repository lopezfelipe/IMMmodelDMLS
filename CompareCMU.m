% Code snippet CompareCMU
%
% Purpose: Compare melt pool widhts predicted for L-PBF of Alloy 625 in an
% EOS M270 machine with experimental data provided by Colt Montgomery
% (cmontogom@andrew.cmu.edu).
% There are two modes available:
% 'optimize': Find the optimal absorption coefficient (A) for the measured
%             data.
% 'predict':  Make predictions for other P-v combinations in the region of
%             validation of the model.

function CompareCMU(mode)

% Read P-v combinations used
P = [50 100 150 195]; % Laser power (W)
v = [0.2 0.4 0.6 0.8 1.0 1.2]; % Scan speed (m/s)
% Experimental error, reported by Montgomery (10.4 microns, 95% C.I.)
e = 10.4*ones(3,1); % Only 3 speeds
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
pred_widths = zeros(length(P),length(v)); % Array of predicted widths

% Find optimal absorption coefficient A
abs = 0.4:0.02:0.8;
errors = zeros(size(abs));
for k =1:length(abs)
    num_pred = 0; % Number of valid comparisons
    for i=3:length(P) % Do only 150W and 195W
        for j=3:length(v)-1 % Do only 600mm/s, 800mm/s, and 1000mm/s
            if isnan(P(i)) || isnan(v(j))
               break; 
            end 
            [T_ss,y_ss,d_T,width] = LPBFsteady(P(i),v(j),abs(k),1.0,5,false);
            pred_widths(i,j) = width.max;
            if (~isnan(widths(i,j))) % Only is data is available, compare
                num_pred = num_pred+1;
                errors(k) = errors(k) + (pred_widths(i,j)-widths(i,j))^2;
            end
        end
    end
    % Find s.d. of prediction error
    errors(k) = sqrt( errors(k)/num_pred );
end
[~,i] = min(errors); A = abs(i);
% Recompute for the chosen A
for i=3:length(P) % Do only 150W and 195W
    for j=3:length(v)-1 % Do only 600mm/s, 800mm/s, and 1000mm/s
        [T_ss,y_ss,d_T,width] = LPBFsteady(P(i),v(j),A,1.0,5,false);
        pred_widths(i,j) = width.max;
    end
end

if strcmp(mode, 'optimize')
    % Plot results
    v_mm_s = 1.0e+3*v(3:end-1);
    figure(1)
    errorbar(v_mm_s,widths(3,3:end-1),e,'rs','LineWidth',2.0); hold on;
    errorbar(v_mm_s,widths(4,3:end-1),e,'bd','LineWidth',2.0);
    plot(v_mm_s,pred_widths(3,3:end-1),'r');
    plot(v_mm_s,pred_widths(4,3:end-1),'b');
    xlabel('Scan speed (mm/s)'); ylabel('Width (\mum)');
    title(['Comparison for A = ',num2str(A)]);
    legend('150 W', '195 W'); xlim([0 1200]); ylim([0 200]); grid on;
elseif strcmp(mode, 'predict')
    % Model predictions
    pred_widths = zeros(length(P),length(v)); % Array of predicted widths
    for i=3:length(P) % Do only 150W and 195W
        for j=3:length(v)-1 % Do only 600mm/s, 800mm/s, and 1000mm/s
            [~,~,~,width] = LPBFsteady(P(i),v(j),A,1.0,5,false);
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
else
    error('Wrong option');
end

    
