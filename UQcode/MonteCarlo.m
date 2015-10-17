%% MonteCarlo for UQ (Uncertainty in inputs)
%% Define nominal conditions for simulation (std. dev.)
perct.P = 0.025; % 2.5%
perct.v = 0.015; % 1.5%
perct.A = 0.25; % 25.0%
perct.hl = 0.05; % 5.0%
perct.T_m = 0.05; % 5.0%
perct.alpha = 0.10; % 10.0%
%% Run simulations
y_m = zeros(1000,1);
for i=1:1000
    [~,~,y_m(i),~] = SLM_steady_MonteCarlo(195.0,0.8,4.0,4,false,perct);
end
save('MonteCarlo');
%% Find mean and std. dev.
widths = 2.0e+6*y_m;
m = mean(widths); S = std(widths);
[f,x]=hist(widths,20); % create histogram
figure(1)
% Plot normalized histogram
bar(x,f/trapz(x,f),'hist','b');hold on;
x = [85:.1:135];
norm = normpdf(x,m,S);
plot(x,norm,'r--','LineWidth',2.0);hold off;
xlabel('Melt pool width (\mum)'); ylabel('Probability density');
title('Normalized histogram of melt pool widths');
legend('Histogram','Normal distribution'); grid on;