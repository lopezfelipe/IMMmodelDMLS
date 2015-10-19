%% MonteCarlo for UQ (Uncertainty in inputs)
%% Define nominal conditions for simulation (std. dev.)
perct.P = 0.025; % 2.5%
perct.v = 0.015; % 1.5%
perct.A = 0.25; % 25.0%
perct.hl = 0.05; % 5.0%
perct.T_m = 0.05; % 5.0%
perct.alpha = 0.10; % 10.0%
%% Run simulations
for i=1:4000
    [~,~,~,width(i)] = SLM_steady_MonteCarlo(195.0,0.8,0.6,5,false,perct);
end
save('MonteCarlo');
m = mean([width.max]); S = std([width.max]);
[f,x]=hist([width.max],20); % create histogram
figure(1)
% Plot normalized histogram
bar(x,f/trapz(x,f),'hist','b');hold on;
x = [60:.1:190];
norm = normpdf(x,m,S);
plot(x,norm,'r--','LineWidth',2.0);hold off;
xlabel('Melt pool width (\mum)'); ylabel('Probability density');
title('Normalized histogram of melt pool widths');
legend('Histogram','Normal distribution'); grid on;
% Returned m = 131.5729 microns, S = 18.6564 microns for N = 4000