%% Calculating stats for each recon

load('JAGS_out.mat')
load('CPS_results.mat')
data_FJ(:,1) = CPS_FJ(:,3); data_FJ(:,2) = recon_FJ(:,3);
data_M(:,1) = CPS_M(:,3); data_M(:,2) = recon_M(:,3);
data_V(:,1) = CPS_V(:,3); data_V(:,2) = recon_V(:,3);

windowSize = 101;
[correlationTS(:,1)] = movingCorrelation(data_FJ, windowSize, 1);
[correlationTS(:,2)] = movingCorrelation(data_M, windowSize, 1);
[correlationTS(:,3)] = movingCorrelation(data_V, windowSize, 1);

figure(1)
subplot(3,1,1)
axis([1000 2000 0.2 1])
plot(recon_FJ(:,1),correlationTS(:,1),'m','linewidth',2)
subplot(3,1,2)
hold on
plot(recon_M(:,1),correlationTS(:,2),'r','linewidth',2)
subplot(3,1,3)
axis([1000 2000 0.2 1])
hold on
plot(recon_V(:,1),correlationTS(:,3),'b','linewidth',2)

plot2svg('recon_runningcorr.svg')