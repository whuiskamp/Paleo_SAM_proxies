%% Plotting the CPS reconstructions with their Bootstrapped uncertainty
% ranges

clear
load('prox_groups/bootstrap_proxies.mat',...
    'boot_FJ_rng','boot_M_rng','boot_V_rng')
load CPS_results.mat
load('JAGS_out.mat','AbramSAM')

% Smooth recons and standardise
zCPS_FJ = zscore(CPS_FJ(:,3)); zCPS_FJ_71 = smooth(zCPS_FJ,71,'rloess');
zCPS_M = zscore(CPS_M(:,3)); zCPS_M_71 = smooth(zCPS_M,71,'rloess');
zCPS_V = zscore(CPS_V(:,3)); zCPS_V_71 = smooth(zCPS_V,71,'rloess');
AbramSAM(:,2) = zscore(AbramSAM(:,2));

% Put bootstrapped results on correct time-scale
boot_FJ_new = flipud(boot_FJ_rng);
boot_M_new = flipud(boot_M_rng);
boot_V_new = flipud(boot_V_rng);

% CPS recons
xaxis = 1001:1995;

figure(1)
subplot(4,1,1)
axis([1000 2000 -4 4])
jbfill(xaxis,boot_FJ_new(:,1)',boot_FJ_new(:,2)','m','k','add',0.2); hold on
plot(CPS_FJ(:,1),zCPS_FJ,'k','linewidth',1)
plot(CPS_FJ(:,1),zCPS_FJ_71,'m','linewidth',2);

subplot(4,1,2)
axis([1000 2000 -4 4])
jbfill(xaxis,boot_M_new(:,1)',boot_M_rng(:,2)','r','k','add',0.2); hold on
plot(CPS_M(:,1),zCPS_M,'k','linewidth',1)
plot(CPS_M(:,1),zCPS_M_71,'r','linewidth',2);

subplot(4,1,3)
axis([1000 2000 -4 4])
jbfill(xaxis,boot_V_new(:,1)',boot_V_rng(:,2)','b','k','add',0.2); hold on
plot(CPS_V(:,1),zCPS_V,'k','linewidth',1)
plot(CPS_V(:,1),zCPS_V_71,'b','linewidth',2);

subplot(4,1,4)
plot(AbramSAM(:,1),AbramSAM(:,2),'k','linewidth',2);
hold
plot(CPS_FJ(:,1),zCPS_FJ_71,'m','linewidth',2);
plot(CPS_M(:,1),zCPS_M_71,'r','linewidth',2);
plot(CPS_V(:,1),zCPS_V_71,'b','linewidth',2);

print('CPS_bootstrap.pdf','-dpdf','-bestfit')

%% Which recons performed best?
clear
load('prox_groups/bootstrap_proxies.mat','all_CPS_FJ_2','all_CPS_M_2','all_CPS_V_2','smpl_dep')
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load SAM_seasonal.mat; 
load Fogt_Jones.mat; FogtJones_SF = flipud(FogtJones_SF); FogtJones_MA = flipud(FogtJones_MA);

nsim = 10000;

for i = 1:nsim
    r_boot_M(1,i) = corr(all_CPS_M_2(1:39,i),Marshall_SAM(21:59,2));
    r_boot_M(2,i) = corr(all_CPS_M_2(1:91,i),FJ_ann(10:100,2));
    r_boot_M(3,i) = corr(all_CPS_M_2(1:109,i),Visbeck_Ann(11:119,2));
    r_boot_FJ(1,i) = corr(all_CPS_FJ_2(1:39,i),Marshall_SAM(21:59,2));
    r_boot_FJ(2,i) = corr(all_CPS_FJ_2(1:91,i),FJ_ann(10:100,2));
    r_boot_FJ(3,i) = corr(all_CPS_FJ_2(1:109,i),Visbeck_Ann(11:119,2));
    r_boot_V(1,i) = corr(all_CPS_V_2(1:39,i),Marshall_SAM(21:59,2));
    r_boot_V(2,i) = corr(all_CPS_V_2(1:91,i),FJ_ann(10:100,2));
    r_boot_V(3,i) = corr(all_CPS_V_2(1:109,i),Visbeck_Ann(11:119,2)); 
end

[~,inx] = max(r_boot_M(2,:));

% Lets try to create a r-val quantile plot WORK IN PROGRESS, what is below
% already works well.
bin = 0:0.05:1.0;
[sort_corr_M sort_corr_ind_M] = sort(r_boot_M(1,:));
sort_corr_FJ = r_boot_M(2,sort_corr_ind_M);
bin_sizes = histc(squeeze(sort_corr_M),bin);

current_index = 1; corr_M_quan = nan(7,length(bin_sizes)); 
for m=1:length(bin_sizes)
    corr_M_quan(:,m) = quantile(reshape(sort_corr_FJ(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

% Test plot

plot(corr_M_quan(1,:),bin)
hold on
for i = 2:7
    plot(corr_M_quan(i,:),bin)
end

%% Plotting
dima = [.15 .58 .3 .3];
dimb = [.15 .55 .3 .3];
dimc = [.15 .52 .3 .3];

figure(1) % Figure for recons calibrated with Marshall
histogram(r_boot_M(1,:),'FaceAlpha',0.5,'FaceColor','r')
hold
histogram(r_boot_M(2,:),'FaceAlpha',0.5,'FaceColor','m')
histogram(r_boot_M(3,:),'FaceAlpha',0.5,'FaceColor','b')
line([0.8125 0.8125],[0 800],'color','r','linestyle','--') % These values taken from recon_plotting_2
line([0.2605 0.2605],[0 800],'color','m','linestyle','--')
line([0.1808 0.1808],[0 800],'color','b','linestyle','--')
xlabel('r (correlation with SAM indices)');
ylabel('Count (# recons.)');
title('Distribution of boot-strapped reconstruction skill (cal. with Marshall)');
a = annotation('textbox',dima,'String',['r max_{M} = ',num2str(max(r_boot_M(1,:)))],'FitBoxToText','on');
a.LineStyle = 'none';
b = annotation('textbox',dimb,'String',['r max_{FJ} = ',num2str(max(r_boot_M(2,:)))],'FitBoxToText','on');
b.LineStyle = 'none';
c = annotation('textbox',dimc,'String',['r max_{V} = ',num2str(max(r_boot_M(3,:)))],'FitBoxToText','on');
c.LineStyle = 'none';

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperSize = [20 16];
print('boot_hist_M','-dpdf','-bestfit')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = [.15 .58 .3 .3];
figure(2) % Figure for recons calibrated with FJ
subplot(2,2,1)
histogram(skill_noise_FJ(1,:))
line([0.6173 0.6173],[0 800],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with Marshall SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with FJ)');
a = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_FJ(:,1))],'FitBoxToText','on');
a.LineStyle = 'none';

dim = [.6 .58 .3 .3];
subplot(2,2,2)
histogram(skill_noise_FJ(2,:))
line([0.4812 0.4812],[0 800],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with FJ SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with FJ)');
b = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_FJ(:,2))],'FitBoxToText','on');
b.LineStyle = 'none';

dim = [.15 .1 .3 .3];
subplot(2,2,3)
histogram(skill_noise_FJ(3,:))
line([0.3148 0.3148],[0 500],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with V SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with FJ)');
c = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_FJ(:,3))],'FitBoxToText','on');
c.LineStyle = 'none';

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperSize = [18 16];
print('noise_FJ','-dpdf','-bestfit')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = [.15 .58 .3 .3];
figure(3) % Figure for recons calibrated with V
subplot(2,2,1)
histogram(skill_noise_V(1,:))
line([0.4082 0.4082],[0 800],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with Marshall SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with V)');
a = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_V(:,1))],'FitBoxToText','on');
a.LineStyle = 'none';

dim = [.6 .58 .3 .3];
subplot(2,2,2)
histogram(skill_noise_V(2,:))
line([0.2047 0.2047],[0 600],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with FJ SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with V)');
b = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_V(:,2))],'FitBoxToText','on');
b.LineStyle = 'none';

dim = [.15 .1 .3 .3];
subplot(2,2,3)
histogram(skill_noise_V(3,:))
line([0.5221 0.5221],[0 800],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with V SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with V)');
c = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_V(:,3))],'FitBoxToText','on');
c.LineStyle = 'none';

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperSize = [18 16];
print('noise_V','-dpdf','-bestfit')

