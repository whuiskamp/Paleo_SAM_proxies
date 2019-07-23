%% Red noise reconstructions
% This script generates CPS reconstructions using noise time-series based on
% the Ebisuzaki method. Autocorrelated noise records which match the
% spectral properties of the proxies are generated and combined into 
% SAM reconstructions. These are then used to assess the significance of the 
% real proxy reconstructions.
% Willem Huiskamp, 2019
clear
load('JAGS_in.mat','zAll_data_shift')
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load SAM_seasonal.mat; 
load Fogt_Jones.mat; FogtJones_SF = flipud(FogtJones_SF); FogtJones_MA = flipud(FogtJones_MA);

for i = 1:52
    start_p(1,i) = min(find(~isnan(zAll_data_shift(:,i))));
    end_p(1,i) = find(~isnan(zAll_data_shift(:,i)),1,'last');
    End_ma(1,i) = 39;
    End_vbk(1,i) = 109;
    End_FJ(1,i) = 91;
    ma_start(1,i) = 20 + start_p(1,i);
    ma_end(1,i) = 59;
    vbk_start(1,i) = 10 + start_p(1,i);
    vbk_end(1,i) = 119;
    FJ_start(1,i) = 9 + start_p(1,i);
    FJ_end(1,i) = 100;
end
End_vbk(1,5) = 103;
vbk_end(1,5) = 113;

nsim = 10000;
noise_prox = nan(995,52,nsim);
r_noise_prox = nan(3,52,nsim);

zAll_data_shift(isnan(zAll_data_shift)) = 0; % This is required as there can be no NaNs in the data when creating synthetic time-series

% Create 10k synthetic, red noise time-series for each proxy.
% Correlate each synthetic proxy with each SAM index
% This loop takes ~4 minutes
tic
for i = 1:size(zAll_data_shift,2)
    value = rand(1); 
    noise_prox(start_p(1,i):end_p(1,i),i,:) = ebisuzaki(zAll_data_shift(start_p(1,i):end_p(1,i),i),nsim,value); % Create 10,000 synthetic time-series' of each proxy
    for j = 1:nsim
        r_noise_prox(1,i,j) = corr(Marshall_SAM(ma_start(1,i):ma_end(1,i),2),noise_prox(start_p(1,i):End_ma(1,i),i,j)); % calculate their correlation with SAM indices
        r_noise_prox(2,i,j) = corr(FJ_ann(FJ_start(1,i):FJ_end(1,i),2),noise_prox(start_p(1,i):End_FJ(1,i),i,j));
        r_noise_prox(3,i,j) = corr(Visbeck_Ann(vbk_start(1,i):vbk_end(1,i),2),noise_prox(start_p(1,i):End_vbk(1,i),i,j));
    end
    i
end
toc

save('noise_proxies.mat','noise_prox','r_noise_prox','-v7.3') % File version required for large size.

%% Create synthetic SAM reconstructions using a simple proxy * r weighting
% r(1-p) makes no sense here, as we know that by design, they have no
% significance.
clear
load('noise_proxies.mat','noise_prox','r_noise_prox')
nsim = 10000;
sized = size(noise_prox,2); smpl_dep = nan(1,995);
noise_SAM_MA = nan(995,nsim);noise_SAM_FJ = nan(995,nsim);noise_SAM_VBK = nan(995,nsim);

for i = 2:sized
    for j = 1:nsim
        noise_prox(:,i,j) = std_correct(noise_prox(:,i-1,j),noise_prox(:,i,j)); % Correct time series so that std and mean are same
    end
end

for i = 1:length(noise_prox(:,1,1))
    smpl_dep(i) = sum(~isnan(squeeze(noise_prox(i,:,1)))); % The sample depth of each network, for variance correction
end

% Create recons. In r_noise_recon, 1 is Marshall, 2 is FJ and 3 is Visbeck.
noise_prox(isnan(noise_prox)) = 0;
for i = 1:nsim
    noise_SAM_MA(:,i) = squeeze(noise_prox(:,:,i))*squeeze(r_noise_prox(1,:,i))';
    noise_SAM_FJ(:,i) = squeeze(noise_prox(:,:,i))*squeeze(r_noise_prox(2,:,i))';
    noise_SAM_VBK(:,i) = squeeze(noise_prox(:,:,i))*squeeze(r_noise_prox(3,:,i))';
end

% Much like the proper recons, we now need to apply variance correction

for i = 1:994
    index(:,i) = smpl_dep(1,i)-smpl_dep(1,i+1)~=0;
end
[~,inx] = find(index); % create index of nest edges

r_end = 995;
n1_start = 1; % establish the first nest, that all others will be scaled to
n1_end = inx(find(inx >= 10,1,'first')); % This is an imperfect fix, but no real alternative

for data = 1:3
    if data == 1
        recon = noise_SAM_MA; recon2 = recon;
    elseif data == 2
        recon = noise_SAM_FJ; recon2 = recon;
    elseif data == 3
        recon = noise_SAM_VBK; recon2 = recon;
    end
    
    for i = 1:length(noise_SAM_MA)
        for l = find(inx >= 10,1,'first'):length(inx) % we only want to start scaling from the end of the first nest
            if l == length(inx) % If we're on the last nest, iust extend to the end of the record
                n2_start = inx(l)+1; n2_end = r_end;
            elseif abs((inx(l)+1)-inx(l+1)) <=3 && l >= length(inx)-2 % If the next is too small and near end of record, extend it to the end.
                n2_start = inx(l)+1; n2_end = r_end;
            elseif abs((inx(l)+1)-inx(l+1)) <=3  % If the next is too small, extend it with the next one.
                n2_start = inx(l)+1; n2_end = inx(l+2);    
            else
                n2_start = inx(l)+1; n2_end = inx(l+1);
            end
                recon2(n2_start:n2_end,i) = std_correct_CPS(recon(n1_start:n1_end,i),recon(n2_start:n2_end,i));
        end
    end
    if data == 1
       all_CPS_M_2 = recon2;
    elseif data == 2
       all_CPS_FJ_2 = recon2;
    else
       all_CPS_V_2 = recon2;
    end
    clear recon recon2
    data
end

save('noise_proxies.mat','all_CPS_M_2','all_CPS_FJ_2','all_CPS_V_2','-append')

%% Calculate zscores and correlate with SAM indices.
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load Fogt_Jones.mat;
load('SAM_seasonal.mat','Visbeck_Ann')
load('noise_proxies.mat','all_CPS_M_2','all_CPS_FJ_2','all_CPS_V_2')
nsim = 10000;


for i = 1:nsim
    skill_noise_M(1,i) = corr(all_CPS_M_2(1:39,i),Marshall_SAM(21:59,2));
    skill_noise_M(2,i) = corr(all_CPS_M_2(1:91,i),FJ_ann(10:100,2));
    skill_noise_M(3,i) = corr(all_CPS_M_2(1:109,i),Visbeck_Ann(11:119,2));
    skill_noise_FJ(1,i) = corr(all_CPS_FJ_2(1:39,i),Marshall_SAM(21:59,2));
    skill_noise_FJ(2,i) = corr(all_CPS_FJ_2(1:91,i),FJ_ann(10:100,2));
    skill_noise_FJ(3,i) = corr(all_CPS_FJ_2(1:109,i),Visbeck_Ann(11:119,2));
    skill_noise_V(1,i) = corr(all_CPS_V_2(1:39,i),Marshall_SAM(21:59,2));
    skill_noise_V(2,i) = corr(all_CPS_V_2(1:91,i),FJ_ann(10:100,2));
    skill_noise_V(3,i) = corr(all_CPS_V_2(1:109,i),Visbeck_Ann(11:119,2)); 
end

% Calculate significance with respect to red noise sims. Corr values from
% recon_plotting_2.m (lines 112-122)
p_SAM_M(:,1) = 1-(sum(0.8125 > skill_noise_M(1,:)))/nsim; 
p_SAM_M(:,2) = 1-(sum(0.2605 > skill_noise_M(2,:)))/nsim;
p_SAM_M(:,3) = 1-(sum(0.1808 > skill_noise_M(3,:)))/nsim;
p_SAM_FJ(:,1) = 1-(sum(0.6173 > skill_noise_FJ(1,:)))/nsim;
p_SAM_FJ(:,2) = 1-(sum(0.4812 > skill_noise_FJ(2,:)))/nsim;
p_SAM_FJ(:,3) = 1-(sum(0.3148 > skill_noise_FJ(3,:)))/nsim;
p_SAM_V(:,1) = 1-(sum(0.4082 > skill_noise_V(1,:)))/nsim;
p_SAM_V(:,2) = 1-(sum(0.2047 > skill_noise_V(2,:)))/nsim;
p_SAM_V(:,3) = 1-(sum(0.5221 > skill_noise_V(3,:)))/nsim;

%% Plotting
dim = [.15 .58 .3 .3];
figure(1) % Figure for recons calibrated with Marshall
subplot(2,2,1)
histogram(skill_noise_M(1,:))
line([0.8125 0.8125],[0 700],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with Marshall SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with Marshall)');
a = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_M(:,1))],'FitBoxToText','on');
a.LineStyle = 'none';

dim = [.6 .58 .3 .3];
subplot(2,2,2)
histogram(skill_noise_M(2,:))
line([0.2605 0.2605],[0 600],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with FJ SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with Marshall)');
b = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_M(:,2))],'FitBoxToText','on');
b.LineStyle = 'none';

dim = [.15 .1 .3 .3];
subplot(2,2,3)
histogram(skill_noise_M(3,:))
line([0.1808 0.1808],[0 600],'color','r','linestyle','--') % These values taken from recon_plotting_2
xlabel('r (correlation with V SAM)');
ylabel('Count (# recons.)');
title('Distribution of red-noise reconstruction skill (cal. with Marshall)');
c = annotation('textbox',dim,'String',['p = ',num2str(p_SAM_M(:,3))],'FitBoxToText','on');
c.LineStyle = 'none';

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperSize = [18 16];
print('noise_M','-dpdf','-fillpage')

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













