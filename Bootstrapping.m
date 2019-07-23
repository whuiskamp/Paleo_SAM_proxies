% Bootstrapping tests
% This script effectively re-does the CPS reconstructions, but this time
% allows re-sampling of proxies to see how sensitive the result is to the
% use of any one proxy.

%% Create Proxy pool
% Load correlation data for each proxy, the proxies themselves
clear
load JAGS_in.mat
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load Fogt_Jones.mat;
load('SAM_seasonal.mat','Visbeck_Ann')

% Define vars
num_reps = 10000;
proxies = nan(length(num_reps),length(zAll_data_shift),52);
smpl_dep = nan(num_reps,length(zAll_data_shift));
corr_FJ = nan(num_reps,1,52); corr_M = nan(num_reps,1,52); corr_V = nan(num_reps,1,52);
rp_FJ = nan(num_reps,1,52); rp_M = nan(num_reps,1,52); rp_V = nan(num_reps,1,52);
zAll_data_2 = zAll_data_shift;

% Correct time series so that std and mean are same (match our longest
% record)
for i = 2:52
    zAll_data_2(:,i) = std_correct(zAll_data_shift(:,1),zAll_data_shift(:,i));
end

% Generate and save groups of 52 proxies with resampling for later
% reconstruction.

tic
for j = 1:num_reps
       [~,idx] = datasample(zAll_data_2(1,:),52,'Replace',true); % randomly sample a proxy pool
       proxies(j,:,:) = zAll_data_2(:,idx); % 1000 randomly sampled networks
       corr_FJ(j,:,:) = corrs_SAM(3,idx); % Corresponding correlation info
       corr_M(j,:,:) = corrs_SAM(1,idx);
       corr_V(j,:,:) = corrs_SAM(2,idx);
       rp_FJ(j,:,:) = skill_FJ(1,idx); % 'skill' value, as per Tierney et al.
       rp_M(j,:,:) = skill_marsh(1,idx);
       rp_V(j,:,:) = skill_vb(1,idx);
       for k = 1:length(zAll_data_2)
            smpl_dep(j,k) = sum(~isnan(squeeze(proxies(j,k,:)))); % The sample depth of each network, for later variance correction
       end
end
toc % takes about 1.5 hours for 10k, or one minute for 1k.
save('prox_groups/bootstrap_proxies.mat',...
        'proxies','corr_FJ','corr_M','corr_V','rp_FJ','rp_M','rp_V','smpl_dep','-v7.3') % File version required for large size.
    
%% Create recons

clear 
load('prox_groups/bootstrap_proxies.mat')
proxies(isnan(proxies)) = 0;
num_reps = 10000;

for i = 1:num_reps
    for j = 1:size(proxies,2) % Goes year-by-year, multiply each proxy with it's skill, then sum
        CPS_FJ(j,i) = squeeze(proxies(i,j,:))'*squeeze(rp_FJ(i,1,:)); % change these from 'corr' to 'rp' as desired
        CPS_M(j,i) = squeeze(proxies(i,j,:))'*squeeze(rp_M(i,1,:));
        CPS_V(j,i) = squeeze(proxies(i,j,:))'*squeeze(rp_V(i,1,:));
    end
end % takes about 30 seconds for 1k.

save('prox_groups/bootstrap_proxies.mat',...
    'CPS_FJ','CPS_M','CPS_V','-append')
%% Correct for Sample depth changes
clear
load('prox_groups/bootstrap_proxies.mat','CPS_FJ','CPS_M','CPS_V','smpl_dep')
num_reps = length(smpl_dep);

tic
for data = 1:3
    if data == 1
        recon = CPS_M; recon2 = nan(size(CPS_M));
    elseif data == 2
        recon = CPS_FJ; recon2 = nan(size(CPS_FJ));
    elseif data == 3
        recon = CPS_V; recon2 = nan(size(CPS_V));
    end
    
    for j = 1:num_reps
        r_end = find(smpl_dep(j,:) == 0,1,'first')-1; % Find the end-point of this recon
        if isempty(r_end)
           r_end = 995;
        end
        index = nan(1,r_end); % Vector for index of transition points
        for k = 1:r_end-1
            index(:,k) = squeeze(smpl_dep(j,k)-smpl_dep(j,k+1))~=0; % Find transition points in network size
        end
        [~,inx] = find(index); % create index of transition points
        n1_start = 1; % establish the first nest, that all others will be scaled to
        n1_end = inx(find(inx >= 10,1,'first')); % This is a dirty fix, maybe improve this later
        recon2(1:n1_end,j) = zscore(recon(1:n1_end,j));
        for l = find(inx >= 10,1,'first'):length(inx) % we only want to start scaling from the end of the first nest
            if inx(l) == r_end
                continue
            elseif l == length(inx) % If we're on the last nest, just extend to the end of the record
                n2_start = inx(l)+1; n2_end = r_end;
            elseif abs((inx(l)+1)-inx(l+1)) <=3 && l >= length(inx)-2 % If the next is too small, extend it with the next one.
                n2_start = inx(l)+1; n2_end = r_end;
            elseif abs((inx(l)+1)-inx(l+1)) <=3  % If the next is too small, extend it with the next one.
                n2_start = inx(l)+1; n2_end = inx(l+2);    
            else
                n2_start = inx(l)+1; n2_end = inx(l+1);
            end
                recon2(n2_start:n2_end,j) = zscore(std_correct_CPS(recon(n1_start:n1_end,j),recon(n2_start:n2_end,j)));
            if max(max(recon2)) >= 10 % This value is arbitrary for the most part
               break
            end
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
    toc
end
toc % Takes 3.5 hrs for 10k and ~6 mins for 1K

save('prox_groups/bootstrap_proxies.mat',...
    'all_CPS_FJ_2','all_CPS_M_2','all_CPS_V_2','-append')

%% Save range as quantiles for later plotting

clear
load('prox_groups/bootstrap_proxies.mat','all_CPS_FJ_2','all_CPS_M_2','all_CPS_V_2','smpl_dep')

for i = 1:size(smpl_dep,2)
    boot_FJ_qn(i,:) = quantile(squeeze(all_CPS_FJ_2(i,:)),[.05 0.25 .5 0.75 .95]);
    boot_M_qn(i,:) = quantile(squeeze(all_CPS_M_2(i,:)),[.05 0.25 .5 0.75 .95]);
    boot_V_qn(i,:) = quantile(squeeze(all_CPS_V_2(i,:)),[.05 0.25 .5 0.75 .95]);
end

save('prox_groups/bootstrap_proxies.mat',...
    'boot_FJ_qn','boot_M_qn','boot_V_qn','-append')

% As a test, lets also save the max/min values

for i = 1:size(smpl_dep,2)
    boot_FJ_rng(i,1) = max(all_CPS_FJ_2(i,:));
    boot_FJ_rng(i,2) = min(all_CPS_FJ_2(i,:));
    boot_M_rng(i,1) = max(all_CPS_M_2(i,:));
    boot_M_rng(i,2) = min(all_CPS_M_2(i,:));
    boot_V_rng(i,1) = max(all_CPS_V_2(i,:));
    boot_V_rng(i,2) = min(all_CPS_V_2(i,:));
end

save('prox_groups/bootstrap_proxies.mat',...
    'boot_FJ_rng','boot_M_rng','boot_V_rng','-append')










