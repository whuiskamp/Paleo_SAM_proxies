% CPS skill range for paleo-proxies
% This will be done using both Marshall and F&J. Similar to the pseudoproxy
% work, these will represent, effectively, different windowsizes.

clear
load JAGS_in.mat
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load Fogt_Jones.mat;
load('SAM_seasonal.mat','Visbeck_Ann')
num_stns = 10:52;
num_reps = 1:1000;
proxies = nan(length(num_reps),length(zAll_data_shift),max(num_stns));
CPS_FJ = nan(995,1000);
CPS_M = nan(995,1000);
CPS_V = nan(995,1000);
all_CPS_FJ = nan(52,995,1000);
all_CPS_M = nan(52,995,1000);
all_CPS_V = nan(52,995,1000);
skill_FJ_CPS = nan(length(num_reps),1);skill_M_CPS = nan(length(num_reps),1);skill_V_CPS = nan(length(num_reps),1);
skill_FJ_M_CPS = nan(max(num_stns),1);skill_M_FJ_CPS = nan(max(num_stns),1);skill_FJ_V_CPS = nan(max(num_stns),1);
skill_M_V_CPS = nan(max(num_stns),1);skill_V_FJ_CPS = nan(max(num_stns),1);skill_V_M_CPS = nan(max(num_stns),1);

skill_FJ_qn = nan(max(num_stns),3);skill_FJ_M_qn = nan(max(num_stns),3);skill_FJ_V_qn = nan(max(num_stns),3);
skill_M_qn = nan(max(num_stns),3);skill_M_FJ_qn = nan(max(num_stns),3);skill_M_V_qn = nan(max(num_stns),3);
skill_V_qn = nan(max(num_stns),3);skill_V_FJ_qn = nan(max(num_stns),3);skill_V_M_qn = nan(max(num_stns),3);

sized = size(zAll_data_shift,2);

len = nan(2,sized);
len(1,:) = 1:sized;
for i = 1:sized
    len(2,i) = max(find(~isnan(zAll_data_shift(:,i)))); % Finds the last value in each column that isn't NaN
end
[~,inx] = sort(len(2,:),2,'descend');
zAll_data_2 = zAll_data_shift(:,inx);
skill_FJ = skill_FJ(:,inx); skill_marsh = skill_marsh(:,inx); skill_vb = skill_vb(:,inx);
inx(2,:) = len(2,inx);
inx = fliplr(inx); % So we can fix variance later

% Correct time series so that std and mean are same
for i = 2:sized
    zAll_data_2(:,i) = std_correct(zAll_data_2(:,i-1),zAll_data_2(:,i));
end

for i = num_stns
    for j = num_reps
       [~,idx] = datasample(zAll_data_shift(1,:),i,'Replace',false);
       proxies(j,:,1:i) = zAll_data_shift(:,idx);
       corr_FJ(j,:,1:i) = corrs_SAM(3,idx);
       corr_M(j,:,1:i) = corrs_SAM(1,idx);
       corr_V(j,:,1:i) = corrs_SAM(2,idx);
       rp_FJ(j,:,1:i) = skill_FJ(1,idx);
       rp_M(j,:,1:i) = skill_marsh(1,idx);
       rp_V(j,:,1:i) = skill_vb(1,idx);
    end
    save(['prox_groups/',num2str(i),'_proxies.mat'],...
        'proxies','corr_FJ','corr_M','corr_V','rp_FJ','rp_M','rp_V')
end

% If you just want your straight up recons
    zAll_data_2(isnan(zAll_data_2)) = 0;
    CPS_FJ(:,1) = All_data_final(:,1);
    CPS_M(:,1) = All_data_final(:,1);
    CPS_V(:,1) = All_data_final(:,1);
for k = 1:length(zAll_data_2)
    CPS_FJ(k,2) = squeeze(zAll_data_2(k,:))*skill_FJ'; 
    CPS_M(k,2) = squeeze(zAll_data_2(k,:))*skill_marsh';
    CPS_V(k,2) = squeeze(zAll_data_2(k,:))*skill_vb';
end

% Correct variance !!! THIS METHOD IS WRONG !!! REDO.
for a = 1:3
    if a == 1
        dataset = CPS_FJ(:,2);
    elseif a == 2 
        dataset = CPS_M(:,2);
    elseif a == 3
        dataset = CPS_V(:,2);
    end    
        dataset(1:inx(2,1),2) = dataset(1:inx(2,1),1);
        for i = 2:42  % This corrects variance, not mean
            if inx(2,i)-inx(2,i-1) == 0
                dataset(inx(2,i-2):inx(2,i),2) = std_correct_CPS(dataset((1:inx(2,1)),1),dataset(inx(2,i-2):inx(2,i),1));
            else
                dataset(inx(2,i-1):inx(2,i),2) = std_correct_CPS(dataset((1:inx(2,1)),1),dataset(inx(2,i-1):inx(2,i),1));
            end
        end
    if a == 1
        CPS_FJ(:,3) = dataset(:,2);
    elseif a == 2
        CPS_M(:,3) = dataset(:,2);
    elseif a == 3 
        CPS_V(:,3) = dataset(:,2);
    end
end
        save('CPS_results.mat','CPS_FJ','CPS_M','CPS_V','-append')

%% Calculate randomly sampled reconstructions  
smpl_dep = nan(size(proxies));
tic
for i = num_stns % This takes 2 hours to run.
    load(['prox_groups/',num2str(i),'_proxies.mat'],...
        'proxies','corr_FJ','corr_M','corr_V','rp_FJ','rp_M','rp_V')
    proxies2 = proxies;
    proxies(isnan(proxies)) = 0;
    for j = num_reps
        for k = 1:length(zAll_data_shift)
            smpl_dep(j,k,i) = sum(~isnan(squeeze(proxies2(j,k,:))));
            CPS_FJ(k,j) = squeeze(proxies(j,k,1:i))'*squeeze(rp_FJ(j,1,1:i)); % change these from 'corr' to 'rp' as wanted
            CPS_M(k,j) = squeeze(proxies(j,k,1:i))'*squeeze(rp_M(j,1,1:i));
            CPS_V(k,j) = squeeze(proxies(j,k,1:i))'*squeeze(rp_V(j,1,1:i));
            
            % calculate sample depth of each proxy network, for later variance correction
            
        end
        
       CPS_FJ(:,j) = zscore(CPS_FJ(:,j));
       CPS_M(:,j) = zscore(CPS_M(:,j));
       CPS_V(:,j) = zscore(CPS_V(:,j));
%        skill_FJ_CPS(j,:) = corr(CPS_FJ(1:91,j),FJ_ann(10:100,2));
%        skill_FJ_M_CPS(j,:) = corr(CPS_FJ(1:39,j),Marshall_SAM(21:59,2));
%        skill_FJ_V_CPS(j,:) = corr(CPS_FJ(1:109,j),Visbeck_Ann(11:119,2));
%        skill_M_CPS(j,:) = corr(CPS_M(1:39,j),Marshall_SAM(21:59,2));
%        skill_M_FJ_CPS(j,:) = corr(CPS_M(1:91,j),FJ_ann(10:100,2));
%        skill_M_V_CPS(j,:) = corr(CPS_M(1:109,j),Visbeck_Ann(11:119,2));
%        skill_V_CPS(j,:) = corr(CPS_V(1:109,j),Visbeck_Ann(11:119,2));
%        skill_V_FJ_CPS(j,:) = corr(CPS_V(1:91,j),FJ_ann(10:100,2));
%        skill_V_M_CPS(j,:) = corr(CPS_V(1:39,j),Marshall_SAM(21:59,2));
    end
    all_CPS_FJ(i,:,:) = CPS_FJ(:,:); % This saves all the individual recons
    all_CPS_M(i,:,:) = CPS_M(:,:);
    all_CPS_V(i,:,:) = CPS_V(:,:);
%     skill_FJ_qn(i,:) = quantile(skill_FJ_CPS,[.05 .5 .95],1); % Save the skill of each network size, as quantiles
%     skill_FJ_M_qn(i,:) = quantile(skill_FJ_M_CPS,[.05 .5 .95],1);
%     skill_FJ_V_qn(i,:) = quantile(skill_FJ_V_CPS,[.05 .5 .95],1);
%     skill_M_qn(i,:) = quantile(skill_M_CPS,[.05 .5 .95],1);
%     skill_M_FJ_qn(i,:) = quantile(skill_M_FJ_CPS,[.05 .5 .95],1);
%     skill_M_V_qn(i,:) = quantile(skill_M_V_CPS,[.05 .5 .95],1);
%     skill_V_qn(i,:) = quantile(skill_V_CPS,[.05 .5 .95],1);
%     skill_V_FJ_qn(i,:) = quantile(skill_V_FJ_CPS,[.05 .5 .95],1);
%     skill_V_M_qn(i,:) = quantile(skill_V_M_CPS,[.05 .5 .95],1);
    clear CPS_FJ CPS_M CPS_V skill_FJ_CPS skill_FJ_M_CPS skill_FJ_V_CPS ...
        skill_M_CPS skill_M_FJ_CPS skill_M_V_CPS skill_V_CPS skill_V_FJ_CPS ...
        skill_V_M_CPS
    i
    toc
end
save('CPS_results.mat','skill_FJ_qn','skill_FJ_M_qn','skill_FJ_V_qn',...
    'skill_M_qn','skill_M_FJ_qn','skill_M_V_qn',...
    'skill_V_qn','skill_V_FJ_qn','skill_V_M_qn','-append');
save('all_CPS_recons.mat','all_CPS_FJ','all_CPS_M','all_CPS_V')
save('prox_groups/sample_depth.mat','smpl_dep')

%% Calculate regional skill 

AA_skill_FJ = mean(abs(skill_FJ(1,[1:12 15:24])));
AA_skill_M = mean(abs(skill_marsh(1,[1:12 15:24])));
AA_skill_V = mean(abs(skill_vb(1,[1:12 15:24])));

SA_skill_FJ = mean(abs(skill_FJ(1,[14 25:38])));
SA_skill_M = mean(abs(skill_marsh(1,[14 25:38])));
SA_skill_V = mean(abs(skill_vb(1,[14 25:38])));

AuNZ_skill_FJ = mean(abs(skill_FJ(1,[13 39:52])));
AuNZ_skill_M = mean(abs(skill_marsh(1,[13 39:52])));
AuNZ_skill_V = mean(abs(skill_vb(1,[13 39:52])));

