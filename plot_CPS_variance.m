%% Plot range of CPS reconstructions, rather than their skill
% I want a figure showing the rage of values in a SAM recon, over the 995
% years and how this changes depending on which index we calibrate with.
clear
load('all_CPS_recons.mat','all_CPS_FJ','all_CPS_M','all_CPS_V')
load('JAGS_in.mat','zAll_data_shift')
load('prox_groups/sample_depth.mat')
for i = 1:995
    num_prx(i) = sum(~isnan(zAll_data_shift(i,:)));
end

num_stns = 52;
num_reps = 1000;
srs_len = 995;
%% Correct for Sample depth changes

tic
for data = 1:3
    if data == 1
        recon = all_CPS_M; recon2 = nan(size(all_CPS_M));
    elseif data == 2
        recon = all_CPS_FJ; recon2 = nan(size(all_CPS_FJ));
    elseif data == 3
        recon = all_CPS_V; recon2 = nan(size(all_CPS_V));
    end
    for i = 45 %10:num_stns
        for j = 1:num_reps
            r_end = find(smpl_dep(j,:,i) == 0,1,'first')-1; % Find the end-point of this recon
            if isempty(r_end)
                r_end = 995;
            end
            index = nan(1,r_end); % Vector for index of transition points
            for k = 1:r_end-1
                index(:,k) = squeeze(smpl_dep(j,k,i)-smpl_dep(j,k+1,i))~=0; % Find transition points in network size
            end
            [~,inx] = find(index); % create index of transition points
            n1_start = 1; % establish the first nest, that all others will be scaled to
            n1_end = inx(find(inx >= 10,1,'first')); % This is a dirty fix, maybe improve this later
            recon2(i,1:n1_end,j) = zscore(recon(i,1:n1_end,j));
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
                recon2(i,n2_start:n2_end,j) = zscore(std_correct_CPS(recon(i,n1_start:n1_end,j),recon(i,n2_start:n2_end,j)));
            end
            if max(max(max(recon2))) >= 10 % This value is arbitrary for the most part
               break
            end
            i;j
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
end
toc 
 

    


%% Plot Variance
CPS_FJ_var = var(squeeze(all_CPS_FJ(:,50,:)),0,3);
CPS_M_var = var(all_CPS_M,0,3);
CPS_V_var = var(all_CPS_V,0,3);

num = 40;
figure(1); plot(CPS_FJ_var(num,:))
%hold;
plot(CPS_M_var(num,:))
plot(CPS_V_var(num,:))

%% Plot range of values
% perhaps a good idea to plot the range for 25 proxies, as this is how many
% are present in Nerilie's recon. Compare this to 40, perhaps?

% Preprocessing
max_FJ_smth = nan(52,995); min_FJ_smth = nan(52,995);
max_M_smth = nan(52,995); min_M_smth = nan(52,995);
max_V_smth = nan(52,995); min_V_smth = nan(52,995);

% Find range of values over 1000 iterations
max_FJ = max(all_CPS_FJ_2,[],3);min_FJ = min(all_CPS_FJ_2,[],3); ave_FJ = nanmean(all_CPS_FJ_2,3);
max_M = max(all_CPS_M_2,[],3);min_M = min(all_CPS_M_2,[],3); ave_M = nanmean(all_CPS_M_2,3);
max_V = max(all_CPS_V_2,[],3);min_V = min(all_CPS_V_2,[],3); ave_V = nanmean(all_CPS_V_2,3);

% Get a smoothed version, as data is quite noisy
% Remember, data only exists for 10:52
for i = 10:52
    max_FJ_smth(i,:) = smooth(squeeze(max_FJ(i,:)),71,'rloess'); min_FJ_smth(i,:) = smooth(squeeze(min_FJ(i,:)),71,'rloess');
    max_M_smth(i,:) = smooth(squeeze(max_M(i,:)),71,'rloess'); min_M_smth(i,:) = smooth(squeeze(min_M(i,:)),71,'rloess');
    max_V_smth(i,:) = smooth(squeeze(max_V(i,:)),71,'rloess'); min_V_smth(i,:) = smooth(squeeze(min_V(i,:)),71,'rloess');
    ave_FJ_smth(i,:) = smooth(squeeze(ave_FJ(i,:)),71,'rloess'); ave_M_smth(i,:) = smooth(squeeze(ave_M(i,:)),71,'rloess');
    ave_V_smth(i,:) = smooth(squeeze(ave_V(i,:)),71,'rloess');
end

% Plot
prxs = 45; % how many proxies in the recon?
figure(1)
hold on
line([0 1000],[0 0],'color','k','linestyle','--')
% FJ 
l1 = plot(max_FJ(prxs,:));
l2 = plot(min_FJ(prxs,:));
l1.Color = [1,0,1,0.1]; % sets colour and transparency 
l2.Color = [1,0,1,0.1];
plot(max_FJ_smth(prxs,:),'color','m','linewidth',1,'linestyle','--');
plot(min_FJ_smth(prxs,:),'color','m','linewidth',1,'linestyle','--');
plot(ave_FJ_smth(prxs,:),'color','m','linewidth',2);

% Marshall
l3 = plot(max_M(prxs,:));
l4 = plot(min_M(prxs,:));
l3.Color = [1,0,0,0.1]; % sets colour and transparency 
l4.Color = [1,0,0,0.1];
plot(max_M_smth(prxs,:),'color','r','linewidth',1,'linestyle','--');
plot(min_M_smth(prxs,:),'color','r','linewidth',1,'linestyle','--');
plot(ave_M_smth(prxs,:),'color','r','linewidth',2);

%Visbeck
l9 = plot(max_V(prxs,:));
l10 = plot(min_V(prxs,:));
l9.Color = [0,0,1,0.1]; % sets colour and transparency 
l10.Color = [0,0,1,0.1];
plot(max_V_smth(prxs,:),'color','b','linewidth',1,'linestyle','--');
plot(min_V_smth(prxs,:),'color','b','linewidth',1,'linestyle','--');
plot(ave_V_smth(prxs,:),'color','b','linewidth',2);


xlabel('Years before 1995')
ylabel('SAM Index')
title(['Range of reconstructed SAM values for ',num2str(prxs),' of 52 randomly sampled proxies'])
print('Range_45p_SAM','-dpdf','-bestfit')












