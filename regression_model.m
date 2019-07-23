% This script reads in the proxy data matrix, tests for correlations with
% various SAM indices and combines them using a multilinear regression
% model. 
% For this script to work, Matjags, weaclim and ebisuzaki packages are required. For
% the regression model to work, the JAGS model must be installed
% Willem Huiskamp, 2015

clear
load('JAGS_in.mat','All_data_final','zAll_data_shift')
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load SAM_seasonal.mat; 
load Fogt_Jones.mat; FogtJones_SF = flipud(FogtJones_SF); FogtJones_MA = flipud(FogtJones_MA);
% Remove proxies you don't want anymore

% All_data_final_long(:,[4 5 6 7 14 17 19 21 22 23 24]) = []; % Removes Villalba data and duplicates.
% All_data_final_long(1:17,:) = []; % Removes 2012-1996
% meta_data(:,[21 22 23 24]) = [];

X = nan(1013,size(zAll_data_shift,2));
corrs_SAM = nan(5,size(zAll_data_shift,2));
corrs_SAM_detr = nan(5,size(zAll_data_shift,2));
p_SAM = nan(5,size(zAll_data_shift,2));
visbeck_syn = nan(91,size(zAll_data_shift,2),10000);
visbeck_syncor = nan(1,size(zAll_data_shift,2),10000);
visbeck_pc = nan(size(zAll_data_shift,2),100);
FJ_SF_syn = nan(79,size(zAll_data_shift,2),10000);
FJ_MA_syn = nan(79,size(zAll_data_shift,2),10000);
FJ_SF_syncor = nan(1,size(zAll_data_shift,2),10000);
FJ_MA_syncor = nan(1,size(zAll_data_shift,2),10000);
FJ_SF_pc = nan(size(zAll_data_shift,2),100);
FJ_MA_pc = nan(size(zAll_data_shift,2),100);

% All_data_final = All_data_final_long;

% Step 1: Standardise the proxies
for i = 2:size(All_data_final,2)
    xmu=nanmean(All_data_final(:,i));
    xsigma=nanstd(All_data_final(:,i));
    zAll_data(:,i-1)=(All_data_final(:,i)-repmat(xmu,length(All_data_final(:,i)),1))./repmat(xsigma,length(All_data_final(:,i)),1);
end

% Correct Tree ring records for the Schulman shift
% This was corrected on 15/07/19 and a new version of 
% zAll_data_shift has been saved to JAGS_in.mat

zAll_data_shift(:,1:24) = zAll_data(:,1:24); % 
zAll_data_shift(1:995,25:52) = zAll_data(2:996,25:52); 
zAll_data_shift(996,:) = [];  % Final version saved in 'JAGS_in'

%% Calculate correlations with SAM
% Step 2: Make sure the response to SAM is the same

load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load SAM_seasonal.mat; 
load Fogt_Jones.mat; FogtJones_SF = flipud(FogtJones_SF); FogtJones_MA = flipud(FogtJones_MA);
load('JAGS_in.mat','zAll_data_shift')
% load SAM_inxs_detrended.mat % This was done to test detrending, leave
% commented out

for i = 1:52
    start(1,i) = min(find(~isnan(zAll_data_shift(:,i))));
    End_ma(1,i) = 39;
    End_vbk(1,i) = 109;
    End_FJ(1,i) = 91;
    ma_start(1,i) = 20 + start(1,i);
    ma_end(1,i) = 59;
    vbk_start(1,i) = 10 + start(1,i);
    vbk_end(1,i) = 119;
    FJ_start(1,i) = 9 + start(1,i);
    FJ_end(1,i) = 100;
end
End_vbk(1,5) = 103;
vbk_end(1,5) = 113;

for i = 1:52 % Use on shifted proxies
    corrs_SAM(1,i) = corr(Marshall_SAM(ma_start(1,i):ma_end(1,i),2),zAll_data_shift(start(1,i):End_ma(1,i),i));
    corrs_SAM(2,i) = corr(Visbeck_Ann(vbk_start(1,i):vbk_end(1,i),2),zAll_data_shift(start(1,i):End_vbk(1,i),i));
    corrs_SAM(3,i) = corr(FJ_ann(FJ_start(1,i):FJ_end(1,i),2),zAll_data_shift(start(1,i):End_FJ(1,i),i));
%    corrs_SAM_detr(1,i) = corr(Marshall_SAM_detr(ma_start(1,i):ma_end(1,i),2),zAll_data_shift(start(1,i):End_ma(1,i),i));
%    corrs_SAM_detr(2,i) = corr(Visbeck_SAM_detr(vbk_start(1,i):vbk_end(1,i),2),zAll_data_shift(start(1,i):End_vbk(1,i),i));
%    corrs_SAM_detr(3,i) = corr(FJ_SAM_detr(FJ_start(1,i):FJ_end(1,i),2),zAll_data_shift(start(1,i):End_FJ(1,i),i));
%     corrs_SAM(4,i) = corr(FogtJones_SF(16:94,2),zAll_data_shift(30:108,i));
%     corrs_SAM(5,i) = corr(FogtJones_MA(17:95,2),zAll_data_shift(30:108,i));
end
save('JAGS_in.mat','corrs_SAM','-append')

% Calculate a separate running correlation criteria using the Fogt index.

for i = 1:52
    runcorr_SAM(:,i) = movingCorrelation(zAll_data_shift
    
    
end
%% Calculate significance for the correlation using synthetic time-series
nsim = 10000;

% Marshall Annual
tic
for i = 1:size(zAll_data_shift,2)
    value = rand(1); 
    range = start(1,i):End_ma(1,i);
    SAM_marshall = ma_start(1,i):ma_end(1,i);
    
    marshall_syn(:,:) = ebisuzaki(zAll_data_shift(range,i),nsim,value); % Create 10,000 synthetic time-series' of each proxy
    for j = 1:10000
        marshall_syncor(:,i,j) = corr(Marshall_SAM(SAM_marshall,2),marshall_syn(:,j)); % calculate their correlation with SAM
    end
    marshall_pc(i,:) = quantile(squeeze(marshall_syncor(:,i,:)),0.01:0.01:1, 1); % Create the probability distribution of each correlation
    if corrs_SAM(1,i) > 0
        p_SAM(1,i) = 1-(sum(corrs_SAM(1,i) > marshall_pc(i,:)))/100;
    elseif corrs_SAM(1,i) < 0
        p_SAM(1,i) = 1-(sum(corrs_SAM(1,i) < marshall_pc(i,:)))/100;
    end
    clear marshall_syn
    i
end

% Visbeck annual

for i = 1:size(zAll_data_shift,2)
    value = rand(1); 
    range = start(1,i):End_vbk(1,i);
    SAM_range = vbk_start(1,i):vbk_end(1,i);
   
    visbeck_syn(:,:) = ebisuzaki(zAll_data_shift(range,i),nsim,value); % Create 10,000 synthetic time-series' of each proxy
    for j = 1:10000
        visbeck_syncor(:,i,j) = corr(Visbeck_Ann(SAM_range,2),visbeck_syn(:,j)); % calculate their correlation with SAM
    end
    visbeck_pc(i,:) = quantile(squeeze(visbeck_syncor(:,i,:)),0.01:0.01:1, 1); % Create the probability distribution of each correlation
    if corrs_SAM(2,i) > 0
        p_SAM(2,i) = 1-(sum(corrs_SAM(2,i) > visbeck_pc(i,:)))/100;
    elseif corrs_SAM(2,i) < 0
        p_SAM(2,i) = 1-(sum(corrs_SAM(2,i) < visbeck_pc(i,:)))/100;
    end
    clear visbeck_syn
    i
end

% FJ annual

for i = 1:size(zAll_data_shift,2)
    value = rand(1); 
    range = start(1,i):End_FJ(1,i);
    SAM_FJ = FJ_start(1,i):FJ_end(1,i);
   
    FJ_syn(:,:) = ebisuzaki(zAll_data_shift(range,i),nsim,value); % Create 10,000 synthetic time-series' of each proxy
    for j = 1:10000
        FJ_syncor(:,i,j) = corr(FJ_ann(SAM_FJ,2),FJ_syn(:,j)); % calculate their correlation with SAM
    end
    FJ_pc(i,:) = quantile(squeeze(FJ_syncor(:,i,:)),0.01:0.01:1, 1); % Create the probability distribution of each correlation
    if corrs_SAM(3,i) > 0
        p_SAM(3,i) = 1-(sum(corrs_SAM(3,i) > FJ_pc(i,:)))/100;
    elseif corrs_SAM(3,i) < 0
        p_SAM(3,i) = 1-(sum(corrs_SAM(3,i) < FJ_pc(i,:)))/100;
    end
    clear FJ_syn
    i
end
toc
% FJ - SF and MA

for i = 1:size(zAll_data,2)
    value = rand(1); 
    range = 30:108;
    SAM_SF = 16:94;
    SAM_MA = 17:95;
    FJ_SF_syn(:,i,:) = ebisuzaki(zAll_data(range,i),nsim,value);
    FJ_MA_syn(:,i,:) = ebisuzaki(zAll_data(range,i),nsim,value);
    for j = 1:10000
        FJ_SF_syncor(:,i,j) = corr(FogtJones_SF(SAM_SF,2),FJ_SF_syn(:,i,j)); 
        FJ_MA_syncor(:,i,j) = corr(FogtJones_MA(SAM_MA,2),FJ_MA_syn(:,i,j));
    end
    FJ_SF_pc(i,:) = quantile(squeeze(FJ_SF_syncor(:,i,:)),0.01:0.01:1, 1);
    FJ_MA_pc(i,:) = quantile(squeeze(FJ_MA_syncor(:,i,:)),0.01:0.01:1, 1);
    if corrs_SAM(4,i) > 0
        p_SAM(4,i) = 1-(sum(corrs_SAM(4,i) > FJ_SF_pc(i,:)))/100;
    elseif corrs_SAM(4,i) < 0
        p_SAM(4,i) = 1-(sum(corrs_SAM(4,i) < FJ_SF_pc(i,:)))/100;
    end
    if corrs_SAM(5,i) > 0
        p_SAM(5,i) = 1-(sum(corrs_SAM(5,i) > FJ_MA_pc(i,:)))/100;
    elseif corrs_SAM(5,i) < 0
        p_SAM(5,i) = 1-(sum(corrs_SAM(5,i) < FJ_MA_pc(i,:)))/100;
    end
end

clear visbeck_syn visbeck_syncor FJ_syn FJ_SF_syncor 
% You have to eyeball which SAM index works better here, including which
% correlations are significant 

% for i = [1 5 6 7 11 12 13 14 15 16 17 18 19 22 23 25 27 28 31 33 34 36 37 38 39 41 42 46 47 48 49 50 51 51 53 56 58] % Whichever proxies needed to be changed
%    zAll_data(:,i) = zAll_data(:,i)*-1;
% end

% Calibrate only using FJ index now

% Filter out records that are well correlated using a r(1-p) method of
% skill: We want the skill score to be bigger than 0.18 - 0.2(1-0.1)
% This is a pretty awfully low bar, but otherwise we exclude all the data.

skill_FJ = corrs_SAM(3,:).*(1-p_SAM(3,:));
skill_marsh = corrs_SAM(1,:).*(1-p_SAM(1,:));
skill_vb = corrs_SAM(2,:).*(1-p_SAM(2,:));

% Next find proxies with only a significant correlation at 0.1

skillful_FJ_p = find(p_SAM(3,:) <= 0.1); 
proxies_pskill = zAll_data(:,skillful_FJ_p);

% Next find only proxies with a correlation above 0.15

skillful_FJ_r = find(abs(corrs_SAM(3,:)) >= 0.15); 
proxies_rskill = zAll_data_shift(:,skillful_FJ_r);

%% Regression Model
%%%% Which data will you use? %%%%
% You have to run the model 3 separate times, once for each SAM index
load JAGS_in.mat
negs = find(corrs_SAM(3,:) < 0); % Change to 1 for Marshall, 2 for visbeck
for i = negs
   zAll_data_shift(:,i) = zAll_data_shift(:,i)*-1;
end

dataset = zAll_data_shift(:,:);
size    = size(dataset,2);

len = nan(2,size);
len(1,:) = 1:size;
for i = 1:size
    len(2,i) = max(find(~isnan(dataset(:,i)))); % Finds the last value in each column that isn't NaN
end
[~,inx] = sort(len(2,:),2,'descend');
zAll_data_2 = dataset(:,inx);
clear len

% Correct time series so that std and mean are same

for i = 2:size
    zAll_data_2(:,i) = std_correct(zAll_data_2(:,1),zAll_data_2(:,i));
end
    
% Step 3: Make the regression model

% Prepare data for model
% Reshape Matrix into vector
X_1 = zAll_data_2';
X_2 = X_1(:);
% Make a similar vector for time
time = repmat(All_data_final(:,1),1,size);
time_1 = time';
time_2 = flipud(time_1(:)-1000); % now gives years as a value of year 1-995
% Find all non-NaN values in X_2 then extract proxies and years
[idx] = find(~isnan(X_2));
prox = X_2(idx);
year = time_2(idx);
Nyrs = max(year);
N = length(prox);

nchains = 3;
datastruct = struct('prox',prox,'N',N,'Nyrs',Nyrs,'year',year);
for i = 1:nchains
    S.int = randn(1);
    S.sigma = rand(1);
    S.sigma_year = rand(1);
    initStructs(i) = S;
end
% The JAGS model is as follows and is saved as a text file model.txt
%
% mod = 'model {
%     for(i in 1:N) {
%       prox[i] ~ dnorm(mu[i], tau)
%       mu[i] <- int + b.year[year[i]]
%     }
% 
%     for(i in 1:Nyrs) {
%       b.year[i] ~ dnorm(0, tau.year)
%     }
% 
%   # Priors
%     tau <- sigma ^ -2
%     sigma ~ dunif(0, 10)
%     tau.year <- sigma_year ^ -2
%     sigma_year ~ dunif(0, 10)
%     int ~ dnorm(0, 0.0001)
% }'

fprintf( 'Running JAGS...\n' );
tic
[samples, stats, structArray] = matjags( ...
    datastruct, ...                     % Observed data   
    fullfile(pwd, 'model.txt'), ...     % File that contains model definition
    initStructs, ...                    % Initial values for latent variables
    'doparallel' , 1, ...               % Parallelization flag
    'nchains', nchains,...              % Number of MCMC chains
    'nburnin', 5000,...                 % Number of burnin steps
    'nsamples', 10000, ...              % Number of samples to extract
    'thin', 1, ...                      % Thinning parameter
    'dic', 1, ...                       % Do the DIC?
    'monitorparams', {'b.year','sigma','sigma_year'}, ...     % List of latent variables to monitor
    'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
    'verbosity' , 1 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 0 );                    % clean up of temporary files?
toc
    
%save('JAGS_out.mat','samples','stats','structArray')

% Scale variance
% All of these are done having removed column 7 from inx
for a = 1:3
    if a == 1
        dataset = recon_FJ(:,2);
    elseif a == 2
        dataset = recon_V(:,2);
    elseif a ==3
        dataset = recon_M(:,2);
    end
    
    dataset(1:inx(2,1),2) = dataset(1:inx(2,1),1);
    for i = 2:41  % This corrects variance, not mean
        if inx(2,i)-inx(2,i-1) == 0
            dataset(inx(2,i-2):inx(2,i),2) = std_correct_CPS(dataset((1:inx(2,1)),1),dataset(inx(2,i-2):inx(2,i),1));
        else
            dataset(inx(2,i-1):inx(2,i),2) = std_correct_CPS(dataset((1:inx(2,1)),1),dataset(inx(2,i-1):inx(2,i),1));
        end
    end
    if a == 1
        recon_FJ(:,3) = dataset(:,2);
    elseif a == 2
        recon_V(:,3) = dataset(:,2);
    elseif a == 3
        recon_M(:,3) = dataset(:,2);
    end
end
save('JAGS_out.mat','recon_FJ','recon_M','recon_V','-append')


% check skill

recon(:,1) = All_data_final(:,1);
recon_H(:,1) = All_data_final(:,1);
recon_L(:,1) = All_data_final(:,1);
recon(:,2) = stats.mean.b_year';
recon_H(:,2) = stats.ci_high.b_year';
recon_L(:,2) = stats.ci_low.b_year';

corr(recon(1:39,2),Marshall_SAM(21:59,2))
corr(recon(1:109,2),Visbeck_Ann(11:119,2))
corr(recon(1:91,2),FJ_ann(10:100,2))


%% Now to compare this with a weighted CPS approach
CPS_SAM = nan(length(zAll_data_shift),2);
CPS_SAM(:,1) = All_data_final(18:size(zAll_data_shift,1),1);
CPS_corrs = corrs_SAM(3,:)'; CPS_corrs = corrs_SAM(1,:)'; % Use Fogt Jones annual correlations or Marshall...
CPS_corrs = skill_FJ'; CPS_corrs = skill_marsh'; % PICK ONE
CPS_data = zAll_data_shift; CPS_data(isnan(CPS_data)) = 0;
for i = 1:length(zAll_data_shift)
    CPS_SAM(i,2) = squeeze(CPS_data(i,:))*CPS_corrs;
end
zSAM = zscore(CPS_SAM(:,2));
zrecon = zscore(recon_FJ(:,2));

corr(CPS_SAM(1:39,2),Marshall_SAM(21:59,2))
corr(CPS_SAM(1:109,2),Visbeck_Ann(11:119,2))
corr(CPS_SAM(1:91,2),FJ_ann(10:100,2))

corr(recon_FJ(1:39,2),Marshall_SAM(21:59,2))
corr(recon_FJ(1:109,2),Visbeck_Ann(11:119,2))
corr(recon_FJ(1:91,2),FJ_ann(10:100,2))

% Basic plots

recon_FJ_71 = smooth(zrecon(:,1),71,'rloess');
CPS_FJ_71 = smooth(zSAM(:,1),71,'rloess');

plot(recon_FJ(:,1),zrecon(:,1),'color','k');
hold on
plot(recon_FJ(:,1),zSAM(:,1),'color','m');
plot(recon_FJ(:,1),recon_FJ_71(:,1),'color','r','linewidth',2); 
plot(recon_FJ(:,1),CPS_FJ_71(:,1),'color','b','linewidth',2);

print -dpdf -painters 'recons_2.pdf'

%% Try an EOF Approach (not used in final analysis)

load JAGS_in.mat
%take a look at how many proxies are informing the reconstruction in each
%year
for i = 1:size(zAll_data_shift,1)
    resolution(i,:) = squeeze(sum(~isnan(zAll_data_shift(i,:))));
end

[EOF,PC,latent,~,var_exp,est_mean] = pca(zAll_data_shift,... % Calculate EOFs of data using als algorithm to account for NaNs
'algorithm','als');

data_interp = PC*EOF' + repmat(est_mean,size(PC,1),1); % Create a new dataset with NaNs interpolated. Result is pretty bad for variables with few values

[EOF2,PC2,latent2,~,var_exp2,est_mean2] = pca(zAll_data_shift,'rows','complete'); % Run again to compare with original

subspace(EOF,EOF2) 







