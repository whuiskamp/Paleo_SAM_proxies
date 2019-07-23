%% Detrending proxies before calibration
%
% This script seeks to apply Empirical Mode Decomposition detrending
% to the paleo-proxies and instrumental indices before calibration in order
% to reduce the impact of the strong positive SAM trend over the satellite
% era. Method is adapted from Wu et al., 2007
% (www.pnas.org/cgi/doi/10.1073/pnas.0701020104). 
% This script requires the matlab function ceemdan.m (updated 2014 version)
% and associated packages. These can be found at 
% http://bioingenieria.edu.ar/grupos/ldnlys/metorres/re_inter.htm#Codigos
%
% Willem Huiskamp, 2019
% 
%   Syntax

% modes=ceemdan(x,Nstd,NR,MaxIter,SNRFlag)
% [modes its]=ceemdan(x,Nstd,NR,MaxIter,SNRFlag)

% I/O Description

% OUTPUT
% modes: contain the obtained modes in a matrix with the rows being the modes        
% its: contain the sifting iterations needed for each mode for each realization (one row for each realization)

% INPUT
% x: signal to decompose
% Nstd: noise standard deviation
% NR: number of realizations
% MaxIter: maximum number of sifting iterations allowed.
% SNRFlag: if equals 1, then the SNR increases for every stage, as in [1].
%          If equals 2, then the SNR is the same for all stages, as in [2]. 

%% Script start
% load data
load('JAGS_in.mat','zAll_data_shift')
load('Fogt_Jones.mat','FJ_ann');load('marshall_SAM.mat','Marshall_SAM')
load('SAM_seasonal.mat','Visbeck_Ann')
Marshall_SAM = flipud(Marshall_SAM);

Marshall_SAM_detr = Marshall_SAM;
FJ_SAM_detr       = FJ_ann;
Visbeck_SAM_detr  = Visbeck_Ann;

Nstd = 0.1; % standard deviation of the white noise added to series'
NR = 500; % Number of realisations 
MaxIter = 5000; % Max number of iterations

%% Detrend and save new SAM Indices
% First detrend the SAM indices themselves
% We only run this over the period 1950-present, over which we'd expect
% anthropogenic impacts.
[modes_M, its_M]=ceemdan(Marshall_SAM(2:59,2),Nstd,NR,MaxIter,1);
[modes_FJ, its_FJ]=ceemdan(FJ_ann(1:55,2),Nstd,NR,MaxIter,1);
[modes_V, its_V]=ceemdan(Visbeck_Ann(1:56,2),Nstd,NR,MaxIter,1);

% Detrend the annual mean SAM indices
% Based on the figures made with the scripts in the section below, we can
% say that the following modes should be removed from each respective index
% to eliminate an anthropogenic-like signal:
% Marshall - mode 5,
% F-J      - mode 5,
% Visbeck  - mode 4, This last one is questionable.

Marshall_SAM_detr(2:59,2) = Marshall_SAM_detr(2:59,2) - modes_M(5,:)';
FJ_SAM_detr(1:55,2)       = FJ_SAM_detr(1:55,2) - modes_FJ(5,:)';
Visbeck_SAM_detr(1:56,2)  = Visbeck_SAM_detr(1:56,2) - modes_V(4,:)';

save('SAM_inxs_detrended.mat','Marshall_SAM_detr','FJ_SAM_detr','Visbeck_SAM_detr','history')

%% Detrend proxy records
%
% we only need to detrend over the calibration period. Using 200 years fits
% all three indexes comfortably
detr_prox = nan(995,52); % proxies, time-series
modes = nan(52,5,43);

for i = 1:52
    start(1,i) = min(find(~isnan(zAll_data_shift(:,i))));
end
p_end = 46; 

[modes, its]=ceemdan(zAll_data_shift(start(:,1):p_end,1),Nstd,NR,MaxIter,1); % Note: Time-series cannot contain NaNs
[modes, its]=ceemdan(zAll_data_shift(4:100,33),Nstd,NR,MaxIter,1);


%% Plotting IMFs
% Plot your IMFs to see which represent the SAM
% 1= Marshall, 2= FJ, 3= V
data = 3;

if data == 1
    modes = modes_M;
    its = its_M;
    orig_data = Marshall_SAM(2:59,2);
    time = Marshall_SAM(2:59,1);
    ti = 'Marshall index IMFs';
elseif data == 2
    modes = modes_FJ;
    its = its_FJ;
    orig_data = FJ_ann(1:55,2);
    time = FJ_ann(1:55,1);
    ti = 'Fogt-Jones index IMFs';
elseif data == 3
    modes = modes_V;
    its = its_V;
    orig_data = Visbeck_Ann(:,2);%orig_data = Visbeck_Ann(1:56,2);
    time = Visbeck_Ann(:,1);%time = Visbeck_Ann(1:56,1);
    ti = 'Visbeck index IMFs';
end

[a, b]=size(modes);

figure;
subplot(a+1,1,1);
plot(time,orig_data); 
ylabel('Proxy/ SAM Index')
title(ti)
set(gca,'xtick',[])
axis tight;

for i=2:a
    subplot(a+1,1,i);
    plot(time,modes(i-1,:));
    ylabel (['IMF ' num2str(i-1)]);
    set(gca,'xtick',[])
    xlim([min(time) max(time)])
end

subplot(a+1,1,a+1)
plot(time,modes(a,:))
ylabel(['IMF ' num2str(a)])
xlim([min(time) max(time)])
xlabel('Year')
if data == 1
    print('IMFs_Marshall_SAM','-dpdf','-fillpage')
    clear modes a b time
elseif data == 2
    print('IMFs_FJ_SAM','-dpdf','-fillpage')
    clear modes a b time
elseif data == 3
    print('IMFs_V_full_SAM','-dpdf','-fillpage')
    clear modes a b time
end
    
figure;
boxplot(its);
print('ITSs_V_SAM','-dpdf','-fillpage')












