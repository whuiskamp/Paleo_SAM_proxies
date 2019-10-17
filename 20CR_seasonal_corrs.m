%% 20CR seasonal correlations
% Import 20CR data and extract only DJF means
clear
SAT = ncread('air.mon.mean.nc','air')-273.15; % (180,91,24,1704) - lon,lat,levels,time
SAT = squeeze(SAT(:,:,1,:)); % Extract only surface values
SLP = ncread('prmsl.mon.mean.nc','prmsl'); % (180,91,164) - lon,lat,time
time = floor(ncread('air.mon.mean.nc','time')/(24*365))+1800;
time = time(1:12:end,:);
lat = ncread('prmsl.mon.mean.nc','lat'); 
lon = ncread('prmsl.mon.mean.nc','lon');
lmask = ncread('20CR_Lmask.nc','LMASK'); % This land mask is NOT exact. It is extracted from the etopo5
% orography file that comes with ferret, regridded onto the 20CR grid. It
% should be good enough for our purposes.

seas = ["DJF","MAM","JJA","SON"];
dima = [.15 .58 .3 .3];
% Extract only one season
for a = seas
    if a == 'DJF'
        start = 12;
    elseif a == 'MAM'
        start = 3;
    elseif a == 'JJA'
        start = 6;
    elseif a == 'SON'
        start = 9;
    end

    for i = start:12:size(SAT,3)-1 % The dataset ends with Dec, not Jan.
        SAT_seas(:,:,ceil(i/12)) = mean(SAT(:,:,i:i+2),3); % Averaging over season, and time is the 3rd dim
        SLP_seas(:,:,ceil(i/12)) = mean(SLP(:,:,i:i+2),3);
    end

    SLP_40 = squeeze(SLP_seas(:,66,:)); SLP_65 = squeeze(mean(SLP_seas(:,78:79,:),2));
    SAM_seas = zscore(squeeze(mean(SLP_40,1)) - squeeze(mean(SLP_65,1)));

    % Mask out the ocean
    lmask2 = nan(180,91,141);
    for i=1:141
        lmask2(:,:,i) = lmask;
    end
    SAT_seas(lmask2==0) = NaN;

    % Create regions for averaging 
    % 1 = Antarctica, 2 = Southern Aus, 3 = NZ, 4 = Sth. America (N), 5 = Sth. America (S)
    
    for region = 1:5
        if region == 1;
            nth = 16;
            sth = 1;
            wst = 1;
            est = 180;
        elseif region == 2;
            nth = 31;
            sth = 25;
            wst = 58;
            est = 78;
        elseif region == 3;
            nth = 28;
            sth = 23;
            wst = 85;
            est = 90;
        elseif region == 4;
            nth = 31;
            sth = 26;
            wst = 144;
            est = 155;
        elseif region == 5;
            nth = 25;
            sth = 19;
            wst = 144;
            est = 155;
        end

        SAT_region = squeeze(nanmean(nanmean(SAT_seas(wst:est,sth:nth,:),1),2));
        % Next, sort your data into years with positive SAM or Negative SAM
        % then calculate regression values for each using Spearman's rank
        % correlation

        [SAM_sort,inx] = sort(SAM_seas,2,'descend');
        SAT_sort = SAT_region(inx);
        % Find where SAM goes from + to -
        chng = find(SAM_sort > 0,1, 'last');
        % regression for +tive SAM
        P = polyfit(SAM_sort(1:chng)',SAT_sort(1:chng),1);
        yfit = P(1)*SAM_sort(1:chng)+P(2);
        % regression for -tive SAM
        Q = polyfit(SAM_sort(chng+1:end)',SAT_sort(chng+1:end),1);
        yfit_2 = Q(1)*SAM_sort(chng+1:end)+Q(2);

        figure(region)
        subplot(2,2,find(contains(seas,a)))
        scatter(SAM_seas,SAT_region,'filled'); hold on
        line([0 0],[min(SAT_region) max(SAT_region)], 'linestyle','--','color','k')
        plot(SAM_sort(1:chng),yfit,'r-'); % Regression for +-tive SAM
        plot(SAM_sort(chng+1:end),yfit_2,'b-'); % Regression for -tive SAM

        xlabel('SAM index')
        ylabel('SAT')
        if region == 1
            ann = annotation('textbox',dima,'String','Antarctica','FitBoxToText','on');
        elseif region == 2
            ann = annotation('textbox',dima,'String','Southern Aus.','FitBoxToText','on');
        elseif region == 3
            ann = annotation('textbox',dima,'String','New Zealand','FitBoxToText','on');
        elseif region == 4
            ann = annotation('textbox',dima,'String','Sth. American (30S-42S)','FitBoxToText','on');
        elseif region == 5
            ann = annotation('textbox',dima,'String','Sth. American (42S-60S)','FitBoxToText','on');
        end
        ann.LineStyle = 'none';
        if a == "DJF"
            title('DJF')
        elseif a == "MAM"
            title('MAM')
        elseif a == "JJA"
            title('JJA')
        elseif a == "SON"
            title('SON')
        end
    end
    clear SAT_seas SLP_seas
    % Make a S.H. map showing ratios of regression slopes for SAT-SAM
    [SAM_sort,inx] = sort(SAM_seas,2,'descend');
    % Find where SAM goes from + to -
    chng = find(SAM_sort > 0,1, 'last');
    SAT_sort_SH = zeros(180,91,141);
    P_SAM = nan(180,91,2); N_SAM = nan(180,91,2);
    yfit  = nan(180,91,65); yfit_2  = nan(180,91,76);
    Rho_P = nan(180,91); Rho_N = nan(180,91);
    for i = 1:180
        for j = 1:91
            if isnan(squeeze(SAT_seas(i,j,1))) == 0
                %tmp = squeeze(SAT_seas(i,j,inx));
                SAT_sort_SH(i,j,:) = squeeze(SAT_seas(i,j,inx));
                % regression for +tive SAM
                P_SAM(i,j,:) = polyfit(SAM_sort(1:chng)',squeeze(SAT_sort_SH(i,j,1:chng)),1);
                yfit(i,j,:) = P_SAM(i,j,1)*SAM_sort(1:chng)+P_SAM(i,j,2);
                Rho_P(i,j,:) = corr(SAM_sort(1:chng)',squeeze(SAT_sort_SH(i,j,1:chng)),'Type','Spearman');
                % regression for -tive SAM
                N_SAM(i,j,:) = polyfit(SAM_sort(chng+1:end)',squeeze(SAT_sort_SH(i,j,chng+1:end)),1);
                yfit_2(i,j,:) =  N_SAM(i,j,1)*SAM_sort(chng+1:end)+N_SAM(i,j,2);
                Rho_N(i,j,:) = corr(SAM_sort(chng+1:end)',squeeze(SAT_sort_SH(i,j,chng+1:end)),'Type','Spearman');
            end
        end
    end
    i=1;j=1;
%end

% First plotting option: If the SAM-SAT correlation in +tive years is different in sign to
% the correlaton in -tive years, OR the slope of the regression is more
% than twice as big in positive years, plot the value for r in +tive SAM
% years
tmp = nan(180,91);
for i = 1:180
    for j = 1:91
        if Rho_P(i,j) > 0 & Rho_N(i,j) < 0 || P_SAM(i,j,2) > 2*N_SAM(i,j,2)
            tmp(i,j) = Rho_P(i,j);
        end
    end
end

lat_S = double(lat);
lon_S = double(lon);
axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -30])
framem
gridm
load coast

contourfm(lat_S,lon_S,(flipud(tmp')))
colormap(b2r(-1,1));
plotm(lat,long,'k','linewidth',2)
title('Spearmans Rho for DJF SAT and SAM in +tive SAM years')


%% Examining the response of real proxies to SAM
clear
% Load in DJF SAM indices 
load('Fogt_Jones.mat','FogtJones_DJF')
FogtJones_DJF = flipud(FogtJones_DJF);
load('marshall_SAM.mat','Marshall_SAM')
Marshall_DJF(:,[1 2]) = flipud(Marshall_SAM(:,[1 6]));
load('SAM_seasonal.mat','Visbeck_DJF')
Visbeck_DJF = flipud(Visbeck_DJF);
% Load our proxies
load('JAGS_in.mat','zAll_data_shift') % Dataset starts at 1995

for i = 1:52
    start(1,i) = min(find(~isnan(zAll_data_shift(:,i))));
    End_ma(1,i) = 39;
    End_vbk(1,i) = 109;
    End_FJ(1,i) = 91;
    ma_start(1,i) = 20 + start(1,i);
    ma_end(1,i) = 59;
    vbk_start(1,i) = 9 + start(1,i);
    vbk_end(1,i) = 118;
    FJ_start(1,i) = 9 + start(1,i);
    FJ_end(1,i) = 100;
end
End_vbk(1,5) = 103;
vbk_end(1,5) = 112;

M_sort = nan(59,52); inx_M = nan(59,52);
M_sort = nan(59,52); inx_M = nan(59,52);
FJ_sort = nan(100,52); inx_FJ = nan(100,52);
FJ_sort = nan(100,52); inx_FJ = nan(100,52);
V_sort = nan(118,52); inx_V = nan(118,52);
V_sort = nan(118,52); inx_V = nan(118,52);

for i = 1:52
    [M_sort(ma_start(1,i):ma_end(1,i),i),inx_M(ma_start(1,i):ma_end(1,i),i)] = sort(Marshall_DJF(ma_start(1,i):ma_end(1,i),2),1,'descend');
    [FJ_sort(FJ_start(1,i):FJ_end(1,i),i),inx_FJ(FJ_start(1,i):FJ_end(1,i),i)] = sort(FogtJones_DJF(FJ_start(1,i):FJ_end(1,i),2),1,'descend');
    [V_sort(vbk_start(1,i):vbk_end(1,i),i),inx_V(vbk_start(1,i):vbk_end(1,i),i)] = sort(Visbeck_DJF(vbk_start(1,i):vbk_end(1,i),2),1,'descend');
    % Find where SAM goes from + to -
    chng_M(:,i) = find(M_sort(:,i) > 0,1, 'last');
    chng_FJ(:,i) = find(FJ_sort(:,i) > 0,1, 'last');
    chng_V(:,i) = find(V_sort(:,i) > 0,1, 'last');
end

proxies_sort_M = nan(39,52); proxies_sort_FJ = nan(91,52); proxies_sort_V = nan(109,52);
P_M = nan(52,2); P_FJ = nan(52,2); P_V = nan(52,2);
N_M = nan(52,2); N_FJ = nan(52,2); N_V = nan(52,2);
Pfit_M = nan(52,2);Pfit_FJ = nan(52,2);Pfit_V = nan(52,2);
Nfit_M = nan(52,2);Nfit_FJ = nan(52,2);Nfit_V = nan(52,2);
yfit_M = nan(52,39); yfit_FJ= nan(52,91); yfit_V = nan(52,109);
yfitN_M = nan(52,39); yfitN_FJ= nan(52,91); yfitN_V = nan(52,109);

for i = 1:52
    tmp_M = zAll_data_shift(start(1,i):End_ma(1,i),i);
    proxies_sort_M(start(1,i):End_ma(1,i),i) = tmp_M(inx_M(ma_start(1,i):ma_end(1,i),i));
    tmp_FJ = zAll_data_shift(start(1,i):End_FJ(1,i),i);
    proxies_sort_FJ(start(1,i):End_FJ(1,i),i) = tmp_FJ(inx_FJ(FJ_start(1,i):FJ_end(1,i),i));
    tmp_V = zAll_data_shift(start(1,i):End_vbk(1,i),i);
    proxies_sort_V(start(1,i):End_vbk(1,i),i) = tmp_V(inx_V(vbk_start(1,i):vbk_end(1,i),i));
    clear tmp_M tmp_FJ tmp_V

    % regression for +tive SAM
    Pfit_M(i,:) = polyfit(M_sort(ma_start(1,i):chng_M(:,i),i),proxies_sort_M(start(:,i):chng_M(:,i)-20,i),1);
    yfit_M(i,start(:,i):chng_M(:,i)-20) = Pfit_M(i,1)*M_sort(ma_start(1,i):chng_M(:,i),i)+Pfit_M(i,2);
    Rho_P_M(i,:) = corr(M_sort(ma_start(1,i):chng_M(:,i),i),proxies_sort_M(start(:,i):chng_M(:,i)-20,i),'Type','Spearman');
    
    Pfit_FJ(i,:) = polyfit(FJ_sort(FJ_start(1,i):chng_FJ(:,i),i),proxies_sort_FJ(start(:,i):chng_FJ(:,i)-9,i),1);
    yfit_FJ(i,start(:,i):chng_FJ(:,i)-9) = Pfit_FJ(i,1)*FJ_sort(FJ_start(1,i):chng_FJ(:,i),i)+Pfit_FJ(i,2);
    Rho_P_FJ(i,:) = corr(FJ_sort(FJ_start(1,i):chng_FJ(:,i),i),proxies_sort_FJ(start(:,i):chng_FJ(:,i)-9,i),'Type','Spearman');
    
    Pfit_V(i,:) = polyfit(V_sort(vbk_start(1,i):chng_V(:,i),i),proxies_sort_V(start(:,i):chng_V(:,i)-9,i),1);
    yfit_V(i,start(:,i):chng_V(:,i)-9) = Pfit_V(i,1)*V_sort(vbk_start(1,i):chng_V(:,i),i)+Pfit_V(i,2);
    Rho_P_V(i,:) = corr(V_sort(vbk_start(1,i):chng_V(:,i),i),proxies_sort_V(start(:,i):chng_V(:,i)-9,i),'Type','Spearman');
    
    % regression for -tive SAM
    Nfit_M(i,:) = polyfit(M_sort(chng_M(:,i)+1:ma_end(1,i),i),proxies_sort_M(chng_M(:,i)-20+1:End_ma(1,i),i),1);
    yfitN_M(i,chng_M(:,i)-20+1:End_ma(1,i)) = Nfit_M(i,1)*M_sort(chng_M(:,i)+1:ma_end(1,i),i)+Nfit_M(i,2);
    Rho_N_M(i,:) = corr(M_sort(chng_M(:,i)+1:ma_end(1,i),i),proxies_sort_M(chng_M(:,i)-20+1:End_ma(1,i),i),'Type','Spearman');
    
    Nfit_FJ(i,:) = polyfit(FJ_sort(chng_FJ(:,i)+1:FJ_end(1,i),i),proxies_sort_FJ(chng_FJ(:,i)-9+1:End_FJ(1,i),i),1); 
    yfitN_FJ(i,chng_FJ(:,i)-9+1:End_FJ(1,i)) = Nfit_FJ(i,1)*FJ_sort(chng_FJ(:,i)+1:FJ_end(1,i),i)+Nfit_FJ(i,2);
    Rho_N_FJ(i,:) = corr(FJ_sort(chng_FJ(:,i)+1:FJ_end(1,i),i),proxies_sort_FJ(chng_FJ(:,i)-9+1:End_FJ(1,i),i),'Type','Spearman');
    
    Nfit_V(i,:) = polyfit(V_sort(chng_V(:,i)+1:vbk_end(1,i),i),proxies_sort_V(chng_V(:,i)-9+1:End_vbk(1,i),i),1);
    yfitN_V(i,chng_V(:,i)-9+1:End_vbk(1,i)) = Nfit_V(i,1)*V_sort(chng_V(:,i)+1:vbk_end(1,i),i)+Nfit_V(i,2);
    Rho_N_V(i,:) = corr(V_sort(chng_V(:,i)+1:vbk_end(1,i),i),proxies_sort_V(chng_V(:,i)-9+1:End_vbk(1,i),i),'Type','Spearman');
end

% plot results
proxies=zAll_data_shift;
figure(1)
for i=1:52
    subplot(5,11,i)
    s = min(proxies(start(1,i):End_ma(1,i),i))-0.5;
    e = max(proxies(start(1,i):End_ma(1,i),i))+0.5;
    scatter(Marshall_DJF(ma_start(1,i):ma_end(1,i),2),proxies(start(1,i):End_ma(1,i),i),'filled'); hold on
    line([0 0],[s e], 'linestyle','--','color','k')
    plot(M_sort(ma_start(1,i):chng_M(:,i),i),yfit_M(i,start(:,i):chng_M(:,i)-20),'r-'); % Regression for +-tive SAM
    plot(M_sort(chng_M(:,i)+1:ma_end(1,i),i),yfitN_M(i,chng_M(:,i)-20+1:End_ma(1,i)),'b-'); % Regression for -tive SAM
    ylim([s e])
end
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('proxies_Marshall_regressions','-dpdf','-besfit')

figure(2)
for i=1:52
    subplot(5,11,i)
    s = min(proxies(start(1,i):End_FJ(1,i),i))-0.5;
    e = max(proxies(start(1,i):End_FJ(1,i),i))+0.5;
    scatter(FogtJones_DJF(FJ_start(1,i):FJ_end(1,i),2),proxies(start(1,i):End_FJ(1,i),i),'filled'); hold on
    line([0 0],[s e], 'linestyle','--','color','k')
    plot(FJ_sort(FJ_start(1,i):chng_FJ(:,i),i),yfit_FJ(i,start(:,i):chng_FJ(:,i)-9),'r-'); % Regression for +-tive SAM
    plot(FJ_sort(chng_FJ(:,i)+1:FJ_end(1,i),i),yfitN_FJ(i,chng_FJ(:,i)-9+1:End_FJ(1,i)),'b-'); % Regression for -tive SAM
    ylim([s e])
end
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('proxies_FJ_regressions','-dpdf','-besfit')

figure(3)
for i=1:52
    subplot(5,11,i)
    s = min(proxies(start(1,i):End_vbk(1,i),i))-0.5;
    e = max(proxies(start(1,i):End_vbk(1,i),i))+0.5;
    scatter(Visbeck_DJF(vbk_start(1,i):vbk_end(1,i),2),proxies(start(1,i):End_vbk(1,i),i),'filled'); hold on
    line([0 0],[s e], 'linestyle','--','color','k')
    plot(V_sort(vbk_start(1,i):chng_V(:,i),i),yfit_V(i,start(:,i):chng_V(:,i)-9),'r-'); % Regression for +-tive SAM
    plot(V_sort(chng_V(:,i)+1:vbk_end(1,i),i),yfitN_V(i,chng_V(:,i)-9+1:End_vbk(1,i)),'b-'); % Regression for -tive SAM
    ylim([s e])
end
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('proxies_Visbeck_regressions','-dpdf','-besfit')


































