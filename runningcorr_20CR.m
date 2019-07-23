% Running correlations of Proxy data with 20CR zonal winds
% This script needs the following packages: b2r, plot2svg, M_map
% Willem Huiskamp, 2015

clear

% load 20CR annual mean zonal winds (850hPa) and proxy data
load('20CR.mat');
load('data_for_correlations.mat','dataMatrix20CR');
load('SAM_20CR.mat');

% For this part, we only want uwind for 1871-1998, so we remove the last 13
% years - this spans 1871-2011 originally

uwind_SF = uwind_SF(1:128,:,:);
uwind_MA = uwind_MA(1:128,:,:);

% Next we need to create 5 evenly spaced calibration windows in the data,
% of the length of our correlation windows.

NUM_YRS = 129;
NUM_CAL_WDW = 5;
lat = 1:46;
lon = 1:180;

proxy_20CR_corrs = nan(3,6,size(NUM_CAL_WDW,1),size(lat,2),size(lon,2)); %3 different window sizes, 6 proxies

for windowsize = [21 31 51]
  
  if windowsize > 20
      clear CAL_WDW
  end
  if windowsize == 21
      w = 1;
  elseif windowsize == 31
      w = 2;
  elseif windowsize == 51
      w = 3;
  end
  
  for c=0:NUM_CAL_WDW-1
       overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/(NUM_CAL_WDW-1));
       CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
  end
     
% Calculate correlations 

    for k = 1:6
        for c = 1:NUM_CAL_WDW
            for i = lat
                for j = lon
                   proxy_20CR_corrs(w,k,c,i,j) = corr(dataMatrix20CR(CAL_WDW(c,:),k), squeeze(ann20CR_1870_1998(j,i,CAL_WDW(c,:))));
%                     if abs(proxy_20CR_corrs(w,k,c,i,j)) < 0.3
%                        proxy_20CR_corrs(w,k,c,i,j) = NaN;
%                     end
                end
            end
        c
        end
    end
    windowsize
end

% Calculate std deviation across each group of windows for each proxy, for
% each window size

std_20CR_corr = nan(3,6,1,size(lat,2),size(lon,2));

for w = 1:3
    for k = 1:6
        for i = lat
            for j = lon
                std_20CR_corr(w,k,:,i,j) = squeeze(nanstd(proxy_20CR_corrs(w,k,:,i,j)));
            end
        end
    end
end

%% Same thing for the seasonal data

clear

load('20CR.mat');
load('data_for_correlations.mat','dataMatrix20CR');
load('SAM_20CR.mat');

% For this part, we only want uwind for 1871-1998, so we remove the last 13
% years - this spans 1871-2011 originally

uwind_SF = uwind_SF(1:128,:,:);
uwind_MA = uwind_MA(1:128,:,:);

NUM_YRS = 128;
NUM_CAL_WDW = 5;
lat = 1:46;
lon = 1:180;

proxy_20CR_SF_corrs = nan(3,6,size(NUM_CAL_WDW,1),size(lat,2),size(lon,2)); %3 different window sizes, 6 proxies
proxy_20CR_MA_corrs = nan(3,6,size(NUM_CAL_WDW,1),size(lat,2),size(lon,2)); %3 different window sizes, 6 proxies

for windowsize = [21 31 51]  % Each window takes 55 seconds to run
  tic 
  if windowsize > 20
      clear CAL_WDW
  end
  if windowsize == 21
      w = 1;
  elseif windowsize == 31
      w = 2;
  elseif windowsize == 51
      w = 3;
  end
  
  for c=0:NUM_CAL_WDW-1
       overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/(NUM_CAL_WDW-1));
       CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
  end
     
% Calculate correlations - note the 1+CAL_WDW as the data matrix starts at
% 1870 but the 20CR data starts at 1871

    for k = 1:6
        for c = 1:NUM_CAL_WDW
            for i = lat
                for j = lon
                   proxy_20CR_SF_corrs(w,k,c,i,j) = corr(dataMatrix20CR(1+CAL_WDW(c,:),k), squeeze(uwind_SF(CAL_WDW(c,:),i,j)));
                   proxy_20CR_MA_corrs(w,k,c,i,j) = corr(dataMatrix20CR(1+CAL_WDW(c,:),k), squeeze(uwind_MA(CAL_WDW(c,:),i,j)));
                     if abs(proxy_20CR_SF_corrs(w,k,c,i,j)) < 0.3
                        proxy_20CR_SF_corrs(w,k,c,i,j) = NaN;
                     end
                     if abs(proxy_20CR_MA_corrs(w,k,c,i,j)) < 0.3
                        proxy_20CR_MA_corrs(w,k,c,i,j) = NaN;
                     end
                end
            end
        c
        end
    end
    windowsize
    toc
end

% Calculate std deviation across each group of windows for each proxy, for
% each window size

std_20CR_SF_corr = nan(3,6,1,size(lat,2),size(lon,2));
std_20CR_MA_corr = nan(3,6,1,size(lat,2),size(lon,2));

for w = 1:3
    for k = 1:6
        for i = lat
            for j = lon
                std_20CR_SF_corr(w,k,:,i,j) = squeeze(nanstd(proxy_20CR_SF_corrs(w,k,:,i,j)));
                std_20CR_MA_corr(w,k,:,i,j) = squeeze(nanstd(proxy_20CR_MA_corrs(w,k,:,i,j)));
            end
        end
    end
end

save('20CR_seasonal_corrs.mat','proxy_20CR_SF_corrs','std_20CR_SF_corr','proxy_20CR_MA_corrs','std_20CR_MA_corr','-append');
%% Plotting!

% Plotting using M_map - Plotting only works if you keep all matlab windows
% on one monitor in earlier matlab versions (WTF)

lats = [-70 -50 -30 -10];

figure
fig_std = tight_subplot(3,6,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for i = 1:18
    axes(fig_std(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

% m_pcolor(lon,lat,squeeze(std_20CR_corr(1,i,:,:,:))); shading flat;
% hold
% m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','k','linewidth',2);
% m_grid('xtick',12,'tickdir','out','ytick',lats,'linest','--');

for a = 1:6
    axes(fig_std(a));
    m_pcolor(lon,lat,squeeze(std_20CR_corr(1,a,:,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','magenta','linewidth',2);
    m_grid('xtick',12,'tickdir','out','ytick',lats,'linest','--');
    caxis manual
    caxis([0 0.8]);
    axes(fig_std(a+6));
    m_pcolor(lon,lat,squeeze(std_20CR_corr(2,a,:,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','magenta','linewidth',2);
    m_grid('xtick',12,'tickdir','out','ytick',lats,'linest','--');
    caxis manual
    caxis([0 0.8]);
    axes(fig_std(a+12));
    m_pcolor(lon,lat,squeeze(std_20CR_corr(3,a,:,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','magenta','linewidth',2);
    m_grid('xtick',12,'tickdir','out','ytick',lats,'linest','--');
    caxis manual
    caxis([0 0.8]);
end


% Set up map
figure
fig_std = tight_subplot(3,6,[0.05 0.01],[0.10 0.01],[0.1 0.01]);
for i = 1:18
    axes(fig_std(i));
    axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -20])
    framem
    load coast
    plotm(lat,long,'m','linewidth',2)
    gridm %If you want gridlines
end

clear lat lon
load('20CR.mat','lat','lon');
lat = double(lat)';
lon = double(lon)';
levels = [0:0.1:0.8];
% Plot Std deviations of correlations across the 5 windows

for i = 1:6
    axes(fig_std(i));
    contourfm(lat,lon,squeeze(std_20CR_corr(1,i,:,:,:)),levels,'linestyle','none')
    caxis manual
    caxis([0 0.8]);
    axes(fig_std(i+6));
    contourfm(lat,lon,squeeze(std_20CR_corr(2,i,:,:,:)),levels,'linestyle','none')
    caxis manual
    caxis([0 0.8]);
    axes(fig_std(i+12));
    contourfm(lat,lon,squeeze(std_20CR_corr(3,i,:,:,:)),levels,'linestyle','none')
    caxis manual
    caxis([0 0.8]);
end

plot2svg('std_corr_20CR.svg')

% work this into that^
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[0 0 3000 3000]);
ti = get(gca,'TightInset');
set(gca,'Position',[ti(1)-0.015 ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
print(gcf, '-depsc', 'map')

% Plot each Proxy's correlation maps

figure
fig_cor = tight_subplot(5,6,[0.05 0.01],[0.10 0.01],[0.1 0.01]);
for i = 1:30
    axes(fig_cor(i));
    axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -20])
    framem
    gridm %If you want gridlines
end

clear lat lon
load('20CR.mat','lat','lon');
lat = double(lat)';
lon = double(lon)';

for i = 1:6  % This will plot for the 31 year window, change the 2 to a 1 or 3 for 21 or 51 year windows.
    axes(fig_cor(i));
    contourfm(lat,lon,squeeze(proxy_20CR_corrs(2,i,1,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
    axes(fig_cor(i+6));
    contourfm(lat,lon,squeeze(proxy_20CR_corrs(2,i,2,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
    axes(fig_cor(i+12));
    contourfm(lat,lon,squeeze(proxy_20CR_corrs(2,i,3,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
    axes(fig_cor(i+18));
    contourfm(lat,lon,squeeze(proxy_20CR_corrs(2,i,4,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
    axes(fig_cor(i+24));
    contourfm(lat,lon,squeeze(proxy_20CR_corrs(2,i,5,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
end
colormap(b2r(-0.6,0.6))

for i = 1:30
    axes(fig_cor(i));
    load coast
    plotm(lat,long,'k','linewidth',2)
end


axes(fig_cor(1));
ylabel('1871-1901')
axes(fig_cor(7));
ylabel('1895-1925')
axes(fig_cor(13));
ylabel('1919-1949')
axes(fig_cor(19));
ylabel('1943-1973')
axes(fig_cor(25));
ylabel('1967-1997')

plot2svg('31yrWdwcorrs.svg')

%% And now the seasonal plot

figure
fig_cor_seas = tight_subplot(5,6,[0.05 0.01],[0.10 0.01],[0.1 0.01]);
for i = 1:30
    axes(fig_cor_seas(i));
    axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -20])
    framem
    gridm %If you want gridlines
end

clear lat lon
load('20CR.mat','lat','lon');
lat = double(lat)';
lon = double(lon)';

for i = 1:6  % This will plot for the 31 year window, change the 2 to a 1 or 3 for 21 or 51 year windows.
    axes(fig_cor_seas(i));
    contourfm(lat,lon,squeeze(proxy_20CR_SF_corrs(2,i,1,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
    axes(fig_cor_seas(i+6));
    contourfm(lat,lon,squeeze(proxy_20CR_SF_corrs(2,i,2,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
    axes(fig_cor_seas(i+12));
    contourfm(lat,lon,squeeze(proxy_20CR_SF_corrs(2,i,3,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
    axes(fig_cor_seas(i+18));
    contourfm(lat,lon,squeeze(proxy_20CR_SF_corrs(2,i,4,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
    axes(fig_cor_seas(i+24));
    contourfm(lat,lon,squeeze(proxy_20CR_SF_corrs(2,i,5,:,:)),'linestyle','none')
     caxis manual
     caxis([-0.6 0.6]);
end
colormap(b2r(-0.6,0.6))

for i = 1:30
    axes(fig_cor_seas(i));
    load coast
    plotm(lat,long,'k','linewidth',2)
end


axes(fig_cor_seas(1));
ylabel('1871-1901')
axes(fig_cor_seas(7));
ylabel('1895-1925')
axes(fig_cor_seas(13));
ylabel('1919-1949')
axes(fig_cor_seas(19));
ylabel('1943-1973')
axes(fig_cor_seas(25));
ylabel('1967-1997')

plot2svg('31yrWdw_seas_corrs.svg')

%% Or with M-Map (which might actually fucking work)

figure
fig_SF = tight_subplot(5,6,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for i = 1:30
    axes(fig_SF(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

clear lat lon
load('20CR.mat','lat','lon');
lat = double(lat)';
lon = double(lon)';

for a = 1:6
    axes(fig_SF(a));
    m_pcolor(lon,lat,squeeze(proxy_20CR_SF_corrs(2,a,1,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([-0.6 0.6]);
    axes(fig_SF(a+6));
    m_pcolor(lon,lat,squeeze(proxy_20CR_SF_corrs(2,a,2,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([-0.6 0.6]); 
    axes(fig_SF(a+12));
    m_pcolor(lon,lat,squeeze(proxy_20CR_SF_corrs(2,a,3,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([-0.6 0.6]); 
    axes(fig_SF(a+18));
    m_pcolor(lon,lat,squeeze(proxy_20CR_SF_corrs(2,a,4,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([-0.6 0.6]); 
    axes(fig_SF(a+24));
    m_pcolor(lon,lat,squeeze(proxy_20CR_SF_corrs(2,a,5,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([-0.6 0.6]); 
end

colormap(b2r(-0.6,0.6))

axes(fig_SF(1));
ylabel('1871-1901')
title('Araucaria')
axes(fig_SF(7));
ylabel('1895-1925')
axes(fig_SF(13));
ylabel('1919-1949')
axes(fig_SF(19));
ylabel('1943-1973')
axes(fig_SF(25));
ylabel('1967-1997')

axes(fig_SF(2));
title('Austrocedrus')
axes(fig_SF(3));
title('Halocarpus')
axes(fig_SF(4));
title('Ice Dust')
axes(fig_SF(5));
title('Nothofagus')
axes(fig_SF(6));
title('Callitris')

letters = 'abcdefghijklmnopqrstuvwxyz';
for i=1:12
    axes(fig_SF(i));
   % set(gca,'YTickLabel',[],'XTickLabel',[])
   % set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
    hold on; plot([0,70],[0.37 0.37],'Color',[0.5 0.5 0.5],'LineWidth',1.5); hold off;
end

plot2svg('31yrWdw_seas_corrs.svg')

%% Seasonal Std Plot
load 20CR_seasonal_corrs.mat

% We leave out the Callitris for this plot as it is not sensitive to
% Spring/Summer precip

figure
fig_SF_std = tight_subplot(3,6,[0.01 0.01],[0.10 0.1],[0.1 0.01]);

for i = 1:18
    axes(fig_SF_std(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

clear lat lon
load('20CR.mat','lat','lon');
lat = double(lat)';
lon = double(lon)';

for a = 1:5
    axes(fig_SF_std(a));
    m_pcolor(lon,lat,squeeze(std_20CR_SF_corr(1,a,:,:,:))); shading flat;
    hold;
    m_coast('color','m','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([0 0.8]);
    axes(fig_SF_std(a+6));
    m_pcolor(lon,lat,squeeze(std_20CR_SF_corr(2,a,:,:,:))); shading flat;
    hold;
    m_coast('color','m','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([0 0.8]); 
    axes(fig_SF_std(a+12));
    m_pcolor(lon,lat,squeeze(std_20CR_SF_corr(3,a,:,:,:))); shading flat;
    hold;
    m_coast('color','m','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([0 0.8]);
end
for a = 1:3
    axes(fig_SF_std(a*6));
    m_pcolor(lon,lat,squeeze(std_20CR_MA_corr(a,6,:,:,:))); shading flat;
    hold;
    m_coast('color','m','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    caxis manual
    caxis([0 0.8]);
end
    

%colormap(b2r(0,0.8))

axes(fig_SF_std(1));
ylabel('21yr Window')
title('Araucaria')
axes(fig_SF_std(7));
ylabel('31yr Window')
axes(fig_SF_std(13));
ylabel('51yr Window')

axes(fig_SF_std(2));
title('Austrocedrus')
axes(fig_SF_std(3));
title('Halocarpus')
axes(fig_SF_std(4));
title('Ice Dust')
axes(fig_SF_std(5));
title('Nothofagus')
axes(fig_SF_std(6));
title('Callitris')

% letters = 'abcdefghijklmnopqrstuvwxyz';
% for i=1:15
%     axes(fig_SF_std(i));
%     text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
%     hold on; plot([0,70],[0.37 0.37],'Color',[0.5 0.5 0.5],'LineWidth',1.5); hold off;
% end

print -dpdf -painters 31yrWdw_seas_std.pdf
plot2svg('31yrWdw_seas_std.svg')

% ######### TO DO ##########
% Change the years (Tree year is the one in which it starts growing)
% Do coherence between visbeck regional SAMs
% Instead of STD plot, just do a difference with the two most recent
% windows, where data is most reliable.


%% Save to netCDF file
% variables are proxy_20CR_corrs(3,6,5,180,46) and
% std_20CR_corrs(3,6,1,180,46)

nccreate('proxycorr_20CR.nc','latitude','Dimensions',...
    {'latitude',46},'Format','classic')
nccreate('proxycorr_20CR.nc','longitude','Dimensions',...
    {'longitude',180},'Format','classic')

ncwrite('proxycorr_20CR.nc','latitude',flipud(lat(1:46,1)))
ncwrite('proxycorr_20CR.nc','longitude',flipud(lon(1:180,1)))

ncwriteatt('proxycorr_20CR.nc','latitude','axis','Y')
ncwriteatt('proxycorr_20CR.nc','latitude','standard_name','latitude')
ncwriteatt('proxycorr_20CR.nc','latitude','units','degrees_north')
ncwriteatt('proxycorr_20CR.nc','longitude','axis','X')
ncwriteatt('proxycorr_20CR.nc','longitude','standard_name','longitude')
ncwriteatt('proxycorr_20CR.nc','longitude','units','degrees_east')

nccreate('proxycorr_20CR.nc','proxy_20CR_corrs','Datatype','single','Dimensions',...
    {'WdwSize',3,'Proxy',6,'Window',5,'latitude',46,'longitude', 180'},'Format','classic')
nccreate('proxycorr_20CR.nc','std_20CR_corr','Datatype','single','Dimensions',...
    {'WdwSize',3,'Proxy',6,'Window',5,'latitude',46,'longitude', 180'},'Format','classic')
            
ncwrite('proxycorr_20CR.nc','proxy_20CR_corrs',proxy_20CR_corrs(1:3,1:6,1:5,1:46,1:180))
ncwrite('proxycorr_20CR.nc','std_20CR_corr',std_20CR_corr(1:3,1:6,1,1:46,1:180))
            
