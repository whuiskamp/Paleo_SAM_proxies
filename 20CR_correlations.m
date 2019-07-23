%% Looking at SAT/ SAM correlations over the 20CR data
% The annual data was created from monthly means, combined with CDO
% Eg: (cdo yearmean prmsl.mon.mean.nc 20CR_SLP_ann.nc)
% 
% Load data
SAT = ncread('20CR_SAT_ann.nc','air')-273.15; % (180,91,24,142) - lon,lat,levels,time
SLP = ncread('20CR_SLP_ann.nc','prmsl'); % (180,91,164) - lon,lat,time
time_SAT = floor(ncread('20CR_SAT_ann.nc','time')/(24*365))+1800; % Convert to years
time_SLP = floor(ncread('20CR_SLP_ann.nc','time')/(24*365))+1800;
lat = ncread('20CR_SAT_ann.nc','lat'); lat_S = lat(61:91,1);
lon = ncread('20CR_SAT_ann.nc','lon');
SAT = squeeze(SAT(:,61:91,1,:)); % Extracts SAT at 1000hPa south of 30S.

% Get everything on a common time-interval, 1871-2000

SLP = SLP(:,:,21:150); SAT = SAT(:,:,1:130);

% Calculate model SAM index
SLP_40 = squeeze(SLP(:,66,:)); SLP_65 = squeeze(mean(SLP(:,78:79,:),2));
SAM = zscore(squeeze(mean(SLP_40,1)) - squeeze(mean(SLP_65,1)));

corrs_Marshall = nan(180,31); % For years 1957-2000
corrs_preMarsh = nan(180,31); % For years 1871-1956

for i = 1:180
    for j = 1:31
        corrs_Marshall(i,j) = squeeze(corr(SAM(1,87:130)',squeeze(SAT(i,j,87:130))));
        corrs_preMarsh(i,j) = squeeze(corr(SAM(1,1:86)',squeeze(SAT(i,j,1:86))));
        if corrs_Marshall(i,j) > 0 & corrs_preMarsh(i,j) <= 0 || corrs_Marshall(i,j) < 0 & corrs_preMarsh(i,j) >= 0
            change(i,j) =  1;
        else change(i,j) = 0;
        end
    end
end

%% Running correlations

% Find significant periodicities in SAM index

wt(SAM); % significant 4 year periodicity and non-significant 8 year. These are both pretty questionable.
windowsize = [4 8 16 32 64];

tic % takes ~ 6mins
for h = 1:size(windowsize,2)
    for i=1:size(SAT,1)
        for j=1:size(SAT,2)
            SAT_SAM_runcorr(i,j,h,:) = movingCorrelation(([squeeze(SAT(i,j,:))'; SAM(1,:)]'),windowsize(:,h),1);
        end
    end
end
toc

SAT_SAM_std = std(SAT_SAM_runcorr,4,'omitnan'); % calculate the standard deviation of the running correlations


%% Plotting
lat_S = double(lat_S);
lon = double(lon);

figure(1)
for i=1:4
    subplot(2,2,i)
    axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -30])
    framem
    gridm
end
load coast
change(change==0) = nan;

subplot(2,2,1)
contourfm(lat_S,lon,corrs_Marshall')
colormap(b2r(-1,1));
plotm(lat,long,'k','linewidth',2)

subplot(2,2,2)
contourfm(lat_S,lon,corrs_preMarsh')
colormap(b2r(-1,1));
plotm(lat,long,'k','linewidth',2)

subplot(2,2,3)
contourfm(lat_S,lon,(corrs_Marshall-corrs_preMarsh)')
colormap(b2r(-1,1));
plotm(lat,long,'k','linewidth',2)

subplot(2,2,4)
contourfm(lat_S,lon,((corrs_Marshall-corrs_preMarsh).*change)')
colormap(b2r(-1,1));
plotm(lat,long,'k','linewidth',2)

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperSize = [20 20];
print('20CR_corrs_multi','-dpdf','-bestfit')


% Matlab is ass, so save as a netCDF and move to cluster for plotting with
% ferret or something.

nccreate('20CR_SAM_corrs.nc','latitude','Dimensions',...
    {'latitude',31},'Format','classic')
nccreate('20CR_SAM_corrs.nc','longitude','Dimensions',...
    {'longitude',180},'Format','classic')

ncwrite('20CR_SAM_corrs.nc','latitude',lat_S(:,1))
ncwrite('20CR_SAM_corrs.nc','longitude',lon(:,1))

ncwriteatt('20CR_SAM_corrs.nc','latitude','axis','Y')
ncwriteatt('20CR_SAM_corrs.nc','latitude','standard_name','latitude')
ncwriteatt('20CR_SAM_corrs.nc','latitude','units','degrees_north')
ncwriteatt('20CR_SAM_corrs.nc','longitude','axis','X')
ncwriteatt('20CR_SAM_corrs.nc','longitude','standard_name','longitude')
ncwriteatt('20CR_SAM_corrs.nc','longitude','units','degrees_east')

nccreate('20CR_SAM_corrs.nc','corrs_Marshall','Datatype','single','Dimensions',...
    {'longitude',180,'latitude',31},'Format','classic')
nccreate('20CR_SAM_corrs.nc','corrs_preMarsh','Datatype','single','Dimensions',...
    {'longitude',180,'latitude',31},'Format','classic')
nccreate('20CR_SAM_corrs.nc','change','Datatype','single','Dimensions',...
    {'longitude',180,'latitude',31},'Format','classic')
            
ncwrite('20CR_SAM_corrs.nc','corrs_Marshall',corrs_Marshall)
ncwrite('20CR_SAM_corrs.nc','corrs_preMarsh',corrs_preMarsh)
ncwrite('20CR_SAM_corrs.nc','change',change)


