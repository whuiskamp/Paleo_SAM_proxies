% Re-shaping data into decadal blocks

% First change Koffman data into annual means
years = floor(Koffman_2014(:,1));
for i = 993:2002
    ind = find(years == i);
    Koffman_annual1(i,2) = nanmean(Koffman_2014(ind,2));
    Koffman_annual1(i,1) = i;
end
Koffman_annual1 = flipud(Koffman_annual1);
Koffman_annual = Koffman_annual1(1:1010,:);

% Next create 10-year averages

dec_mean = nanmean(reshape(Koffman_annual(:,2),10,[]))';

%% 20CR SLP into SAM index - spans 1871-2012

% Change from monthly mean to annual mean (Jan-Dec)

SLP = nc_varget('20CR_SLP_mon.nc','prmsl');
lat = nc_varget('20CR_SLP_mon.nc','lat');
lon = nc_varget('20CR_SLP_mon.nc','lon');

X = reshape(SLP,12,(1704/12),91,180);
SLP_annual = squeeze(nanmean(X,1));
clear X
% Calculating the SAM index

% defined in Gallant et al. 2013
% as the difference between the normalized zonally
% averaged sea level pressure anomalies at 40S and 60S,

slp_40s = squeeze(SLP_annual(:,66,:)); % calculate SLP for 40S and 60S
slp_60s = squeeze(SLP_annual(:,76,:)); 
slp_40_mean = detrend(nanmean(slp_40s,2)); % take the zonal mean and detrend
slp_60_mean = detrend(nanmean(slp_60s,2));
SAM1 = slp_40_mean - slp_60_mean;

for i = 1:length(SAM1)
   SAM2(i) = (SAM1(i)-min(SAM1))/(max(SAM1)-min(SAM1)); % normalise
end

for i = 1:length(SAM2)
    SAM(i) = SAM2(i) - mean(SAM2); % calculate anomalies
end

%% 20CR Sept-Feb and March - Aug Zonal winds @ 850hPa

uwind = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/20CR/20CR_septfeb.nc','uwnd');     % 852,24,91,180
lat   = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/20CR/20CR_septfeb.nc','lat');
lon   = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/20CR/20CR_septfeb.nc','lon');
time = (1871:1:2011);

% the 4th pressure level is 850
uwind = squeeze(uwind(3:848,4,:,:)); % This eliminates the Jan-Feb from 1871 as well as the last 4 months of 2012

X = reshape(uwind,6,(846/6),91,180);
uwind_SF = squeeze(nanmean(X,1));
clear X

save('SAM_20CR.mat','uwind_SF','lat','lon','time','SAM') % Use append flag if you want to add more later

uwind = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/20CR/20CR_MA.nc','uwnd');     % 852,24,91,180
lat   = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/20CR/20CR_MA.nc','lat');
lon   = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/20CR/20CR_MA.nc','lon');
time = (1871:1:2011);

% the 4th pressure level is 850
uwind = squeeze(uwind(1:852,4,:,:)); % Includes 2012 - not really important though.

X = reshape(uwind,6,(852/6),91,180);
uwind_MA = squeeze(nanmean(X,1));
clear X

save('SAM_20CR.mat','uwind_MA','-append')

%% JRA55 Sept-Feb Zonal winds @ 850hPa

clear
uwind_SF = nc_varget('JRA55/JRA55_SF.nc','UGRD_GDS0_ISBL_S123');
time_SF  = floor(nc_varget('JRA55/JRA55_SF.nc','initial_time0_encoded')/1e6);
uwind_MA = nc_varget('JRA55/JRA55_MA.nc','UGRD_GDS0_ISBL_S123');
time_MA  = floor(nc_varget('JRA55/JRA55_MA.nc','initial_time0_encoded')/1e6);

lat = nc_varget('JRA55/JRA55_SF.nc','g0_lat_1');
lon = nc_varget('JRA55/JRA55_SF.nc','g0_lon_2');

uwind_SF = uwind_SF(3:254,:,:); % Eliminates the Jan-Feb from 1958 and last 4 months of 2000
uwind_MA = uwind_MA(1:252,:,:); % Removes the year 2000

X = reshape(uwind_SF,6,(252/6),73,288);
uwind_SF = squeeze(nanmean(X,1));
clear X

X = reshape(uwind_MA,6,(252/6),73,288);
uwind_MA = squeeze(nanmean(X,1));
clear X

save('JR55.mat','uwind_SF','uwind_MA','lat','lon','time_SF','time_MA') % Use append flag if you want to add more later

%% ERA-Interim Sept-Feb Zonal winds @ 850hPa

clear
uwind_SF = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/ERA-Interim/ERA_SF.nc','u');
time_SF  = floor(nc_varget('../../../../../../../srv/ccrc/data40/z3215716/ERA-Interim/ERA_SF.nc','time')/8760 + 1900);
uwind_MA = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/ERA-Interim/ERA_MA.nc','u');
time_MA  = floor(nc_varget('../../../../../../../srv/ccrc/data40/z3215716/ERA-Interim/ERA_MA.nc','time')/8760 + 1900);

lat = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/ERA-Interim/era_int_uwind.nc','latitude');
lon = nc_varget('../../../../../../../srv/ccrc/data40/z3215716/ERA-Interim/era_int_uwind.nc','longitude');

uwind_SF = uwind_SF(3:128,:,:); % Eliminates the Jan-Feb from 1979 and last 4 months of 2000 + everything after
uwind_MA = uwind_MA(1:126,:,:); % Removes the year 2000-2015

X = reshape(uwind_SF,6,(126/6),721,1440);
uwind_SF = squeeze(nanmean(X,1));
clear X

X = reshape(uwind_MA,6,(126/6),721,1440);
uwind_MA = squeeze(nanmean(X,1));
clear X

save('ERA_INT_seas.mat','uwind_SF','uwind_MA','lat','lon','time_SF','time_MA') % Use append flag if you want to add more later
