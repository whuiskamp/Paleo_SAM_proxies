%% Spatial correlations
clear
ncid = netcdf.open('JRA55_uwind_850hPa.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

load('EOF1_2.mat')
lat = g0_lat_1;
lon = g0_lon_2;
uwind = UGRD_GDS0_ISBL_S123;

% To simply plot wind field;
%contourf(g0_lon_2,g0_lat_1,squeeze(UGRD_GDS0_ISBL_S123(:,:,1)'))

%Compute annual means

X = reshape(uwind,144,37,12,(660/12));
uwind_annual = squeeze(nanmean(X,3));

%% JRA55

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(EOF1(343:397,2),squeeze(uwind_annual(j,i,1:55)));
     EOF1_corr(i,j)=corr(1,2); 
     EOF1_p(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
%      [corr,pval]=corrcoef(EOF1(378:397,2),squeeze(uwind_annual(j,i,22:41))); % For 1979-1998
%      EOF1_corr_79(i,j)=corr(1,2);                                            % For 1979-1998
%      EOF1_p_79(i,j)=pval(1,2);                                               % For 1979-1998
  end
end

[row,col] = find(EOF1_p>0.05);
EOF1_corr2 = EOF1_corr;

for i = 1:length(col)
    EOF1_corr2(row(i),col(i)) = NaN;
end
% 
% [row,col] = find(EOF1_p_79>0.05);
% EOF1_corr2_79 = EOF1_corr_79;
% 
% for i = 1:length(col)
%     EOF1_corr2_79(row(i),col(i)) = NaN;
% end

% figure(1)
% contourf(lon,lat,EOF1_corr,80,'linestyle','none')
% figure(2)
% contourf(lon,lat,EOF1_p,80,'linestyle','none')
% figure(3)
% contourf(lon,lat,EOF1_corr2,80,'linestyle','none')

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(EOF2(343:397,2),squeeze(uwind_annual(j,i,1:55)));
     EOF2_corr(i,j)=corr(1,2); 
     EOF2_p(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
%      [corr,pval]=corrcoef(EOF2(378:397,2),squeeze(uwind_annual(j,i,22:41))); % For 1979-1998
%      EOF2_corr_79(i,j)=corr(1,2);                                            % For 1979-1998
%      EOF2_p_79(i,j)=pval(1,2);                                               % For 1979-1998
  end
end

[row,col] = find(EOF2_p>0.05);
EOF2_corr2 = EOF2_corr;

for i = 1:length(col)
    EOF2_corr2(row(i),col(i)) = NaN;
end
% 
% 
% [row,col] = find(EOF2_p_79>0.05);
% EOF2_corr2_79 = EOF2_corr_79;
% 
% for i = 1:length(col)
%     EOF2_corr2_79(row(i),col(i)) = NaN;
% end

% figure(4)
% contourf(lon,lat,EOF2_corr,80,'linestyle','none')
% figure(5)
% contourf(lon,lat,EOF2_p,80,'linestyle','none')
% figure(6)
% contourf(lon,lat,EOF2_corr2,80,'linestyle','none')

nccreate('EOF_JRA_corr.nc','latitude','Dimensions',...
    {'latitude',37},'Format','classic')
nccreate('EOF_JRA_corr.nc','longitude','Dimensions',...
    {'longitude',144},'Format','classic')
ncwrite('EOF_JRA_corr.nc','latitude',flipud(lat(1:37,1)))
ncwrite('EOF_JRA_corr.nc','longitude',flipud(lon(1:144,1)))

ncwriteatt('EOF_JRA_corr.nc','latitude','axis','Y')
ncwriteatt('EOF_JRA_corr.nc','latitude','standard_name','latitude')
ncwriteatt('EOF_JRA_corr.nc','latitude','units','degrees_north')
ncwriteatt('EOF_JRA_corr.nc','longitude','axis','X')
ncwriteatt('EOF_JRA_corr.nc','longitude','standard_name','longitude')
ncwriteatt('EOF_JRA_corr.nc','longitude','units','degrees_east')

nccreate('EOF_JRA_corr.nc','EOF1_corr','Datatype','single','Dimensions',...
    {'latitude',37','longitude', 144},'Format','classic')
nccreate('EOF_JRA_corr.nc','EOF2_corr','Datatype','single','Dimensions',...
    {'latitude',37','longitude', 144},'Format','classic')

% nccreate('EOF_JRA_corr.nc','EOF1_corr_79','Datatype','single','Dimensions',...
%     {'latitude',37','longitude', 144},'Format','classic')
% nccreate('EOF_JRA_corr.nc','EOF2_corr_79','Datatype','single','Dimensions',...
%     {'latitude',37','longitude', 144},'Format','classic')

ncwrite('EOF_JRA_corr.nc','EOF1_corr',EOF1_corr2(1:37,1:144))
ncwrite('EOF_JRA_corr.nc','EOF2_corr',EOF2_corr2(1:37,1:144))

% ncwrite('EOF_JRA_corr.nc','EOF1_corr_79',EOF1_corr2_79(1:37,1:144))
% ncwrite('EOF_JRA_corr.nc','EOF2_corr_79',EOF2_corr2_79(1:37,1:144))

ncwriteatt('EOF_JRA_corr.nc','EOF1_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_JRA_corr.nc','EOF1_corr','long_name','Correlation of EOF1 with JRA-55 at 95% significance level, 1955-1998')
ncwriteatt('EOF_JRA_corr.nc','EOF2_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_JRA_corr.nc','EOF2_corr','long_name','Correlation of EOF2 with JRA-55 at 95% significance level, 1955-1998')

% ncwriteatt('EOF_JRA_corr.nc','EOF1_corr_79','standard_name','Correlation at 95% sig.')
% ncwriteatt('EOF_JRA_corr.nc','EOF1_corr_79','long_name','Correlation of EOF1 with JRA-55 at 95% significance level, 1979-1998')
% ncwriteatt('EOF_JRA_corr.nc','EOF2_corr_79','standard_name','Correlation at 95% sig.')
% ncwriteatt('EOF_JRA_corr.nc','EOF2_corr_79','long_name','Correlation of EOF2 with JRA-55 at 95% significance level, 1979-1998')

%% NCEP2

load('NCEP2.mat')
load('EOF_1_2_3.mat')

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(EOF1(378:397,2),squeeze(NCEP2_79_98(j,i,1:20)));
     EOF1_corr_NCEP2(i,j)=corr(1,2); 
     EOF1_p_NCEP2(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF1_p_NCEP2>0.05);
EOF1_corr2_NCEP2 = EOF1_corr_NCEP2;

for i = 1:length(col)
    EOF1_corr2_NCEP2(row(i),col(i)) = NaN;
end

% figure(10)
% contourf(lon,lat,EOF1_corr_NCEP2,80,'linestyle','none')
% figure(11)
% contourf(lon,lat,EOF1_p_NCEP2,80,'linestyle','none')
% figure(12)
% contourf(lon,lat,EOF1_corr2_NCEP2,80,'linestyle','none')

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(EOF2(378:397,2),squeeze(NCEP2_79_98(j,i,1:20)));
     EOF2_corr_NCEP2(i,j)=corr(1,2); 
     EOF2_p_NCEP2(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF2_p_NCEP2>0.05);
EOF2_corr2_NCEP2 = EOF2_corr_NCEP2;

for i = 1:length(col)
    EOF2_corr2_NCEP2(row(i),col(i)) = NaN;
end

% figure(13)
% contourf(lon,lat,EOF2_corr_NCEP2,80,'linestyle','none')
% figure(14)
% contourf(lon,lat,EOF2_p_NCEP2,80,'linestyle','none')
% figure(15)
% contourf(lon,lat,EOF2_corr2_NCEP2,80,'linestyle','none')

% figure(16)
% contourf(lon,lat,EOF3_corr_NCEP2,80,'linestyle','none')
% figure(17)
% contourf(lon,lat,EOF3_p_NCEP2,80,'linestyle','none')
% figure(18)
% contourf(lon,lat,EOF3_corr2_NCEP2,80,'linestyle','none')

nccreate('EOF_NCEP2_corr.nc','latitude','Dimensions',...
    {'latitude',37},'Format','classic')
nccreate('EOF_NCEP2_corr.nc','longitude','Dimensions',...
    {'longitude',144},'Format','classic')
ncwrite('EOF_NCEP2_corr.nc','latitude',flipud(lat(1:37,1)))
ncwrite('EOF_NCEP2_corr.nc','longitude',flipud(lon(1:144,1)))

ncwriteatt('EOF_NCEP2_corr.nc','latitude','axis','Y')
ncwriteatt('EOF_NCEP2_corr.nc','latitude','standard_name','latitude')
ncwriteatt('EOF_NCEP2_corr.nc','latitude','units','degrees_north')
ncwriteatt('EOF_NCEP2_corr.nc','longitude','axis','X')
ncwriteatt('EOF_NCEP2_corr.nc','longitude','standard_name','longitude')
ncwriteatt('EOF_NCEP2_corr.nc','longitude','units','degrees_east')

nccreate('EOF_NCEP2_corr.nc','EOF1_corr','Datatype','single','Dimensions',...
    {'latitude',37,'longitude', 144},'Format','classic')
nccreate('EOF_NCEP2_corr.nc','EOF2_corr','Datatype','single','Dimensions',...
    {'latitude',37,'longitude', 144},'Format','classic')

ncwrite('EOF_NCEP2_corr.nc','EOF1_corr',EOF1_corr2_NCEP2(1:37,1:144))
ncwrite('EOF_NCEP2_corr.nc','EOF2_corr',EOF2_corr2_NCEP2(1:37,1:144))

ncwriteatt('EOF_NCEP2_corr.nc','EOF1_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_NCEP2_corr.nc','EOF1_corr','long_name','Correlation of EOF1 with NCEP2 at 95% significance level')
ncwriteatt('EOF_NCEP2_corr.nc','EOF2_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_NCEP2_corr.nc','EOF2_corr','long_name','Correlation of EOF2 with NCEP2 at 95% significance level')

%% ERA-Interim
clear
ncid = netcdf.open('ERA-Interim/Uwind850_ERA_Int_79_98.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

load('EOF1_2.mat')
load('All_proxies_commonT.mat')

%Compute annual means

X = reshape(u,480,121,12,(240/12));
uwind_annual = squeeze(nanmean(X,3));

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(EOF1(378:397,2),squeeze(uwind_annual(j,i,1:20)));
     EOF1_corr_ERA(i,j)=corr(1,2); 
     EOF1_p_ERA(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF1_p_ERA>0.05);
EOF1_corr2_ERA = EOF1_corr_ERA;

for i = 1:length(col)
    EOF1_corr2_ERA(row(i),col(i)) = NaN;
end

% figure(19)
% contourf(longitude,latitude,EOF1_corr_ERA,80,'linestyle','none')
% figure(20)
% contourf(longitude,latitude,EOF1_p_ERA,80,'linestyle','none')
% figure(21)
% contourf(longitude,latitude,EOF1_corr2_ERA,80,'linestyle','none')

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(EOF2(378:397,2),squeeze(uwind_annual(j,i,1:20)));
     EOF2_corr_ERA(i,j)=corr(1,2); 
     EOF2_p_ERA(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF2_p_ERA>0.05);
EOF2_corr2_ERA = EOF2_corr_ERA;

for i = 1:length(col)
    EOF2_corr2_ERA(row(i),col(i)) = NaN;
end

% figure(22)
% contourf(longitude,latitude,EOF2_corr_ERA,80,'linestyle','none')
% figure(23)
% contourf(longitude,latitude,EOF2_p_ERA,80,'linestyle','none')
% figure(24)
% contourf(longitude,latitude,EOF2_corr2_ERA,80,'linestyle','none')

% Correlations for individual series which show coherence in wavelets

% Araucaria (South America)

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(NArauca(378:397,2),squeeze(uwind_annual(j,i,1:20)));
     Arauca_corr_ERA(i,j)=corr(1,2); 
     Arauca_p_ERA(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Arauca_p_ERA>0.05);
Arauca_corr2_ERA = Arauca_corr_ERA;

for i = 1:length(col)
    Arauca_corr2_ERA(row(i),col(i)) = NaN;
end

% Austrocedrus (South America)

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(NAustro(378:397,2),squeeze(uwind_annual(j,i,1:20)));
     Austro_corr_ERA(i,j)=corr(1,2); 
     Austro_p_ERA(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Austro_p_ERA>0.05);
Austro_corr2_ERA = Austro_corr_ERA;

for i = 1:length(col)
    Austro_corr2_ERA(row(i),col(i)) = NaN;
end

% Nothofagus (South America)

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(NNothof(378:397,2),squeeze(uwind_annual(j,i,1:20)));
     Nothof_corr_ERA(i,j)=corr(1,2); 
     Nothof_p_ERA(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Nothof_p_ERA>0.05);
Nothof_corr2_ERA = Nothof_corr_ERA;

for i = 1:length(col)
    Nothof_corr2_ERA(row(i),col(i)) = NaN;
end

% Halocarpus (New Zealand)

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(NHaloca(378:397,2),squeeze(uwind_annual(j,i,1:20)));
     Haloca_corr_ERA(i,j)=corr(1,2); 
     Haloca_p_ERA(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Haloca_p_ERA>0.05);
Haloca_corr2_ERA = Haloca_corr_ERA;

for i = 1:length(col)
    Haloca_corr2_ERA(row(i),col(i)) = NaN;
end

% Lake Tay (Western Australia)

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(NTay(378:397,2),squeeze(uwind_annual(j,i,1:20)));
     Tay_corr_ERA(i,j)=corr(1,2); 
     Tay_p_ERA(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Tay_p_ERA>0.05);
Tay_corr2_ERA = Tay_corr_ERA;

for i = 1:length(col)
    Tay_corr2_ERA(row(i),col(i)) = NaN;
end

% Dust (West Antarctica)

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(NKoff(378:397,2),squeeze(uwind_annual(j,i,1:20)));
     Koff_corr_ERA(i,j)=corr(1,2); 
     Koff_p_ERA(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Koff_p_ERA>0.05);
Koff_corr2_ERA = Koff_corr_ERA;

for i = 1:length(col)
    Koff_corr2_ERA(row(i),col(i)) = NaN;
end

% figure(25)
% contourf(longitude,latitude,Koff_corr_ERA,80,'linestyle','none')
% figure(26)
% contourf(longitude,latitude,Koff_p_ERA,80,'linestyle','none')
% figure(27)
% contourf(longitude,latitude,Koff_corr2_ERA,80,'linestyle','none')

nccreate('EOF_ERA_corr.nc','latitude','Dimensions',...
    {'latitude',121},'Format','classic')
nccreate('EOF_ERA_corr.nc','longitude','Dimensions',...
    {'longitude',480},'Format','classic')
ncwrite('EOF_ERA_corr.nc','latitude',flipud(latitude(1:121,1)))
ncwrite('EOF_ERA_corr.nc','longitude',flipud(longitude(1:480,1)))

ncwriteatt('EOF_ERA_corr.nc','latitude','axis','Y')
ncwriteatt('EOF_ERA_corr.nc','latitude','standard_name','latitude')
ncwriteatt('EOF_ERA_corr.nc','latitude','units','degrees_north')
ncwriteatt('EOF_ERA_corr.nc','longitude','axis','X')
ncwriteatt('EOF_ERA_corr.nc','longitude','standard_name','longitude')
ncwriteatt('EOF_ERA_corr.nc','longitude','units','degrees_east')

nccreate('EOF_ERA_corr.nc','EOF1_corr','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
nccreate('EOF_ERA_corr.nc','EOF2_corr','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')

nccreate('EOF_ERA_corr.nc','Arau_corr','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
nccreate('EOF_ERA_corr.nc','Aust_corr','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
nccreate('EOF_ERA_corr.nc','Halo_corr','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
nccreate('EOF_ERA_corr.nc','Noth_corr','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
nccreate('EOF_ERA_corr.nc','Koff_corr','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
nccreate('EOF_ERA_corr.nc','Tay_corr','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')

 nccreate('EOF_ERA_corr.nc','Koff_corr_1','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
 nccreate('EOF_ERA_corr.nc','Koff_p','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
ncwrite('EOF_ERA_corr.nc','Koff_corr_1',Koff_corr_ERA(1:121,1:480))
ncwrite('EOF_ERA_corr.nc','Koff_p',Koff_p_ERA(1:121,1:480))

ncwrite('EOF_ERA_corr.nc','EOF1_corr',EOF1_corr2_ERA(1:121,1:480))
ncwrite('EOF_ERA_corr.nc','EOF2_corr',EOF2_corr2_ERA(1:121,1:480))

ncwrite('EOF_ERA_corr.nc','Arau_corr',Arauca_corr2_ERA(1:121,1:480))
ncwrite('EOF_ERA_corr.nc','Aust_corr',Austro_corr2_ERA(1:121,1:480))
ncwrite('EOF_ERA_corr.nc','Halo_corr',Haloca_corr2_ERA(1:121,1:480))
ncwrite('EOF_ERA_corr.nc','Noth_corr',Nothof_corr2_ERA(1:121,1:480))
ncwrite('EOF_ERA_corr.nc','Koff_corr',Koff_corr2_ERA(1:121,1:480))
ncwrite('EOF_ERA_corr.nc','Tay_corr',Tay_corr2_ERA(1:121,1:480))

ncwriteatt('EOF_ERA_corr.nc','EOF1_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','EOF1_corr','long_name','Correlation of EOF1 with ERA-Interim at 95% significance level, 1979-1998')
ncwriteatt('EOF_ERA_corr.nc','EOF2_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','EOF2_corr','long_name','Correlation of EOF2 with ERA-Interim at 95% significance level, 1979-1998')

ncwriteatt('EOF_ERA_corr.nc','Arau_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','Arau_corr','long_name','Correlation of Araucaria with ERA at 95% significance level, 1979-1998')
ncwriteatt('EOF_ERA_corr.nc','Aust_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','Aust_corr','long_name','Correlation of Austrocedrus with ERA at 95% significance level, 1979-1998')
ncwriteatt('EOF_ERA_corr.nc','Halo_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','Halo_corr','long_name','Correlation of Halocarpus with ERA at 95% significance level, 1979-1998')
ncwriteatt('EOF_ERA_corr.nc','Noth_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','Noth_corr','long_name','Correlation of Nothofagus with ERA at 95% significance level, 1979-1998')
ncwriteatt('EOF_ERA_corr.nc','Koff_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','Koff_corr','long_name','Correlation of Coarse particle dust % with ERA at 95% significance level, 1979-1998')
ncwriteatt('EOF_ERA_corr.nc','Tay_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','Tay_corr','long_name','Correlation of Callitris with ERA at 95% significance level, 1979-1998')

%% ERA40: Coarse and fine res.
% Coarse (2.5x2.5) resoltion
clear
ncid = netcdf.open('ERA40_uwind_850hPa_coarse.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

load('EOF1_2.mat')

%Compute annual means

X = reshape(u,144,37,12,(492/12));
uwind_annual = squeeze(nanmean(X,3));

for i=1:length(latitude),
  for j=1:length(longitude),
%      [corr,pval]=corrcoef(EOF1(357:397,2),squeeze(uwind_annual(j,i,1:41)));
%      EOF1_corr_ERA40c(i,j)=corr(1,2); 
%      EOF1_p_ERA40c(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
     [corr,pval]=corrcoef(EOF1(378:397,2),squeeze(uwind_annual(j,i,22:41)));
     EOF1_corr_ERA40_79(i,j)=corr(1,2); 
     EOF1_p_ERA40_79(i,j)=pval(1,2);
  end
end

% [row,col] = find(EOF1_p_ERA40c>0.05);
% EOF1_corr2_ERA40c = EOF1_corr_ERA40c;
% 
% for i = 1:length(col)
%     EOF1_corr2_ERA40c(row(i),col(i)) = NaN;
% end
% 
[row,col] = find(EOF1_p_ERA40_79>0.05);
EOF1_corr2_ERA40_79 = EOF1_corr_ERA40_79;

for i = 1:length(col)
    EOF1_corr2_ERA40_79(row(i),col(i)) = NaN;
end

% figure(19)
% contourf(longitude,latitude,EOF1_corr_ERA40c,80,'linestyle','none')
% figure(20)
% contourf(longitude,latitude,EOF1_p_ERA40c,80,'linestyle','none')
% figure(21)
% contourf(longitude,latitude,EOF1_corr2_ERA40c,80,'linestyle','none')

for i=1:length(latitude),
  for j=1:length(longitude),
%      [corr,pval]=corrcoef(EOF2(357:397,2),squeeze(uwind_annual(j,i,1:41)));
%      EOF2_corr_ERA40c(i,j)=corr(1,2); 
%      EOF2_p_ERA40c(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
     [corr,pval]=corrcoef(EOF2(378:397,2),squeeze(uwind_annual(j,i,22:41)));
     EOF2_corr_ERA40_79(i,j)=corr(1,2); 
     EOF2_p_ERA40_79(i,j)=pval(1,2);
  end
end

% [row,col] = find(EOF2_p_ERA40c>0.05);
% EOF2_corr2_ERA40c = EOF2_corr_ERA40c;
% 
% for i = 1:length(col)
%     EOF2_corr2_ERA40c(row(i),col(i)) = NaN;
% end

[row,col] = find(EOF2_p_ERA40_79>0.05);
EOF2_corr2_ERA40_79 = EOF2_corr_ERA40_79;

for i = 1:length(col)
    EOF2_corr2_ERA40_79(row(i),col(i)) = NaN;
end

% figure(22)
% contourf(longitude,latitude,EOF2_corr_ERA40c,80,'linestyle','none')
% figure(23)
% contourf(longitude,latitude,EOF2_p_ERA40c,80,'linestyle','none')
% figure(24)
% contourf(longitude,latitude,EOF2_corr2_ERA40c,80,'linestyle','none')


% figure(25)
% contourf(longitude,latitude,EOF3_corr_ERA40c,80,'linestyle','none')
% figure(26)
% contourf(longitude,latitude,EOF3_p_ERA40c,80,'linestyle','none')
% figure(27)
% contourf(longitude,latitude,EOF3_corr2_ERAc,80,'linestyle','none')

nccreate('EOF_ERA40c_corr.nc','latitude','Dimensions',...
    {'latitude',37},'Format','classic')
nccreate('EOF_ERA40c_corr.nc','longitude','Dimensions',...
    {'longitude',144},'Format','classic')
ncwrite('EOF_ERA40c_corr.nc','latitude',flipud(latitude(1:37,1)))
ncwrite('EOF_ERA40c_corr.nc','longitude',flipud(longitude(1:144,1)))

ncwriteatt('EOF_ERA40c_corr.nc','latitude','axis','Y')
ncwriteatt('EOF_ERA40c_corr.nc','latitude','standard_name','latitude')
ncwriteatt('EOF_ERA40c_corr.nc','latitude','units','degrees_north')
ncwriteatt('EOF_ERA40c_corr.nc','longitude','axis','X')
ncwriteatt('EOF_ERA40c_corr.nc','longitude','standard_name','longitude')
ncwriteatt('EOF_ERA40c_corr.nc','longitude','units','degrees_east')

% nccreate('EOF_ERA40c_corr.nc','EOF1_corr_ERA40c','Datatype','single','Dimensions',...
%     {'latitude',37,'longitude', 144},'Format','classic')
% nccreate('EOF_ERA40c_corr.nc','EOF2_corr_ERA40c','Datatype','single','Dimensions',...
%     {'latitude',37,'longitude', 144},'Format','classic')

nccreate('EOF_ERA40c_corr.nc','EOF1_corr_ERA40_79','Datatype','single','Dimensions',...
    {'latitude',37,'longitude', 144},'Format','classic')
nccreate('EOF_ERA40c_corr.nc','EOF2_corr_ERA40_79','Datatype','single','Dimensions',...
    {'latitude',37,'longitude', 144},'Format','classic')

% ncwrite('EOF_ERA40c_corr.nc','EOF1_corr_ERA40c',EOF1_corr2_ERA40c(1:37,1:144))
% ncwrite('EOF_ERA40c_corr.nc','EOF2_corr_ERA40c',EOF2_corr2_ERA40c(1:37,1:144))

ncwrite('EOF_ERA40c_corr.nc','EOF1_corr_ERA40_79',EOF1_corr2_ERA40_79(1:37,1:144))
ncwrite('EOF_ERA40c_corr.nc','EOF2_corr_ERA40_79',EOF2_corr2_ERA40_79(1:37,1:144))

ncwriteatt('EOF_ERA40c_corr.nc','EOF1_corr_ERA40_79','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA40c_corr.nc','EOF1_corr_ERA40_79','long_name','Correlation of EOF1 with ERA40 at 95% significance level, 1958-1998')
ncwriteatt('EOF_ERA40c_corr.nc','EOF2_corr_ERA40_79','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA40c_corr.nc','EOF2_corr_ERA40_79','long_name','Correlation of EOF2 with ERA40 at 95% significance level, 1958-1998')

% ncwriteatt('EOF_ERA40c_corr.nc','EOF1_corr_ERA40c','standard_name','Correlation at 95% sig.')
% ncwriteatt('EOF_ERA40c_corr.nc','EOF1_corr_ERA40c','long_name','Correlation of EOF1 with ERA40 at 95% significance level, 1958-1998')
% ncwriteatt('EOF_ERA40c_corr.nc','EOF2_corr_ERA40c','standard_name','Correlation at 95% sig.')
% ncwriteatt('EOF_ERA40c_corr.nc','EOF2_corr_ERA40c','long_name','Correlation of EOF2 with ERA40 at 95% significance level, 1958-1998')

% Fine (0.75x0.75) resolution

clear
ncid = netcdf.open('ERA40_uwind_850hPa_fine.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

load('EOF1_2.mat')

%Compute annual means

X = reshape(u,480,121,12,(492/12));
uwind_annual = squeeze(nanmean(X,3));

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(EOF1(357:397,2),squeeze(uwind_annual(j,i,1:41)));
     EOF1_corr_ERA40f(i,j)=corr(1,2); 
     EOF1_p_ERA40f(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF1_p_ERA40f>0.05);
EOF1_corr2_ERA40f = EOF1_corr_ERA40f;

for i = 1:length(col)
    EOF1_corr2_ERA40f(row(i),col(i)) = NaN;
end

% figure(19)
% contourf(longitude,latitude,EOF1_corr_ERA40f,80,'linestyle','none')
% figure(20)
% contourf(longitude,latitude,EOF1_p_ERA40f,80,'linestyle','none')
% figure(21)
% contourf(longitude,latitude,EOF1_corr2_ERA40f,80,'linestyle','none')

for i=1:length(latitude),
  for j=1:length(longitude),
     [corr,pval]=corrcoef(EOF2(357:397,2),squeeze(uwind_annual(j,i,1:41)));
     EOF2_corr_ERA40f(i,j)=corr(1,2); 
     EOF2_p_ERA40f(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF2_p_ERA40f>0.05);
EOF2_corr2_ERA40f = EOF2_corr_ERA40f;

for i = 1:length(col)
    EOF2_corr2_ERA40f(row(i),col(i)) = NaN;
end

% figure(22)
% contourf(longitude,latitude,EOF2_corr_ERA40f,80,'linestyle','none')
% figure(23)
% contourf(longitude,latitude,EOF2_p_ERA40f,80,'linestyle','none')
% figure(24)
% contourf(longitude,latitude,EOF2_corr2_ERA40f,80,'linestyle','none')

% figure(25)
% contourf(longitude,latitude,EOF3_corr_ERA40f,80,'linestyle','none')
% figure(26)
% contourf(longitude,latitude,EOF3_p_ERA40f,80,'linestyle','none')
% figure(27)
% contourf(longitude,latitude,EOF3_corr2_ERAc,80,'linestyle','none')

nccreate('EOF_ERA_corr.nc','EOF1_corr_ERA40f','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')
nccreate('EOF_ERA_corr.nc','EOF2_corr_ERA40f','Datatype','single','Dimensions',...
    {'latitude',121,'longitude', 480},'Format','classic')

ncwrite('EOF_ERA_corr.nc','EOF1_corr_ERA40f',EOF1_corr2_ERA40f(1:121,1:480))
ncwrite('EOF_ERA_corr.nc','EOF2_corr_ERA40f',EOF2_corr2_ERA40f(1:121,1:480))

ncwriteatt('EOF_ERA_corr.nc','EOF1_corr_ERA40f','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','EOF1_corr_ERA40f','long_name','Correlation of EOF1 with ERA40 at 95% significance level, 1958-1998')
ncwriteatt('EOF_ERA_corr.nc','EOF2_corr_ERA40f','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_ERA_corr.nc','EOF2_corr_ERA40f','long_name','Correlation of EOF2 with ERA40 at 95% significance level, 1958-1998')

%% 20CR
clear
load('20CR.mat')
load('EOF1_2.mat')
load('All_proxies_commonT.mat')
uwind = single(ann20CR_1870_1998);

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(EOF1(269:397,2),squeeze(uwind(j,i,1:129)));
     EOF1_corr_20CR(i,j)=corr(1,2); 
     EOF1_p_20CR(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF1_p_20CR>0.05);
EOF1_corr2_20CR = EOF1_corr_20CR;

for i = 1:length(col)
    EOF1_corr2_20CR(row(i),col(i)) = NaN;
end

% figure(19)
% contourf(lon,lat,EOF1_corr_20CR,80,'linestyle','none')
% figure(20)
% contourf(lon,lat,EOF1_p_20CR,80,'linestyle','none')
% figure(21)
% contourf(lon,lat,EOF1_corr2_20CR,80,'linestyle','none')

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(EOF2(269:397,2),squeeze(uwind(j,i,1:129)));
     EOF2_corr_20CR(i,j)=corr(1,2); 
     EOF2_p_20CR(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF2_p_20CR>0.05);
EOF2_corr2_20CR = EOF2_corr_20CR;

for i = 1:length(col)
    EOF2_corr2_20CR(row(i),col(i)) = NaN;
end

% figure(22)
% contourf(lon,lat,EOF2_corr_20CR,80,'linestyle','none')
% figure(23)
% contourf(lon,lat,EOF2_p_20CR,80,'linestyle','none')
% figure(24)
% contourf(lon,lat,EOF2_corr2_20CR,80,'linestyle','none')

% figure(25)
% contourf(lon,lat,EOF3_corr_20CR,80,'linestyle','none')
% figure(26)
% contourf(lon,lat,EOF3_p_20CR,80,'linestyle','none')
% figure(27)
% contourf(lon,lat,EOF3_corr2_20CR,80,'linestyle','none')

% Correlations for individual series which show coherence in wavelets

% Araucaria (South America)

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(NArauca(269:397,2),squeeze(uwind(j,i,1:129)));
     Arauca_corr_20CR(i,j)=corr(1,2); 
     Arauca_p_20CR(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Arauca_p_20CR>0.05);
Arauca_corr2_20CR = Arauca_corr_20CR;

for i = 1:length(col)
    Arauca_corr2_20CR(row(i),col(i)) = NaN;
end

% Austrocedrus (South America)

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(NAustro(269:397,2),squeeze(uwind(j,i,1:129)));
     Austro_corr_20CR(i,j)=corr(1,2); 
     Austro_p_20CR(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Austro_p_20CR>0.05);
Austro_corr2_20CR = Austro_corr_20CR;

for i = 1:length(col)
    Austro_corr2_20CR(row(i),col(i)) = NaN;
end

% Nothofagus (South America)

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(NNothof(269:397,2),squeeze(uwind(j,i,1:129)));
     Nothof_corr_20CR(i,j)=corr(1,2); 
     Nothof_p_20CR(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Nothof_p_20CR>0.05);
Nothof_corr2_20CR = Nothof_corr_20CR;

for i = 1:length(col)
    Nothof_corr2_20CR(row(i),col(i)) = NaN;
end

% Halocarpus (New Zealand)

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(NHaloca(269:397,2),squeeze(uwind(j,i,1:129)));
     Haloca_corr_20CR(i,j)=corr(1,2); 
     Haloca_p_20CR(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Haloca_p_20CR>0.05);
Haloca_corr2_20CR = Haloca_corr_20CR;

for i = 1:length(col)
    Haloca_corr2_20CR(row(i),col(i)) = NaN;
end

% Lake Tay (Western Australia)

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(NTay(269:397,2),squeeze(uwind(j,i,1:129)));
     Tay_corr_20CR(i,j)=corr(1,2); 
     Tay_p_20CR(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Tay_p_20CR>0.05);
Tay_corr2_20CR = Tay_corr_20CR;

for i = 1:length(col)
    Tay_corr2_20CR(row(i),col(i)) = NaN;
end

% Dust (West Antarctica)

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(NKoff(269:397,2),squeeze(uwind(j,i,1:129)));
     Koff_corr_20CR(i,j)=corr(1,2); 
     Koff_p_20CR(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(Koff_p_20CR>0.05);
Koff_corr2_20CR = Koff_corr_20CR;

for i = 1:length(col)
    Koff_corr2_20CR(row(i),col(i)) = NaN;
end

nccreate('EOF_20CR_corr.nc','latitude','Dimensions',...
    {'latitude',46},'Format','classic')
nccreate('EOF_20CR_corr.nc','longitude','Dimensions',...
    {'longitude',180},'Format','classic')
ncwrite('EOF_20CR_corr.nc','latitude',flipud(lat(1:46,1)))
ncwrite('EOF_20CR_corr.nc','longitude',flipud(lon(1:180,1)))

ncwriteatt('EOF_20CR_corr.nc','latitude','axis','Y')
ncwriteatt('EOF_20CR_corr.nc','latitude','standard_name','latitude')
ncwriteatt('EOF_20CR_corr.nc','latitude','units','degrees_north')
ncwriteatt('EOF_20CR_corr.nc','longitude','axis','X')
ncwriteatt('EOF_20CR_corr.nc','longitude','standard_name','longitude')
ncwriteatt('EOF_20CR_corr.nc','longitude','units','degrees_east')

nccreate('EOF_20CR_corr.nc','EOF1_corr','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','EOF2_corr','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')

nccreate('EOF_20CR_corr.nc','Arau_corr','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','Aust_corr','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','Halo_corr','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','Noth_corr','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','Koff_corr','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','Tay_corr','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')

ncwrite('EOF_20CR_corr.nc','EOF1_corr',EOF1_corr2_20CR(1:46,1:180))
ncwrite('EOF_20CR_corr.nc','EOF2_corr',EOF2_corr2_20CR(1:46,1:180))

ncwrite('EOF_20CR_corr.nc','Arau_corr',Arauca_corr2_20CR(1:46,1:180))
ncwrite('EOF_20CR_corr.nc','Aust_corr',Austro_corr2_20CR(1:46,1:180))
ncwrite('EOF_20CR_corr.nc','Halo_corr',Haloca_corr2_20CR(1:46,1:180))
ncwrite('EOF_20CR_corr.nc','Noth_corr',Nothof_corr2_20CR(1:46,1:180))
ncwrite('EOF_20CR_corr.nc','Koff_corr',Koff_corr2_20CR(1:46,1:180))
ncwrite('EOF_20CR_corr.nc','Tay_corr',Tay_corr2_20CR(1:46,1:180))

ncwriteatt('EOF_20CR_corr.nc','EOF1_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','EOF1_corr','long_name','Correlation of EOF1 with 20CR at 95% significance level, 1870-1998')
ncwriteatt('EOF_20CR_corr.nc','EOF2_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','EOF2_corr','long_name','Correlation of EOF2 with 20CR at 95% significance level, 1870-1998')

ncwriteatt('EOF_20CR_corr.nc','Arau_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','Arau_corr','long_name','Correlation of Araucaria with 20CR at 95% significance level, 1870-1998')
ncwriteatt('EOF_20CR_corr.nc','Aust_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','Aust_corr','long_name','Correlation of Austrocedrus with 20CR at 95% significance level, 1870-1998')
ncwriteatt('EOF_20CR_corr.nc','Halo_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','Halo_corr','long_name','Correlation of Halocarpus with 20CR at 95% significance level, 1870-1998')
ncwriteatt('EOF_20CR_corr.nc','Noth_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','Noth_corr','long_name','Correlation of Nothofagus with 20CR at 95% significance level, 1870-1998')
ncwriteatt('EOF_20CR_corr.nc','Koff_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','Koff_corr','long_name','Correlation of Coarse particle dust % with 20CR at 95% significance level, 1870-1998')
ncwriteatt('EOF_20CR_corr.nc','Tay_corr','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','Tay_corr','long_name','Correlation of Callitris with 20CR at 95% significance level, 1870-1998')

%% 20CR - split: 1870-1970; 1971-1998
clear
load('20CR.mat')
load('EOF1_2.mat')
uwind = single(ann20CR_1870_1998);

% Early - pre 1970

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(EOF1(269:369,2),squeeze(uwind(j,i,1:101)));
     EOF1_corr_20CR_early(i,j)=corr(1,2); 
     EOF1_p_20CR_early(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF1_p_20CR_early>0.05);
EOF1_corr2_20CR_early = EOF1_corr_20CR_early;

for i = 1:length(col)
    EOF1_corr2_20CR_early(row(i),col(i)) = NaN;
end


for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(EOF2(269:369,2),squeeze(uwind(j,i,1:101)));
     EOF2_corr_20CR_early(i,j)=corr(1,2); 
     EOF2_p_20CR_early(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
  end
end

[row,col] = find(EOF2_p_20CR_early>0.05);
EOF2_corr2_20CR_early = EOF2_corr_20CR_early;

for i = 1:length(col)
    EOF2_corr2_20CR_early(row(i),col(i)) = NaN;
end

% Late - post 1970

for i=1:length(lat),
  for j=1:length(lon),
     %[corr,pval]=corrcoef(EOF1(370:397,2),squeeze(uwind(j,i,102:129)));
     [corr,pval]=corrcoef(EOF1(378:397,2),squeeze(uwind(j,i,110:129)));    % For 1979-1998
     %EOF1_corr_20CR_late(i,j)=corr(1,2); 
     %EOF1_p_20CR_late(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
     EOF1_corr_20CR_79(i,j)=corr(1,2);                                     % For 1979-1998
     EOF1_p_20CR_79(i,j)=pval(1,2);                                        % For 1979-1998
  end
end

% [row,col] = find(EOF1_p_20CR_late>0.05);
% EOF1_corr2_20CR_late = EOF1_corr_20CR_late;
% 
% for i = 1:length(col)
%     EOF1_corr2_20CR_late(row(i),col(i)) = NaN;
% end

[row,col] = find(EOF1_p_20CR_79>0.05);
EOF1_corr2_20CR_79 = EOF1_corr_20CR_79;

for i = 1:length(col)
    EOF1_corr2_20CR_79(row(i),col(i)) = NaN;
end


for i=1:length(lat),
  for j=1:length(lon),
     %[corr,pval]=corrcoef(EOF2(370:397,2),squeeze(uwind(j,i,102:129)));
     [corr,pval]=corrcoef(EOF2(378:397,2),squeeze(uwind(j,i,110:129)));    % For 1979-1998
     %EOF2_corr_20CR_late(i,j)=corr(1,2); 
     %EOF2_p_20CR_late(i,j)=pval(1,2);  % although as I said pval here can be evaluated better for most geophysical timeseries. 2/sqrt(N)
     EOF2_corr_20CR_79(i,j)=corr(1,2);                                     % For 1979-1998
     EOF2_p_20CR_79(i,j)=pval(1,2);                                        % For 1979-1998
  end
end

% [row,col] = find(EOF2_p_20CR_late>0.05);
% EOF2_corr2_20CR_late = EOF2_corr_20CR_late;
% 
% for i = 1:length(col)
%     EOF2_corr2_20CR_late(row(i),col(i)) = NaN;
% end

[row,col] = find(EOF2_p_20CR_79>0.05);
EOF2_corr2_20CR_79 = EOF2_corr_20CR_79;

for i = 1:length(col)
    EOF2_corr2_20CR_79(row(i),col(i)) = NaN;
end

% nccreate('EOF_20CR_corr.nc','EOF1_corr_early','Datatype','single','Dimensions',...
%     {'latitude',46,'longitude', 180},'Format','classic')
% nccreate('EOF_20CR_corr.nc','EOF2_corr_early','Datatype','single','Dimensions',...
%     {'latitude',46,'longitude', 180},'Format','classic')
% nccreate('EOF_20CR_corr.nc','EOF1_corr_late','Datatype','single','Dimensions',...
%     {'latitude',46,'longitude', 180},'Format','classic')
% nccreate('EOF_20CR_corr.nc','EOF2_corr_late','Datatype','single','Dimensions',...
%     {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','EOF1_corr_79','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','EOF2_corr_79','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')
nccreate('EOF_20CR_corr.nc','EOF3_corr_79','Datatype','single','Dimensions',...
    {'latitude',46,'longitude', 180},'Format','classic')

% ncwrite('EOF_20CR_corr.nc','EOF1_corr_early',EOF1_corr2_20CR_early(1:46,1:180))
% ncwrite('EOF_20CR_corr.nc','EOF2_corr_early',EOF2_corr2_20CR_early(1:46,1:180))
% ncwrite('EOF_20CR_corr.nc','EOF1_corr_late',EOF1_corr2_20CR_late(1:46,1:180))
% ncwrite('EOF_20CR_corr.nc','EOF2_corr_late',EOF2_corr2_20CR_late(1:46,1:180))
ncwrite('EOF_20CR_corr.nc','EOF1_corr_79',EOF1_corr2_20CR_79(1:46,1:180))
ncwrite('EOF_20CR_corr.nc','EOF2_corr_79',EOF2_corr2_20CR_79(1:46,1:180))

% ncwriteatt('EOF_20CR_corr.nc','EOF1_corr_early','standard_name','Correlation at 95% sig.')
% ncwriteatt('EOF_20CR_corr.nc','EOF1_corr_early','long_name','Correlation of EOF1 with 20CR at 95% significance level, 1870-1970')
% ncwriteatt('EOF_20CR_corr.nc','EOF2_corr_early','standard_name','Correlation at 95% sig.')
% ncwriteatt('EOF_20CR_corr.nc','EOF2_corr_early','long_name','Correlation of EOF2 with 20CR at 95% significance level, 1870-1970')
% ncwriteatt('EOF_20CR_corr.nc','EOF1_corr_late','standard_name','Correlation at 95% sig.')
% ncwriteatt('EOF_20CR_corr.nc','EOF1_corr_late','long_name','Correlation of EOF1 with 20CR at 95% significance level, 1971-1998')
% ncwriteatt('EOF_20CR_corr.nc','EOF2_corr_late','standard_name','Correlation at 95% sig.')
% ncwriteatt('EOF_20CR_corr.nc','EOF2_corr_late','long_name','Correlation of EOF2 with 20CR at 95% significance level, 1971-1998')

ncwriteatt('EOF_20CR_corr.nc','EOF1_corr_79','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','EOF1_corr_79','long_name','Correlation of EOF1 with 20CR at 95% significance level, 1979-1998')
ncwriteatt('EOF_20CR_corr.nc','EOF2_corr_79','standard_name','Correlation at 95% sig.')
ncwriteatt('EOF_20CR_corr.nc','EOF2_corr_79','long_name','Correlation of EOF2 with 20CR at 95% significance level, 1979-1998')

%% 20CR timeseries of max zonal wind speed latitude

% 20CR correlation/ linear regression stuff

% clear
% load('20CR.mat')
% load('EOF_1_2_3.mat')
% uwind = single(ann20CR_1870_1998);
% uwind_mean = nanmean(uwind,1);
% 
% [max,lat_m] = max(uwind_mean,[],2);
% latm2 = squeeze(lat_m);
% 
% for i=1:length(latm2)
%    uwind_max(i)=lat(latm2(i));
% end
% 
% [corr,pval]=corrcoef(EOF3(269:397,2),uwind_max(1:129));
% %x = EOF3(269:397,2);         %rsq = 0.0167, pretty bad
% %x = NArauca(269:397,2);      %rsq = 0.0113, pretty bad
% %x = NAustro(269:397,2);      %rsq = 0.0085, pretty bad
% %x = NHaloca(269:397,2);      %rsq = 0.0169, pretty bad
% y = uwind_max(1:129)';
% p = polyfit(x,y,1);
% yfit = polyval(p,x);
% yresid = y - yfit;
% SSresid = sum(yresid.^2);    
% SStotal = (length(y)-1) * var(y);     
% rsq = 1 - SSresid/SStotal      
      
