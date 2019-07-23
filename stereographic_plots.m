% Polar Stereographic Plots 

%% Set up map
figure
axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -30])
framem
load coast
plotm(lat,long,'k','linewidth',2)
gridm %If you want gridlines
clear lat long

%% Load and plot data

ncid = netcdf.open('ERA-Interim/Uwind850_ERA_Int_79_98.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

lat = im2double(latitude);
lon = im2double(longitude);
lon(480,1) = 179.9; %Eliminates white line, as lon only goes from -180:179.25
u1 = double(u);


contourfm(lat,lon,squeeze(nanmean(u1,3)*10e-4)',30,'linestyle','none')
%contourfm(lat,lon,squeeze(nanmean(u1,3))',30,'linestyle','none')

ice = ([-79.6 45.7; -75 0.08; -75 -84; -81 -148; -77.8 -102.9; -78.1 -95.7; -77.01 -89.1;...
        -86.5 -108; -78 106; -64 -57; -90 0; -79 -112; -72.48 159.06; -66.77 112.81; -84 43; ...
        -70.86 11.54; -75 0; -79.38 -111.24; -77.68 -124; -64.2 -57.7; -75.92 -84.25; -79.46 -112.09]); %lat,lon of ice proxy sites
trees = ([-33 120; -38.08 -71.08; -40.1 -72.05; -43 -72.08; -39.06 -71.04; -42.08 -71.06; -39 -71.08; -37.14 -71;...
          -39.03 -71.05; -41.2 -71.9; -54.2 -68.7; -41.02 -71.13; -41.03 -72.04; -41.02 -71.13; -35.1 -71.01;...
          -42 147; -43 170; -42 148; -42 172; -39 177; -42 146; -39 175; -42 146; -39 175; -36 174; -39 174;...
          -47 168; -41 173; -36 148]);
lake = ([-33.14 -70.15]);

plotm(ice,'Marker','o','MarkerFaceColor','white','linestyle','none','color','k','markersize',7)
plotm(trees,'Marker','o','MarkerFaceColor','m','linestyle','none','color','k','markersize',7)
plotm(lake,'Marker','o','MarkerFaceColor','blue','linestyle','none','color','k','markersize',7)

plot2svg('prox_loc_revised.svg')

%% Draw a SH Map of just the westerlies

figure
axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 -0],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
framem
load coast
plotm(lat,long,'k','linewidth',2)
gridm %If you want gridlines
clear lat long

ncid = netcdf.open('ERA-Interim/Uwind850_ERA_Int_79_98.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

lat = im2double(latitude);
lon = im2double(longitude);
lon(480,1) = 180; %Eliminates white line, as lon only goes from -180:179.25
u1 = double(u);

contourfm(lat,lon,squeeze(nanmean(u1,3)*10e-4)','linestyle','none')
colormap(b2r(-20,20))

plot2svg('SHW.svg')


