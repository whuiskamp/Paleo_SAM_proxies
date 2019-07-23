% Steig data
%% Spatial Correlations
clear
load Steig.mat NSteig_common
load marshall_SAM.mat
NSteig = NSteig_common;

Marshall_SAM = flipud(Marshall_SAM);
correlations = nan(2,16);

for i = 2:17
    [r, p] = corrcoef(Marshall_SAM(33:59,2),NSteig(1:27,i));
    correlations(1,i-1) = r(2);
    correlations(2,i-1) = p(2);
end

load ('SAM_20CR.mat','uwind_SF','uwind_MA','time');
steig_20CR_SF = nan(16,91,180);
steig_20CR_MA = nan(16,91,180);
NSteig = flipud(NSteig);

for i = 2:17
    for j = 1:91
        for k = 1:180
            steig_20CR_SF(i-1,j,k) = corr(NSteig(1:91,i),uwind_SF(23:113,j,k));
            steig_20CR_MA(i-1,j,k) = corr(NSteig(1:91,i),uwind_MA(23:113,j,k));
        end
    end
    i
end

save('steig_corrs.mat','steig_20CR_MA','steig_20CR_SF','lat_20CR','lon_20CR','time_20CR');

% JRA55

clear uwind_MA uwind_SF time_20CR steig_20CR_MA steig_20CR_SF lat_20CR lon_20CR
load ('JRA55_season.mat','uwind_SF','uwind_MA','time','lat','lon');
load ('Steig.mat', 'NSteig_common')
NSteig = flipud(NSteig_common);

steig_JRA55_SF = nan(16,73,288);
steig_JRA55_MA = nan(16,73,288);
% NSteig = flipud(NSteig);

for i = 2:17
    for j = 1:73
        for k = 1:288
            steig_JRA55_SF(i-1,j,k) = corr(NSteig(66:91,i),uwind_SF(1:26,j,k));
            steig_JRA55_MA(i-1,j,k) = corr(NSteig(66:91,i),uwind_MA(1:26,j,k));
        end
    end
    i
end

save('steig_corrs.mat','steig_JRA55_MA','steig_JRA55_SF','lat_JRA55','lon_JRA55','time_JRA55','-append');

% ERA-Interim

clear
load ERA_INT_seas.mat
load Steig.mat
NSteig = NSteig_common;
NSteig = flipud(NSteig);

steig_JRA55_SF = nan(16,721,1440);
steig_JRA55_MA = nan(16,721,1440);

for i = 2:17
    for j = 1:721
        for k = 1:1440
            steig_JRA55_SF(i-1,j,k) = corr(NSteig(87:91,i),uwind_SF(1:5,j,k));
            steig_JRA55_MA(i-1,j,k) = corr(NSteig(87:91,i),uwind_MA(1:5,j,k));
        end
    end
    i
end

save('steig_corrs.mat','steig_ERA_MA','steig_ERA_SF','lat_ERA','lon_ERA','time_ERA','-append');

%% Plotting
% 20CR

figure
fig_20CR_SF = tight_subplot(4,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);

for i = 1:16
    axes(fig_20CR_SF(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

lat = double(lat)';
lon = double(lon)';

good_SF = steig_20CR_SF; 
good_SF(abs(steig_20CR_SF) < 0.3) = nan;

for a = 1:16
    axes(fig_20CR_SF(a));
    m_pcolor(lon,lat,squeeze(steig_20CR_SF(a,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    colormap(b2r(-1,1));
    caxis manual
    caxis([-1 1]);
    m_contour(lon,lat,squeeze(good_SF(a,:,:)),'linewidth',1.5,'color','k');
end

figure(2)
fig_20CR_MA = tight_subplot(4,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);

for i = 1:16
    axes(fig_20CR_MA(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

lat = double(lat_20CR)';
lon = double(lon_20CR)';

good_MA = steig_20CR_MA; 
good_MA(abs(steig_20CR_MA) < 0.3) = nan;

for a = 1:16
    axes(fig_20CR_MA(a));
    m_pcolor(lon_20CR,lat_20CR,squeeze(steig_20CR_MA(a,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    colormap(b2r(-1,1));
    caxis manual
    caxis([-1 1]);
    m_contour(lon_20CR,lat_20CR,squeeze(good_MA(a,:,:)),'linewidth',1.5,'color','k');
end

% JRA55

figure(3)
fig_JRA55_SF = tight_subplot(4,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);

for i = 1:16
    axes(fig_JRA55_SF(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

lat = double(lat_JRA55)';
lon = double(lon_JRA55)';

good_SF = steig_JRA55_SF; 
good_SF(abs(steig_JRA55_SF) < 0.3) = nan;

for a = 1:16
    axes(fig_JRA55_SF(a));
    m_pcolor(lon,lat,squeeze(steig_JRA55_SF(a,:,:))); shading flat;
    hold
    m_coast('color','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    colormap(b2r(-1,1));
    caxis manual
    caxis([-1 1]);
    m_contour(lon,lat,squeeze(good_SF(a,:,:)),'linewidth',1.5,'color','k');
end

print -dpdf -painters steig_JRA55_SF.pdf

figure(4)
fig_JRA55_MA = tight_subplot(4,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);

for i = 1:16
    axes(fig_JRA55_MA(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

lat = double(lat_JRA55)';
lon = double(lon_JRA55)';

good_MA = steig_JRA55_MA; 
good_MA(abs(steig_JRA55_MA) < 0.3) = nan;

for a = 1:16
    axes(fig_JRA55_MA(a));
    m_pcolor(lon_JRA55,lat_JRA55,squeeze(steig_JRA55_MA(a,:,:))); shading flat;
    hold
    m_coast('color','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    colormap(b2r(-1,1));
    caxis manual
    caxis([-1 1]);
    m_contour(lon_JRA55,lat_JRA55,squeeze(good_MA(a,:,:)),'linewidth',1.5,'color','k');
end

print -dpdf -painters steig_JRA55_MA.pdf

% ERA-Int
clear
load('steig_corrs.mat','steig_ERA_MA','steig_ERA_SF','lat_ERA','lon_ERA');
figure
fig_ERA_SF = tight_subplot(4,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);

for i = 1:16
    axes(fig_ERA_SF(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

lat = double(lat_ERA)';
lon = double(lon_ERA)';

good_SF = steig_ERA_SF; 
good_SF(abs(steig_ERA_SF) < 0.3) = nan;

for a = 1:16
    axes(fig_ERA_SF(a));
    m_pcolor(lon,lat,squeeze(steig_ERA_SF(a,:,:))); shading flat;
    hold
    m_contour(lon,lat,squeeze(good_SF(a,:,:)),'linewidth',1.5,'color','k');
    m_coast('color','m','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    colormap(b2r(-1,1));
    caxis manual
    caxis([-1 1]);
end

%Do NOT save this thing as an SVG, its over 1/2 GB (I never did find out
%how big it would get...)

figure(2)
fig_ERA_MA = tight_subplot(4,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);

for i = 1:16
    axes(fig_ERA_MA(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

lat = double(lat_ERA)';
lon = double(lon_ERA)';

good_MA = steig_ERA_MA; 
good_MA(abs(steig_ERA_MA) < 0.3) = nan;

for a = 1:16
    axes(fig_ERA_MA(a));
    m_pcolor(lon_ERA,lat_ERA,squeeze(steig_ERA_MA(a,:,:))); shading flat;
    hold
    m_coast('patch',[0.7 0.7 0.7],'facealpha',0,'edgecolor','m','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    colormap(b2r(-1,1));
    caxis manual
    caxis([-1 1]);
    m_contour(lon_ERA,lat_ERA,squeeze(good_MA(a,:,:)),'linewidth',1.5,'color','k');
end

%% Wavelets
clear
load Steig.mat

for i = 2:17
    test = NSteig_common(:,[1 i]);
    figure
    clf
    wt(test)
    plot2svg(['Steig_common',num2str(i),'.svg'])
end

%% SAM index crosswavelets

clear
load Steig.mat
load SAM_seasonal.mat

fig_FJ_cwc_SF = tight_subplot(4,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);

for i = 2:17
    axes(fig_FJ_cwc_SF(i-1));
    test = NSteig_common(:,[1 i]);
    wtc(test,FogtJones_SF)
end
plot2svg(['Steig_FJ_cwc_SF.svg'])

fig_FJ_cwc_MA = tight_subplot(4,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);
for i = 2:17
    axes(fig_FJ_cwc_MA(i-1));
    test = NSteig_common(:,[1 i]);
    wtc(test,FogtJones_MA)
end
plot2svg(['Steig_FJ_cwc_MA.svg'])
    
%% All ice data EOF calculation

clear
load('All_Ice.mat','NAll_ice_common');

[COEFF,SCORE,latent,tsquare,var_explained] = EOF_calc(NAll_ice_common(:,2:20));

EOF1 = SCORE(:,2);

% correlate with JRA55

load ('JRA55_season.mat','uwind_SF','uwind_MA','time','lat','lon');

EOF1_JRA55_SF = nan(73,288);
EOF1_JRA55_MA = nan(73,288);
EOF1 = flipud(EOF1);

for j = 1:73
    for k = 1:288
        EOF1_JRA55_SF(j,k) = corr(EOF1(66:91),uwind_SF(1:26,j,k));
        EOF1_JRA55_MA(j,k) = corr(EOF1(66:91),uwind_MA(1:26,j,k));
    end
end

figure
m_proj('stereographic','lat',-90,'long',0,'radius',90);

lat = double(lat)';
lon = double(lon)';

good_SF = EOF1_JRA55_SF; 
good_SF(abs(EOF1_JRA55_SF) < 0.3) = nan;

m_pcolor(lon,lat,squeeze(EOF1_JRA55_SF(:,:))); shading flat;
hold
m_coast('color','k','linewidth',2);
m_grid('xticklabels',[],'yticklabels',[],'linest','--');
colormap(b2r(-1,1));
caxis manual
caxis([-1 1]);
m_contour(lon,lat,squeeze(good_SF(:,:)),'linewidth',1.5,'color','k');

print -dpdf -painters EOF1_JRA55_SF.pdf

figure(2)

m_proj('stereographic','lat',-90,'long',0,'radius',90);

good_MA = EOF1_JRA55_MA; 
good_MA(abs(EOF1_JRA55_MA) < 0.3) = nan;

m_pcolor(lon,lat,squeeze(EOF1_JRA55_MA(:,:))); shading flat;
hold
m_coast('color','k','linewidth',2);
m_grid('xticklabels',[],'yticklabels',[],'linest','--');
colormap(b2r(-1,1));
caxis manual
caxis([-1 1]);
m_contour(lon,lat,squeeze(good_MA(:,:)),'linewidth',1.5,'color','k');


print -dpdf -painters EOF1_JRA55_MA.pdf

% Correlate with 20CR

clear
load ('20CR');


EOF1_20CR_SF = nan(73,288);
EOF1_20CR_MA = nan(73,288);
EOF1_20CR_ANN = nan(180,46);
EOF1 = SCORE(:,4);
EOF1 = flipud(EOF1);

for j = 1:180
    for k = 1:46
        [EOF1_20CR_ANN(j,k) P_20CR_ANN(j,k)] = corr(EOF1,squeeze(ann20CR_1870_1998(j,k,24:114)));
    end
end

[row,col] = find(P_20CR_ANN>0.05);

for i = 1:length(col)
    EOF1_20CR_ANN(row(i),col(i)) = NaN;
end

load('SAM_20CR.mat')
for j = 1:91
    for k = 1:180
        [EOF1_20CR_SF(j,k) P_20CR_SF(j,k)] = corr(EOF1,uwind_SF(23:113,j,k));
        [EOF1_20CR_MA(j,k) P_20CR_MA(j,k)]= corr(EOF1,uwind_MA(23:113,j,k));
    end
end

[row,col] = find(P_20CR_SF>0.05);

for i = 1:length(col)
    EOF1_20CR_SF(row(i),col(i)) = NaN;
end

figure
m_proj('stereographic','lat',-90,'long',0,'radius',90);

lat = double(lat)';
lon = double(lon)';

good_ANN = EOF1_20CR_ANN; 
good_ANN(abs(EOF1_20CR_ANN) < 0.3) = nan;

m_pcolor(lon,lat,squeeze(EOF1_20CR_ANN(:,:)')); shading flat;
hold
m_coast('color','k','linewidth',2);
m_grid('xticklabels',[],'yticklabels',[],'linest','--');
colormap(b2r(-1,1));
caxis manual
caxis([-1 1]);
m_contour(lon,lat,squeeze(good_ANN(:,:)'),'linewidth',1.5,'color','k');

print -dpdf -painters EOF3_20CR_ANN_sig.pdf

figure(2)
m_proj('stereographic','lat',-90,'long',0,'radius',90);

lat = double(lat)';
lon = double(lon)';

good_SF = EOF1_20CR_SF; 
good_SF(abs(EOF1_20CR_SF) < 0.3) = nan;

m_pcolor(lon,lat,squeeze(EOF1_20CR_SF(:,:))); shading flat;
hold
m_coast('color','k','linewidth',2);
m_grid('xticklabels',[],'yticklabels',[],'linest','--');
colormap(b2r(-1,1));
caxis manual
caxis([-1 1]);
m_contour(lon,lat,squeeze(good_SF(:,:)),'linewidth',1.5,'color','k');

print -dpdf -painters EOF2_20CR_SF.pdf

figure(3)

[row,col] = find(P_20CR_MA>0.05);

for i = 1:length(col)
    EOF1_20CR_MA(row(i),col(i)) = NaN;
end

m_proj('stereographic','lat',-90,'long',0,'radius',90);

good_MA = EOF1_20CR_MA; 
good_MA(abs(EOF1_20CR_MA) < 0.3) = nan;

m_pcolor(lon,lat,squeeze(EOF1_20CR_MA(:,:))); shading flat;
hold
m_coast('color','k','linewidth',2);
m_grid('xticklabels',[],'yticklabels',[],'linest','--');
colormap(b2r(-1,1));
caxis manual
caxis([-1 1]);
m_contour(lon,lat,squeeze(good_MA(:,:)),'linewidth',1.5,'color','k');


print -dpdf -painters EOF2_20CR_MA.pdf
