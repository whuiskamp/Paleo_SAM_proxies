%% Scripts for ALL data over common time period (1893 - 1983)
% Run the EOF
clear
load NAll_data.mat
[COEFF,SCORE,latent,tsquare,var_explained] = EOF_calc(NAll_data_common(:,2:25));

save('NAll_data.mat','COEFF','SCORE','latent','tsquare','var_explained','-append');

[COEFF_S,SCORE_S,latent_S,tsquare_S,var_explained_S] = EOF_calc(NAll_data_shifted(:,2:25));

save('NAll_data.mat','COEFF_S','SCORE_S','latent_S','tsquare_S','var_explained_S','-append');

clear
load('NAll_data.mat','SCORE_S')
load ('SAM_20CR.mat','uwind_SF','uwind_MA','time');
load('20CR.mat')
SCORE_S = flipud(SCORE_S);

EOF1_20CR_SF = nan(91,180); P_EOF1_20CR_SF = nan(91,180);
EOF1_20CR_MA = nan(91,180); P_EOF1_20CR_MA = nan(91,180);
EOF1_20CR_AN = nan(46,180); P_EOF2_20CR_SF = nan(91,180);
EOF2_20CR_SF = nan(91,180); P_EOF2_20CR_MA = nan(91,180);
EOF2_20CR_MA = nan(91,180); P_EOF3_20CR_SF = nan(91,180);
EOF2_20CR_AN = nan(46,180); P_EOF3_20CR_MA = nan(91,180);
EOF3_20CR_SF = nan(91,180); P_EOF4_20CR_SF = nan(91,180);
EOF3_20CR_MA = nan(91,180); P_EOF4_20CR_MA = nan(91,180);
EOF3_20CR_AN = nan(46,180); P_EOF1_20CR_AN = nan(46,180);
EOF4_20CR_SF = nan(91,180); P_EOF2_20CR_AN = nan(46,180);
EOF4_20CR_MA = nan(91,180); P_EOF3_20CR_AN = nan(46,180);
EOF4_20CR_AN = nan(46,180); P_EOF4_20CR_AN = nan(46,180);

for j = 1:91
    for k = 1:180
        [EOF1_20CR_SF(j,k) P_EOF1_20CR_SF(j,k) ]=  corr(SCORE_S(:,1),uwind_SF(23:113,j,k));
        [EOF1_20CR_MA(j,k) P_EOF1_20CR_MA(j,k) ] =  corr(SCORE_S(:,1),uwind_MA(23:113,j,k));
        [EOF2_20CR_SF(j,k) P_EOF2_20CR_SF(j,k) ] =  corr(SCORE_S(:,2),uwind_SF(23:113,j,k));
        [EOF2_20CR_MA(j,k) P_EOF2_20CR_MA(j,k) ] =  corr(SCORE_S(:,2),uwind_MA(23:113,j,k));
        [EOF3_20CR_SF(j,k) P_EOF3_20CR_SF(j,k) ] =  corr(SCORE_S(:,3),uwind_SF(23:113,j,k));
        [EOF3_20CR_MA(j,k) P_EOF3_20CR_MA(j,k) ] =  corr(SCORE_S(:,3),uwind_MA(23:113,j,k));
        [EOF4_20CR_SF(j,k) P_EOF4_20CR_SF(j,k) ] =  corr(SCORE_S(:,4),uwind_SF(23:113,j,k));
        [EOF4_20CR_MA(j,k) P_EOF4_20CR_MA(j,k) ] =  corr(SCORE_S(:,4),uwind_MA(23:113,j,k));
    end
end
for i = 1:46
    for j = 1:180
        [EOF1_20CR_AN(i,j) P_EOF1_20CR_AN(i,j)] =  corr(SCORE_S(:,1),squeeze(ann20CR_1870_1998(j,i,24:114)));
        [EOF2_20CR_AN(i,j) P_EOF2_20CR_AN(i,j)] =  corr(SCORE_S(:,2),squeeze(ann20CR_1870_1998(j,i,24:114)));
        [EOF3_20CR_AN(i,j) P_EOF3_20CR_AN(i,j)] =  corr(SCORE_S(:,3),squeeze(ann20CR_1870_1998(j,i,24:114)));
        [EOF4_20CR_AN(i,j) P_EOF4_20CR_AN(i,j)] =  corr(SCORE_S(:,4),squeeze(ann20CR_1870_1998(j,i,24:114)));
    end
end

save('Alldata_corrs.mat','EOF1_20CR_SF','EOF1_20CR_MA','EOF1_20CR_AN','EOF2_20CR_SF','EOF2_20CR_MA','EOF2_20CR_AN',...
     'EOF3_20CR_SF','EOF3_20CR_MA','EOF3_20CR_AN','EOF4_20CR_SF','EOF4_20CR_MA','EOF4_20CR_AN',...
     'P_EOF1_20CR_SF','P_EOF1_20CR_MA','P_EOF1_20CR_AN','P_EOF2_20CR_SF','P_EOF2_20CR_MA','P_EOF2_20CR_AN',...
     'P_EOF3_20CR_SF','P_EOF3_20CR_MA','P_EOF3_20CR_AN','P_EOF4_20CR_SF','P_EOF4_20CR_MA','P_EOF4_20CR_AN');


%% Plotting - Timeseries plots

figure(1)
plot(var_explained_S,'+','linestyle','-') 
line([0 25],[50 50],'linestyle','--','color','black')
ylabel('Cumulative Variance Explained')
xlabel('No. of PCs')

% From this we want to retain the first 4 PCs. This is arbitrary, the
% threshold is that we can explain 50% of the variance.

rng = nan(91,2);
for i = 1:91
    rng(i,1) = min(NAll_data_shifted(i,2:25));
    rng(i,2) = max(NAll_data_shifted(i,2:25));
end

for i = 2:25
    run11yrmean(:,i) = tsmovavg(NAll_data_shifted(:,i),'s',11,1);
end

figure(2)
subplot(4,1,1)
jbfill(NAll_data_shifted(:,1)',rng(:,2)',rng(:,1)','b','k',[],0.3);
hold on
for i = 2:25
    plot(NAll_data_shifted(6:86,1),run11yrmean(11:91,i),...  % plot so that running mean is centered
        'color',rand(1,3),...
        'linewidth',2)
end

EOF1_11 = tsmovavg(SCORE_S(:,1),'s',11,1);
EOF2_11 = tsmovavg(SCORE_S(:,2),'s',11,1);
EOF3_11 = tsmovavg(SCORE_S(:,3),'s',11,1);

subplot(4,1,2)
hold on
plot(NAll_data_shifted(1:91,1),SCORE_S(:,1),'color','m','linewidth',2)
plot(NAll_data_shifted(6:86,1),EOF1_11(11:91),'color','k','linewidth',2)
subplot(4,1,3)
hold on
plot(NAll_data_shifted(1:91,1),SCORE_S(:,2),'color','r','linewidth',2)
plot(NAll_data_shifted(6:86,1),EOF2_11(11:91),'color','k','linewidth',2)
subplot(4,1,4)
hold on
plot(NAll_data_shifted(1:91,1),SCORE_S(:,3),'color','b','linewidth',2)
plot(NAll_data_shifted(6:86,1),EOF3_11(11:91),'color','k','linewidth',2)


%% Spatial correlations

% 20CR
clear
load Alldata_corrs.mat
load('20CR.mat','lon','lat')
lat_AN = lat;
lon_AN = lon;
load('SAM_20CR.mat','lon','lat')
lat = double(lat)';
lon = double(lon)';
lat_AN = double(lat_AN)';
lon_AN = double(lon_AN)';

EOF_seas = nan(8,91,180);
EOF_seas(1,:,:) = EOF1_20CR_SF;
EOF_seas(2,:,:) = EOF2_20CR_SF;
EOF_seas(3,:,:) = EOF3_20CR_SF;
EOF_seas(4,:,:) = EOF4_20CR_SF;
EOF_seas(5,:,:) = EOF1_20CR_MA;
EOF_seas(6,:,:) = EOF2_20CR_MA;
EOF_seas(7,:,:) = EOF3_20CR_MA;
EOF_seas(8,:,:) = EOF4_20CR_MA;
EOF_AN = nan(4,46,180);
EOF_AN(1,:,:) = EOF1_20CR_AN;
EOF_AN(2,:,:) = EOF2_20CR_AN;
EOF_AN(3,:,:) = EOF3_20CR_AN;
EOF_AN(4,:,:) = EOF4_20CR_AN;
P_EOF_seas = nan(8,91,180);
P_EOF_seas(1,:,:) = P_EOF1_20CR_SF;
P_EOF_seas(2,:,:) = P_EOF2_20CR_SF;
P_EOF_seas(3,:,:) = P_EOF3_20CR_SF;
P_EOF_seas(4,:,:) = P_EOF4_20CR_SF;
P_EOF_seas(5,:,:) = P_EOF1_20CR_MA;
P_EOF_seas(6,:,:) = P_EOF2_20CR_MA;
P_EOF_seas(7,:,:) = P_EOF3_20CR_MA;
P_EOF_seas(8,:,:) = P_EOF4_20CR_MA;
P_EOF_AN = nan(4,46,180);
P_EOF_AN(1,:,:) = P_EOF1_20CR_AN;
P_EOF_AN(2,:,:) = P_EOF2_20CR_AN;
P_EOF_AN(3,:,:) = P_EOF3_20CR_AN;
P_EOF_AN(4,:,:) = P_EOF4_20CR_AN;

figure
fig_20CR_SF = tight_subplot(3,4,[0.01 0.01],[0.10 0.1],[0.1 0.01]);
for i = 1:12
    axes(fig_20CR_SF(i));
    m_proj('stereographic','lat',-90,'long',0,'radius',90);
end

for a = 1:8
    clear row col
    [row,col] = find(squeeze(P_EOF_seas(a,:,:)) > 0.05);
    for i = 1:length(col)
        EOF_seas(a,row(i),col(i)) = NaN;
    end
    good_EOF = squeeze(EOF_seas(a,:,:)); 
    good_EOF(abs(good_EOF) < 0.3) = nan;
    axes(fig_20CR_SF(a));
    m_pcolor(lon,lat,squeeze(EOF_seas(a,:,:))); shading flat;
    hold
    m_coast('color','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    colormap(b2r(-1,1));
    caxis manual
    caxis([-1 1]);
    m_contour(lon,lat,squeeze(good_EOF(:,:)),'linewidth',1.5,'color','k');
end

for a = 1:4
    clear row col
    [row,col] = find(squeeze(P_EOF_AN(a,:,:)) > 0.05);
    for i = 1:length(col)
        EOF_AN(a,row(i),col(i)) = NaN;
    end
    good_EOF = squeeze(EOF_AN(a,:,:)); 
    good_EOF(abs(good_EOF) < 0.3) = nan;
    axes(fig_20CR_SF(a+8));
    m_pcolor(lon_AN,lat_AN,squeeze(EOF_AN(a,:,:))); shading flat;
    hold
    m_coast('color','k','linewidth',2);
    m_grid('xticklabels',[],'yticklabels',[],'linest','--');
    colormap(b2r(-1,1));
    caxis manual
    caxis([-1 1]);
    m_contour(lon_AN,lat_AN,squeeze(good_EOF(:,:)),'linewidth',1.5,'color','k');
end

print -dpdf -painters 'Alldata_20CR_EOFs.pdf'








