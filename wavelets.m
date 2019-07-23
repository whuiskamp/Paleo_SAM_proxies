%% Wavelet plots

load All_proxies_commonT
load marshall_SAM.mat
load Fogt_Jones.mat
load('SAM_seasonal.mat','Visbeck_Ann')
load('Villalba2012.mat')

Araucaria = flipud(Araucaria);
Austrocedrus = flipud(Austrocedrus);
Halocarpus = flipud(Halocarpus);
FJ_ann = flipud(FJ_ann); Visbeck_Ann = flipud(Visbeck_Ann);

% Check that our null hypothesis is correct (power decay should match red
% noise spectrum - broadly)

X=NNothof(:,2);
[P,freq]=pburg(zscore(X),7,[],1);
aa=ar1(X);
Ptheoretical=(1-aa.^2)./(abs(1-aa.*exp(-2*pi*i*freq))).^2;
semilogy(freq,P/sum(P),freq,Ptheoretical/sum(Ptheoretical),'r');
legend('observed',sprintf('Theoretical AR1=%.2f',aa),'location','best')

%% Create normalised data for each pair of time series
% The overlap function outputs zscores of the two time-series over a common
% interval
[zArau_Aust zAustro_Arau] = overlap(Araucaria,Austrocedrus);
[zArau_Halo zHalo_Arau] = overlap(Araucaria,Halocarpus);
[zAustro_Halo zHalo_Aust] = overlap(Austrocedrus, Halocarpus);

[zAustro_FJ zFJ_Austro] = overlap(Austrocedrus,FJ_ann);
[zAustro_Visbeck zVisbeck_Aust] = overlap(Austrocedrus,Visbeck_Ann);
[zArau_FJ zFJ_Arau] = overlap(Araucaria,FJ_ann);
[zArau_Visbeck zVisbeck_Arau] = overlap(Araucaria,Visbeck_Ann);
[zHalo_FJ zFJ_Halo] = overlap(Halocarpus,FJ_ann);
[zHalo_Visbeck zVisbeck_Halo] = overlap(Halocarpus,Visbeck_Ann);


%% Plot the continuous wavelet transforms
time = squeeze(NTay(:,1));
marshallSAM(:,:) = Marshall_SAM(:,1:2); % can now make wavelets with Marshall SAM

Fig1 = tight_subplot(3,3,[0.05 0.03],[0.10 0.01],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
tlim=[min(time(:)) max(time(:))];
axes(Fig1(1));
wt(marshallSAM)
axes(Fig1(2))
wt(FogtJones_DJF)
axes(Fig1(3))
wt(FogtJones_MAM)
axes(Fig1(4));
wt(NAustro);
set(gca,'xlim',tlim);
% caxis manual
% caxis([1/16 4]);
axes(Fig1(5));
wt(NArauca);
set(gca,'xlim',tlim)
axes(Fig1(6));
wt(NNothof);
set(gca,'xlim',tlim)
axes(Fig1(7));
wt(NHaloca);
set(gca,'xlim',tlim)
axes(Fig1(8));
wt(NTay);
set(gca,'xlim',tlim)
axes(Fig1(9));
wt(NKoff);
set(gca,'xlim',tlim)
% Be warned, using plot2svg with these figures outputs the wrong values on
% colourbars
plot2svg('wavelet_paleo.svg')

%% Cross wavelet transform or coherence - just find and replace wtc with xwt

Fig2 = tight_subplot(3,3,[0.05 0.03],[0.10 0.01],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig2(1));
wtc(zAustro_Arau,zArau_Aust)
colorbar('off')
axes(Fig2(2));
wtc(zAustro_Halo,zHalo_Aust)
colorbar('off')
axes(Fig2(3));
wtc(zArau_Halo,zHalo_Arau)
colorbar('off')
axes(Fig2(4));
wtc(zFJ_Arau,zArau_FJ)
axes(Fig2(5));
wtc(zFJ_Austro,zAustro_FJ)
axes(Fig2(6));
wtc(zFJ_Halo,zHalo_FJ)
axes(Fig2(7));
wtc(zVisbeck_Arau,zArau_Visbeck)
axes(Fig2(8));
wtc(zVisbeck_Aust,zAustro_Visbeck)
axes(Fig2(9));
wtc(zVisbeck_Halo,zHalo_Visbeck)

plot2svg('wave_coher_multi_2.svg')

%print(gcf,'wave_coher_multi','-depsc2','-painters') % using painters makes a proper vector graphic

Fig3 = tight_subplot(2,3,[0.05 0.03],[0.10 0.01],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig3(1));
wtc(NAustro,FogtJones_DJF)
colorbar('off')
axes(Fig3(2));
wtc(FogtJones_DJF,NNothof)
colorbar('off')
axes(Fig3(3));
wtc(FogtJones_DJF,NHaloca)
colorbar('off')
axes(Fig3(4));
wtc(FogtJones_DJF,NTay)
colorbar('off')
axes(Fig3(5));
wtc(FogtJones_DJF,NKoff)
colorbar('off')
axes(Fig3(6));
wtc(FogtJones_DJF,NNothof)
colorbar('off')
plot2svg('wave_coher_Fogt_DJF.svg')

Fig4 = tight_subplot(2,3,[0.05 0.03],[0.10 0.01],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig4(1));
wtc(NAustro,FogtJones_MAM)
colorbar('off')
axes(Fig4(2));
wtc(FogtJones_MAM,NNothof)
colorbar('off')
axes(Fig4(3));
wtc(FogtJones_MAM,NHaloca)
colorbar('off')
axes(Fig4(4));
wtc(FogtJones_MAM,NTay)
colorbar('off')
axes(Fig4(5));
wtc(FogtJones_MAM,NKoff)
colorbar('off')
axes(Fig4(6));
wtc(FogtJones_MAM,NNothof)
colorbar('off')
plot2svg('wave_coher_Fogt_MAM.svg')

% All the SAM index wavelets
load SAM_seasonal.mat
Fig5 = tight_subplot(2,3,[0.05 0.03],[0.10 0.01],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig5(1));
wt(marshall_MA)
axes(Fig5(2));
wt(FogtJones_MA)
axes(Fig5(3));
wt(visbeck_MA)
axes(Fig5(4));
wt(marshall_SF)
axes(Fig5(5));
wt(visbeck_SF)
axes(Fig5(6));
wt(FogtJones_SF)
plot2svg('wt_SAMs_SF.svg')

%% CWC for the SAM Sept-Feb records and Paleo

load All_proxies_commonT
load SAM_seasonal.mat

Fig6 = tight_subplot(2,3,[0.05 0.03],[0.10 0.1],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig6(1));
wtc(marshall_SF,NArauca)
title('Araucaria')
xlabel([])
axes(Fig6(2));
wtc(marshall_SF,NAustro)
title('Austrocedrus')
xlabel([])
axes(Fig6(3));
wtc(marshall_SF,NHaloca)
title('Halocarpus')
xlabel([])
axes(Fig6(4));
wtc(marshall_SF,NKoff)
title('Ice dust')
axes(Fig6(5));
wtc(marshall_SF,NNothof)
title('Nothofagus')
axes(Fig6(6));
wtc(marshall_MA,NTay)
title('Callitris')

plot2svg('cwc_marshall_seas.svg')

Fig7 = tight_subplot(3,3,[0.05 0.03],[0.10 0.1],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig7(1));
wtc(FogtJones_SF,NArauca)
title('Araucaria')
xlabel([])
axes(Fig7(2));
wtc(FogtJones_SF,NAustro)
title('Austrocedrus')
xlabel([])
axes(Fig7(3));
wtc(FogtJones_SF,NHaloca)
title('Halocarpus')
xlabel([])
axes(Fig7(4));
wtc(FogtJones_SF,NKoff)
title('Ice dust')
xlabel([])
axes(Fig7(5));
wtc(FogtJones_SF,NNothof)
title('Nothofagus')
axes(Fig7(6));
wtc(FogtJones_MA,NTay)
title('Callitris')
axes(Fig7(7));
wtc(FogtJones_MA,NMulva)
title('\delta D')


plot2svg('cwc_fogt_seas.svg')

Fig8 = tight_subplot(2,3,[0.05 0.03],[0.10 0.1],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig8(1));
wtc(visbeck_SF,NArauca)
title('Araucaria')
xlabel([])
axes(Fig8(2));
wtc(visbeck_SF,NAustro)
title('Austrocedrus')
xlabel([])
axes(Fig8(3));
wtc(visbeck_SF,NHaloca)
title('Halocarpus')
xlabel([])
axes(Fig8(4));
wtc(visbeck_SF,NKoff)
title('Ice dust')
axes(Fig8(5));
wtc(visbeck_SF,NNothof)
title('Nothofagus')
axes(Fig8(6));
wtc(visbeck_MA,NTay)
title('Callitris')

plot2svg('cwc_visbeck_glob_seas.svg')

% Plot the coherence of each proxy with each regional Visbeck SAM (for the right season)
% Antarctic record uses Aus as it is in this sector, broadly speaking (Ehhhh...).

Fig9 = tight_subplot(3,3,[0.05 0.03],[0.10 0.1],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig9(1));
wtc(visbeck_SF_SA,NArauca)
title('Araucaria')
xlabel([])
axes(Fig9(2));
wtc(visbeck_SF_SA,NAustro)
title('Austrocedrus')
xlabel([])
axes(Fig9(3));
wtc(visbeck_SF_AU,NHaloca)
title('Halocarpus')
xlabel([])
axes(Fig9(4));
wtc(visbeck_SF_AU,NKoff)
title('Ice dust')
xlabel([])
axes(Fig9(5));
wtc(visbeck_SF_SA,NNothof)
title('Nothofagus')
axes(Fig9(6));
wtc(visbeck_MA_AU,NTay)
title('Callitris')
axes(Fig9(7));
wtc(visbeck_SF_SA,NMulva)
title('\delta D')

plot2svg('cwc_visbeck_seas.svg')

% Coherence between regional Visbeck SAM indices.

Fig10 = tight_subplot(2,3,[0.05 0.03],[0.10 0.1],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig10(1));
wtc(visbeck_SF_SA,visbeck_SF_AU)
title('S. America - Aus/NZ')
xlabel([])
axes(Fig10(2));
wtc(visbeck_SF_SA,visbeck_SF_AF)
title('S. America - Africa')
xlabel([])
axes(Fig10(3));
wtc(visbeck_SF_AU,visbeck_SF_AF)
title('Aus/NZ - Africa')
xlabel([])
axes(Fig10(4));
wtc(visbeck_SF_SA,visbeck_SF)
title('S. America - S. Hemisphere')
axes(Fig10(5));
wtc(visbeck_SF_AU,visbeck_SF)
title('Aus/NZ - S. Hemisphere')
axes(Fig10(6));
wtc(visbeck_SF_AF,visbeck_SF)
title('Africa - S. Hemisphere')

plot2svg('cwc_visbeck_regions_SF.svg')

%% CWC for SAM indices and recons

Fig1 = tight_subplot(6,3,[0.05 0.03],[0.10 0.01],[0.1 0.01]); %First set of square brackets lets you adjust gaps between subplots
axes(Fig1(1));
wtc(Marshall_SAM(:,[1 2]),CPS_M(:,[1 3]))
axes(Fig1(2));
wtc(FJ_ann,CPS_M(:,[1 3]))
axes(Fig1(3));
wtc(Visbeck_Ann,CPS_M(:,[1 3]))

axes(Fig1(4));
wtc(Marshall_SAM(:,[1 2]),recon_M(:,[1 3]))
axes(Fig1(5));
wtc(FJ_ann,recon_M(:,[1 3]))
axes(Fig1(6));
wtc(Visbeck_Ann,recon_M(:,[1 3]))

axes(Fig1(7));
wtc(Marshall_SAM(:,[1 2]),CPS_FJ(:,[1 3]))
axes(Fig1(8));
wtc(FJ_ann,CPS_FJ(:,[1 3]))
axes(Fig1(9));
wtc(Visbeck_Ann,CPS_FJ(:,[1 3]))

axes(Fig1(10));
wtc(Marshall_SAM(:,[1 2]),recon_FJ(:,[1 3]))
axes(Fig1(11));
wtc(FJ_ann,recon_FJ(:,[1 3]))
axes(Fig1(12));
wtc(Visbeck_Ann,recon_FJ(:,[1 3]))

axes(Fig1(13));
wtc(Marshall_SAM(:,[1 2]),CPS_V(:,[1 3]))
axes(Fig1(14));
wtc(FJ_ann,CPS_V(:,[1 3]))
axes(Fig1(15));
wtc(Visbeck_Ann,CPS_V(:,[1 3]))

axes(Fig1(16));
wtc(Marshall_SAM(:,[1 2]),recon_V(:,[1 3]))
axes(Fig1(17));
wtc(FJ_ann,recon_V(:,[1 3]))
axes(Fig1(18));
wtc(Visbeck_Ann,recon_V(:,[1 3]))

plot2svg('recon_SAM_wtc.svg')



