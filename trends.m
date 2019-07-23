% Trends

load('CPS_results.mat')
load('JAGS_out.mat')
load('EOF_SAM.mat')
load('CO2.mat')

wdw = 101;

for i = ceil(wdw/2):995-floor(wdw/2)
    p(i,:) = polyfit(CPS_FJ(i-floor(wdw/2):i+floor(wdw/2),1),CPS_FJ(i-floor(wdw/2):i+floor(wdw/2),3),1);
    %y() = polyval(p,CPS_FJ(1:40,1));
end

p(946:995,1) = NaN;

p(find(p==0)) = nan;
plot(CPS_FJ(:,1),p(:,1))

%EOF of recons

matrix(:,1) = (zscore(recon_FJ(:,3)));
matrix(:,2) = (zscore(recon_M(:,3)));
matrix(:,3) = (zscore(recon_V(:,3)));
matrix(:,4) = (zscore(CPS_FJ(:,3)));
matrix(:,5) = (zscore(CPS_M(:,3)));
matrix(:,6) = (zscore(CPS_V(:,3)));

[COEFF,SCORE,rotscores,latent,tsquare,var_explained] = EOF_calc(matrix);

%% Plotting
final_recon_71 = smooth(rotscores(:,1),71,'rloess');
zZhang_71 = smooth(zscore(Zhang_SAM(:,2)),71,'rloess');
zVillalba_71 = smooth(zscore(Villalba_SAM(:,2)),71,'rloess');
zAbram_71 = zscore(AbramSAM(:,2));

subplot(3,1,1)
plot(CPS_FJ(:,1),final_recon_71(:,1),'k','linewidth',2);
hold on
plot(AbramSAM(:,1),zAbram_71(:,1),'b','linewidth',2);
plot(Zhang_SAM(:,1),zZhang_71(:,1),'r','linewidth',2);
plot(Villalba_SAM(:,1),zVillalba_71(:,1),'m','linewidth',2);
axis([1000 2010 -3 3])

subplot(3,2,3)
scatter(CO2_LAWD(:,1),CO2_LAWD(:,2),'r','linewidth',2)
hold on
plot(CO2_LawD_splinefit(:,1),CO2_LawD_splinefit(:,2),'r','linewidth',2)
scatter(CO2_WAIS(:,1),CO2_WAIS(:,2),'b','linewidth',2)
axis([1000 1860 270 290])

subplot(3,2,4)
scatter(CO2_LAWD(:,1),CO2_LAWD(:,2),'r','linewidth',2)
hold on
plot(CO2_LawD_splinefit(:,1),CO2_LawD_splinefit(:,2),'r','linewidth',2)
scatter(CO2_WAIS(:,1),CO2_WAIS(:,2),'b','linewidth',2)
axis([1860 2010 290 380])
%set(gca,'XAxisLocation','bottom')
%set(gca,'XAxisLocation','top')
line([1860 1860],[290 380],'linestyle','--','color','k')
set(gca,'YAxisLocation','right')
plot2svg('SAM_CO2_1.svg')
clf
subplot(3,1,3)
%Plot this on a separate figure - matlab fucks up.
wtc(EOF_SAM,CO2_LawD_splinefit)
plot2svg('SAM_CO2.svg')

%Stats for the EOF recon

corr(EOF_SAM(1:39,2),Marshall_SAM(21:59,2)) % 0.7088 with ci of 0.5114 - 0.8249
corr(EOF_SAM(1:91,2),FJ_ann(10:100,2)) % 0.4333 with ci of 0.2714 - 0.5708
corr(EOF_SAM(1:109,2),Visbeck_Ann(11:119,2)) % 0.4656 with ci of 0.2933 - 0.6072

corr(EOF_SAM(:,2),AbramSAM(:,3)) % 0.4795 with ci of 0.4293 - 0.5258
corr(EOF_SAM(1:587,2),Villalba_SAM(12:598,2)) %  -0.0613
corr(EOF_SAM(1:496,2),Zhang_SAM(13:508,2)) % 0.0667

corr(CPS_M(:,2),AbramSAM(:,3)) % 0.4943
corr(CPS_M(1:587,2),Villalba_SAM(12:598,2)) % 0.1305
corr(CPS_M(1:496,2),Zhang_SAM(13:508,2)) % 0.0760

corr(CPS_FJ(:,2),AbramSAM(:,3)) % 0.4918
corr(CPS_FJ(1:587,2),Villalba_SAM(12:598,2)) % 0.1494
corr(CPS_FJ(1:496,2),Zhang_SAM(13:508,2)) % 0.1575

corr(CPS_V(:,2),AbramSAM(:,3)) % 0.0791
corr(CPS_V(1:587,2),Villalba_SAM(12:598,2)) % -0.2208
corr(CPS_V(1:496,2),Zhang_SAM(13:508,2)) % 0.0030

corr(recon_M(:,3),AbramSAM(:,3)) % 0.3422
corr(recon_M(1:587,3),Villalba_SAM(12:598,2)) % -0.1459
corr(recon_M(1:496,3),Zhang_SAM(13:508,2)) % 0.0705

corr(recon_FJ(:,3),AbramSAM(:,3)) % 0.4491
corr(recon_FJ(1:587,3),Villalba_SAM(12:598,2)) % 0.1744
corr(recon_FJ(1:496,3),Zhang_SAM(13:508,2)) % 0.0386

corr(recon_V(:,3),AbramSAM(:,3)) % 0.0365
corr(recon_V(1:587,3),Villalba_SAM(12:598,2)) % -0.054 
corr(recon_V(1:496,3),Zhang_SAM(13:508,2)) % -0.0059

rhos1000 = bootstrp(1000,'corr',EOF_SAM(:,2),AbramSAM(:,3));
ci = bootci(5000,@corr,EOF_SAM(:,2),AbramSAM(:,3))
ci = bootci(5000,@corr,EOF_SAM(1:39,2),Marshall_SAM(21:59,2))
ci = bootci(5000,@corr,EOF_SAM(1:91,2),FJ_ann(10:100,2))
ci = bootci(5000,@corr,EOF_SAM(1:109,2),Visbeck_Ann(11:119,2))


