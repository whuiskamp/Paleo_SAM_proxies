% Plotting the different SAM reconstructions

clear
load JAGS_out.mat
load CPS_results.mat
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load Fogt_Jones.mat;
load('SAM_seasonal.mat','Visbeck_Ann')
recon_FJ_71 = smooth(recon_FJ(:,2),71,'rloess'); zrecon_FJ_71 = smooth(zscore(recon_FJ(:,2)),71,'rloess'); zrecon_FJ_71_scaled = smooth(zscore(recon_FJ(:,3)),71,'rloess');
recon_M_71 = smooth(recon_M(:,2),71,'rloess'); zrecon_M_71 = smooth(zscore(recon_M(:,2)),71,'rloess'); zrecon_M_71_scaled = smooth(zscore(recon_M(:,3)),71,'rloess');
recon_V_71 = smooth(recon_V(:,2),71,'rloess'); zrecon_V_71 = smooth(zscore(recon_V(:,2)),71,'rloess'); zrecon_V_71_scaled = smooth(zscore(recon_V(:,3)),71,'rloess');
%% Plotting

zCPS_FJ = zscore(CPS_FJ(:,3)); zCPS_FJ_71 = smooth(zCPS_FJ,71,'rloess');
zCPS_M = zscore(CPS_M(:,3)); zCPS_M_71 = smooth(zCPS_M,71,'rloess');
zCPS_V = zscore(CPS_V(:,3)); zCPS_V_71 = smooth(zCPS_V,71,'rloess');

AbramSAM(:,2) = zscore(AbramSAM(:,2));

subplot(4,1,1)
jbfill(recon_FJ(:,1)',recon_H_FJ(:,2)',recon_L_FJ(:,2)','m','m','add',0.3);
hold on
plot(recon_FJ(:,1),recon_FJ(:,2),'k');
plot(recon_FJ(:,1),recon_FJ_71,'m','linewidth',2);
subplot(4,1,2)
jbfill(recon_M(:,1)',recon_H_M(:,2)',recon_L_M(:,2)','r','r','add',0.3);
hold on
plot(recon_M(:,1),recon_M(:,2),'k');
plot(recon_M(:,1),recon_M_71,'r','linewidth',2);
subplot(4,1,3)
jbfill(recon_V(:,1)',recon_H_V(:,2)',recon_L_V(:,2)','b','b','add',0.3);
hold on
plot(recon_V(:,1),recon_V(:,2),'k');
plot(recon_V(:,1),recon_V_71,'b','linewidth',2);
subplot(4,1,4)
plot(recon_FJ(:,1),zrecon_FJ_71_scaled,'m','linewidth',2);
hold on
plot(recon_M(:,1),zrecon_M_71_scaled,'r','linewidth',2);
plot(recon_V(:,1),zrecon_V_71_scaled,'b','linewidth',2);
plot(AbramSAM(:,1),AbramSAM(:,2),'k','linewidth',2);

plot2svg('JAGS_recons.svg')

% CPS Skill
xaxis = 10:52;
figure(2)
jbfill(xaxis,skill_FJ_qn(10:52,3)',skill_FJ_qn(10:52,1)','m','k','add',0.3);
hold on
jbfill(xaxis,skill_M_qn(10:52,3)',skill_M_qn(10:52,1)','r','k','add',0.3);
jbfill(xaxis,skill_V_qn(10:52,3)',skill_V_qn(10:52,1)','b','k','add',0.3);
hold on
plot(skill_FJ_qn(:,2),'m','linewidth',2)
plot(skill_M_qn(:,2),'r','linewidth',2)
plot(skill_V_qn(:,2),'b','linewidth',2)
line([10 10],[0.1 0.9],'linestyle','--','color','k')
line([20 20],[0.1 0.9],'linestyle','--','color','k')
line([30 30],[0.1 0.9],'linestyle','--','color','k')
line([40 40],[0.1 0.9],'linestyle','--','color','k')
line([50 50],[0.1 0.9],'linestyle','--','color','k')

line([0 60],[0.2 0.2],'linestyle','--','color','k')
line([0 60],[0.3 0.3],'linestyle','--','color','k')
line([0 60],[0.4 0.4],'linestyle','--','color','k')
line([0 60],[0.5 0.5],'linestyle','--','color','k')
line([0 60],[0.6 0.6],'linestyle','--','color','k')
line([0 60],[0.7 0.7],'linestyle','--','color','k')
line([0 60],[0.8 0.8],'linestyle','--','color','k')

plot2svg('CPS_skill.svg')

% CPS recons

figure(3)
subplot(4,1,1)
axis([1000 2000 -4 4])
plot(CPS_FJ(:,1),zCPS_FJ,'k','linewidth',1)
hold
plot(CPS_FJ(:,1),zCPS_FJ_71,'m','linewidth',2);
subplot(4,1,2)
axis([1000 2000 -4 4])
plot(CPS_M(:,1),zCPS_M,'k','linewidth',1)
hold
plot(CPS_M(:,1),zCPS_M_71,'r','linewidth',2);
axis([1000 2000 -4 4]);
subplot(4,1,3)
axis([1000 2000 -4 4])
plot(CPS_V(:,1),zCPS_V,'k','linewidth',1)
hold
plot(CPS_V(:,1),zCPS_V_71,'b','linewidth',2);
subplot(4,1,4)
plot(AbramSAM(:,1),AbramSAM(:,2),'k','linewidth',2);
hold
plot(CPS_FJ(:,1),zCPS_FJ_71,'m','linewidth',2);
plot(CPS_M(:,1),zCPS_M_71,'r','linewidth',2);
plot(CPS_V(:,1),zCPS_V_71,'b','linewidth',2);

plot2svg('CPS_recons_2.svg')

%% Skill of each recon.
corr(recon_M(1:39,2),Marshall_SAM(21:59,2)) % 0.7394 
corr(recon_M(1:109,2),Visbeck_Ann(11:119,2)) % 0.2579
corr(recon_M(1:91,2),FJ_ann(10:100,2)) % 0.2885

corr(recon_FJ(1:39,2),Marshall_SAM(21:59,2)) % 0.5705
corr(recon_FJ(1:109,2),Visbeck_Ann(11:119,2)) % 0.2821
corr(recon_FJ(1:91,2),FJ_ann(10:100,2)) % 0.4174

corr(recon_V(1:39,3),Marshall_SAM(21:59,2)) % 0.4644
corr(recon_V(1:109,3),Visbeck_Ann(11:119,2)) % 0.4771
corr(recon_V(1:91,3),FJ_ann(10:100,2)) % 0.1909

corr(CPS_M(1:39,2),Marshall_SAM(21:59,2)) % 0.8125
corr(CPS_M(1:109,2),Visbeck_Ann(11:119,2)) % 0.1808
corr(CPS_M(1:91,2),FJ_ann(10:100,2)) % 0.2605

corr(CPS_FJ(1:39,2),Marshall_SAM(21:59,2)) % 0.6173
corr(CPS_FJ(1:109,2),Visbeck_Ann(11:119,2)) % 0.3148
corr(CPS_FJ(1:91,2),FJ_ann(10:100,2)) % 0.4812

corr(CPS_V(1:39,2),Marshall_SAM(21:59,2)) % 0.4082
corr(CPS_V(1:109,2),Visbeck_Ann(11:119,2)) % 0.5221
corr(CPS_V(1:91,2),FJ_ann(10:100,2)) % 0.2047

corr(CPS_M(:,2),AbramSAM(:,3)) % 0.4943
corr(CPS_FJ(:,2),AbramSAM(:,3)) % 0.4918
corr(CPS_V(:,2),AbramSAM(:,3)) % 0.0791

corr(recon_M(:,3),AbramSAM(:,3)) % 0.3422
corr(recon_FJ(:,3),AbramSAM(:,3)) % 0.4491
corr(recon_V(:,3),AbramSAM(:,3)) % 0.0365

corr(recon_M(:,3),CPS_M(:,3)) % 0.6788
corr(recon_FJ(:,3),CPS_FJ(:,3)) % 0.8127
corr(recon_V(:,3),CPS_V(:,3)) % 0.8080

%% Wavelets of previous SAM recons

wt(AbramSAM(:,[1 3]))
line([1500 1500],[0 256],'linestyle','--','color','white','linewidth',2)
plot2svg('wt_Abram.svg')
figure(2)
wt(Villalba_SAM)
line([1500 1500],[0 256],'linestyle','--','color','white','linewidth',2)
plot2svg('wt_Villalba.svg')
figure(3)
wt(Zhang_SAM)
plot2svg('wt_zhang.svg')

%% Wavelets of Villalba data

clear all
load('All_proxies_commonT')
wtc(NArauca,NAustro)
plot2svg('wtc_Arau_Aust.svg')
wtc(NArauca,NHaloca)
plot2svg('wtc_Arau_Halo')
wtc(NAustro,NHaloca)
plot2svg('wtc_Austro_Halo')

wt(NArauca)
wt(NAustro)
wt(NHaloca)

% SAM timeseries

AbramSAM(:,2) = zscore(AbramSAM(:,2));
Villalba_SAM(:,3) = zscore(Villalba_SAM(:,2)); Villalba_SAM(:,3) =  smooth(Villalba_SAM(:,3),71,'rloess');
Zhang_SAM(:,3) = zscore(Zhang_SAM(:,2)); Zhang_SAM(:,3) =  smooth(Zhang_SAM(:,3),71,'rloess');

plot(AbramSAM(:,1),AbramSAM(:,2),'k','linewidth',2);
hold
plot(Villalba_SAM(:,1),Villalba_SAM(:,3),'m','linewidth',2);
plot(Zhang_SAM(:,1),Zhang_SAM(:,3),'r','linewidth',2);



