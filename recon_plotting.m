% Plotting the different SAM reconstructions

clear
load JAGS_out.mat
load CPS_results.mat
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load Fogt_Jones.mat;
load('SAM_seasonal.mat','Visbeck_Ann')

%% Scale variance of recons
% zrecon_FJ = zscore(recon_FJ); zrecon_FJ(:,1) = recon_FJ(:,1);
% zrecon_H_FJ = zscore(recon_H_FJ); zrecon_H_FJ(:,1) = recon_H_FJ(:,1);
% zrecon_L_FJ = zscore(recon_L_FJ); zrecon_L_FJ(:,1) = recon_L_FJ(:,1);
% zrecon_M = zscore(recon_M); zrecon_M(:,1) = recon_M(:,1);
% zrecon_H_M = zscore(recon_H_M); zrecon_H_M(:,1) = recon_H_M(:,1);
% zrecon_L_M = zscore(recon_L_M); zrecon_L_M(:,1) = recon_L_M(:,1);
% zrecon_V = zscore(recon_V); zrecon_V(:,1) = recon_V(:,1);
% zrecon_H_V = zscore(recon_H_V); zrecon_H_V(:,1) = recon_H_V(:,1);
% zrecon_L_V = zscore(recon_L_V); zrecon_L_V(:,1) = recon_L_V(:,1);
zCPS_FJ = zscore(CPS_FJ(:,3));
zCPS_M = zscore(CPS_M(:,3));
zCPS_V = zscore(CPS_V(:,3));

% recon_FJ_2(:,2) = std_correct(FJ_ann(10:100,2),recon_FJ(:,2));
% recon_H_FJ_2 = std_correct_CPS(FJ_ann(10:100,2),recon_H_FJ(:,2));
% recon_L_FJ_2 = std_correct_CPS(FJ_ann(10:100,2),recon_L_FJ(:,2));
% recon_M_2(:,2) = std_correct(Marshall_SAM(21:59,2),recon_M(:,2));
% recon_H_M_2 = std_correct_CPS(Marshall_SAM(21:59,2),recon_H_M(:,2));
% recon_L_M_2 = std_correct_CPS(Marshall_SAM(21:59,2),recon_L_M(:,2));
% recon_V_2(:,2) = std_correct(Visbeck_Ann(11:119,2),recon_V(:,2));
% recon_H_V_2 = std_correct_CPS(Visbeck_Ann(11:119,2),recon_H_V(:,2));
% recon_L_V_2 = std_correct_CPS(Visbeck_Ann(11:119,2),recon_L_V(:,2));

% zCPS_FJ(:,2) = std_correct(FJ_ann(10:100,2),zCPS_FJ);
% zCPS_M(:,2) = std_correct(Marshall_SAM(21:59,2),zCPS_M);
% zCPS_V(:,2) = std_correct(Visbeck_Ann(11:119,2),zCPS_V);

recon_FJ_7 = smooth(recon_FJ(:,2),7,'rloess'); recon_FJ_71 = smooth(recon_FJ(:,2),71,'rloess'); recon_FJ_71_scaled = smooth(recon_FJ(:,3),71,'rloess');
recon_M_7 = smooth(recon_M(:,2),7,'rloess'); recon_M_71 = smooth(recon_M(:,2),71,'rloess'); recon_M_71_scaled = smooth(recon_M(:,3),71,'rloess'); 
recon_V_7 = smooth(recon_V(:,2),7,'rloess');  recon_V_71 = smooth(recon_V(:,2),71,'rloess'); recon_V_71_scaled = smooth(recon_V(:,3),71,'rloess'); 


zCPS_FJ_71 = smooth(zCPS_FJ(:,2),71,'rloess');
zCPS_M_71 = smooth(zCPS_M(:,2),71,'rloess');
zCPS_V_71 = smooth(zCPS_V(:,2),71,'rloess');

%%

AbramSAM(:,2) = zscore(AbramSAM(:,2));

subplot(4,1,1)
jbfill(recon_FJ(:,1)',recon_H_FJ(:,2)',recon_L_FJ(:,2)','m','k','add',0.3);
hold on
plot(recon_FJ(:,1),recon_FJ_71,'m','linewidth',2);
subplot(4,1,2)
jbfill(recon_M(:,1)',recon_H_M(:,2)',recon_L_M(:,2)','r','k','add',0.3);
hold on
plot(recon_M(:,1),recon_M_71,'r','linewidth',2);
subplot(4,1,3)
jbfill(recon_V(:,1)',recon_H_V(:,2)',recon_L_V(:,2)','b','k','add',0.3);
hold on
plot(recon_V(:,1),recon_V_71,'b','linewidth',2);
subplot(4,1,4)
plot(recon_FJ(:,1),recon_FJ_71_scaled,'m','linewidth',2);
hold on
plot(recon_M(:,1),recon_M_71_scaled,'r','linewidth',2);
plot(recon_V(:,1),recon_V_71_scaled,'b','linewidth',2);
plot(AbramSAM(:,1),AbramSAM(:,2),'k','linewidth',2);

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
plot(CPS_FJ(:,1),zCPS_FJ(:,2),'k','linewidth',1)
hold
plot(CPS_FJ(:,1),zCPS_FJ_71,'m','linewidth',2);
subplot(4,1,2)
plot(CPS_M(:,1),zCPS_M(:,2),'k','linewidth',1)
hold
plot(CPS_M(:,1),zCPS_M_71,'r','linewidth',2);
axis([1000 2000 -5 5]);
subplot(4,1,3)
plot(CPS_V(:,1),zCPS_V(:,2),'k','linewidth',1)
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
corr(zCPS_M(1:39,2),Marshall_SAM(21:59,2))
corr(zCPS_M(1:109,2),Visbeck_Ann(11:119,2))
corr(zCPS_M(1:91,2),FJ_ann(10:100,2))

corr(zCPS_FJ(1:39,2),Marshall_SAM(21:59,2))
corr(zCPS_FJ(1:109,2),Visbeck_Ann(11:119,2))
corr(zCPS_FJ(1:91,2),FJ_ann(10:100,2))

corr(zCPS_V(1:39,2),Marshall_SAM(21:59,2))
corr(zCPS_V(1:109,2),Visbeck_Ann(11:119,2))
corr(zCPS_V(1:91,2),FJ_ann(10:100,2))

corr(zCPS_M(:,2),AbramSAM(:,3))
corr(zCPS_FJ(:,2),AbramSAM(:,3))
corr(zCPS_V(:,2),AbramSAM(:,3))

corr(zrecon_M(:,2),AbramSAM(:,3)) 
corr(zrecon_FJ(:,2),AbramSAM(:,3))
corr(zrecon_V(:,2),AbramSAM(:,3)) 

corr(zrecon_M(:,2),CPS_M(:,2))
corr(zrecon_FJ(:,2),CPS_FJ(:,2))
corr(zrecon_V(:,2),CPS_V(:,2))