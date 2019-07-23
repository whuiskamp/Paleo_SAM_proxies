%Tree Ring data preliminary plots
% ave files are created using a 10 year running mean

import 'Tree rings.mat'
figure(1)
plot(Hihitahi(1:452,1),H_avg(1:452,2))
hold
plot(Oroko(1:1000,1),O_avg(1:1000,2),'Color',[1 0 1])
plot(Takapari(1:628,1),T_avg(1:628,2),'Color',[0 0 0])

line([1750 1750],[0.4 1.8],'linestyle','--','color','black','linewidth',1);
line([1850 1850],[0.4 1.8],'linestyle','--','color','black','linewidth',1);

line([1500 1500],[0.4 1.8],'linestyle','--','color','black','linewidth',1);
line([1600 1600],[0.4 1.8],'linestyle','--','color','black','linewidth',1);

line([1350 1350],[0.4 1.8],'linestyle','--','color','black','linewidth',1);
line([1450 1450],[0.4 1.8],'linestyle','--','color','black','linewidth',1);

axis([1200 2000 0.5 1.6])

figure(2)
plot(Hihitahi(1:452,1),H_avg_20(1:452,2))
hold
plot(Oroko(1:1000,1),O_avg_20(1:1000,2),'Color',[1 0 1])
plot(Takapari(1:628,1),T_avg_20(1:628,2),'Color',[0 0 0])

line([1750 1750],[0.4 1.8],'linestyle','--','color','black','linewidth',1);
line([1850 1850],[0.4 1.8],'linestyle','--','color','black','linewidth',1);

line([1500 1500],[0.4 1.8],'linestyle','--','color','black','linewidth',1);
line([1600 1600],[0.4 1.8],'linestyle','--','color','black','linewidth',1);

line([1350 1350],[0.4 1.8],'linestyle','--','color','black','linewidth',1);
line([1450 1450],[0.4 1.8],'linestyle','--','color','black','linewidth',1);

axis([1200 2000 0.7 1.45])
