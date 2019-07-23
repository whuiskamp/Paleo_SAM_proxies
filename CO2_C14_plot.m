% CO2 / 14C relationship plot
data(1:4940,:) = []; 
YRBP(:,2) = perMil1(:);
C14(:,1) = C14(:,1)-1950; C14(:,1) = C14(:,1)*-1;

CO2_LAWD = flipud(CO2_LAWD);

C14(:,1) = data(:,1); C14(:,2) = data(:,4);
C14_spline(:,1) = 1000:1:1950;
C14_spline(:,2) = interp1(C14(:,1),C14(:,2),1000:1:1950);
CO2_spline = CO2_LawD_splinefit(1000:1950,:);

plot(zscore(C14_spline(:,2)),zscore(CO2_spline(:,2)));
line([0 0],[-2 5],'linestyle','--','color','k')
line([-2 2],[0 0],'linestyle','--','color','k')

plotyy(C14_spline(:,1),zscore(C14_spline(:,2)),CO2_LAWD(29:78,1),zscore(CO2_LAWD(29:78,2)));

%% 

load marshall_SAM

marsh_S = smooth(Marshall_SAM(:,2),21,'rloess');

plot(Marshall_SAM(1:59,1),Marshall_SAM(1:59,2))
hold
plot(Marshall_SAM(1:59,1),marsh_S(1:59,1),'linewidth',2,'color','r')