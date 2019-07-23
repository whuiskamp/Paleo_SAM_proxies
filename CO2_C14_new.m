% Plotting 14C/CO2 for the last 1000 years
% New C14 data is from the SHCal13 record by Hogg et al 2013
% New high res CO2 data from Ahn et al 2012

load('CO2_1000_new.mat')
import SHCal13  %This won't work, but hey ...

[AX,H1,H2] = plotyy(data(4941:5141,1),data(4941:5141,4),untitled(1:109,1),untitled(1:109,2));
set(get(AX(1),'Ylabel'),'String','\Delta ^{14}C (Permill)','color','blue')
set(get(AX(2),'Ylabel'),'String','CO2','color','red')  
set(AX(2), 'YColor', 'red') %DC14 is already blue...

xlabel('Years BP')

line([100 100],[-270 320],'linestyle','--','color','black','linewidth',1);
line([200 200],[-270 320],'linestyle','--','color','black','linewidth',1);

line([350 350],[-270 320],'linestyle','--','color','black','linewidth',1);
line([450 450],[-270 320],'linestyle','--','color','black','linewidth',1);

line([500 500],[-270 320],'linestyle','--','color','black','linewidth',1);
line([600 600],[-270 320],'linestyle','--','color','black','linewidth',1);