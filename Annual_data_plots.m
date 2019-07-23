% Annual Data plots
%% Normalising Data
clear
load ('Koffman_2014.mat')
load ('Villalba2012.mat')
load ('Cullen_Grierson09.mat')
load ('EOF1_2.mat')
load ('Law_Dome_CO2.mat')
load ('Steig.mat')

Koffman_annual = flipud(Koffman_annual);
Koffman_annual_noNaNs = Koffman_annual(~any(isnan(Koffman_annual),2),:);

%% Normalise Data - Changed! Now normalised only over period 1600 - 1998

[rows,~]=size(Koffman_annual(601:999,2));
colMax=max(abs(Koffman_annual(601:999,2)),[],1); %# take max absolute value to account for negative numbers
NKoff=Koffman_annual(601:999,2)./repmat(colMax,rows,1);

min = min(Koffman_annual(601:999,2));
max = max(Koffman_annual(601:999,2));
NKoff2 = ((Koffman_annual(601:999,2) - repmat(min,rows,1))./repmat(max-min,rows,1));

[rows,~]=size(Lake_Tay(1:397,2)); 
colMax=max(abs(Lake_Tay(1:397,2)),[],1); %# take max absolute value to account for negative numbers
NTay=Lake_Tay(1:397,2)./repmat(colMax,rows,1);

[rows,~]=size(Villalba2012(192:590,2)); 
colMax=max(abs(Villalba2012(192:590,2)),[],1); %# take max absolute value to account for negative numbers
NAustro=Villalba2012(192:590,2)./repmat(colMax,rows,1);

[rows,~]=size(Villalba2012(192:590,3)); 
colMax=max(abs(Villalba2012(192:590,3)),[],1); %# take max absolute value to account for negative numbers
NArauca=Villalba2012(192:590,3)./repmat(colMax,rows,1);

[rows,~]=size(Villalba2012(192:590,4));
colMax=max(abs(Villalba2012(192:590,4)),[],1); %# take max absolute value to account for negative numbers
NHaloca=Villalba2012(192:590,4)./repmat(colMax,rows,1);

[rows,~]=size(Villalba2012(601:999,7)); 
colMax=max(abs(Villalba2012(601:999,7)),[],1); %# take max absolute value to account for negative numbers
NNothof=Villalba2012(601:999,7)./repmat(colMax,rows,1);

[rows,~]=size(CO2_LawD(1600:1998,2)); 
colMax=max(abs(CO2_LawD(1600:1998,2)),[],1); %# take max absolute value to account for negative numbers
NCO2=CO2_LawD(1600:1998,2)./repmat(colMax,rows,1);

[rows,~]=size(EOF1(1:397,2)); 
colMax=max(abs(EOF1(1:397,2)),[],1); %# take max absolute value to account for negative numbers
NEOF1=EOF1(1:397,2)./repmat(colMax,rows,1);

[rows,~]=size(EOF2(1:397,2)); 
colMax=max(abs(EOF2(1:397,2)),[],1); %# take max absolute value to account for negative numbers
NEOF2=EOF2(1:397,2)./repmat(colMax,rows,1);

[rows,~]=size(Mulvaney(:,2));
min = min(Mulvaney(:,2));
max = max(Mulvaney(:,2));
NMulva(:,1) = Mulvaney(:,1);
NMulva(:,2) = ((Mulvaney(:,2) - repmat(min,rows,1))./repmat(max-min,rows,1));

% Normalise the new Steig dataset (Figure out which ones to use later)

NSteig = nan(1002,17);
NSteig(:,1) = Steig(:,1);
for i = 2:17
    [rows,~]=size(Steig(:,i));
    min = min(Steig(:,i));
    max = max(Steig(:,i));
    NSteig(:,i) = ((Steig(:,i) - repmat(min,rows,1))./repmat(max-min,rows,1));
    clear min max
end

% Alternatively, normalise them over a common interval

NSteig_common = nan(91,17);
NSteig_common(:,1) = Steig(19:109,1);
for i = 2:17
    [rows,~]=size(Steig(19:109,i));
    min = min(Steig(19:109,i));
    max = max(Steig(19:109,i));
    NSteig_common(1:91,i) = ((Steig(19:109,i) - repmat(min,rows,1))./repmat(max-min,rows,1));
    clear min max
end

save('Steig.mat','NSteig_common','-append')

% All ice core data, normalised over common interval 1893-1983

NAll_ice_common = nan(91,20);
NAll_ice_common(:,1) = All_ice_common(1:91,1);
for i = 2:20
    [rows,~]=size(All_ice_common(1:91,i));
    min = min(All_ice_common(1:91,i));
    max = max(All_ice_common(1:91,i));
    NAll_ice_common(1:91,i) = ((All_ice_common(1:91,i) - repmat(min,rows,1))./repmat(max-min,rows,1));
    clear min max
end


% Calculate running means for different temporal windows
  
Koff_11 = tsmovavg(NKoff(1:399,1),'s',11,1);
Tay_11 = tsmovavg(NTay(1:397,1),'s',11,1);
Austro_11 = tsmovavg(NAustro(1:399,1),'s',11,1);
Aurauc_11 = tsmovavg(NArauca(1:399,1),'s',11,1);
Haloca_11 = tsmovavg(NHaloca(1:399,1),'s',11,1);
Nothof_11 = tsmovavg(NNothof(1:399,1),'s',11,1);
CO2_11    = tsmovavg(NCO2(1:399,1),'s',11,1);
EOF1_11   = tsmovavg(NEOF1(1:397,1),'s',11,1);
EOF2_11   = tsmovavg(NEOF2(1:397,1),'s',11,1);

%% Plotting

figure(1)
subplot(9,1,1)
plot(Koffman_annual(601:999,1),NKoff(1:399,1))
hold
plot(Koffman_annual(606:994,1),Koff_11(11:399,1),'color','black','linewidth',2)
axis([1600 2010 0 1])
line([1750 1750],[0 1],'linestyle','--','color','black')
line([1850 1850],[0 1],'linestyle','--','color','black')
title('Ice dust coarse particle %')

subplot(9,1,2)
plot(Lake_Tay(1:397,1),NTay(1:397,1),'Color',[0 0.498039215686275 0])
hold
plot(Lake_Tay(6:392,1),Tay_11(11:397,1),'Color','black','linewidth',2)
axis([1600 2010 0 1])
line([1750 1750],[0 1],'linestyle','--','color','black')
line([1850 1850],[0 1],'linestyle','--','color','black')
title('Lake Tay tree ring index')

subplot(9,1,3)
plot(Villalba2012(192:590,1),NAustro(1:399,1),'Color',[1 1 0])
hold
plot(Villalba2012(197:585,1),Austro_11(11:399,1),'Color','black','linewidth',2)
axis([1600 2010 0 1])
line([1750 1750],[0 1],'linestyle','--','color','black')
line([1850 1850],[0 1],'linestyle','--','color','black')
title('Austrocedrus tree ring index - SA')

subplot(9,1,4)
plot(Villalba2012(192:590,1),NArauca(1:399,1),'Color',[1 0 0])
hold
plot(Villalba2012(197:585,1),Aurauc_11(11:399,1),'Color','black','linewidth',2)
axis([1600 2010 0 1])
line([1750 1750],[0 1],'linestyle','--','color','black')
line([1850 1850],[0 1],'linestyle','--','color','black')
title('Araucaria tree ring index - SA')

subplot(9,1,5)
plot(Villalba2012(192:590,1),NHaloca(1:399,1),'Color',[1 0.694117647058824 0.392156862745098])
hold
plot(Villalba2012(197:585,1),Haloca_11(11:399,1),'Color','black','linewidth',2)
axis([1600 2010 0 1])
line([1750 1750],[0 1],'linestyle','--','color','black')
line([1850 1850],[0 1],'linestyle','--','color','black')
title('Halocarpus tree ring index - NZ')

subplot(9,1,6)
plot(Villalba2012(601:999,6),NNothof(1:399,1),'Color',[0.47843137254902 0.0627450980392157 0.894117647058824])
hold
plot(Villalba2012(606:994,6),Nothof_11(11:399,1),'Color','black','linewidth',2)
axis([1600 2010 0 1])
line([1750 1750],[0 1],'linestyle','--','color','black')
line([1850 1850],[0 1],'linestyle','--','color','black')
title('Nothofagus tree ring index - SA')

subplot(9,2,13)
plot(CO2_LawD(1600:1859,1),NCO2(1:260,1))
%hold
%plot(CO2_LawD(1:2003,1),CO2_10(1:2003,1),'Color','black','linewidth',2)
axis([1600 1860 0.74 0.8])
line([1750 1750],[0 1],'linestyle','--','color','black')
line([1850 1850],[0 1],'linestyle','--','color','black')
title('CO2 Law Dome')

subplot(9,2,14)
plot(CO2_LawD(1860:1998,1),NCO2(261:399,1))
axis([1860 2010 0.76 1])
set(gca,'YAxisLocation','right')

subplot(9,1,8)
plot(EOF1(1:397,1),NEOF1(1:397,1),'color','cyan')
hold
plot(EOF1(6:392,1),EOF1_11(11:397,1),'color','black','linewidth',2)
axis([1600 2010 -1 1])
line([1750 1750],[-1 1],'linestyle','--','color','black')
line([1850 1850],[-1 1],'linestyle','--','color','black')
title('PCA - Axis 1')

subplot(9,1,9)
plot(EOF2(1:397,1),NEOF2(1:397,1),'color','cyan')
hold
plot(EOF2(6:392,1),EOF2_11(11:397,1),'color','black','linewidth',2)
axis([1600 2010 -1 1])
line([1750 1750],[-1 1],'linestyle','--','color','black')
line([1850 1850],[-1 1],'linestyle','--','color','black')
title('PCA - Axis 2')
xlabel('Years')






























