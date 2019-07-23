% Data nomalisation 
%% Common time period 1893-1983

clear
load ('Koffman_2014.mat')
load ('Villalba2012.mat')
Villalba2012 = flipud(Villalba2012);
load ('Cullen_Grierson09.mat')
Lake_Tay = flipud(Lake_Tay);
load ('All_Ice.mat','NAll_ice_common') % This includes all the Steig records as well asl JRI (Mulvaney) and Koffman dust record.

NAll_data_common = nan(91,25);
NAll_data_common(:,1) = NAll_ice_common(1:91,1);

NAll_data_common(:,2:20) = NAll_ice_common(:,2:20);

% Villalba data
[rows,~]=size(NAll_ice_common(:,1));
min = min(Villalba2012(435:525,2));
max = max(Villalba2012(435:525,2));
NAll_data_common(:,21) = ((Villalba2012(435:525,2) - repmat(min,rows,1))./repmat(max-min,rows,1));

clear min max
min = min(Villalba2012(435:525,3));
max = max(Villalba2012(435:525,3));
NAll_data_common(:,22) = ((Villalba2012(435:525,3) - repmat(min,rows,1))./repmat(max-min,rows,1));

clear min max
min = min(Villalba2012(435:525,4));
max = max(Villalba2012(435:525,4));
NAll_data_common(:,23) = ((Villalba2012(435:525,4) - repmat(min,rows,1))./repmat(max-min,rows,1));

clear min max
min = min(Villalba2012(26:116,7));
max = max(Villalba2012(26:116,7));
NAll_data_common(:,24) = ((Villalba2012(26:116,7) - repmat(min,rows,1))./repmat(max-min,rows,1));

% Cullen and Grierson
clear min max
min = min(Lake_Tay(23:113,2));
max = max(Lake_Tay(23:113,2));
NAll_data_common(:,25) = ((Lake_Tay(23:113,2) - repmat(min,rows,1))./repmat(max-min,rows,1));

save('NAll_data.mat','NAll_data_common')

%% Same again but with tree rings shifted back one year (But not lake Tay, their record is sensitive to autumn/ winter)

NAll_data_shifted = NAll_data_common;

% Villalba data
min = min(Villalba2012(434:524,2));
max = max(Villalba2012(434:524,2));
NAll_data_shifted(:,21) = ((Villalba2012(434:524,2) - repmat(min,rows,1))./repmat(max-min,rows,1));

clear min max
min = min(Villalba2012(434:524,3));
max = max(Villalba2012(434:524,3));
NAll_data_shifted(:,22) = ((Villalba2012(434:524,3) - repmat(min,rows,1))./repmat(max-min,rows,1));

clear min max
min = min(Villalba2012(434:524,4));
max = max(Villalba2012(434:524,4));
NAll_data_shifted(:,23) = ((Villalba2012(434:524,4) - repmat(min,rows,1))./repmat(max-min,rows,1));

clear min max
min = min(Villalba2012(25:115,7));
max = max(Villalba2012(25:115,7));
NAll_data_shifted(:,24) = ((Villalba2012(25:115,7) - repmat(min,rows,1))./repmat(max-min,rows,1));

save('NAll_data.mat','NAll_data_shifted','-append')
