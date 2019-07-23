% Outputs the indices for the start and end of an overlapping period of
% time for two time-series'
% Requires the wavelet toolbox by Grinsted '04
% Time series' must be formatted the same way, from long ago to present.

function [zts1, zts2] = overlap(ts1,ts2)

overlap_1(2,1) = nan;
overlap_2(2,1) = nan;

[x,dt]=formatts(ts1);
[y,~]=formatts(ts2);
t=(max(x(1,1),y(1,1)):dt:min(x(end,1),y(end,1)))'; %common time period

[overlap_1(1,1), ~] = find(ts1==(min(t)));
[overlap_2(1,1), ~] = find(ts2==(min(t)));

[overlap_1(2,1), ~] = find(ts1==(max(t)));
[overlap_2(2,1), ~] = find(ts2==(max(t)));

zts1(:,1) = ts1(overlap_1(1,1):overlap_1(2,1),1);
zts1(:,2) = zscore(ts1(overlap_1(1,1):overlap_1(2,1),2));
zts2(:,1) = ts2(overlap_2(1,1):overlap_2(2,1),1);
zts2(:,2) = zscore(ts2(overlap_2(1,1):overlap_2(2,1),2));

end
