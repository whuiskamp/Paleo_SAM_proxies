% This function makes it so that the standard deviation and mean of ts2, a
% longer time series, is equal to the mean and standard deviation of ts1.

function [ts_correct] = std_correct(ts1, ts2)
    len = max(find(~isnan(ts2)));
    len2 = max(find(~isnan(ts1)));
    if len < len2
        ts2(len+1:len2) = NaN;
    elseif len > len2
        ts1(len2+1:len) = NaN;
    end
    ts(:,1) = ts1(:); ts(:,2) = ts2(:);
    diffs = nanstd(ts(1:len,1)) / nanstd(ts(1:len,2));
    ts(:,2) = ts(:,2) * diffs;
     diffm = nanmean(ts(1:len,1)) - nanmean(ts(1:len,2));    
     ts(:,2) = ts(:,2) + diffm;                              
    ts_correct(:,1) = ts(:,2);
    clear ts ts1 ts2
end

