%% 20CR seasonal correlations
% Import 20CR data and extract only DJF means
clear
SAT = ncread('air.mon.mean.nc','air')-273.15; % (180,91,24,142) - lon,lat,levels,time
SAT = squeeze(SAT(:,:,1,:)); % Extract only surface values
SLP = ncread('prmsl.mon.mean.nc','prmsl'); % (180,91,164) - lon,lat,time
time = floor(ncread('air.mon.mean.nc','time')/(24*365))+1800;
time = time(1:12:end,:);
lat = ncread('prmsl.mon.mean.nc','lat'); 
lon = ncread('prmsl.mon.mean.nc','lon');
lmask = ncread('20CR_Lmask.nc','LMASK'); % This land mask is NOT exact. It is extracted from the etopo5
% orography file that comes with ferret, regridded onto the 20CR grid. It
% should be good enough for our purposes.

seas = ["DJF","MAM","JJA","SON"]
dima = [.15 .58 .3 .3];
% Extract only one season
for a = seas
    if a == 'DJF'
        start = 12;
    elseif a == 'MAM'
        start = 3;
    elseif a == 'JJA'
        start = 6;
    elseif a == 'SON'
        start = 9;
    end

    for i = start:12:size(SAT,3)-1 % The dataset ends with Dec, not Jan.
        SAT_seas(:,:,ceil(i/12)) = mean(SAT(:,:,i:i+2),3); % Averaging over season, and time is the 3rd dim
        SLP_seas(:,:,ceil(i/12)) = mean(SLP(:,:,i:i+2),3);
    end

    SLP_40 = squeeze(SLP_seas(:,66,:)); SLP_65 = squeeze(mean(SLP_seas(:,78:79,:),2));
    SAM_seas = zscore(squeeze(mean(SLP_40,1)) - squeeze(mean(SLP_65,1)));

    % Mask out the ocean
    lmask2 = nan(180,91,141);
    for i=1:141
        lmask2(:,:,i) = lmask;
    end
    SAT_seas(lmask2==0) = NaN;

    % Create regions for averaging 
    % 1 = Antarctica, 2 = Southern Aus, 3 = NZ, 4 = Sth. America (N), 5 = Sth. America (S)
    
    for region = 1:5
        if region == 1;
            nth = 16;
            sth = 1;
            wst = 1;
            est = 180;
        elseif region == 2;
            nth = 31;
            sth = 25;
            wst = 58;
            est = 78;
        elseif region == 3;
            nth = 28;
            sth = 23;
            wst = 85;
            est = 90;
        elseif region == 4;
            nth = 31;
            sth = 26;
            wst = 144;
            est = 155;
        elseif region == 5;
            nth = 25;
            sth = 19;
            wst = 144;
            est = 155;
        end

        SAT_region = squeeze(nanmean(nanmean(SAT_seas(wst:est,sth:nth,:),1),2));
        % Next, sort your data into years with positive SAM or Negative SAM
        % then calculate regression values for each using Spearman's rank
        % correlation

        [SAM_sort,inx] = sort(SAM_seas,2,'descend');
        SAT_sort = SAT_region(inx);
        % Find where SAM goes from + to -
        chng = find(SAM_sort > 0,1, 'last');
        % regression for +-tive SAM
        P = polyfit(SAM_sort(1:chng)',SAT_sort(1:chng),1);
        yfit = P(1)*SAM_sort(1:chng)+P(2);
        % regression for -tive SAM
        Q = polyfit(SAM_sort(chng+1:end)',SAT_sort(chng+1:end),1);
        yfit_2 = Q(1)*SAM_sort(chng+1:end)+Q(2);

        figure(region)
        subplot(2,2,find(contains(seas,a)))
        scatter(SAM_seas,SAT_region,'filled'); hold on
        line([0 0],[min(SAT_region) max(SAT_region)], 'linestyle','--','color','k')
        plot(SAM_sort(1:chng),yfit,'r-'); % Regression for +-tive SAM
        plot(SAM_sort(chng+1:end),yfit_2,'b-'); % Regression for -tive SAM

        xlabel('SAM index')
        ylabel('SAT')
        if region == 1
            ann = annotation('textbox',dima,'String','Antarctica','FitBoxToText','on');
        elseif region == 2
            ann = annotation('textbox',dima,'String','Southern Aus.','FitBoxToText','on');
        elseif region == 3
            ann = annotation('textbox',dima,'String','New Zealand','FitBoxToText','on');
        elseif region == 4
            ann = annotation('textbox',dima,'String','Sth. American (30S-42S)','FitBoxToText','on');
        elseif region == 5
            ann = annotation('textbox',dima,'String','Sth. American (42S-60S)','FitBoxToText','on');
        end
        ann.LineStyle = 'none';
        if a == "DJF"
            title('DJF')
        elseif a == "MAM"
            title('MAM')
        elseif a == "JJA"
            title('JJA')
        elseif a == "SON"
            title('SON')
        end
    end
    clear SAT_seas SLP_seas
end


% What if we want a scatter plot for each grid cell?
% It makes a big mess...
figure(1); hold on
for i = sth:nth
    for j= wst:est
        if isnan(squeeze(SAT_seas(j,i,1))) == 0
            scatter(SAM_seas,squeeze(SAT_seas(j,i,:)),'filled')
        end
    end
end




