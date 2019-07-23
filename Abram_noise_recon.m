%% Abram recons
%
% This script uses all of the proxies employed in Abram et al. 2014 to 
% asses its significance with respect to red noise (using % red noise 
% proxy reconstructions)

clear
load marshall_SAM.mat; Marshall_SAM = flipud(Marshall_SAM);
load SAM_seasonal.mat; 
load Fogt_Jones.mat; FogtJones_SF = flipud(FogtJones_SF); FogtJones_MA = flipud(FogtJones_MA);
load Abrm_prox.mat;

for i = 1:52
    start_p(1,i) = min(find(~isnan(zAll_data_shift(:,i))));
    end_p(1,i) = find(~isnan(zAll_data_shift(:,i)),1,'last');
    End_ma(1,i) = 39;
    End_vbk(1,i) = 109;
    End_FJ(1,i) = 91;
    ma_start(1,i) = 20 + start_p(1,i);
    ma_end(1,i) = 59;
    vbk_start(1,i) = 10 + start_p(1,i);
    vbk_end(1,i) = 119;
    FJ_start(1,i) = 9 + start_p(1,i);
    FJ_end(1,i) = 100;
end
End_vbk(1,5) = 103;
vbk_end(1,5) = 113;

nsim = 10000;
noise_prox = nan(995,52,nsim);
r_noise_prox = nan(3,52,nsim);