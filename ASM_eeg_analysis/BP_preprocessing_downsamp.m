%% pre-process the bandpower data to then perform spectral analysis via wavelet transform and FOOOF

% load the bandpower data 
load([ptID '_bandpower.mat'])


% down sample the data to fit the ASM load
down_samp_factor = 6;
delta_downsamp = nanmean(reshape([pt_data.delta(:); nan(mod(-numel(pt_data.delta),down_samp_factor),1)],down_samp_factor,[]));
theta_downsamp = nanmean(reshape([pt_data.theta(:); nan(mod(-numel(pt_data.theta),down_samp_factor),1)],down_samp_factor,[]));
alpha_downsamp = nanmean(reshape([pt_data.alpha(:); nan(mod(-numel(pt_data.alpha),down_samp_factor),1)],down_samp_factor,[]));
beta_downsamp = nanmean(reshape([pt_data.beta(:); nan(mod(-numel(pt_data.beta),down_samp_factor),1)],down_samp_factor,[]));
gamma_downsamp = nanmean(reshape([pt_data.gamma(:); nan(mod(-numel(pt_data.gamma),down_samp_factor),1)],down_samp_factor,[]));
time_downsamp = nanmean(reshape([pt_data.times(:); nan(mod(-numel(pt_data.times),down_samp_factor),1)],down_samp_factor,[]));