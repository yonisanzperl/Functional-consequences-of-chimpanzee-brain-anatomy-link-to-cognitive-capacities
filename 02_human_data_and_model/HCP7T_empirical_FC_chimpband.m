% HCP7T_empirical_FC_chimpband.m
% Empirical target for the human working-point sweep run with the chimp settings
% (model_Humans_HCP7T_asChimp.m): per-subject FC and metastability of the HCP 7T
% resting state filtered in the chimp model band (0.02-0.03 Hz).
% Preprocessing identical to HCP_7T_empiricalmeasures.m (Schaefer 1000 -> 100, TR = 1 s,
% Tmax = 900, cut = 50, 2nd-order Butterworth, filtfilt, Hilbert on the 100 parcels);
% only the band changes.
% Output: HCP7T_empirical_FC_chimpband.mat
%   FCemp        (182 x 100 x 100) subject FCs in 0.02-0.03 Hz
%   FCemp_global (100 x 100)       group-mean FC (fit target of the sweep)
%   Metaemp      (182 x 1)         metastability (std over time of the Kuramoto order parameter)
%   flp, fhi
% Needs hcp7t_rfMRI_REST1_PA_schaefer1000.mat (1.1 GB, in Monkey_MVH/HumanBehaviour_7T) and
% schaefer1000to100.mat. Takes a few minutes.

clear all; close all
addpath('/Users/cbc/Documents/sanz/Monkey_MVH/HumanBehaviour_7T')
load('/Users/cbc/Documents/sanz/Monkey_MVH/HumanBehaviour_7T/hcp7t_rfMRI_REST1_PA_schaefer1000.mat');
load('schaefer1000to100.mat');

NSUB = numel(subject); N = 1000; NP = 100;
TR = 1; Tmax = 900; cut = 50;
flp = 0.02; fhi = 0.03;                                   % chimp model band
[bfilt, afilt] = butter(2, [flp fhi]/(1/(2*TR)));
for i = 1:NP, partition{i} = find(schaefer1000to100 == i); end

FCemp = nan(NSUB, NP, NP); Metaemp = nan(NSUB,1);
for nsub = 1:NSUB
    ts = subject{nsub}.schaeferts(:,1:Tmax);
    sf = nan(N, Tmax);
    for seed = 1:N
        x = detrend(ts(seed,:) - nanmean(ts(seed,:)));
        if ~any(isnan(x)), sf(seed,:) = filtfilt(bfilt, afilt, x); end
    end
    sf = sf(:, cut:end-cut);
    ts2 = zeros(NP, size(sf,2));
    for i = 1:NP, ts2(i,:) = nanmean(sf(partition{i},:), 1); end
    Ph = zeros(size(ts2));
    for i = 1:NP, Ph(i,:) = angle(hilbert(ts2(i,:) - mean(ts2(i,:)))); end
    kop = abs(nansum(complex(cos(Ph), sin(Ph)))) / NP;
    FCemp(nsub,:,:) = corrcoef(ts2');
    Metaemp(nsub) = nanstd(kop);
    if mod(nsub,20) == 0, fprintf('subject %d/%d\n', nsub, NSUB); end
end
FCemp_global = squeeze(nanmean(FCemp, 1));
save('HCP7T_empirical_FC_chimpband.mat', 'FCemp', 'FCemp_global', 'Metaemp', 'flp', 'fhi');
fprintf('Saved HCP7T_empirical_FC_chimpband.mat: metastability %.3f +/- %.3f (n = %d)\n', ...
    nanmean(Metaemp), nanstd(Metaemp), sum(~isnan(Metaemp)));
