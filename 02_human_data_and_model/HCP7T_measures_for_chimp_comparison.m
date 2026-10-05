% HCP7T_measures_for_chimp_comparison.m
% Human empirical counterpart of the chimpanzee Fig 2 / Fig 3 measures.
% Computes, for every HCP 7T subject, EXACTLY the measures used for the chimps
% (same definitions as functional_measures_CHIMP_MVH.m and Fig3_CHIMP_MVH.m),
% in TWO bands:
%   band 1 = 0.02-0.03 Hz  (the band of the chimpanzee model; main comparison)
%   band 2 = 0.008-0.08 Hz (standard empirical band; robustness)
% and exports them with the g-factor and demographics in one small file:
%   HCP7T_measures_for_chimp_comparison.mat
%
% Preprocessing follows HCP_7T_empiricalmeasures.m (Schaefer 1000 -> 100, TR = 1 s,
% Tmax = 900, cut = 50). g-factor as in Behavior_Ceff_7T_CHimPaper.m (factoran, 1 factor,
% 10 NIH toolbox tests).
% Needs: hcp7t_rfMRI_REST1_PA_schaefer1000.mat, schaefer1000to100.mat,
%        hcpbehaviouraldata.mat, BCT (modularity_und, clustering_coef_wu, betweenness_wei),
%        Signal Processing + Statistics Toolboxes. Run from Monkey_MVH/HumanBehaviour_7T/.

clear all; close all
addpath(genpath('/Users/cbc/Documents/matlab_stuffs/BCT/'));   % as in Fig3_CHIMP_MVH.m
addpath '/Users/cbc/Documents/sanz/Monkey_MVH/HumanBehaviour_7T'

load('hcp7t_rfMRI_REST1_PA_schaefer1000.mat');   % subject{}.schaeferts, subject{}.id
load('schaefer1000to100.mat');
NSUB = numel(subject);  N = 1000;  NP = 100;
TR = 1;  Tmax = 900;  cut = 50;  ignLag = 2;     % as in the chimp code: corr2(Amplitude(1:end-2), kop(3:end)) (was 3 by mistake; r = 0.998 between the two, results unchanged)
for i = 1:NP, partition{i} = find(schaefer1000to100 == i); end
bands = [0.02 0.03; 0.008 0.08];
bandNames = {'0.02-0.03 Hz (chimp model band)', '0.008-0.08 Hz (standard)'};
nB = size(bands,1);

names = {'meanFC','Ignition','EdgeTurbulence','EntropyeFCD','FCmodularity', ...
         'FCclustering','FCbetweenness','Metastability','Synchrony'};
M = nan(NSUB, numel(names), nB);
ids = nan(NSUB,1);

for b = 1:nB
    [bf, af] = butter(2, bands(b,:)/(1/(2*TR)));
    for s = 1:NSUB
        ids(s) = subject{s}.id;
        ts = subject{s}.schaeferts(:,1:Tmax);
        sf = nan(N, Tmax);
        for i = 1:N
            x = detrend(ts(i,:) - nanmean(ts(i,:)));
            if ~any(isnan(x)), sf(i,:) = filtfilt(bf, af, x); end
        end
        sf = sf(:, cut:end-cut);
        ts2 = zeros(NP, size(sf,2));
        for i = 1:NP, ts2(i,:) = nanmean(sf(partition{i},:), 1); end

        % --- same definitions as the chimp measures ---
        Amp = zeros(size(ts2)); Ph = zeros(size(ts2));
        for i = 1:NP
            h = hilbert(ts2(i,:) - mean(ts2(i,:)));
            Amp(i,:) = abs(h); Ph(i,:) = angle(h);
        end
        kop = abs(nansum(complex(cos(Ph), sin(Ph)))) / NP;
        FC = corrcoef(ts2');
        ig = zeros(NP,1);
        for i = 1:NP, ig(i) = corr2(Amp(i,1:end-ignLag), kop(1+ignLag:end)); end
        zt = zscore(ts2');                                   % T x NP
        edgeTurbNode = std(zt .* mean(zt,2), [], 1);         % = std over time of GBCEdge
        Isub = find(tril(ones(NP),-1));
        E = zeros(numel(Isub), size(zt,1));
        for t = 1:size(zt,1), tmp = zt(t,:)'*zt(t,:); E(:,t) = tmp(Isub); end
        FCD = (E'*E) ./ (vecnorm(E)'*vecnorm(E));
        IsubT = find(tril(ones(size(FCD)),-1));
        [~, Q] = modularity_und(FC);
        M(s,:,b) = [ nanmean(FC(:)), mean(ig), mean(edgeTurbNode), ...
                     0.5*log(2*pi*var(FCD(IsubT))) + 0.5, Q, ...
                     mean(abs(clustering_coef_wu(FC))), mean(betweenness_wei(FC)), ...
                     std(kop), mean(kop) ];
        if mod(s,20) == 0, fprintf('band %d: subject %d/%d\n', b, s, NSUB); end
    end
end

%% ---- g-factor and demographics (as in Behavior_Ceff_7T_CHimPaper.m) ----
load('hcpbehaviouraldata.mat');  behav = hcpbehaviouraldata;
try, fprintf('Behaviour columns 1-6: %s\n', strjoin(behav.Properties.VariableNames(1:6), ', ')); catch, end
behavsub = nan(size(behav,1),1);
for jj = 1:size(behav,1), behavsub(jj) = behav{jj,1}; end
cols = [122 129 118 158 120 145 156 116 125 127];   % the 10 tests used for the g-factor
tests = nan(NSUB, numel(cols));  gender = cell(NSUB,1);  age = cell(NSUB,1);
for s = 1:NSUB
    idx = find(behavsub == ids(s));
    for c = 1:numel(cols), tests(s,c) = behav{idx, cols(c)}; end
    gender{s} = behav{idx,4};  age{s} = behav{idx,5};    % unrestricted HCP: col 4 Gender, col 5 Age (bins)
end
ok = all(~isnan(tests), 2);
g = nan(NSUB,1);
[~,~,~,~,F] = factoran(tests(ok,:), 1);
g(ok) = F;
if corr(g(ok), tests(ok,1)) < 0, g = -g; end          % orient: higher g = better (PMAT24)

save('HCP7T_measures_for_chimp_comparison.mat', 'ids','names','bands','bandNames','M', ...
     'g','tests','gender','age','ok');
fprintf('Saved HCP7T_measures_for_chimp_comparison.mat (%d subjects, %d with g-factor)\n', NSUB, sum(ok));
