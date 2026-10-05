% Fig2_surrogate_connectomes.m
% Surrogate connectomes for Fig 2 (panel i): the same model and the same measures as the
% real connectomes (functional_measures_CHIMP_MVH.m, Sep 2026 rerun), run on randomised
% versions of each animal's structural connectome, at the G* of that animal's REAL connectome.
%
%   P.family = 'Degree'   : randmio_und(C, P.degreeIter)      keeps each region's degree
%   P.family = 'Strength' : null_model_und_sign(C)            keeps degree, strength and
%                                                             weight distribution
% One surrogate SET = one randomised connectome per animal, simulated with P.nTrials trials,
% saved as one file with the same name pattern and RH fields as before, so that
% Fig2_REAL0926_full.m and Fig2c_surrogate_rank_test.m read it unchanged:
%   <P.outDir>/functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_Degree_001.mat
%   <P.outDir>/functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_Strength__001.mat
% RH fields: Integration, Ignition, IgStd, EdgeMetaNode, EdgesMeta, EntroEdgesMeta,
% Metasim_sub2, Synchro_sub2, Qs, FCsim (59 animals; NaN for animals not simulated).
%
% Identical to functional_measures_CHIMP_MVH.m: C = C/mean(C(:)) BEFORE randomising,
% a = -0.02, noise 0.04, TR = 0.72 s, dt = 0.1*TR/2, 1000 s transient, Tmax = 1200 samples,
% band 0.02-0.03 Hz (2nd-order Butterworth, filtfilt), cut = 10 samples, heterogeneous node
% frequencies f = 0.025*(1 + 0.1*rand) drawn once per run (here once per surrogate set,
% with a fixed seed), G* = max of the smoothed metastability curve of the real connectome
% (meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat, Grange = (0:0.1:10)/1000).
% Same measure definitions (Ignition lag 2, edge turbulence = std over time of the node
% edge time series, entropy of eFCD, metastability/synchrony = std/mean of the Kuramoto
% order parameter, integration = nanmean of the trial-averaged FC incl. diagonal).
% Not computed (not used in Fig 2): LZ complexity, Granger, event-based ignition.
% Only speed-ups: trials integrated in batches (same equations), animals in parallel.
%
% Speed: one trial ~1-2 s (integration + eFCD). One set = 43 animals x 200 trials
% ~ 3-4 h serial, ~35-45 min with 6 workers. 20 sets per family ~ 12-15 h with 6 workers.
% Per-animal files are written to P.outDir/tmp_<family>_<set>/, so an interrupted run resumes.
% Needs BCT (randmio_und, modularity_und) and null_model_und_sign.m / randmio_und_signed.m
% (in this folder), Signal Processing Toolbox, Parallel Computing Toolbox (P.nWorkers > 0).
% Run from Fig2_measures_vs_cognition. Then in Fig2_REAL0926_full.m set
%   P.nullDir = 'surrogates_sweepfreq';  P.nullNoise = 0.04;

clear all; close all
addpath(genpath('/Users/cbc/Documents/matlab_stuffs/BCT/'));

%% =================== PARAMETERS (edit only here) ===================
P.family     = 'Degree';          % 'Degree' or 'Strength'
P.sets       = 21:50;              % surrogate sets to run (e.g. 1:20, or 21:100 later)
P.nTrials    = 200;               % trials per animal, as for the real connectomes
P.batch      = 50;                % trials integrated together (memory: N x 2 x batch)
P.nWorkers   = 6;                 % 0 = serial
P.degreeIter = 1000;              % randmio_und iterations (as in the original comment)
P.onlyAnalysed = true;            % true: the 43 animals of Fig 2; false: all with SC
P.outDir     = 'surrogates_sweepfreq';
P.seed       = 2026;
% ---- model (identical to functional_measures_CHIMP_MVH.m) ----
P.N = 114; P.a = -0.02; P.sig = 0.04; P.TR = 0.72; P.dt = 0.1*P.TR/2;
P.Ttrans = 1000; P.Tmax = 1200; P.cut = 10; P.flp = 0.02; P.fhi = 0.03; P.f0 = 0.025;
P.Grange = (0:0.1:10)/1000;
%% ===================================================================

[P.bfilt, P.afilt] = butter(2, [P.flp P.fhi]/(1/(2*P.TR)));
load('meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat', 'Metasim_sub2');
load('CHIMP_CWAS_DATA_FORDECO.mat', 'CJ_CHIMP', 'AapCog_normalized');
nA = size(CJ_CHIMP,3);
Ms = zeros(nA, size(Metasim_sub2,2));
for s = 1:nA, Ms(s,:) = smooth(squeeze(mean(Metasim_sub2(s,:,:),3))); end
[Meta, Gidx] = max(Ms, [], 2);
hasSC = squeeze(~all(all(isnan(CJ_CHIMP),1),2));
if P.onlyAnalysed
    animals = find(Meta > 0 & hasSC); animals = setdiff(animals, [23 27]);
    animals = animals(~isnan(AapCog_normalized(animals)));
else
    animals = find(hasSC);
end
animals = animals(:)';
fprintf('%s surrogates: %d animals, sets %d-%d, %d trials each\n', P.family, numel(animals), P.sets(1), P.sets(end), P.nTrials);
if ~exist(P.outDir, 'dir'), mkdir(P.outDir); end
if P.nWorkers > 0 && isempty(gcp('nocreate')), parpool(P.nWorkers); end
switch P.family
    case 'Degree',   tag = 'Degree_';
    case 'Strength', tag = 'Strength__';
    otherwise, error('P.family must be ''Degree'' or ''Strength''');
end

for set = P.sets
    fOut = fullfile(P.outDir, sprintf('functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_%s%03d.mat', tag, set));
    if exist(fOut, 'file'), fprintf('set %d exists, skipped\n', set); continue; end
    tSet = tic;
    rng(P.seed + set);                                        % node frequencies of this set
    Pset = P; Pset.w = 2*pi*P.f0*(rand(P.N,1)*0.1 + 1);
    tmpDir = fullfile(P.outDir, sprintf('tmp_%s_%03d', P.family, set));
    if ~exist(tmpDir, 'dir'), mkdir(tmpDir); end
    Call = CJ_CHIMP(:,:,animals); Gall = P.Grange(Gidx(animals));
    parfor (k = 1:numel(animals), P.nWorkers)
        f = fullfile(tmpDir, sprintf('animal%02d.mat', animals(k)));
        if exist(f, 'file'), continue; end
        rng(P.seed + 1000*set + animals(k), 'twister');       % surrogate + noise stream
        C = Call(:,:,k); C(isnan(C)) = 0;
        C = C / mean(C(:));                                    % as in the real run, before randomising
        if strcmp(Pset.family, 'Degree'), C = randmio_und(C, Pset.degreeIter);
        else, C = null_model_und_sign(C); end
        out = run_animal(C, Gall(k), Pset);
        parsave(f, out);
    end
    % ---- assemble in the RH format ----
    RH = struct();
    flds = {'Integration','Ignition','IgStd','EdgesMeta','EntroEdgesMeta','Metasim_sub2','Synchro_sub2','Qs'};
    for q = 1:numel(flds), RH.(flds{q}) = nan(1, nA); end
    RH.EdgeMetaNode = nan(nA, P.N); RH.FCsim = nan(nA, P.N, P.N);
    for k = 1:numel(animals)
        o = load(fullfile(tmpDir, sprintf('animal%02d.mat', animals(k)))); o = o.out;
        for q = 1:numel(flds), RH.(flds{q})(animals(k)) = o.(flds{q}); end
        RH.EdgeMetaNode(animals(k),:) = o.EdgeMetaNode; RH.FCsim(animals(k),:,:) = o.FCsim;
    end
    RH.sigma = P.sig; RH.family = P.family; RH.set = set; RH.Gstar_from = 'real connectome'; RH.P = rmfield(Pset, {'bfilt','afilt'});
    save(fOut, 'RH');
    fprintf('set %d saved (%.1f min)\n', set, toc(tSet)/60);
end

%% =================== local functions ===================
function out = run_animal(C, G, P)
    N = size(C,1); nT = P.nTrials;
    Isub = find(tril(ones(N),-1));
    FCacc = zeros(N); v = zeros(nT, 6); etn = zeros(nT, N); igsd = zeros(nT,1);
    for c0 = 1:P.batch:nT
        c1 = min(nT, c0 + P.batch - 1);
        X = sim_batch(C, G, c1 - c0 + 1, P);
        for t = c0:c1
            m = trial_measures(X(:,:,t-c0+1), P, Isub);
            FCacc = FCacc + m.FC;
            v(t,:) = [m.meta, m.sync, m.ign, m.edgesMeta, m.entropy, 0];
            etn(t,:) = m.edgeTurbNode; igsd(t) = m.igStd;
        end
    end
    FCm = FCacc / nT;
    [~, Q] = modularity_und(FCm);
    out.Integration    = nanmean(nanmean(FCm));
    out.Qs             = Q;
    out.Metasim_sub2   = mean(v(:,1));
    out.Synchro_sub2   = mean(v(:,2));
    out.Ignition       = mean(v(:,3));
    out.EdgesMeta      = mean(v(:,4));
    out.EntroEdgesMeta = mean(v(:,5));
    out.IgStd          = std(igsd);
    out.EdgeMetaNode   = mean(etn, 1);
    out.FCsim          = FCm;
end

function X = sim_batch(C, G, B, P)
    % Hopf network, B trials at once; same equations, dt, transient and sampling rule as
    % functional_measures_CHIMP_MVH.m (z = [x y], zz = [y x], omega = [-w w]).
    N = size(C,1); wC = G*C; deg = sum(wC,2); w = P.w; a = P.a; dt = P.dt; ds = sqrt(dt)*P.sig;
    x = 0.1*ones(N,B); y = 0.1*ones(N,B);
    for t = 0:dt:P.Ttrans
        cx = wC*x - deg.*x; cy = wC*y - deg.*y; r2 = x.*x + y.*y;
        xn = x + dt*(a*x - w.*y - x.*r2 + cx) + ds*randn(N,B);
        y  = y + dt*(a*y + w.*x - y.*r2 + cy) + ds*randn(N,B);
        x  = xn;
    end
    X = zeros(N, P.Tmax, B); nn = 0;
    for t = 0:dt:((P.Tmax-1)*P.TR)
        cx = wC*x - deg.*x; cy = wC*y - deg.*y; r2 = x.*x + y.*y;
        xn = x + dt*(a*x - w.*y - x.*r2 + cx) + ds*randn(N,B);
        y  = y + dt*(a*y + w.*x - y.*r2 + cy) + ds*randn(N,B);
        x  = xn;
        if abs(mod(t, P.TR)) < 0.01 && nn < P.Tmax, nn = nn + 1; X(:,nn,:) = reshape(x, N, 1, B); end
    end
end

function m = trial_measures(x, P, Isub)
    % x: N x Tmax; same definitions as functional_measures_CHIMP_MVH.m
    N = size(x,1);
    sf = zeros(size(x));
    for i = 1:N, sf(i,:) = filtfilt(P.bfilt, P.afilt, detrend(x(i,:) - nanmean(x(i,:)))); end
    ts2 = sf(:, P.cut:end-P.cut);
    Amp = zeros(size(ts2)); Ph = zeros(size(ts2));
    for i = 1:N
        h = hilbert(ts2(i,:) - mean(ts2(i,:))); Amp(i,:) = abs(h); Ph(i,:) = angle(h);
    end
    m.FC = corrcoef(ts2');
    kop = abs(nansum(complex(cos(Ph), sin(Ph)))) / N;
    m.meta = nanstd(kop); m.sync = nanmean(kop);
    ig = zeros(N,1);
    for i = 1:N, ig(i) = corr2(Amp(i,1:end-2), kop(3:end)); end
    m.ign = mean(ig); m.igStd = std(ig);
    z = zscore(ts2');                                   % T x N
    T = size(z,1); E = zeros(numel(Isub), T); GBC = zeros(N, T);
    for t = 1:T
        tmp = z(t,:)' * z(t,:); E(:,t) = tmp(Isub); GBC(:,t) = mean(tmp)';
    end
    m.edgesMeta = std(E(:));
    m.edgeTurbNode = std(GBC, [], 2)';
    FCD = (E'*E) ./ (vecnorm(E)'*vecnorm(E));
    IsubT = find(tril(ones(T),-1));
    m.entropy = 0.5*log(2*pi*var(FCD(IsubT))) + 0.5;
end

function parsave(f, out)
    save(f, 'out');
end
