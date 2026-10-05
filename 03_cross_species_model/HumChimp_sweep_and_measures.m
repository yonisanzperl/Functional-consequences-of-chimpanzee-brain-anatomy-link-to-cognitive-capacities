% HumChimp_sweep_and_measures.m
% Human vs chimpanzee comparison with the Fig 2 model, run IDENTICALLY for both species
% on the harmonised connectomes (V1 = density-matched, V2 = density + weight-matched).
%
% Model and measures follow the Fig 2 pipeline
% (NOSUSEFUL_Fig2_2026.09.23/CHIMP_CWAS_metastability_empirica.m and
%  functional_measures_CHIMP_MVH.m):
%   Hopf network, C = C/mean(C(:)), a = -0.02, noise sig = 0.04, TR = 0.72 s,
%   dt = 0.1*TR/2, transient 1000 s, 1200 samples, band-pass 0.02-0.03 Hz (2nd-order
%   Butterworth, filtfilt), 10 samples cut at each end, node frequency 0.025 Hz.
%   G* = max of the smoothed (trial-averaged) metastability curve.
%   Measures at G*: trial-averaged FC -> meanFC (nanmean incl. diagonal, as before),
%   FC modularity (modularity_und), hemispheric FC, FC-SC correlation; per-trial ->
%   metastability, synchrony, Ignition (lag 2), edge turbulence, entropy eFCD.
%
% Only differences from the Fig 2 scripts, to make it exactly the same for both species:
%   * the SAME parameters are used for the G search and for the measures (Fig 2 used
%     heterogeneous node frequencies in the search and homogeneous ones in the measures;
%     here P.omega applies to both);
%   * node frequencies are drawn once (fixed seed) and shared by all subjects of both
%     species, so homologous regions get the same frequency;
%   * one shared G range for both species, searched coarse-to-fine (the model is not
%     changed, only which G values are simulated);
%   * the measures use new trials, independent of those that chose G*.
%
% Input : SC_HumChimp_harmonized.mat (from HumChimp_SC_harmonize.m)
% Output: HumChimp_<V>_sweep_measures.mat for each version, with, per species (H / Ch):
%         Gstar, atEdge, Meta curve (Sweep), and R* structs with the measures.
%         The last section compares the species (rank-sum, rank-biserial r) for V1 and V2.
% Per-subject files are written to a temporary folder, so an interrupted run resumes.
% Needs: Signal Processing Toolbox, modularity_und (BCT), Parallel Computing Toolbox
% for P.nWorkers > 0. Run from Monkey_MVH/HumChimp_compare/.

clear all; close all
addpath(genpath('/Users/cbc/Documents/matlab_stuffs/'));   % BCT
addpath(fullfile('..','DataGustavo'));                       % modularity_und copy

%% =================== PARAMETERS (edit only here) ===================
%P.versions  = {'V2','V1'};        % V2 first (main analysis), then V1 (sensitivity)
P.versions  = {'V0'};        % V2 first (main analysis), then V1 (sensitivity)
P.a         = -0.02;
P.sig       = 0.04;
P.TR        = 0.72;
P.dt        = 0.1*P.TR/2;
P.flp       = 0.02;  P.fhi = 0.03;
P.f0        = 0.025;
P.omega     = 'heterogeneous';    % 'heterogeneous': f0*(1 + 0.1*rand), as in the Fig 2 G search
                                  % 'homogeneous'  : f0 for every node, as in the Fig 2 measures
P.Tmax      = 1200;
P.Ttrans    = 1000;               % s
P.cut       = 10;
P.ignLag    = 2;
P.Grange    = (0:0.1:40)/1000;    % shared grid, step 0.0001 as in Fig 2, up to 0.04
P.coarseStep= 20;                 % coarse grid step 0.002
P.fineHalf  = 20;                 % fine window: coarse peak +/- 0.002
P.maxShift  = 3;
P.nSel      = 40;                 % trials per G for the G search (Fig 2: 20 x 2)
P.nMeas     = 200;                % trials for the measures at G* (Fig 2: 200)
P.batch     = 100;                % systems integrated together
P.seed      = 2026;
P.nWorkers  = 6;                  % 0 = serial
%% ===================================================================

S = load('SC_HumChimp_harmonized.mat');
N = size(S.SChum_V2,1);
[P.bfilt, P.afilt] = butter(2, [P.flp P.fhi]/(1/(2*P.TR)));
rng(P.seed);
switch P.omega
    case 'heterogeneous', P.fdiff = P.f0*(rand(1,N)*0.1 + 1);
    case 'homogeneous',   P.fdiff = P.f0*ones(1,N);
end
P.w = 2*pi*P.fdiff(:);
lh = startsWith(S.regions, 'ctx-lh'); rh = startsWith(S.regions, 'ctx-rh');

for iv = 1:numel(P.versions)
    V = P.versions{iv};
    SCh = S.(['SChum_' V]);  SCc = S.(['SCchimp_' V]);
    nH = size(SCh,3); nC = size(SCc,3);
    jobs = [ones(nH,1) (1:nH)'; 2*ones(nC,1) (1:nC)'];      % [species subject]
    tmpDir = sprintf('HumChimp_%s_per_subject', V);
    if ~exist(tmpDir,'dir'), mkdir(tmpDir); end
    fprintf('\n===== %s: %d humans + %d chimps, N = %d =====\n', V, nH, nC, N);

    tAll = tic;
    parfor (j = 1:size(jobs,1), P.nWorkers)
        sp = jobs(j,1); s = jobs(j,2);
        f = fullfile(tmpDir, sprintf('sp%d_sub%02d.mat', sp, s));
        if exist(f,'file'), continue; end
        t0 = tic;
        if sp == 1, C = SCh(:,:,s); else, C = SCc(:,:,s); end
        C = C / mean(C(:));
        rng(P.seed + 1000*sp + s);
        out = run_subject(C, P, lh, rh);
        out.seconds = toc(t0);
        parsave(f, out);
        fprintf('%s %s %2d: G* = %.4f%s  meanFC %.3f  Ign %.3f  (%.0f min)\n', V, ...
            char('H'*(sp==1) + 'C'*(sp==2)), s, out.Gstar, repmat(' (AT EDGE)',1,out.atEdge), ...
            out.R.Integration, out.R.Ignition, out.seconds/60);
    end
    fprintf('%s done in %.1f h\n', V, toc(tAll)/3600);

    % ---- assemble ----
    for sp = 1:2
        if sp == 1, n = nH; tag = 'H'; else, n = nC; tag = 'Ch'; end
        R = struct(); Sweep = nan(n, numel(P.Grange)); Gstar = nan(n,1); atEdge = false(n,1);
        for s = 1:n
            o = load(fullfile(tmpDir, sprintf('sp%d_sub%02d.mat', sp, s))); o = o.out;
            Sweep(s,:) = o.metaCurve; Gstar(s) = o.Gstar; atEdge(s) = o.atEdge;
            fn = fieldnames(o.R);
            for k = 1:numel(fn)
                v = o.R.(fn{k});
                if isscalar(v), R.(fn{k})(s,1) = v; else, R.(fn{k})(s,:) = v(:)'; end
            end
        end
        OUT.(['R' tag]) = R; OUT.(['Sweep' tag]) = Sweep;
        OUT.(['Gstar' tag]) = Gstar; OUT.(['atEdge' tag]) = atEdge;
    end
    OUT.P = P; OUT.Grange = P.Grange; OUT.regions = S.regions;
    OUT.subjH = S.subjH; OUT.subjC = S.subjC; OUT.version = V;
    save(sprintf('HumChimp_%s_sweep_measures.mat', V), '-struct', 'OUT', '-v7');
    fprintf('Saved HumChimp_%s_sweep_measures.mat (G* at edge: %d humans, %d chimps)\n', ...
        V, sum(OUT.atEdgeH), sum(OUT.atEdgeCh));
    clear OUT
end

%% =================== species comparison (V1 and V2) ===================
mNames = {'Integration','Ignition','EdgeTurb','EntroEdgesMeta','Qs','Metasim','Syncsim', ...
          'FCintehe','FCSCcorr','Gstar'};
mLabels = {'Mean FC','Ignition','Edge turbulence','Entropy eFCD','FC modularity', ...
           'Metastability','Synchrony','Interhemispheric FC','FC-SC correlation','G*'};
cH = [0.17 0.37 0.54]; cC = [0.78 0.44 0.16];
figure('Color','w','Position',[40 40 1500 560]);
for iv = 1:numel(P.versions)
    V = P.versions{iv};
    f = sprintf('HumChimp_%s_sweep_measures.mat', V);
    if ~exist(f,'file'), continue; end
    D = load(f);
    fprintf('\n===== %s: human (n = %d) vs chimp (n = %d) =====\n', V, numel(D.GstarH), numel(D.GstarCh));
    fprintf('%-20s %11s %11s %9s %7s\n', 'measure', 'human', 'chimp', 'p', 'r_rb');
    for m = 1:numel(mNames)
        if strcmp(mNames{m},'Gstar'), xh = D.GstarH; xc = D.GstarCh;
        else, xh = D.RH.(mNames{m}); xc = D.RCh.(mNames{m}); end
        [p, ~, st] = ranksum(xh, xc);
        rb = 2*(st.ranksum - numel(xh)*(numel(xh)+1)/2)/(numel(xh)*numel(xc)) - 1;
        fprintf('%-20s %11.4g %11.4g %9.2g %+7.2f\n', mLabels{m}, median(xh), median(xc), p, rb);
        subplot(numel(P.versions), numel(mNames), (iv-1)*numel(mNames) + m); hold on
        scatter(1 + 0.1*randn(numel(xh),1), xh, 10, cH, 'filled', 'MarkerFaceAlpha', 0.5);
        scatter(2 + 0.1*randn(numel(xc),1), xc, 10, cC, 'filled', 'MarkerFaceAlpha', 0.5);
        plot([0.7 1.3], median(xh)*[1 1], 'k-', 'LineWidth', 2);
        plot([1.7 2.3], median(xc)*[1 1], 'k-', 'LineWidth', 2);
        set(gca,'XTick',[1 2],'XTickLabel',{'Hum','Chimp'},'FontSize',8,'TickDir','out'); xlim([0.4 2.6]); box off
        if p < 0.05, fw = 'bold'; else, fw = 'normal'; end
        title({mLabels{m}, sprintf('r_{rb} = %+.2f, p = %.2g', rb, p)}, 'FontSize', 8, 'FontWeight', fw);
        if m == 1, ylabel(V, 'FontWeight', 'bold'); end
    end
end
fprintf('\nr_rb > 0: human higher. Report differences that hold in BOTH V1 and V2.\n');

%% =================== local functions ===================
function out = run_subject(C, P, lh, rh)
    nG = numel(P.Grange);
    meta = nan(nG, 1);
    % G search, coarse
    cIdx = 1:P.coarseStep:nG;
    meta = sweep_points(C, P, cIdx, meta);
    curve = smooth(meta(cIdx)); [~, k] = max(curve); center = cIdx(k);
    % G search, full resolution around the coarse peak
    for it = 1:P.maxShift+1
        win = max(1, center-P.fineHalf) : min(nG, center+P.fineHalf);
        meta = sweep_points(C, P, win(isnan(meta(win))), meta);
        curve = smooth(meta(win)); [~, k] = max(curve); gidx = win(k);
        if (gidx == win(1) && win(1) > 1) || (gidx == win(end) && win(end) < nG), center = gidx; else, break; end
    end
    out.metaCurve = meta'; out.Gstar = P.Grange(gidx);
    out.atEdge = gidx <= 2 || gidx >= nG - 1;

    % measures at G*, new trials
    N = size(C,1); FCacc = zeros(N); v = zeros(P.nMeas, 6); etn = zeros(P.nMeas, N);
    for c0 = 1:P.batch:P.nMeas
        c1 = min(P.nMeas, c0 + P.batch - 1);
        X = sim_batch(C, repmat(out.Gstar, 1, c1-c0+1), P);
        for t = c0:c1
            m = measures(X(:,:,t-c0+1), P, true);
            FCacc = FCacc + m.FC;
            v(t,:) = [m.meta, m.sync, m.ignition, mean(m.edgeTurbNode), m.entropy, 0];
            etn(t,:) = m.edgeTurbNode;
        end
    end
    FCm = FCacc / P.nMeas;
    [~, Q] = modularity_und(FCm);
    R.Integration    = nanmean(FCm(:));
    R.Qs             = Q;
    R.Metasim        = mean(v(:,1));
    R.Syncsim        = mean(v(:,2));
    R.Ignition       = mean(v(:,3));
    R.EdgeTurb       = mean(v(:,4));
    R.EntroEdgesMeta = mean(v(:,5));
    R.EdgeMetaNode   = mean(etn, 1);
    R.FClh           = mean(mean(FCm(lh,lh)));
    R.FCrh           = mean(mean(FCm(rh,rh)));
    R.FCintehe       = mean(mean(FCm(lh,rh)));
    R.FCSCcorr       = corr(FCm(:), C(:));
    R.FCsim          = FCm(:)';
    out.R = R;
end

function meta = sweep_points(C, P, idx, meta)
    if isempty(idx), return; end
    gOf = repelem(idx(:)', P.nSel);
    val = nan(numel(gOf), 1);
    for c0 = 1:P.batch:numel(gOf)
        c1 = min(numel(gOf), c0 + P.batch - 1);
        X = sim_batch(C, P.Grange(gOf(c0:c1)), P);
        for b = c0:c1
            m = measures(X(:,:,b-c0+1), P, false);
            val(b) = m.meta;
        end
    end
    for g = idx, meta(g) = mean(val(gOf == g)); end
end

function X = sim_batch(C, Gvec, P)
    % Hopf network, several systems at once; same equations, dt, transient and sampling
    % as the Fig 2 scripts (z = [x y], zz = [y x], omega = [-w w]).
    N = size(C,1); B = numel(Gvec); G = Gvec(:)';
    deg = sum(C,2); w = P.w; a = P.a; dt = P.dt; ds = sqrt(dt)*P.sig;
    x = 0.1*ones(N,B); y = 0.1*ones(N,B);
    nTrans = numel(0:dt:P.Ttrans);
    tRun = 0:dt:((P.Tmax-1)*P.TR);
    keep = abs(mod(tRun, P.TR)) < 0.01;                 % same sampling rule as Fig 2
    for k = 1:nTrans
        cx = C*x; cy = C*y; r2 = x.*x + y.*y;
        xn = x + dt*(a*x - w.*y - x.*r2 + G.*(cx - deg.*x)) + ds*randn(N,B);
        y  = y + dt*(a*y + w.*x - y.*r2 + G.*(cy - deg.*y)) + ds*randn(N,B);
        x  = xn;
    end
    X = zeros(N, P.Tmax, B); nn = 0;
    for k = 1:numel(tRun)
        cx = C*x; cy = C*y; r2 = x.*x + y.*y;
        xn = x + dt*(a*x - w.*y - x.*r2 + G.*(cx - deg.*x)) + ds*randn(N,B);
        y  = y + dt*(a*y + w.*x - y.*r2 + G.*(cy - deg.*y)) + ds*randn(N,B);
        x  = xn;
        if keep(k) && nn < P.Tmax, nn = nn + 1; X(:,nn,:) = reshape(x, N, 1, B); end
    end
end

function m = measures(x, P, full)
    % x: N x Tmax. Same definitions as functional_measures_CHIMP_MVH.m
    N = size(x,1);
    Xd = detrend(x' - mean(x'));
    sf = filtfilt(P.bfilt, P.afilt, Xd);
    ts2 = sf(P.cut:end-P.cut, :);                        % T' x N
    H = hilbert(ts2 - mean(ts2));
    Amp = abs(H); Ph = angle(H);
    kop = abs(sum(exp(1i*Ph), 2)) / N;
    m.meta = std(kop); m.sync = mean(kop);
    m.FC = corrcoef(ts2);
    L = P.ignLag;
    A = Amp(1:end-L,:); A = A - mean(A); kk = kop(1+L:end); kk = kk - mean(kk);
    m.ignition = mean((A'*kk) ./ sqrt(sum(A.^2)' * sum(kk.^2)));
    zt = zscore(ts2);
    m.edgeTurbNode = std(zt .* mean(zt,2), [], 1);
    m.entropy = NaN;
    if full
        Isub = find(tril(ones(N),-1));
        E = zeros(numel(Isub), size(zt,1));
        for t = 1:size(zt,1), tmp = zt(t,:)'*zt(t,:); E(:,t) = tmp(Isub); end
        FCD = (E'*E) ./ (vecnorm(E)'*vecnorm(E));
        IsubT = find(tril(ones(size(FCD)),-1));
        m.entropy = 0.5*log(2*pi*var(FCD(IsubT))) + 0.5;
    end
end

function parsave(f, out)
    save(f, 'out');
end
