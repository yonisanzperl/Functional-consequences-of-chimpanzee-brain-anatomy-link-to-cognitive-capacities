% HumChimp_SC_harmonize.m
% Harmonise human and chimpanzee structural connectomes (same BB50 homologous
% parcellation, same DSI Studio GQI pipeline, FA 0.05) before running the Hopf model,
% so that model-derived dynamical measures can be compared across species.
%
% Same parcellation and tractography remove the biggest differences, but the species
% can still differ in (i) density (number of non-zero edges), (ii) overall streamline
% counts / strength, and (iii) the shape of the weight distribution. Each of these
% moves the model measures (meanFC, Ignition, ...) independently of topology.
%
% Three versions are built, from least to most harmonised:
%   V0  raw                     : clean, then C/mean(C)  (= current pipeline)
%   V1  density-matched         : keep the strongest edges so every connectome has the
%                                 same density (the minimum over all subjects of both
%                                 species), then C/mean(C)
%   V2  density + weight-matched: V1, then quantile-normalise the edge weights of every
%                                 subject to one reference distribution pooled over
%                                 both species (the rank order of edges, i.e. the
%                                 topology, is kept; only the weight values are
%                                 equalised), then C/mean(C).   <-- recommended
%
% For each version it reports, per species, density, isolated nodes, strength
% heterogeneity, weight-distribution shape and (if BCT is on the path) SC graph
% measures, with the human-vs-chimp rank-sum test and the rank-biserial effect size.
% A good harmonisation removes the differences in density/strength/weights while the
% topological differences (clustering, modularity, path length) remain what they are.
%
% Before harmonising: regions that are not homologous or not reliably tracked are removed
% in both species, and subjects that fail QC (empty, density < 0.15, or any disconnected
% region, which would leave a Hopf node uncoupled and lower meanFC) are excluded.
%
% Output: SC_HumChimp_harmonized.mat (-v7) with SChum_V0/V1/V2, SCchimp_V0/V1/V2
% (N x N x subjects), ready for the G sweep + measures scripts (use exactly the same
% model parameters and G range in both species).
%
% Run from Monkey_MVH/HumChimp_compare/. Needs connectivity_HUMAN_BB50_v7.mat and
% connectivity_CHIMP_BB50_v7.mat. BCT optional (for the graph measures).
% Checked in Python on 23 Sep 2026: 74 regions, 58 humans, 22 chimpanzees, density 0.190.

clear all; close all
addpath(genpath('/Users/cbc/Documents/matlab_stuffs/'));   % BCT, as in your scripts

P.metric       = 1;          % 1 = number of streamlines (see weightDescriptions)
P.dropRegions  = {'unknown','-FH'};   % 'unknown' = unlabelled cortex (not a homologous area);
                                      % FH is disconnected in several chimps -> removed bilaterally
P.minDensity   = 0.15;       % QC: subjects below this density = failed tractography
P.density      = [];         % [] = minimum density over all subjects that pass QC
rng(2026);

%% ---- load ----
H = load('connectivity_HUMAN_BB50_v7.mat');      % v7 copies of the DSI Studio v7.3 files
C = load('connectivity_CHIMP_BB50_v7.mat');
assert(isequal(cellstr(H.regionDescriptions), cellstr(C.regionDescriptions)), 'the two species must share the parcellation');
wd = cellstr(H.weightDescriptions); fprintf('weight used: %s\n', strtrim(wd{P.metric}));
SCh = double(squeeze(H.connectivity(:,:,P.metric,:)));
SCc = double(squeeze(C.connectivity(:,:,P.metric,:)));
for s = 1:size(SCh,3), SCh(:,:,s) = clean_sc(SCh(:,:,s)); end
for s = 1:size(SCc,3), SCc(:,:,s) = clean_sc(SCc(:,:,s)); end

% ---- regions: drop non-homologous / unreliable ones in BOTH species ----
reg = strtrim(cellstr(H.regionDescriptions)); reg = reg(:);
drop = false(numel(reg),1);
for k = 1:numel(P.dropRegions), drop = drop | contains(reg, P.dropRegions{k}); end
fprintf('dropped regions: %s\n', strjoin(reg(drop), ', '));
SCh = SCh(~drop,~drop,:); SCc = SCc(~drop,~drop,:); regions = reg(~drop);
N = size(SCh,1);

% ---- subject QC: empty, failed (low density) or with a disconnected region ----
qc = @(W) nnz(triu(W,1))/(size(W,1)*(size(W,1)-1)/2) >= P.minDensity && all(sum(W,2) > 0);
okH = arrayfun(@(s) qc(SCh(:,:,s)), 1:size(SCh,3));
okC = arrayfun(@(s) qc(SCc(:,:,s)), 1:size(SCc,3));
subjH = strtrim(cellstr(H.subjects)); subjC = strtrim(cellstr(C.subjects));
fprintf('QC excluded humans: %s\n', strjoin(subjH(~okH), ', '));
fprintf('QC excluded chimps: %s\n', strjoin(subjC(~okC), ', '));
SCh = SCh(:,:,okH); SCc = SCc(:,:,okC); subjH = subjH(okH); subjC = subjC(okC);
nH = size(SCh,3); nC = size(SCc,3);
fprintf('N = %d regions; %d humans, %d chimpanzees\n', N, nH, nC);

%% ---- target density ----
dH = arrayfun(@(s) density_und(SCh(:,:,s)), 1:nH);
dC = arrayfun(@(s) density_und(SCc(:,:,s)), 1:nC);
fprintf('raw density: human %.3f [%.3f-%.3f], chimp %.3f [%.3f-%.3f]\n', ...
        median(dH), min(dH), max(dH), median(dC), min(dC), max(dC));
if isempty(P.density), P.density = min([dH dC]); end
fprintf('target density = %.3f\n', P.density);

%% ---- build V0, V1, V2 ----
V = struct();
V(1).name = 'V0 raw';                 V(1).h = SCh;  V(1).c = SCc;
V(2).name = 'V1 density-matched';     V(2).h = SCh;  V(2).c = SCc;
for s = 1:nH, V(2).h(:,:,s) = threshold_prop(SCh(:,:,s), P.density); end
for s = 1:nC, V(2).c(:,:,s) = threshold_prop(SCc(:,:,s), P.density); end
V(3).name = 'V2 density+weights';
[V(3).h, V(3).c] = quantile_norm(V(2).h, V(2).c);

% keep the un-normalised versions for the diagnostics, then C/mean(C) as in the model
for v = 1:3
    V(v).hraw = V(v).h; V(v).craw = V(v).c;
    for s = 1:nH, W = V(v).h(:,:,s); V(v).h(:,:,s) = W / mean(W(:)); end
    for s = 1:nC, W = V(v).c(:,:,s); V(v).c(:,:,s) = W / mean(W(:)); end
end

%% ---- diagnostics ----
hasBCT = exist('clustering_coef_wu','file') == 2 && exist('modularity_und','file') == 2 ...
         && exist('efficiency_wei','file') == 2;
dnames = {'density','isolated nodes','mean weight (raw)','strength CV','weight CV','weight skewness'};
if hasBCT, dnames = [dnames, {'clustering','modularity Q','global efficiency'}]; end
nD = numel(dnames);
D = cell(3,2);
for v = 1:3
    D{v,1} = sc_diag(V(v).hraw, V(v).h, hasBCT);
    D{v,2} = sc_diag(V(v).craw, V(v).c, hasBCT);
end

for v = 1:3
    fprintf('\n===== %s =====\n', V(v).name);
    fprintf('%-20s %12s %12s %9s %9s\n', 'measure', 'human (med)', 'chimp (med)', 'p', 'r_rb');
    for m = 1:nD
        xh = D{v,1}(:,m); xc = D{v,2}(:,m);
        [p, rb] = rank_test(xh, xc);
        fprintf('%-20s %12.4g %12.4g %9.2g %+9.2f\n', dnames{m}, median(xh), median(xc), p, rb);
    end
end
fprintf('\nr_rb = rank-biserial effect size (positive = human > chimp).\n');
fprintf('Goal: in V2, density / strength / weight rows show no species difference;\n');
fprintf('topology rows (clustering, modularity, efficiency) show what is left.\n');

%% ---- figure ----
cH = [0.17 0.37 0.54]; cC = [0.78 0.44 0.16];
figure('Color','w','Position',[40 40 190*nD 560]);
for v = 1:3
    for m = 1:nD
        subplot(3, nD, (v-1)*nD + m); hold on
        xh = D{v,1}(:,m); xc = D{v,2}(:,m);
        scatter(1 + 0.12*randn(nH,1), xh, 12, cH, 'filled', 'MarkerFaceAlpha', 0.5);
        scatter(2 + 0.12*randn(nC,1), xc, 12, cC, 'filled', 'MarkerFaceAlpha', 0.5);
        plot([0.7 1.3], median(xh)*[1 1], 'k-', 'LineWidth', 2);
        plot([1.7 2.3], median(xc)*[1 1], 'k-', 'LineWidth', 2);
        [p, rb] = rank_test(xh, xc);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Hum','Chimp'}, 'FontSize', 8, 'TickDir','out');
        xlim([0.4 2.6]); box off
        if p < 0.05, fw = 'bold'; else, fw = 'normal'; end
        title({dnames{m}, sprintf('r_{rb} = %+.2f, p = %.2g', rb, p)}, 'FontSize', 8, 'FontWeight', fw);
        if m == 1, ylabel(V(v).name, 'FontWeight','bold'); end
    end
end
sgtitle(sprintf('Human (n = %d) vs chimpanzee (n = %d) SC, BB50 parcellation, target density %.3f', ...
        nH, nC, P.density), 'FontSize', 10);

%% ---- save ----
SChum_V0 = V(1).h; SChum_V1 = V(2).h; SChum_V2 = V(3).h;
SCchimp_V0 = V(1).c; SCchimp_V1 = V(2).c; SCchimp_V2 = V(3).c;
save('SC_HumChimp_harmonized.mat', 'SChum_V0','SChum_V1','SChum_V2', ...
     'SCchimp_V0','SCchimp_V1','SCchimp_V2', 'regions','subjH','subjC', 'P', 'D', 'dnames', '-v7');
fprintf('\nSaved SC_HumChimp_harmonized.mat\n');

%% ===================== local functions =====================
function W = clean_sc(W)
    W(isnan(W)) = 0; W(W < 0) = 0;
    W = (W + W') / 2;
    W(1:size(W,1)+1:end) = 0;
end

function d = density_und(W)
    N = size(W,1);
    d = nnz(triu(W,1)) / (N*(N-1)/2);
end

function W = threshold_prop(W, d)
    % keep the strongest fraction d of all possible edges (ties broken by order)
    N = size(W,1);
    iu = find(triu(true(N),1));
    w = W(iu);
    k = round(d * numel(iu));
    [~, o] = sort(w, 'descend');
    keep = false(size(w)); keep(o(1:min(k, nnz(w)))) = true;
    w(~keep) = 0;
    W = zeros(N); W(iu) = w; W = W + W';
end

function [Hn, Cn] = quantile_norm(Hs, Cs)
    % every subject (both species) gets the same edge-weight distribution: the mean of
    % the sorted non-zero weights over all subjects. Needs equal edge counts (V1).
    N = size(Hs,1); iu = find(triu(true(N),1));
    all_ = cat(3, Hs, Cs); S = size(all_,3);
    k = arrayfun(@(s) nnz(all_(iu + (s-1)*N*N)), 1:S);
    k = min(k);
    ref = zeros(k,1);
    for s = 1:S
        w = all_(iu + (s-1)*N*N); w = sort(w(w > 0), 'descend');
        ref = ref + w(1:k);
    end
    ref = ref / S;
    out = zeros(size(all_));
    for s = 1:S
        w = all_(iu + (s-1)*N*N);
        [~, o] = sort(w, 'descend');
        wn = zeros(size(w)); wn(o(1:k)) = ref;
        W = zeros(N); W(iu) = wn; out(:,:,s) = W + W';
    end
    Hn = out(:,:,1:size(Hs,3)); Cn = out(:,:,size(Hs,3)+1:end);
end

function X = sc_diag(Wraw, Wnorm, hasBCT)
    S = size(Wraw,3); X = [];
    for s = 1:S
        W = Wraw(:,:,s); A = Wnorm(:,:,s); N = size(W,1);
        w = W(triu(true(N),1)); w = w(w > 0);
        st = sum(A,2);
        row = [density_und(W), sum(sum(W,2) == 0), mean(W(:)), std(st)/mean(st), ...
               std(w)/mean(w), skewness(w)];
        if hasBCT
            [~, Q] = modularity_und(A);
            row = [row, mean(clustering_coef_wu(A / max(A(:)))), Q, efficiency_wei(A)];   % efficiency_wei converts weights to lengths itself
        end
        X(s,:) = row; %#ok<AGROW>
    end
end

function [p, rb] = rank_test(xh, xc)
    [p, ~, st] = ranksum(xh, xc);
    n1 = numel(xh); n2 = numel(xc);
    U = st.ranksum - n1*(n1+1)/2;
    rb = 2*U/(n1*n2) - 1;                      % positive = human > chimp
end
