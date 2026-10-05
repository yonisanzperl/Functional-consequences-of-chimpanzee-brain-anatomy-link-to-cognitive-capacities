% Fig4_heldout_explanation.m
% Fig 4: which measures explain cognition in held-out animals?
% Each predictor is used ALONE in a linear model; every animal's score is predicted from
% a model fitted to the other n-1 animals (exact leave-one-out via the hat matrix).
%   held-out r   correlation between observed and held-out predicted scores
%   held-out R2  1 - PRESS / total sum of squares (0 = no better than the mean score;
%                negative = worse than the mean)
%   p            permutation test (P.nPerm shuffles of the scores), one-sided on held-out r
%   q            Benjamini-Hochberg FDR within each family (structural, dynamical); across all printed too
% Predictors: 13 SC measures (density, strength, 11 network metrics from Fig 3) and the
% model DYNAMICAL measures (REAL_0926): mean FC, Ignition, edge turbulence, entropy eFCD, plus
% synchrony and metastability automatically when the file has RH.Synchro_sub2 (Sep 2026 rerun).
% FC modularity moved to Fig 3 (FC topology), 25 Sep 2026; previous version kept as
% Fig4_heldout_explanation_withModularity_OLD.m.
% Secondary: does Ignition add to a 5-measure SC model (and SC to Ignition)?
%
% Needs cmocean.m (in this folder), results_Fig3_network_metrics.mat (Fig3_network_metrics_SC_FC.m) and the REAL_0926
% measures file. Run from its own figure folder (all inputs are copied there).
% NOTE for the text: Ignition was the strongest measure in Fig 2 on the same animals, so
% this is held-out explanation, not an independent prediction.

clear all; close all
P.fileFC = 'functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_REAL_0926_001.mat';
P.fileData = 'CHIMP_CWAS_DATA_FORDECO.mat';
P.nPerm  = 5000;
rng(2026);

F3 = load('results_Fig3_network_metrics.mat');      % indx, y, SCm, names
R = load(P.fileFC); RH = R.RH;
load(P.fileData, 'CJ_CHIMP');
indx = F3.indx(:); y = F3.y(:); n = numel(y);
if iscell(F3.names), F3.names = F3.names(:)'; else, F3.names = cellstr(F3.names)'; end
dens = nan(n,1); stren = nan(n,1);
for i = 1:n
    C = squeeze(CJ_CHIMP(:,:,indx(i))); C(isnan(C)) = 0;
    dens(i) = mean(C(:) > 0); stren(i) = mean(C(:));
end

% 26 Sep 2026: mean FC removed (r = 0.98 with synchrony); 5 dynamical measures
Xd   = [RH.Synchro_sub2(indx)', RH.Ignition(indx)', mean(RH.EdgeMetaNode(indx,:),2), RH.EntroEdgesMeta(indx)', RH.Metasim_sub2(indx)'];
labD = {'Synchrony','Ignition','Edge turbulence','Entropy eFCD','Metastability'};
colD = [38 16 22 30 10];                           % Fig 2 colour positions (thermal, 44)
X = [dens, stren, F3.SCm, Xd];
lab = [{'SC density','SC strength'}, strcat({'SC '}, F3.names), labD];
isDyn = [false(1, 2 + numel(F3.names)), true(1, numel(labD))];
nP = size(X,2);

%% ---- single-predictor held-out fits ----
perms_ = zeros(n, P.nPerm);
for k = 1:P.nPerm, perms_(:,k) = randperm(n)'; end
[r_ho, R2_ho, p] = deal(nan(nP,1));
for j = 1:nP
    [r_ho(j), R2_ho(j)] = loo_fit(X(:,j), y);
    nul = zeros(P.nPerm,1);
    for k = 1:P.nPerm, nul(k) = loo_fit(X(:,j), y(perms_(:,k))); end
    p(j) = (1 + sum(nul >= r_ho(j))) / (P.nPerm + 1);
end
q = nan(nP,1);                                     % FDR within each family (13 structural, 5 dynamical), as in Fig 3
q(~isDyn) = fdr_bh(p(~isDyn)); q(isDyn) = fdr_bh(p(isDyn));
qAll = fdr_bh(p);                                  % FDR across all predictors (printed for reference)
[~, o] = sort(R2_ho, 'descend');
fprintf('\n%-26s %10s %10s %8s %8s %9s\n', 'predictor', 'held-out r', 'held-out R2', 'p', 'q family', 'q all');
for j = o', fprintf('%-26s %+10.2f %+10.3f %8.4f %8.3f %9.3f\n', lab{j}, r_ho(j), R2_ho(j), p(j), q(j), qAll(j)); end

%% ---- secondary: Ignition on top of SC (and SC on top of Ignition) ----
% NOT REPORTED in the ms (28 Sep 2026): the combined SC model has negative held-out R2, so adding Ignition only
% brings a bad model close to zero. Kept (printed only) in case a reviewer asks.
% 28 Sep 2026: combined structural model = density, strength + one metric per family of Fig 3
% (segregation: modularity; integration: global efficiency; centrality: betweenness; resilience: assortativity).
% Variable name sc5 kept for continuity; it now holds 6 measures. Previous set: density, strength, modularity,
% transitivity, strength CV.
sc5 = [dens, stren, F3.SCm(:,3), F3.SCm(:,5), F3.SCm(:,7), F3.SCm(:,9)];
ign = RH.Ignition(indx)';
[~, R2sc] = loo_fit(sc5, y); [~, R2full] = loo_fit([sc5 ign], y);
nul = zeros(P.nPerm,1);
for k = 1:P.nPerm, [~, t] = loo_fit([sc5 ign(randperm(n))], y); nul(k) = t - R2sc; end
pGain = (1 + sum(nul >= R2full - R2sc)) / (P.nPerm + 1);
[Fa, pa] = nested_F(sc5, [sc5 ign], y); [Fb, pb] = nested_F(ign, [sc5 ign], y);
fprintf('\nSC-6 held-out R2 = %+.3f; SC-6 + Ignition = %+.3f; gain = %+.3f (permutation p = %.4f)\n', R2sc, R2full, R2full - R2sc, pGain);
fprintf('In-sample nested F: Ignition added to SC-6 p = %.4f; SC-6 added to Ignition p = %.3f\n', pa, pb);

%% ---- figure ----
% colours: the model measures keep their Fig 2 colours (cmocean 'thermal', 44: integration 4,
% Ignition 16, edge turbulence 22, entropy 30, synchrony 38, metastability 10); structural in grey.
cSC = [0.65 0.68 0.72];
Col = cmocean('thermal', 44);
dynCol = Col(colD,:);
lab{find(strcmp(lab,'Entropy eFCD'))} = 'Entropy (eFCD)';
dynIdx = find(isDyn);
cIgn = dynCol(2,:);

figure('Color','w','Position',[60 60 1300 520]);
ax = axes('position', [.13 .11 .3347 .815]); hold on   % panel a (same place as the old subplot(1,2,1))
for t = 1:nP
    j = o(t);
    if isDyn(j), c = dynCol(dynIdx == j,:); else, c = cSC; end
    barh(nP - t + 1, R2_ho(j), 0.7, 'FaceColor', c, 'EdgeColor', 'none');
    if p(j) < 0.001, txt = '  p < 0.001'; else, txt = sprintf('  p = %.3f', p(j)); end
    if q(j) < 0.05, txt = [txt ' *']; end
    text(max(R2_ho(j),0) + 0.003, nP - t + 1, txt, 'FontSize', 7.5);
end
xline(0,'k-'); set(gca, 'YTick', 1:nP, 'YTickLabel', fliplr(lab(o)), 'TickDir','out', 'FontSize', 8); box off
xlabel('held-out R^2 (variance explained in held-out animals)'); xlim([min(R2_ho)-0.02, max(R2_ho)+0.09]);
title({'a  Each measure alone (leave-one-out)', 'coloured: model dynamics (as in Fig 2); grey: structural measures; * FDR q < 0.05 within family'}, 'FontWeight','normal', 'FontSize', 9);

% b (28 Sep 2026): four scatters, observed vs held-out predicted score: combined SC model (6 measures; grey) and each
% dynamical measure that explains cognition in held-out animals (synchrony, Ignition, edge turbulence),
% in their Fig 2 colours; same axes in the four panels. Previous panel b: Fig4_heldout_explanation_panelb_OLD.m
[~, R2sc5, yhS] = loo_fit(sc5, y);
nulS = zeros(P.nPerm,1);                                 % permutation p for the combined SC model (as in panel a)
for k = 1:P.nPerm, nulS(k) = loo_fit(sc5, y(perms_(:,k))); end
pS5 = (1 + sum(nulS >= corr(y, yhS))) / (P.nPerm + 1);
bNames = {'SC, 6 measures', 'Synchrony', 'Ignition', 'Edge turbulence'};
YH = zeros(n, 4); YH(:,1) = yhS; bCol = zeros(4,3); bCol(1,:) = cSC * 0.85;
bP = nan(1,4); bQ = nan(1,4); bP(1) = pS5;
for b = 2:4
    jd = find(strcmp(labD, bNames{b})); [~, ~, YH(:,b)] = loo_fit(Xd(:,jd), y);
    bCol(b,:) = dynCol(jd,:); jj = find(strcmp(lab, bNames{b})); bP(b) = p(jj); bQ(b) = q(jj);
end
xl = [min(YH(:)) max(YH(:))] + [-.05 .05]*range(YH(:)); yl = [min(y) - .04, max(y) + .08];
bx = [.535 .765]; by = [.56 .10]; bw = .17; bh = .36;     % 2 x 2 grid right of panel a
for b = 1:4
    ax = axes('position', [bx(1 + mod(b-1,2)), by(1 + (b > 2)), bw, bh]); hold on
    yh = YH(:,b); c = bCol(b,:);
    mdl = fitlm(yh, y); xx = linspace(min(yh), max(yh), 100)'; [yp, yci] = predict(mdl, xx);
    fill([xx; flipud(xx)], [yci(:,1); flipud(yci(:,2))], c, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(xx, yp, '-', 'Color', c, 'LineWidth', 1.8);
    plot(xl, [1 1]*mean(y), ':', 'Color', [.5 .5 .5]);   % mean score (held-out R2 = 0 reference)
    plot(yh, y, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
    R2b = 1 - sum((y - yh).^2) / sum((y - mean(y)).^2);
    if bP(b) < 0.001, ptxt = 'p < 0.001'; else, ptxt = sprintf('p = %.3f', bP(b)); end
    if isnan(bQ(b)), stxt = ptxt; else, stxt = sprintf('%s, q = %.3f', ptxt, bQ(b)); end
    text(.04, .98, {sprintf('held-out r = %.2f, R^2 = %.2f', corr(y, yh), R2b), stxt}, 'units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 7.5);
    xlim(xl); ylim(yl); pbaspect([1 1 1]); set(ax, 'TickDir', 'out', 'FontSize', 7.5, 'Box', 'off')
    title(bNames{b}, 'FontWeight', 'bold', 'FontSize', 9, 'Color', c)
    if b > 2, xlabel('held-out predicted score'); end
    if mod(b,2) == 1, ylabel('observed score (PCTB)'); end
    if b == 1, text(-.35, 1.12, 'b', 'units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold'); end
end
print('-dsvg','Fig4_heldout_explanation.svg'); %print('-dpng','-r300','Fig4_heldout_explanation.png');

save('results_Fig4_heldout.mat', 'lab','X','y','r_ho','R2_ho','p','q','R2sc','R2full','pGain','pa','pb','YH','pS5','P');

%% ===================== local functions =====================
function [r, R2, yh] = loo_fit(X, y)
    A = [ones(numel(y),1) X];
    H = A * pinv(A'*A) * A';
    e = y - H*y;
    yh = y - e ./ (1 - diag(H));
    r = corr(y, yh);
    R2 = 1 - sum((y - yh).^2) / sum((y - mean(y)).^2);
end

function [F, p] = nested_F(Xr, Xf, y)
    n = numel(y);
    A0 = [ones(n,1) Xr]; A1 = [ones(n,1) Xf];
    r0 = sum((y - A0*(A0\y)).^2); r1 = sum((y - A1*(A1\y)).^2);
    d1 = size(A1,2) - size(A0,2); d2 = n - size(A1,2);
    F = ((r0 - r1)/d1) / (r1/d2); p = 1 - fcdf(F, d1, d2);
end

function q = fdr_bh(p)
    p = p(:); m = numel(p); [ps_, o] = sort(p);
    qs_ = ps_ .* m ./ (1:m)'; qs_ = flipud(cummin(flipud(qs_)));
    q = nan(m,1); q(o) = min(qs_, 1);
end
