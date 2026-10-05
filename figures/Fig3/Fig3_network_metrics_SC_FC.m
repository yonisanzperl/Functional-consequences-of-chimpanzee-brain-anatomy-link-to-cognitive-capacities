% Fig3_network_metrics_SC_FC.m
% Fig 3: full set of network metrics on the structural connectome (SC) and on the
% model functional connectivity (FC, REAL_0926: noise 0.04, same model as Fig 2),
% related to cognition in the same 43 chimpanzees as Fig 2.
%
% 11 metrics, computed identically on SC and FC (weights >= 0, diagonal 0, divided by
% the maximum weight; FC dense = all positive correlations kept):
%   segregation : clustering, transitivity, modularity, local efficiency
%   integration : global efficiency, characteristic path length
%   hubs        : betweenness, participation, assortativity, strength heterogeneity (CV)
%   level       : mean strength
% Path metrics use connection lengths 1/w. Local efficiency = global efficiency of the
% subgraph of each node's neighbours (Latora & Marchiori), averaged over nodes.
%
% MAIN FIG 3 (framing: FC network metrics outperform SC network metrics)
%   A  r with cognition of each metric on SC and on FC; FDR (Benjamini-Hochberg) applied
%      within each connectome type (11 SC tests, 11 FC tests), as stated in the Methods
%   B  |r| on FC vs |r| on SC for each metric; sign test (binomial) over the 11 metrics;
%      Williams tests per metric are printed (not significant individually: report the
%      consistent pattern, not single-metric differences)
%   C  the FC effects controlling for the working point G* (G* itself vs cognition printed)
% SUPPLEMENTARY FIG S3 (for reviewers)
%   A  relation of each FC metric to mean FC (weighted metrics covary with overall strength)
%   B  real FC vs strength-preserving rewired FC
%   C  correlation with cognition of (real - surrogate mean)
%
% Surrogate validation (FC): P.nSurr strength-preserving rewirings of each animal's FC.
% Each step picks four distinct regions a,b,c,d and moves weight
%   w_ab, w_cd  +delta   and   w_ad, w_cb  -delta
% so every region's strength (and therefore mean FC) is unchanged, and weights stay in
% [0,1]. P.nSwaps x (number of edges) steps per surrogate. Then:
%   rank p   : is the real r with cognition stronger than the r of surrogate k computed
%              across all animals, for k = 1..nSurr?  p = (1 + #|r_null| >= |r_real|)/(nSurr+1)
%   topology part : real metric minus the mean of that animal's surrogates, correlated
%              with cognition (what the arrangement of connections adds beyond strengths).
%
% Needs BCT (clustering_coef_wu, transitivity_wu, modularity_und, distance_wei,
% betweenness_wei, participation_coef, assortativity_wei), Statistics Toolbox, and for
% speed the Parallel Computing Toolbox (P.nWorkers > 0). Run from its own figure folder (all inputs are copied there).
% Same algorithms as the Python version used for the 24 Sep 2026 results (random streams differ).

clear all; close all
addpath(genpath('/Users/cbc/Documents/matlab_stuffs/BCT/'));

%% ---- parameters ----
P.fileFC   = 'functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_REAL_0926_001.mat';
P.fileMeta = 'meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat';
P.fileData = 'CHIMP_CWAS_DATA_FORDECO.mat';
P.Grange   = (0:0.1:10)/1000;       % grid of the coupling search (101 values)
P.nSurr    = 100;                   % surrogates per animal, put 100
P.nSwaps   = 10;                    % rewiring steps per edge
P.nWorkers = 6;
P.seed     = 2026;
P.useSavedMetrics = true;            % reuse SCm/FCm from results_Fig3_network_metrics*.mat if the animals match (seconds instead of minutes)
P.doSurrogates = false;              % FC rewiring for the supplementary figure (slow; not in the current ms version)

%% ---- data and subject selection (as in Fig2_CHIMP_MVH.m) ----
load(P.fileMeta, 'Metasim_sub2');
load(P.fileData, 'CJ_CHIMP', 'AapCog_normalized');
R = load(P.fileFC); RH = R.RH;
nA = size(CJ_CHIMP,3);
for s = 1:nA, Ms(s,:) = smooth(squeeze(mean(Metasim_sub2(s,:,:),3))); end
[Meta, Gidx] = max(Ms, [], 2);
indx = find(Meta > 0);
indx = setdiff(indx, [23 27]);
indx = indx(~isnan(AapCog_normalized(indx)));
y = AapCog_normalized(indx); y = y(:); n = numel(indx);
Gstar = P.Grange(Gidx(indx))'; meanFC = RH.Integration(indx); meanFC = meanFC(:);
N = size(CJ_CHIMP,1);
fprintf('n = %d animals\n', n);

names = {'Clustering','Transitivity','Modularity','Local efficiency','Global efficiency', ...
         'Char. path length','Betweenness','Participation','Assortativity', ...
         'Strength heterogeneity','Mean strength'};
nM = numel(names);

%% ---- metrics on SC and dense FC ----
SCm = nan(n,nM); FCm = nan(n,nM); loaded = false;
if P.useSavedMetrics
    for fr = {'results_Fig3_network_metrics_main.mat', 'results_Fig3_network_metrics.mat'}
        if exist(fr{1}, 'file')
            Sv = load(fr{1}, 'SCm', 'FCm', 'indx', 'P');
            sameFile = isfield(Sv,'P') && isfield(Sv.P,'fileFC') && strcmp(Sv.P.fileFC, P.fileFC);
            if isequal(Sv.indx(:), indx(:)) && sameFile && ~isempty(Sv.FCm)
                SCm = Sv.SCm; FCm = Sv.FCm; loaded = true;
                d_ = dir(fr{1}); fprintf('Network metrics loaded from %s (saved %s)\n', fr{1}, d_.date); break
            end
        end
    end
end
if ~loaded
    tic
    for i = 1:n
        s = indx(i);
        SCm(i,:) = net_metrics(prep(squeeze(CJ_CHIMP(:,:,s))));
        FCm(i,:) = net_metrics(prep(squeeze(RH.FCsim(s,:,:))));
        if mod(i,10) == 0, fprintf('metrics: %d/%d animals (%.1f min)\n', i, n, toc/60); end
    end
end

%% ---- statistics ----
[rS,pS,rF,pF,rFg,pFg,rFm,pFm,looLo,looHi] = deal(nan(nM,1));
for j = 1:nM
    [rS(j),pS(j)] = corr(SCm(:,j), y);
    [rF(j),pF(j)] = corr(FCm(:,j), y);
    [rFg(j),pFg(j)] = partialcorr(FCm(:,j), y, Gstar);
    [rFm(j),pFm(j)] = partialcorr(FCm(:,j), y, meanFC);
    loo = arrayfun(@(k) corr(FCm(setdiff(1:n,k),j), y(setdiff(1:n,k))), 1:n);
    looLo(j) = min(loo); looHi(j) = max(loo);
end
qS = fdr_bh(pS); qF = fdr_bh(pF); qFg = fdr_bh(pFg);           % FDR within each connectome type
[rGc, pGc] = corr(Gstar, y);
fprintf('\nG*: r with cognition %+.2f (p = %.3f), with mean FC %+.2f\n', rGc, pGc, corr(Gstar,meanFC));
% FC vs SC: sign test over metrics and Williams test per metric
nUp = sum(abs(rF) > abs(rS)); pSign = 2*binocdf(min(nUp, nM-nUp), nM, 0.5); pSign = min(pSign,1);
[tW, pW] = deal(nan(nM,1));
for j = 1:nM
    r23 = corr(FCm(:,j), SCm(:,j)) * sign(rF(j)) * sign(rS(j));
    [tW(j), pW(j)] = williams(abs(rF(j)), abs(rS(j)), r23, n);
end
fprintf('|r| FC > |r| SC in %d of %d metrics (sign test p = %.3f)\n', nUp, nM, pSign);
fprintf('\n%-24s | %-16s | %-16s | %-14s | %-15s | %-15s | %s\n', 'metric', 'SC r (q)', 'FC r (q)', 'FC>SC Williams p', 'FC r | G* (p)', 'FC r | mFC (p)', 'FC LOO range');
for j = 1:nM
    fprintf('%-24s | %+.2f (%.3f)    | %+.2f (%.3f)    | %14.3f | %+.2f (%.3f)   | %+.2f (%.3f)   | %+.2f to %+.2f\n', ...
        names{j}, rS(j), qS(j), rF(j), qF(j), pW(j), rFg(j), pFg(j), rFm(j), pFm(j), looLo(j), looHi(j));
end

if P.doSurrogates
%% ---- surrogate validation (strength-preserving FC rewiring) ----
Msur = nan(n, P.nSurr, nM);
FCall = cell(n,1);
for i = 1:n, FCall{i} = prep(squeeze(RH.FCsim(indx(i),:,:))); end
tic
parfor (i = 1:n, P.nWorkers)
    rs = RandStream('mt19937ar', 'Seed', P.seed + indx(i));
    tmp = nan(P.nSurr, nM);
    for k = 1:P.nSurr
        S = rewire_strength(FCall{i}, P.nSwaps, rs);
        tmp(k,:) = net_metrics(S);          % no re-normalisation: strengths stay identical
    end
    Msur(i,:,:) = tmp;
end
fprintf('\nSurrogates done in %.1f min\n', toc/60);

[rNullMean, rNullSD, pRank, rTopo, pTopo] = deal(nan(nM,1));
for j = 1:nM
    rn = arrayfun(@(k) corr(squeeze(Msur(:,k,j)), y), 1:P.nSurr);
    rNullMean(j) = mean(rn); rNullSD(j) = std(rn);
    pRank(j) = (1 + sum(abs(rn) >= abs(rF(j)))) / (P.nSurr + 1);
    topo = FCm(:,j) - mean(squeeze(Msur(:,:,j)), 2);
    [rTopo(j), pTopo(j)] = corr(topo, y);
end
qTopo = fdr_bh(pTopo);
fprintf('\n%-24s | %7s | %-18s | %7s | %-22s\n', 'FC metric', 'real r', 'surrogate r', 'rank p', 'topology part r (q)');
for j = 1:nM
    fprintf('%-24s | %+7.2f | %+.2f +- %.2f      | %7.3f | %+.2f (p %.3f, q %.3f)\n', names{j}, rF(j), ...
        rNullMean(j), rNullSD(j), pRank(j), rTopo(j), pTopo(j), qTopo(j));
end

end

%% ---- MAIN FIG 3 (28 Sep 2026: SC | model FC columns grouped by metric family; previous: Fig3_network_metrics_SC_FC_summaryonly_OLD.m) ----
%   a  left column: cognition vs each SC metric (grey); right column: same metric on the model FC (blue);
%      rows grouped by family (segregation, integration, centrality, resilience & strength); square panels,
%      metric values on the y axis; linear fit with shaded 95% CI; r and q (FDR within connectome type)
%      above each panel (bold = q < 0.05); p printed in the Command Window
%      left margin: network scheme of each family (highlighted in black: triangles / shortest path /
%      connector hub / strong, mutually linked nodes)
%   b  r with cognition of each metric, same family order: SC (grey), model FC (blue), filled = q < 0.05;
%      FC controlling for G* (open blue diamonds)
%   c  |r| on model FC vs |r| on SC; sign test over the 11 metrics
cSC = [0.55 0.58 0.62]; cFC = [0.17 0.37 0.54]; band = [0.94 0.95 0.96];
short = {'Clustering','Transitivity','Modularity','Local eff.','Global eff.','Path length', ...
         'Betweenness','Participation','Assortativity','Strength het.','Mean strength'};
famName = {'Segregation','Integration','Centrality',{'Resilience &','strength'}};
famIdx  = {[1 2 3 4], [5 6], [7 8], [9 10 11]};          % metrics are already stored in this order
tcr = tinv(0.975, n-2); rcrit = tcr / sqrt(n - 2 + tcr^2);
figure
xSize = 18; ySize = 24;
set(gcf,'PaperUnits','centimeters','PaperPosition',[(21-xSize)/2 (30-ySize)/2 xSize ySize])
set(gcf,'Position',[40 40 xSize*40 ySize*40],'Color','w')
sq = 1.35;  w = sq/xSize; h = sq/ySize;                   % square panels of 1.35 cm
rowStep = (sq + .42)/ySize; famGap = .40/ySize; top = 1 - 1.2/ySize;
colX = [.235, .235 + w + .075];                           % SC | model FC (left margin holds family schemes)
ypos = zeros(1,nM); famTop = zeros(1,4); famBot = zeros(1,4); yc = top - h;
for f = 1:4
    famTop(f) = yc + h;
    for j = famIdx{f}, ypos(j) = yc; yc = yc - rowStep; end
    famBot(f) = ypos(famIdx{f}(end)); yc = yc - famGap;
end
Mcol = {SCm, FCm}; rCol = {rS, rF}; qCol = {qS, qF}; cCol = {cSC, cFC}; colName = {'SC','model FC'};
for cI = 1:2
    annotation('textbox', [colX(cI) - .02, top + .1/ySize, w + .04, .5/ySize], 'String', colName{cI}, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','EdgeColor','none', ...
        'FontSize', 8.5, 'FontWeight','bold', 'Color', cCol{cI});
    for f = 1:4
        for j = famIdx{f}
            ax = axes('position', [colX(cI) ypos(j) w h]); hold on
            v = Mcol{cI}(:,j); c = cCol{cI};
            mdl = fitlm(y, v); xx = linspace(min(y), max(y), 100)'; [yp, yci] = predict(mdl, xx);
            fill([xx; flipud(xx)], [yci(:,1); flipud(yci(:,2))], c, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
            plot(xx, yp, '-', 'Color', [.45 .45 .45], 'LineWidth', 1);
            plot(y, v, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'MarkerSize', 3, 'LineWidth', 0.3);
            yl = [min([v; yci(:)]) max([v; yci(:)])]; ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
            xlim([min(y) max(y)] + [-.03 .03]*range(y))
            if qCol{cI}(j) < 0.05, fw = 'bold'; else, fw = 'normal'; end
            title(sprintf('r = %.2f, q = %.3f', rCol{cI}(j), qCol{cI}(j)), 'FontSize', 5.5, 'FontWeight', fw)
            set(ax, 'TickDir','out', 'FontSize', 5.5, 'LineWidth', .6, 'Box','off', 'TickLength', [.03 .03])
            ax.YAxis.TickLabelFormat = '%.3g'; ax.YAxis.Exponent = 0;
            yticks(ax, linspace(ax.YLim(1) + .1*diff(ax.YLim), ax.YLim(2) - .1*diff(ax.YLim), 3));   % metric values
            if j ~= famIdx{f}(end), set(ax, 'XTickLabel', []); end
            if j == nM, xlabel('PCTB', 'FontSize', 7); end
            if cI == 1, ylabel(short{j}, 'FontSize', 6.5); end
        end
    end
end
axT = axes('position', [0 0 1 1], 'Visible', 'off'); xlim(axT, [0 1]); ylim(axT, [0 1]); hold(axT, 'on');
iSz = 1.75;                                                % scheme size (cm)
for f = 1:4                                                % family brackets + network scheme + name
    xb = .155; ym = (famTop(f) + famBot(f))/2;
    annotation('line', [xb+.008 xb], [famTop(f) famTop(f)], 'Color', [.3 .3 .3], 'LineWidth', .8);
    annotation('line', [xb xb], [famBot(f) famTop(f)], 'Color', [.3 .3 .3], 'LineWidth', .8);
    annotation('line', [xb xb+.008], [famBot(f) famBot(f)], 'Color', [.3 .3 .3], 'LineWidth', .8);
    axI = axes('position', [xb - .012 - iSz/xSize, ym - iSz/ySize/2, iSz/xSize, iSz/ySize]);
    family_scheme(axI, f);
    text(axT, .008, ym, famName{f}, 'Rotation', 90, 'HorizontalAlignment','center', ...
        'VerticalAlignment','top', 'FontSize', 7, 'FontWeight','bold', 'Color', [.3 .3 .3]);
end
text(axT, .01, top + .5/ySize, 'a', 'fontsize', 12, 'fontweight', 'bold')

% b: r per metric (same family order, separators between families)
yy = nM:-1:1;
ax = axes('position', [.62 .53 .345 .40]); hold on
patch([-rcrit rcrit rcrit -rcrit], [0.4 0.4 nM+0.6 nM+0.6], band, 'EdgeColor','none');
k = 0;
for f = 1:4
    idx = famIdx{f}; k = k + numel(idx);
    if f < 4, yline(nM - k + 0.5, 'Color', [.75 .75 .75], 'LineWidth', .6); end
    fn = famName{f}; if iscell(fn), fn = strjoin(fn, ' '); end
    text(0.63, mean(yy(idx)), fn, 'Rotation', -90, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize', 6, 'Color', [.3 .3 .3]);
end
for j = 1:nM
    if qS(j) < 0.05, fS = cSC; else, fS = 'w'; end
    if qF(j) < 0.05, fF = cFC; else, fF = 'w'; end
    plot(rS(j), yy(j)+0.2, 'o', 'MarkerFaceColor', fS, 'MarkerEdgeColor', cSC, 'MarkerSize', 5, 'LineWidth', 1.2);
    plot(rF(j), yy(j), 'o', 'MarkerFaceColor', fF, 'MarkerEdgeColor', cFC, 'MarkerSize', 5, 'LineWidth', 1.2);
    plot(rFg(j), yy(j)-0.2, 'd', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cFC, 'MarkerSize', 4.5, 'LineWidth', 1);
end
h1 = plot(nan,nan,'o','MarkerFaceColor',cSC,'MarkerEdgeColor',cSC,'MarkerSize',5);
h2 = plot(nan,nan,'o','MarkerFaceColor',cFC,'MarkerEdgeColor',cFC,'MarkerSize',5);
h3 = plot(nan,nan,'d','MarkerFaceColor','w','MarkerEdgeColor',cFC,'MarkerSize',4.5);
xline(0,'k:'); xlim([-0.6 0.6]); ylim([0.4 nM+0.6]);
set(ax,'YTick',1:nM,'YTickLabel',fliplr(short),'TickDir','out','FontSize',6.5,'Box','off')
xlabel('r with cognition (filled: FDR q < 0.05)', 'FontSize', 7)
lg = legend([h1 h2 h3], {'SC','model FC','model FC, controlling G*'}, 'Orientation','horizontal','Box','off','FontSize',6);
lg.Position(1:2) = [.64 .455];
text(-.36, 1.02, 'b', 'units','normalized','fontsize',12,'fontweight','bold')

% c: |r| FC vs |r| SC
ax = axes('position', [.62 .06 .34 .34*xSize/ySize]); hold on
plot([0 .5], [0 .5], 'k:');
dxy = repmat([.008 .006], nM, 1);                          % label offsets to avoid overlaps
dxy(3,:) = [-.075 .004]; dxy(4,:) = [.012 -.013]; dxy(6,:) = [-.02 -.03]; dxy(8,:) = [-.085 .006]; dxy(10,:) = [.008 -.016];
for j = 1:nM
    if abs(rF(j)) > abs(rS(j)), c = cFC; else, c = cSC; end
    plot(abs(rS(j)), abs(rF(j)), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
    text(abs(rS(j)) + dxy(j,1), abs(rF(j)) + dxy(j,2), short{j}, 'FontSize', 5.5);
end
xlim([0 .5]); ylim([0 .5]); axis square
set(ax,'TickDir','out','FontSize',6.5,'Box','off')
xlabel('|r| with cognition, SC', 'FontSize', 7); ylabel('|r| with cognition, model FC', 'FontSize', 7)
title(sprintf('FC > SC in %d of %d metrics (sign test p = %.3f)', nUp, nM, pSign), 'FontSize', 7, 'FontWeight','normal')
text(-.24, 1.04, 'c', 'units','normalized','fontsize',12,'fontweight','bold')
print('-dpng','-r300','Fig3_main_SC_vs_FC.png'); print('-dsvg','Fig3_main_SC_vs_FC.svg');

if P.doSurrogates
%% ---- SUPPLEMENTARY FIG S3 ----
figure('Color','w','Position',[60 60 1500 520]);
subplot(1,3,1); hold on
rm = arrayfun(@(j) corr(FCm(:,j), meanFC), 1:nM);
plot(rm, rF, 'o', 'MarkerFaceColor', cFC, 'MarkerEdgeColor','w', 'MarkerSize', 7);
plot([-1 1], [-1 1]*corr(meanFC,y), '--', 'Color', cSur);
for j = 1:nM, text(rm(j)+0.03, rF(j), names{j}, 'FontSize', 8); end
xline(0,'k:'); yline(0,'k:'); xlim([-1 1]); ylim([-0.55 0.55]); box off; set(gca,'TickDir','out');
xlabel('r(FC metric, mean FC) across animals'); ylabel('r with cognition');
title('A  FC metrics and overall FC strength');

subplot(1,3,2); hold on
patch([-0.30 0.30 0.30 -0.30], [0.4 0.4 nM+0.6 nM+0.6], band, 'EdgeColor','none');
for j = 1:nM
    rn = arrayfun(@(k) corr(squeeze(Msur(:,k,j)), y), 1:P.nSurr);
    plot(rn, yy(j) + 0.08*randn(1,P.nSurr), '.', 'Color', cSur, 'MarkerSize', 6);
    plot(rF(j), yy(j), 'o', 'MarkerFaceColor', cFC, 'MarkerEdgeColor','w', 'MarkerSize', 8);
end
xline(0,'k:'); set(gca,'YTick',fliplr(yy),'YTickLabel',fliplr(names),'TickDir','out'); xlim([-0.6 0.6]); ylim([0.4 nM+0.6]); box off
xlabel('r with cognition'); title(sprintf('B  Real FC (blue) vs %d strength-preserving rewirings', P.nSurr));

subplot(1,3,3); hold on
patch([-0.30 0.30 0.30 -0.30], [0.4 0.4 nM+0.6 nM+0.6], band, 'EdgeColor','none');
for j = 1:nM
    if pTopo(j) < 0.05, c = cFC; else, c = [0.55 0.69 0.83]; end
    barh(yy(j), rTopo(j), 0.6, 'FaceColor', c, 'EdgeColor', 'none');
    text(rTopo(j) + 0.02*sign(rTopo(j)+eps), yy(j), sprintf('p = %.3f', pTopo(j)), 'FontSize', 7.5, ...
        'HorizontalAlignment', ternary(rTopo(j) >= 0, 'left', 'right'));
end
xline(0,'k-'); set(gca,'YTick',[],'TickDir','out'); xlim([-0.6 0.6]); ylim([0.4 nM+0.6]); box off
xlabel('r with cognition of (real - surrogate mean)'); title('C  Beyond each region''s strength');

end   % P.doSurrogates (supplementary figure)
% without surrogates, save to a separate file so the existing results with surrogates are not overwritten
if P.doSurrogates, fRes = 'results_Fig3_network_metrics.mat';
else, fRes = 'results_Fig3_network_metrics_main.mat'; [Msur, rNullMean, rNullSD, pRank, rTopo, pTopo, qTopo] = deal([]); end
save(fRes, 'names','indx','y','Gstar','meanFC','SCm','FCm','Msur', ...
     'rS','pS','qS','rF','pF','qF','rFg','pFg','qFg','rFm','pFm','looLo','looHi', ...
     'nUp','pSign','tW','pW','rGc','pGc', ...
     'rNullMean','rNullSD','pRank','rTopo','pTopo','qTopo','P');

%% ===================== local functions =====================
function W = prep(W)
    W(isnan(W)) = 0; W = (W + W')/2; W(1:size(W,1)+1:end) = 0; W(W < 0) = 0;
    W = W / max(W(:));
end

function v = net_metrics(W)
    n = size(W,1); k = sum(W > 0, 2); s = sum(W, 2);
    cc = clustering_coef_wu(W);
    tr = transitivity_wu(W);
    [Ci, Q] = modularity_und(W);
    L = zeros(n); L(W > 0) = 1 ./ W(W > 0);
    D = dist_fast(L); off = ~eye(n);
    Eg = mean(1 ./ D(off));
    fin = isfinite(D) & off; cpl = mean(D(fin));
    el = zeros(n,1);
    for u = 1:n
        nb = find(W(u,:) > 0);
        if numel(nb) > 1
            Ls = L(nb,nb); Ds = dist_fast(Ls); o2 = ~eye(numel(nb));
            el(u) = mean(1 ./ Ds(o2));
        end
    end
    bc = betweenness_wei(L);          % scale differs from networkx by a constant (irrelevant for r)
    pc = participation_coef(W, Ci);
    ast = assortativity_wei(W, 0);
    v = [mean(cc), tr, Q, mean(el), Eg, cpl, mean(bc), mean(pc), ast, std(s)/mean(s), mean(s)];
end

function D = dist_fast(L)
    % shortest-path lengths on connection lengths L (0 = no connection); same as BCT distance_wei
    % but with MATLAB's compiled graph routines (much faster for the 114-node dense FC)
    D = distances(graph(L, 'upper'));
    D(1:size(D,1)+1:end) = 0;
end

function W = rewire_strength(W, nSwaps, rs)
    n = size(W,1); nE = n*(n-1)/2; done = 0;
    while done < nSwaps*nE
        q = randi(rs, n, 1, 4);
        a = q(1); b = q(2); c = q(3); d = q(4);
        if a==b || a==c || a==d || b==c || b==d || c==d, continue; end
        lo = max(-min(W(a,b), W(c,d)), max(W(a,d), W(c,b)) - 1);
        hi = min(min(W(a,d), W(c,b)), 1 - max(W(a,b), W(c,d)));
        de = lo + rand(rs)*(hi - lo);
        W(a,b) = W(a,b) + de; W(b,a) = W(a,b);
        W(c,d) = W(c,d) + de; W(d,c) = W(c,d);
        W(a,d) = W(a,d) - de; W(d,a) = W(a,d);
        W(c,b) = W(c,b) - de; W(b,c) = W(c,b);
        done = done + 1;
    end
end

function [t, p] = williams(r12, r13, r23, n)
    % Williams (1959) test for two dependent correlations sharing variable 1 (cognition)
    detR = 1 - r12^2 - r13^2 - r23^2 + 2*r12*r13*r23; rbar = (r12 + r13)/2;
    t = (r12 - r13) * sqrt((n-1)*(1+r23) / (2*(n-1)/(n-3)*detR + rbar^2*(1-r23)^3));
    p = 2*tcdf(-abs(t), n-3);
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

function q = fdr_bh(p)
    p = p(:); m = numel(p); [ps_, o] = sort(p);
    qs_ = ps_ .* m ./ (1:m)'; qs_ = flipud(cummin(flipud(qs_)));
    q = nan(m,1); q(o) = min(qs_, 1);
end

function family_scheme(ax, f)
% small network cartoon for each family of graph measures (Fig 3a)
g = [.62 .62 .62]; k = [.12 .12 .12]; hold(ax, 'on');
switch f
    case 1   % segregation: two dense modules, triangles highlighted
        P = [.12 .70; .38 .88; .40 .55; .15 .40; .62 .45; .88 .60; .85 .22; .60 .12];
        scheme_edges(ax, P, [1 2; 2 3; 3 1; 1 4; 3 4; 5 6; 6 7; 7 8; 8 5; 5 7; 3 5], g, .7);
        scheme_edges(ax, P, [1 2; 2 3; 3 1; 5 6; 6 7; 7 5], k, 1.2);
        scheme_nodes(ax, P, [4 8], 14, g); scheme_nodes(ax, P, [1 2 3 5 6 7], 14, k);
    case 2   % integration: shortest path between distant nodes
        P = [.08 .50; .28 .80; .30 .25; .52 .55; .72 .85; .72 .20; .92 .50; .50 .05];
        scheme_edges(ax, P, [1 2; 1 3; 2 4; 3 4; 4 5; 4 6; 5 7; 6 7; 3 8; 8 6; 2 5], g, .7);
        scheme_edges(ax, P, [1 2; 2 5; 5 7], k, 1.5);
        scheme_nodes(ax, P, [3 4 6 8], 14, g); scheme_nodes(ax, P, [2 5], 14, k); scheme_nodes(ax, P, [1 7], 24, k);
    case 3   % centrality: connector hub between modules
        P = [.50 .50; .12 .80; .25 .92; .10 .52; .88 .80; .75 .92; .90 .52; .30 .12; .70 .12];
        scheme_edges(ax, P, [2 3; 3 4; 2 4; 5 6; 6 7; 5 7; 8 9], g, .7);
        scheme_edges(ax, P, [1 2; 1 4; 1 5; 1 7; 1 8; 1 9], k, 1.1);
        scheme_nodes(ax, P, 2:9, 14, g); scheme_nodes(ax, P, 1, 50, k);
    case 4   % resilience & strength: node strength and assortative mixing
        P = [.20 .75; .50 .88; .80 .72; .25 .25; .55 .42; .85 .20; .08 .45]; S = [42 58 38 10 14 8 8];
        scheme_edges(ax, P, [1 2; 2 3; 1 3], k, 1.7);
        scheme_edges(ax, P, [4 5; 5 6; 4 7; 1 7; 3 5], g, .5);
        scheme_nodes(ax, P, 4:7, S(4:7), g); scheme_nodes(ax, P, 1:3, S(1:3), k);
end
xlim(ax, [-.05 1.05]); ylim(ax, [-.05 1.05]); axis(ax, 'equal', 'off');
end

function scheme_edges(ax, P, pr, c, lw)
for ii = 1:size(pr,1), plot(ax, P(pr(ii,:),1), P(pr(ii,:),2), '-', 'Color', c, 'LineWidth', lw); end
end

function scheme_nodes(ax, P, idx, sz, c)
scatter(ax, P(idx,1), P(idx,2), sz, c, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', .4);
end
