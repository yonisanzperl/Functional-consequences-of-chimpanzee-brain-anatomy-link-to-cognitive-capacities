% Fig2_REAL0926_full.m
% Full Fig 2, in the style of the original Fig2_CHIMP_MVH.m, with the noise-matched
% measures (REAL_0926) and the 43 chimpanzees with SC and cognition.
%
%   a-e  scatter: cognitive score (PCTB) vs each model DYNAMICAL measure: synchrony, Ignition,
%        edge turbulence, entropy of eFCD, metastability (Sep 2026 rerun with the sweep node
%        frequencies). Mean FC (integration) is printed only (r = 0.98 with synchrony). FC modularity moved to Fig 3 (FC topology), 25 Sep 2026;
%        the previous version is kept as Fig2_REAL0926_full_withModularity_OLD.m.
%        cmocean 'thermal' colours, fitlm line (grey) with the 95% confidence band shaded in
%        the measure colour (26 Sep 2026), R and p top-left, n bottom-left
%   next robustness: r (filled; line = leave-one-out range) and partial r controlling for
%        SC density, SC strength and G* (open); grey band = not significant at n = 43
%   next working point: optimal coupling G* vs cognition
%   last real connectomes vs degree-preserving and degree + strength-preserving surrogate
%        connectomes (r of each surrogate run across animals). The surrogate runs in
%        P.nullDir were simulated with noise 0.005; set P.nullNoise = 0.04 once they are
%        rerun with the REAL_0926 settings and the "preliminary" label disappears.
%
% Also prints a table: r, p, FDR q (across the dynamical measures), partial r (p, q), leave-one-out range,
% surrogate mean r and rank p. Needs cmocean.m (in this folder) and the Statistics Toolbox.
% Run from this folder.

clear all; close all

P.fileMeta  = 'meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat';
P.fileData  = 'CHIMP_CWAS_DATA_FORDECO.mat';
P.fileFC    = 'functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_REAL_0926_001.mat';
P.nullDir   = 'surrogates_sweepfreq';    % Fig2_surrogate_connectomes.m (Sep 2026); old: 'validation_nulls_noise0.005_TO_RERUN'
P.nullNoise = 0.04;                       % noise of the surrogate runs (0.04 = matched to the real run)
P.Grange    = (0:0.1:10)/1000;            % coupling grid of the G sweep

load(P.fileMeta, 'Metasim_sub2');
load(P.fileData, 'CJ_CHIMP', 'AapCog_normalized');
load(P.fileFC, 'RH');

%% ---- subjects (as in Fig2_CHIMP_MVH.m) ----
nA = 59;
for ii = 1:nA, Metasim_sub(ii,:) = smooth(squeeze(mean(Metasim_sub2(ii,:,:),3))); end
[Meta, Gidx] = max(Metasim_sub, [], 2);
indx = find(Meta > 0);
indx = setdiff(indx, [23,27]);                       % SC outliers
indx = indx(~isnan(AapCog_normalized(indx)));        % animal without cognition score
n = numel(indx);
cog = AapCog_normalized(indx); cog = cog(:);
Gstar = P.Grange(Gidx(indx))'; Gstar = Gstar(:);
dens = nan(n,1); stren = nan(n,1);
for i = 1:n
    C = squeeze(CJ_CHIMP(:,:,indx(i))); C(isnan(C)) = 0;
    dens(i) = mean(C(:) > 0); stren(i) = mean(C(:));
end

%% ---- measures ----
% dynamical measures; colour positions in cmocean('thermal',44) are fixed per measure and shared
% with Figs 4-6: Integration 4, Ignition 16, edge turbulence 22, entropy 30, synchrony 38, metastability 10
% 26 Sep 2026: integration (mean FC) is no longer a separate measure (r = 0.98 with synchrony across
% animals); it is printed below for the supplementary table.
assert(isfield(RH, 'Synchro_sub2'), 'measures file without synchrony: use the Sep 2026 rerun');
D = { 'Synchrony',             'Synchrony',       38, @(R) R.Synchro_sub2(:);
      'Ignition',              'Ignition',        16, @(R) R.Ignition(:);
      'Edge turbulence',       'Edge turbulence', 22, @(R) mean(R.EdgeMetaNode,2);
      'Entropy (eFCD)',        'Entropy',         30, @(R) R.EntroEdgesMeta(:);
      'Metastability',         'Metastability',   10, @(R) R.Metasim_sub2(:) };
names = D(:,1)'; short = D(:,2)'; nM = numel(names);
getX = @(R) cell2mat(cellfun(@(g) safeget(g, R), D(:,4)', 'UniformOutput', false));
X = getX(RH); X = X(indx,:);
Col = cmocean('thermal', 44);
cols = Col([D{:,3}],:);
letters = 'abcdefghij';
tcr = tinv(0.975, n-2); rcrit = tcr / sqrt(n - 2 + tcr^2);

%% ---- statistics ----
[r, p, rp, pp, looLo, looHi] = deal(nan(nM,1));
for k = 1:nM
    [r(k), p(k)] = corr(X(:,k), cog);
    [rp(k), pp(k)] = partialcorr(X(:,k), cog, [dens stren Gstar]);
    loo = arrayfun(@(i) corr(X(setdiff(1:n,i),k), cog(setdiff(1:n,i))), 1:n);
    looLo(k) = min(loo); looHi(k) = max(loo);
end
q = fdr_bh(p); qp = fdr_bh(pp);
[rG, pGs] = corr(Gstar, cog);

% surrogate connectomes
fams = {'Degree',   'functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_Degree_*.mat';
        'Strength', 'functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_Strength__*.mat'};
rNull = cell(2,1); pRank = nan(nM,2);
isSurr = ~strcmp(short, 'Metastability');     % surrogate analysis without metastability (not related to cognition in e-f), 26 Sep 2026
for f = 1:2
    d = dir(fullfile(P.nullDir, fams{f,2}));
    d = d(~cellfun(@isempty, regexp({d.name}, '_\d{3}\.mat$', 'once')));
    rn = nan(numel(d), nM);
    for s = 1:numel(d)
        S = load(fullfile(d(s).folder, d(s).name), 'RH'); Xs = getX(S.RH); Xs = Xs(indx,:);
        for k = find(isSurr), rn(s,k) = corr(Xs(:,k), cog); end
    end
    rNull{f} = rn;
    for k = 1:nM
        if ~isSurr(k) || all(isnan(rn(:,k))), continue; end      % metastability left out / measure absent
        pRank(k,f) = (1 + sum(abs(rn(:,k)) >= abs(r(k)))) / (sum(~isnan(rn(:,k))) + 1);
    end
    nNull(f) = numel(d);
end
qRank = nan(nM,2);                               % FDR of the rank p across the surrogate-tested measures, per null family
for f = 1:2, qRank(isSurr,f) = fdr_bh(pRank(isSurr,f)); end

fprintf('\nn = %d; G* vs cognition r = %+.2f (p = %.2f)\n', n, rG, pGs);
mFC = RH.Integration(indx); mFC = mFC(:); [rI, pI] = corr(mFC, cog); [rIp, pIp] = partialcorr(mFC, cog, [dens stren Gstar]);
fprintf('Supplementary: mean FC vs synchrony r = %.2f; mean FC vs cognition r = %+.2f (p = %.4f); partial r = %+.2f (p = %.4f)\n', ...
    corr(mFC, X(:,1)), rI, pI, rIp, pIp);
fprintf('Surrogates: %d degree-preserving and %d degree + strength-preserving sets (rank p, FDR q across %d measures)\n', nNull(1), nNull(2), sum(isSurr));
fprintf('%-22s %6s %7s %6s | %8s %7s %6s | %-14s | %-24s %-24s\n', 'measure','r','p','q','partial','p','q','LOO range','degree null r (p, q)','strength null r (p, q)');
for k = 1:nM
    fprintf('%-22s %+6.2f %7.4f %6.3f | %+8.2f %7.4f %6.3f | %+.2f to %+.2f | %+.2f (%.3f, %.3f)      %+.2f (%.3f, %.3f)\n', names{k}, ...
        r(k), p(k), q(k), rp(k), pp(k), qp(k), looLo(k), looHi(k), mean(rNull{1}(:,k)), pRank(k,1), qRank(k,1), mean(rNull{2}(:,k)), pRank(k,2), qRank(k,2));
end

%% ---- figure ----
figure
xSize = 20; ySize = 10;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters','PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')
W = min(.13, .68/nM); xs = linspace(.065, .985 - W, nM);

for k = 1:nM
    ax = axes('position', [xs(k) .60 W .33]); hold on
    y = X(:,k); c = cols(k,:);
    mdl = fitlm(cog, y);                                          % linear fit, 95% confidence band (shaded)
    xx = linspace(min(cog), max(cog), 100)'; [yp, yci] = predict(mdl, xx);
    fill([xx; flipud(xx)], [yci(:,1); flipud(yci(:,2))], c, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(xx, yp, '-', 'Color', [.45 .45 .45], 'LineWidth', 1.3);
    plot(cog, y, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'MarkerSize', 7, 'LineWidth', 0.5);
    yl = [min([y; yci(:)]) max([y; yci(:)])]; ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.25*diff(yl)]);   % room for the text
    if p(k) >= 0.001, ptxt = sprintf('p = %.3f', p(k)); else, ptxt = 'p < 0.001'; end
    text(.05, .98, {sprintf('R = %.2f', r(k)), sprintf('%s, q = %.3f', ptxt, q(k))}, 'units','normalized','fontsize',7.5,'VerticalAlignment','top')   % q = FDR across the measures
    text(.05, .06, sprintf('\\itn\\rm = %d', n), 'units','normalized','fontsize',7.5)
    text(-.25, 1.05, letters(k), 'units','normalized','fontsize',12,'fontweight','bold')
    xlabel('Cognitive score (PCTB)'); ylabel(names{k})
    set(ax, 'TickDir','out', 'FontSize', 8, 'LineWidth', 0.8, 'Box', 'off')
    axis('square')
end

% f: robustness
yy = nM:-1:1;
ax = axes('position', [.10 .09 .22 .36]); hold on
patch([-rcrit rcrit rcrit -rcrit], [0.4 0.4 nM+0.9 nM+0.9], [.94 .95 .96], 'EdgeColor','none');
for k = 1:nM
    plot([looLo(k) looHi(k)], [yy(k)+.12 yy(k)+.12], '-', 'Color', cols(k,:), 'LineWidth', 1.5);
    plot(r(k), yy(k)+.12, 'o', 'MarkerFaceColor', cols(k,:), 'MarkerEdgeColor', cols(k,:), 'MarkerSize', 6);
    plot(rp(k), yy(k)-.12, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cols(k,:), 'MarkerSize', 6, 'LineWidth', 1.5);
end
h1 = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k'); h2 = plot(nan, nan, 'ko', 'MarkerFaceColor', 'w');
xline(0, 'k:'); xlim([-0.7 0.7]); ylim([0.4 nM+0.9]);
set(ax, 'YTick', fliplr(yy), 'YTickLabel', fliplr(short), 'TickDir','out', 'FontSize', 8, 'Box','off')
xlabel('r with cognition'); title('Robustness', 'FontWeight', 'normal', 'FontSize', 8.5)
legend([h1 h2], {'r (line: leave-one-out range)', 'controlling SC density, strength, G*'}, 'Location','northeast', 'Box','off', 'FontSize', 6.5)
text(-.35, 1.05, letters(nM+1), 'units','normalized','fontsize',12,'fontweight','bold')

% g: working point
ax = axes('position', [.42 .09 .20 .36]); hold on
plot(cog, Gstar*1000, 'o', 'MarkerFaceColor', [.54 .58 .63], 'MarkerEdgeColor', 'w', 'MarkerSize', 7);
text(.05, .98, sprintf('R = %.2f, p = %.2f', rG, pGs), 'units','normalized','fontsize',7.5,'VerticalAlignment','top')
xlabel('Cognitive score (PCTB)'); ylabel('Optimal coupling G^* (\times10^{-3})')
set(ax, 'TickDir','out', 'FontSize', 8, 'Box','off'); title('Working point', 'FontWeight', 'normal', 'FontSize', 8.5)
text(-.3, 1.05, letters(nM+2), 'units','normalized','fontsize',12,'fontweight','bold')

% h: surrogate connectomes -- measures with an association with cognition (raw or controlled);
% metastability is left out (not related to cognition in panels e and f), 26 Sep 2026
kS = find(isSurr); nS = numel(kS); yS = nS:-1:1;
ax = axes('position', [.76 .09 .22 .36]); hold on
patch([-rcrit rcrit rcrit -rcrit], [0.4 0.4 nS+0.9 nS+0.9], [.94 .95 .96], 'EdgeColor','none');
cN = [.72 .75 .78; .49 .53 .58];
for i = 1:nS
    k = kS(i);
    for f = 1:2
        rn = rNull{f}(:,k);
        plot(rn, yS(i) + (1.5-f)*0.3 + 0.03*randn(size(rn)), '.', 'Color', cN(f,:), 'MarkerSize', 7);
    end
    plot(r(k), yS(i), 'o', 'MarkerFaceColor', cols(k,:), 'MarkerEdgeColor', 'w', 'MarkerSize', 7);
end
h1 = plot(nan, nan, '.', 'Color', cN(1,:), 'MarkerSize', 10); h2 = plot(nan, nan, '.', 'Color', cN(2,:), 'MarkerSize', 10);
xline(0, 'k:'); xlim([-0.7 0.7]); ylim([0.4 nS+0.9]);
set(ax, 'YTick', 1:nS, 'YTickLabel', fliplr(short(kS)), 'TickDir','out', 'FontSize', 8, 'Box','off')
xlabel('r with cognition'); legend([h1 h2], {'degree-preserving SC', 'degree + strength-preserving SC'}, 'Location','northeast','Box','off','FontSize',6.5)
qMax = max(qRank(isSurr,:),[],'all');          % largest FDR-corrected rank p over measures and null models
text(.03, .03, sprintf('FDR-corrected rank p \\leq %.3f (%d + %d surrogates)', qMax, nNull(1), nNull(2)), ...
    'units','normalized', 'FontSize', 6.5, 'Color', [.35 .35 .35])
ttl = 'Real vs surrogate connectomes';
if P.nullNoise ~= 0.04, ttl = [ttl ' (preliminary)']; end
title(ttl, 'FontWeight', 'normal', 'FontSize', 8.5)
text(-.35, 1.05, letters(nM+3), 'units','normalized','fontsize',12,'fontweight','bold')

save('results_Fig2_REAL0926.mat', 'names','indx','cog','X','Gstar','dens','stren', ...
     'r','p','q','rp','pp','qp','looLo','looHi','rG','pGs','rNull','pRank','qRank','nNull','P');
 print('-dsvg','Fig2_REAL0926_full.svg'); 
 
 %print('-dpng','-r300','Fig2_REAL0926_full.png');

function x = safeget(g, R)
    try, x = g(R); catch, x = nan(numel(R.Integration),1); end   % field missing (e.g. surrogate files)
end

function q = fdr_bh(p)
    p = p(:); m = numel(p); [ps_, o] = sort(p);
    qs_ = ps_ .* m ./ (1:m)'; qs_ = flipud(cummin(flipud(qs_)));
    q = nan(m,1); q(o) = min(qs_, 1);
end
