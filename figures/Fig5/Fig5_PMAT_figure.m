% Fig5_PMAT_figure.m
% Fig 5 (human validation, part 2): the chimpanzee measures computed on empirical
% HCP 7T resting-state fMRI (n = 181, chimp band 0.02-0.03 Hz) and fluid reasoning
% (PMAT24), in the style and colours of Fig 2.
% Panel letters start at 'c' (a-b = working point, Fig5a_working_point.m).
%   c-g  PMAT24 vs the DYNAMICAL measures (synchrony, Ignition, edge turbulence, entropy eFCD,
%        metastability): plain Pearson r of raw PMAT24 vs raw measure, as in Fig 2; grey fit line,
%        95% confidence band in the measure colour; r, p and FDR q (across the 5 measures).
%   h    robustness, as Fig 2f: r (filled; line = leave-one-out range, r recomputed leaving out
%        each participant) and partial r controlling for age (HCP bins) and sex (open circles);
%        grey band = |r| not significant at n. Replaces the supplementary partial-correlation
%        figure (FigS_PMAT_partial.m, removed 29 Sep 2026).
%   i    same measures in both species: plain Pearson r with 95% CI (chimpanzee: model,
%        PCTB, n = 43; human: empirical, PMAT24, n = 181); "~=" marks a significant species
%        difference (Fisher z, p < 0.05). Chimp values from the Sep 2026 rerun (sweep frequencies).
%   j    held-out prediction, as the chimp Fig 4a: held-out R2 (bars; raw PMAT24 and measure, no
%        covariates; exact leave-one-out; p from 5000 permutations of the scores; * FDR q < 0.05)
%   k-m  observed vs held-out predicted PMAT24 for synchrony, entropy and metastability (as Fig 4b)
% Mean FC printed only (r = 0.93 with synchrony).
% The figure is saved as Fig5_PMAT.png (300 dpi) and Fig5_PMAT.svg.
% Needs: HCP7T_measures_for_chimp_comparison.mat, HCP7T_sex.mat, cmocean.m (all in this
% folder); Statistics Toolbox. Previous version (with MODE = 'partial' supplement):
% Fig5_PMAT_figure_withMODE_OLD.m.

clear all; close all
load('HCP7T_measures_for_chimp_comparison.mat');   % M (subjects x 9 measures x 2 bands), tests, age, ok
load('HCP7T_sex.mat');                             % sex (1 = M)

%% ---- data (band 1 = 0.02-0.03 Hz) ----
if iscell(age)
    ageN = cellfun(@(c) str2double(regexp(char(string(c)),'\d+','match','once')), age);
else
    ageN = double(age);
end
ageN = ageN(:); sex = double(sex(:));
pmat = tests(:,1);                                  % PMAT24_A_CR
X = M(:,[9 2 3 4 8],1);                             % Synchrony, Ignition, EdgeTurbulence, EntropyeFCD, Metastability (mean FC: printed only)
mFC = M(:,1,1);
keep = logical(ok(:)) & ~isnan(pmat) & ~isnan(ageN) & ~isnan(sex) & all(~isnan(X),2);
y = pmat(keep); X = X(keep,:); mFC = mFC(keep); ageN = ageN(keep); sex = sex(keep); n = numel(y);

names = {'Synchrony','Ignition','Edge turbulence','Entropy (eFCD)','Metastability'};
short = {'Synchrony','Ignition','Edge turbulence','Entropy','Metastability'};
nM = numel(names);
Col = cmocean('thermal', 44); cols = Col([38 16 22 30 10],:);  % Fig 2 colour positions
chimpR = [-0.43 -0.43 -0.41 -0.14 -0.17]; nC = 43;             % Fig 2, Sep 2026 rerun: synchrony, Ignition, turbulence, entropy, metastability
L0 = 'c';                                                        % first panel letter (a-b = working point)

%% ---- statistics ----
Z = [ones(n,1) ageN sex]; res = @(v) v - Z*(Z\v);
ry = res(y);
[pr, pp, rH, pH] = deal(nan(nM,1));
for j = 1:nM
    [pr(j), pp(j)] = partialcorr(X(:,j), y, [ageN sex]);
    [rH(j), pH(j)] = corr(X(:,j), y);
end
rS = rH; pS = pH; rLab = 'r';
q = fdr_bh(pS); qp = fdr_bh(pp);
[looLo, looHi] = deal(nan(nM,1));                  % leave-one-out range of r (as Fig 2f)
for j = 1:nM
    loo = arrayfun(@(i) corr(X(setdiff(1:n,i),j), y(setdiff(1:n,i))), 1:n);
    looLo(j) = min(loo); looHi(j) = max(loo);
end
tcr = tinv(0.975, n-2); rcrit = tcr / sqrt(n - 2 + tcr^2);
% held-out prediction, as the chimp Fig 4 (raw PMAT24 and raw measure, no covariates; 29 Sep 2026), exact
% leave-one-out linear prediction (hat matrix), held-out r and R2; p from 5000 permutations of the scores
% (one-sided on held-out r, same shuffles for all measures); FDR across the 5 measures
nPerm = 5000; rng(2026);
perms_ = zeros(n, nPerm); for k = 1:nPerm, perms_(:,k) = randperm(n)'; end
[rHO, R2HO, pHO] = deal(nan(nM,1)); YHO = nan(n, nM);
for j = 1:nM
    xr = X(:,j);
    [rHO(j), R2HO(j), YHO(:,j)] = loo_fit(xr, y);
    nul = zeros(nPerm,1);
    for k = 1:nPerm, nul(k) = loo_fit(xr, y(perms_(:,k))); end
    pHO(j) = (1 + sum(nul >= rHO(j))) / (nPerm + 1);
end
qHO = fdr_bh(pHO);
ciH = tanh(atanh(rH) + [-1 1]*1.96/sqrt(n-3));
ciC = tanh(atanh(chimpR(:)) + [-1 1]*1.96/sqrt(nC-3));
pDiff = 2*normcdf(-abs((atanh(chimpR(:)) - atanh(rH)) ./ sqrt(1/(nC-3) + 1/(n-3))));
fprintf('q = FDR across the %d measures\n', nM);
[rI, pI] = corr(mFC, y); [rIp, pIp] = partialcorr(mFC, y, [ageN sex]);
fprintf('Supplementary: mean FC vs synchrony r = %.2f; mean FC vs PMAT24 r = %+.2f (p = %.4f); partial r = %+.2f (p = %.4f)\n', corr(mFC, X(:,1)), rI, pI, rIp, pIp);
fprintf('%-22s %10s %10s %8s %8s\n','measure','held-out r','held-out R2','p','q');
for j = 1:nM, fprintf('%-22s %+10.2f %+10.3f %8.4f %8.3f\n', names{j}, rHO(j), R2HO(j), pHO(j), qHO(j)); end
fprintf('%-22s %7s %7s %6s | %-14s | %9s %7s %6s | %8s %8s\n','measure','r','p','q','LOO range','partial r','p','q','chimp r','p diff');
for j = 1:nM, fprintf('%-22s %+7.2f %7.4f %6.3f | %+.2f to %+.2f | %+9.2f %7.4f %6.3f | %+8.2f %8.3f\n', names{j}, rH(j), pH(j), q(j), looLo(j), looHi(j), pr(j), pp(j), qp(j), chimpR(j), pDiff(j)); end

%% ---- figure ----
figure
xSize = 18; ySize = 18; y0 = .72; hA = .22;     % 3 rows: scatters c-g; h, i; held-out j-m
set(gcf,'PaperUnits','centimeters','PaperPosition',[(21-xSize)/2 (30-ySize)/2 xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')
W = .11; xs = linspace(.065, .985 - W, nM);
for k = 1:nM
    ax = axes('position', [xs(k) y0 W hA]); hold on
    c = cols(k,:);
    xk = y; yk = X(:,k);
    mdl = fitlm(xk, yk);
    xx = linspace(min(xk), max(xk), 100)'; [yp, yci] = predict(mdl, xx);
    fill([xx; flipud(xx)], [yci(:,1); flipud(yci(:,2))], c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(xx, yp, '-', 'Color', [.45 .45 .45], 'LineWidth', 1.3);
    plot(xk, yk, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'MarkerSize', 3, 'LineWidth', 0.3);
    yl = [min(yk) max(yk)]; ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
    if pS(k) < 0.001, ptxt = sprintf('p < 0.001, q = %.3f', q(k)); else, ptxt = sprintf('p = %.3f, q = %.3f', pS(k), q(k)); end
    if q(k) < 0.05, fw = 'bold'; else, fw = 'normal'; end
    title({sprintf('%s = %.2f', rLab, rS(k)), ptxt}, 'FontSize', 6.5, 'FontWeight', fw)   % stats above the panel
    lt = char(L0+k-1);
    text(-.3, 1.35, lt, 'units','normalized','fontsize',12,'fontweight','bold')
    xlabel('PMAT24'); ylabel(short{k}, 'FontSize', 7.5);
    set(ax, 'TickDir','out', 'FontSize', 8, 'LineWidth', 0.8, 'Box', 'off')
    axis('square')
end

% h: robustness (as Fig 2f)
yy = nM:-1:1;
ax = axes('position', [.15 .40 .26 .22]); hold on
patch([-rcrit rcrit rcrit -rcrit], [0.4 0.4 nM+0.9 nM+0.9], [.94 .95 .96], 'EdgeColor','none');
for k = 1:nM
    plot([looLo(k) looHi(k)], [yy(k)+.12 yy(k)+.12], '-', 'Color', cols(k,:), 'LineWidth', 1.5);
    plot(rH(k), yy(k)+.12, 'o', 'MarkerFaceColor', cols(k,:), 'MarkerEdgeColor', cols(k,:), 'MarkerSize', 6);
    plot(pr(k), yy(k)-.12, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cols(k,:), 'MarkerSize', 6, 'LineWidth', 1.5);
end
h1 = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k'); h2 = plot(nan, nan, 'ko', 'MarkerFaceColor', 'w');
xline(0, 'k:'); xlim([-0.5 0.5]); ylim([0.4 nM+0.9]);
set(ax, 'YTick', fliplr(yy), 'YTickLabel', fliplr(short), 'TickDir','out', 'FontSize', 8, 'Box','off')
xlabel('r with PMAT24')
legend([h1 h2], {'r (line: leave-one-out range)', 'controlling age and sex'}, 'Location','northoutside', 'Box','off', 'FontSize', 6.5)
text(-.5, 1.05, char(L0+nM), 'units','normalized','fontsize',12,'fontweight','bold')

% i: both species (plain r in both, as in Fig 2)
yy = nM:-1:1;
ax = axes('position', [.58 .40 .38 .22]); hold on
for j = 1:nM
    plot(ciC(j,:), [yy(j)+.14 yy(j)+.14], '-', 'Color', cols(j,:), 'LineWidth', 1.6);
    hC = plot(chimpR(j), yy(j)+.14, 'd', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cols(j,:), 'LineWidth', 1.6, 'MarkerSize', 6);
    plot(ciH(j,:), [yy(j)-.14 yy(j)-.14], '-', 'Color', cols(j,:), 'LineWidth', 1.6);
    hH = plot(rH(j), yy(j)-.14, 'o', 'MarkerFaceColor', cols(j,:), 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
    if pDiff(j) < 0.05, text(0.66, yy(j), '\neq', 'FontSize', 11, 'Color', [.54 .58 .63]); end
end
h1 = plot(nan, nan, 'kd', 'MarkerFaceColor', 'w'); h2 = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k');
xline(0, 'k:'); xlim([-0.72 0.72]); ylim([0.4 nM+0.9]);
set(ax, 'YTick', 1:nM, 'YTickLabel', fliplr(short), 'TickDir','out', 'FontSize', 8, 'Box','off')

legend([h1 h2], {sprintf('Chimpanzees (model, PCTB, n = %d)', nC), sprintf('Humans (fMRI, PMAT24, n = %d)', n)}, ...
    'Location','northoutside','Box','off','FontSize',6.8)
xlabel({'r with cognition (95% CI)', '\neq: species differ, p < 0.05'})
text(-.16, 1.05, char(L0+nM+1), 'units','normalized','fontsize',12,'fontweight','bold')

% j: held-out prediction, as the chimp Fig 4a (bars) ...
ax = axes('position', [.15 .06 .18 .22]); hold on
for k = 1:nM
    barh(yy(k), R2HO(k), 0.6, 'FaceColor', cols(k,:), 'EdgeColor', 'none');
    if pHO(k) < 0.001, txt = '  p < 0.001'; else, txt = sprintf('  p = %.3f', pHO(k)); end
    if qHO(k) < 0.05, txt = [txt ' *']; end
    text(max(R2HO(k),0) + 0.002, yy(k), txt, 'FontSize', 6.5);
end
xline(0, 'k-'); ylim([0.4 nM+0.9]); xlim([min([R2HO; 0]) - 0.01, max(R2HO) + 0.06]);
set(ax, 'YTick', fliplr(yy), 'YTickLabel', fliplr(short), 'TickDir','out', 'FontSize', 8, 'Box','off')
xlabel('held-out R^2'); title({'Held-out prediction', '* FDR q < 0.05'}, 'FontWeight', 'normal', 'FontSize', 8.5)
text(-.7, 1.05, char(L0+nM+2), 'units','normalized','fontsize',12,'fontweight','bold')

% k-m: ... and observed vs held-out predicted PMAT24 (as the chimp Fig 4b) for synchrony and entropy
%      (significant held-out predictors) and metastability; same axes in the three panels
idxS = [1 4 5];                                     % Synchrony, Entropy (eFCD), Metastability
YS = YHO(:,idxS); xl = [min(YS(:)) max(YS(:))] + [-.05 .05]*range(YS(:)); yl = [min(y) - 1, max(y) + 1];
for t = 1:numel(idxS)
    k = idxS(t); c = cols(k,:); yh = YHO(:,k);
    ax = axes('position', [.42 + (t-1)*.195, .06, .16, .20]); hold on
    mdl = fitlm(yh, y); xx = linspace(min(yh), max(yh), 100)'; [yp, yci] = predict(mdl, xx);
    fill([xx; flipud(xx)], [yci(:,1); flipud(yci(:,2))], c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(xx, yp, '-', 'Color', c, 'LineWidth', 1.5);
    plot(yh, y, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'MarkerSize', 3, 'LineWidth', 0.3);
    if pHO(k) < 0.001, ptxt = sprintf('p < 0.001, q = %.3f', qHO(k)); else, ptxt = sprintf('p = %.3f, q = %.3f', pHO(k), qHO(k)); end
    if qHO(k) < 0.05, fw = 'bold'; else, fw = 'normal'; end
    xlim(xl); ylim(yl); pbaspect([1 1 1]); set(ax, 'TickDir','out', 'FontSize', 7, 'Box','off')
    title({['\color[rgb]{' sprintf('%.3f,%.3f,%.3f', c) '}\bf ' short{k}], ['\color{black}\rm ' sprintf('held-out r = %.2f, R^2 = %.2f', rHO(k), R2HO(k))], ptxt}, 'FontSize', 6.5, 'FontWeight', fw)
    xlabel('held-out predicted PMAT24'); if t == 1, ylabel('PMAT24'); end
    text(-.3, 1.45, char(L0+nM+2+t), 'units','normalized','fontsize',12,'fontweight','bold')
end

print('-dpng','-r300','Fig5_PMAT.png'); print('-dsvg','Fig5_PMAT.svg');
save('results_Fig5_PMAT.mat', 'names','n','rH','pH','q','pr','pp','qp','looLo','looHi','rHO','R2HO','pHO','qHO','YHO','chimpR','pDiff');

function [r, R2, yh] = loo_fit(x, y)
    A = [ones(numel(y),1) x]; H = A * pinv(A'*A) * A';
    e = y - H*y; yh = y - e ./ (1 - diag(H));
    r = corr(y, yh); R2 = 1 - sum((y - yh).^2) / sum((y - mean(y)).^2);
end

function q = fdr_bh(p)
    p = p(:); m = numel(p); [ps_, o] = sort(p);
    qs_ = ps_ .* m ./ (1:m)'; qs_ = flipud(cummin(flipud(qs_)));
    q = nan(m,1); q(o) = min(qs_, 1);
end
