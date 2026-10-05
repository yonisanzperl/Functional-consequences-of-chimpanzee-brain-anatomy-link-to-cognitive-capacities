% Fig6_alt_model_level.m
% ALTERNATIVE Fig 6 (30 Sep 2026): methodological framing, model level only.
% The model infers functional dynamics from each connectome, which allows a fair comparison between
% species even where no fMRI is available. Only the model's own quantities are compared: the
% metastability landscape across the global coupling G, the working point G* (maximum of metastability)
% and metastability at G*. No link to cognition, no selection among dynamical or topological measures.
% The full version (dynamical measures, FC modularity, linking panel) is kept in Fig6_cross_species.m.
%   a  metastability vs G (mean +- s.d. across individuals, V2 sweep; each individual's sweep is
%      interpolated on a common grid) with each individual's G* marked above the curves
%   b  G* humans vs chimpanzees (V2), raincloud (box, points, half violin): r_rb (> 0 = human higher) and
%      Mann-Whitney p in V2, V0 and V1
%   c  metastability at G* (200 simulations at G*), same statistics
% The other dynamical measures (synchrony, Ignition, edge turbulence, entropy eFCD) are printed only: no
% species difference in V0, V1 or V2 (reported in the text). Panel d (harmonisation incl. FC-SC corr) removed.
% Colours: paper palette (cmocean thermal, 44 levels): metastability #10 (as in Figs 2 and 5), G* #33/#37;
% humans full tone, chimpanzees 50% lighter tint. Needs cmocean.m (in this folder).
% Data: HumChimp_V0/V1/V2_sweep_measures.mat (58 humans, 22 chimpanzees, BB50, 74 regions).
% Saves Fig6_alt_model_level.png (300 dpi), .svg and results_Fig6_alt.mat.

clear all; close all
Vs = {'V0','V1','V2'}; iV = 3;                         % V2 shown in a-c
Col = cmocean('thermal', 44); tint = @(c, f) c*(1-f) + f;          % paper palette (Figs 2-5)
mH = Col(10,:); mC = tint(mH, .5);      % metastability colour (as in Figs 2 and 5): humans dark, chimpanzees light
gH = Col(33,:); gC = Col(37,:);         % G*: two close orange/yellow tones of the same palette (humans #33, chimpanzees #37)
cH = mH; cC = mC;
for v = 1:3, Dv{v} = load(sprintf('HumChimp_%s_sweep_measures.mat', Vs{v})); end
D = Dv{iV}; nH = numel(D.GstarH); nC = numel(D.GstarCh);

%% ---- statistics ----
lab = {'Metastability','G*'};
[rb, p] = deal(nan(2,3));
for v = 1:3
    [rb(1,v), p(1,v)] = rbiserial(Dv{v}.RH.Metasim(:), Dv{v}.RCh.Metasim(:));
    [rb(2,v), p(2,v)] = rbiserial(Dv{v}.GstarH(:), Dv{v}.GstarCh(:));
end
fprintf('\nSpecies differences (r_rb > 0 = human higher; Mann-Whitney p)\n%-20s %16s %16s %16s\n', '', 'V0', 'V1', 'V2');
for k = 1:2, fprintf('%-20s %+6.2f (p %6.2g) %+6.2f (p %6.2g) %+6.2f (p %6.2g)\n', lab{k}, rb(k,1), p(k,1), rb(k,2), p(k,2), rb(k,3), p(k,3)); end
fprintf('Median G*: humans %.4f, chimpanzees %.4f (V2)\n', median(D.GstarH), median(D.GstarCh));
om = {'Syncsim','Ignition','EdgeTurb','EntroEdgesMeta','Metasim'}; ol = {'Synchrony','Ignition','Edge turbulence','Entropy eFCD','Metastability'};
fprintf('\nAll dynamical measures (FDR across the 5 within each version)\n');
for v = 1:3
    [ro, po] = deal(nan(5,1));
    for k = 1:5, [ro(k), po(k)] = rbiserial(Dv{v}.RH.(om{k})(:), Dv{v}.RCh.(om{k})(:)); end
    qo = fdr_bh(po);
    for k = 1:5, fprintf('%s %-16s r_rb %+5.2f  p %.3g  p_FDR %.3g\n', Vs{v}, ol{k}, ro(k), po(k), qo(k)); end
end

gg = 0:0.0005:0.02;                                    % common grid for the sweep curves
curv = @(S) cell2mat(arrayfun(@(i) interp1(D.Grange(~isnan(S(i,:))), S(i,~isnan(S(i,:))), gg), (1:size(S,1))', 'UniformOutput', false));
CH = curv(D.SweepH); CC = curv(D.SweepCh);

%% ---- figure ----
figure
xSize = 18; ySize = 7.5;
set(gcf,'PaperUnits','centimeters','PaperPosition',[(21-xSize)/2 (30-ySize)/2 xSize ySize])
set(gcf,'Position',[40 40 xSize*50 ySize*50],'Color','w')

% a: metastability landscape
ax = axes('position', [.08 .18 .42 .60]); hold on
M_ = {CH, CC}; C_ = {cH, cC};
for s = 1:2
    m = mean(M_{s}, 1); sd = std(M_{s}, [], 1);
    fill([gg fliplr(gg)]*1e3, [m-sd fliplr(m+sd)], C_{s}, 'FaceAlpha', .2, 'EdgeColor', 'none');
    h(s) = plot(gg*1e3, m, '-', 'Color', C_{s}, 'LineWidth', 1.6);
end
yl = ylim; top = yl(2);
scatter(D.GstarH*1e3, top + .004 + .001*randn(nH,1), 6, gH, 'filled', 'MarkerFaceAlpha', .6, 'Clipping', 'off');
scatter(D.GstarCh*1e3, top + .010 + .001*randn(nC,1), 7, gC, 'filled', 'd', 'MarkerFaceAlpha', .8, 'Clipping', 'off');
xline(median(D.GstarH)*1e3, '--', 'Color', gH, 'LineWidth', .8); xline(median(D.GstarCh)*1e3, '--', 'Color', gC, 'LineWidth', .9);   % median G*
ylim(yl); xlabel('global coupling G (\times10^{-3})'); ylabel('metastability (model)')
legend(h, {sprintf('humans (n = %d)', nH), sprintf('chimpanzees (n = %d)', nC)}, 'Location', 'southeast', 'Box', 'off', 'FontSize', 6.5)
text(0, 1.13, 'individual G*', 'units','normalized', 'FontSize', 6, 'Color', [.4 .4 .4])
set(ax, 'TickDir','out', 'FontSize', 7, 'LineWidth', .7, 'Box','off')
text(-.13, 1.2, 'a', 'units','normalized','fontsize',12,'fontweight','bold')

% b, c: G* and metastability at G*
vals = {D.GstarH*1e3, D.GstarCh*1e3, 'working point G* (\times10^{-3})', 2, {gH, gC}; D.RH.Metasim, D.RCh.Metasim, 'metastability at G*', 1, {mH, mC}};
for i = 1:2
    ax = axes('position', [.61 + (i-1)*.215, .18, .15, .58]); hold on
    raincloud({vals{i,1}(:), vals{i,2}(:)}, vals{i,5}, {'o', 'd'});
    k = vals{i,4}; r = rb(k,iV); pp = p(k,iV);
    if pp < 0.05, st = '*'; fw = 'bold'; else, st = ''; fw = 'normal'; end
    pf = @(x) ternary(x < 0.001, 'p < 0.001', sprintf('p = %.3f', x));
    title({sprintf('V2: r_{rb} = %+.2f%s, %s', r, st, pf(pp)), sprintf('V0: r_{rb} = %+.2f, %s', rb(k,1), pf(p(k,1))), ...
           sprintf('V1: r_{rb} = %+.2f, %s', rb(k,2), pf(p(k,2)))}, 'FontSize', 6, 'FontWeight', fw)
    ylabel(vals{i,3}, 'FontSize', 7); set(ax, 'TickDir','out', 'FontSize', 6.5, 'LineWidth', .7, 'Box','off')
    text(-.8, 1.3, char('a' + i), 'units','normalized','fontsize',12,'fontweight','bold')
end

print('-dpng','-r300','Fig6_alt_model_level.png'); print('-dsvg','Fig6_alt_model_level.svg');
save('results_Fig6_alt.mat', 'lab', 'rb', 'p', 'gg', 'CH', 'CC');

%% ===================== local functions =====================
function raincloud(vals, cols, marks)
% raincloud per class: box (left), jittered points (middle), half violin (right)
    for x0 = 1:numel(vals)
        v = vals{x0}(:); c = cols{x0};
        q = prctile(v, [25 50 75]); iq = q(3) - q(1);
        lo = min(v(v >= q(1) - 1.5*iq)); hi = max(v(v <= q(3) + 1.5*iq));
        bx = x0 - .28; w = .12;
        plot([bx bx], [lo q(1)], 'k-', 'LineWidth', .7); plot([bx bx], [q(3) hi], 'k-', 'LineWidth', .7);
        rectangle('Position', [bx - w/2, q(1), w, q(3) - q(1)], 'FaceColor', c, 'EdgeColor', 'k', 'LineWidth', .7);
        plot([bx - w/2, bx + w/2], [q(2) q(2)], 'k-', 'LineWidth', .9);
        scatter(x0 - .08 + .12*rand(size(v)), v, 6, c, 'filled', 'Marker', marks{x0}, 'MarkerFaceAlpha', .8);
        yy = linspace(min(v), max(v), 200); dd = ksdensity(v, yy); dd = dd / max(dd) * .32;
        fill([x0 + .1 + dd, x0 + .1*ones(1,200)], [yy, fliplr(yy)], c, 'EdgeColor', 'none', 'FaceAlpha', .85);
    end
    xlim([.55 2.5]); set(gca, 'XTick', 1:numel(vals), 'XTickLabel', {'Hum','Chimp'});
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function q = fdr_bh(p)
    p = p(:); m = numel(p); [ps_, o] = sort(p);
    qs_ = ps_ .* m ./ (1:m)'; qs_ = flipud(cummin(flipud(qs_)));
    q = nan(m,1); q(o) = min(qs_, 1);
end

function [rb, p] = rbiserial(h, c)
    [p, ~, st] = ranksum(h, c);
    rb = 2*(st.ranksum - numel(h)*(numel(h)+1)/2)/(numel(h)*numel(c)) - 1;
end
