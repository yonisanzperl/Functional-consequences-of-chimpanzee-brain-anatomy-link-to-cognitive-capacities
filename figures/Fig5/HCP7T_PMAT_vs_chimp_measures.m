% HCP7T_PMAT_vs_chimp_measures.m
% Do the dynamical measures that relate to cognition in chimpanzees (Fig 2 / Fig 3)
% behave the same way in humans?
%
% Humans:  HCP 7T resting state (REST1 PA), empirical FC, Schaefer 100,
%          band 0.02-0.03 Hz (the band of the chimpanzee model).
% Score:   PMAT24_A_CR (Penn matrices: non-verbal fluid reasoning), the human test
%          closest in content to the chimpanzee PCTB physical-cognition battery.
% Control: age and sex (partial correlations).
%
% For each measure:
%   r, rho          Pearson and Spearman correlation with PMAT24 (no covariates)
%   partial r       controlling for age + sex, with 95% CI (Fisher z, n-3-k) and p
%                   (t test, df = n-2-k)
%   q               Benjamini-Hochberg FDR across the 5 pre-specified Fig 2 measures
%                   (meanFC, Ignition, Edge turbulence, Entropy eFCD, FC modularity)
%   held-out r      PMAT24 and the measure are both residualised on age + sex; each
%                   subject's score is then predicted from a line fitted to the other
%                   n-1 subjects (exact leave-one-out via the hat matrix), and
%                   correlated with the observed score. p = permutation test (5000
%                   shuffles of the score, one-sided: held-out r larger than chance).
%
% Figures:
%   1  PMAT24 vs each measure (partial-residual scatter, both adjusted for age + sex)
%   2  Same measures in both species: plain Pearson r with 95% CI in both (no covariates,
%      since chimpanzee age/sex are not available), with a Fisher z test of the species difference
%   3  Chimpanzee r vs human r (plain Pearson in both)
%
% Needs: HCP7T_measures_for_chimp_comparison.mat (HCP7T_measures_for_chimp_comparison.m),
%        HCP7T_sex.mat. Statistics Toolbox. MATLAB R2018b+.
% Run from Monkey_MVH/HumanBehaviour_7T/.
%
% Chimpanzee values: Pearson r with PCTB cognition, n = 43, noise-matched simulations
% (functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_REAL_0926_001.mat; Fig 2 and Fig 3 of v5).
%
% Caveat: PMAT24 was chosen after looking at all 10 NIH/Penn tests (see the
% supplementary loop at the end); HCP family structure and head motion are not modelled.

clear all; close all
rng(2026);

load('HCP7T_measures_for_chimp_comparison.mat');   % ids, names, bands, M, tests, gender, age, ok
load('HCP7T_sex.mat');                             % sex (1 = M, 0 = F)
b = 1;                                             % band 1 = 0.02-0.03 Hz (chimp band)
fprintf('Band: %s\n', bandNames{b});

%% ---- data ----
if iscell(age)
    ageN = cellfun(@(c) str2double(regexp(char(string(c)),'\d+','match','once')), age);
else
    ageN = double(age);
end
ageN = ageN(:); sex = double(sex(:));
pmat = tests(:,1);                                 % column 1 = PMAT24_A_CR
X    = M(:,:,b);
keep = logical(ok(:)) & ~isnan(pmat) & ~isnan(ageN) & ~isnan(sex) & all(~isnan(X),2);
y = pmat(keep); X = X(keep,:); ageN = ageN(keep); sex = sex(keep);
n = numel(y); covs = [ageN sex]; k = size(covs,2);
fprintf('n = %d (%d males); PMAT24 range %d-%d; r(PMAT,sex) = %+.2f, r(PMAT,age) = %+.2f\n', ...
        n, sum(sex), min(y), max(y), corr(y,sex), corr(y,ageN));

labels = {'Mean FC','Ignition','Edge turbulence','Entropy eFCD','FC modularity', ...
          'FC clustering','FC betweenness','Metastability','Synchrony'};
nM = numel(names);
pre = 1:5;                                         % pre-specified Fig 2 measures
chimpR = [-0.42 -0.43 -0.41 -0.14 +0.38 -0.41 +0.32 -0.17 -0.43];   % Sep 2026 rerun (sweep frequencies), n = 43; clustering, betweenness = dense model FC (Fig 3)
chimpN = 43;                                       % 95% CI of the chimp r (Fisher z, no covariates)
chimpLo = tanh(atanh(chimpR) - 1.96/sqrt(chimpN-3));
chimpHi = tanh(atanh(chimpR) + 1.96/sqrt(chimpN-3));

%% ---- statistics ----
Z  = [ones(n,1) covs];
res = @(v) v - Z*(Z\v);
ry = res(y);
nPerm = 5000;
perms_ = zeros(n, nPerm);
for p = 1:nPerm, perms_(:,p) = randperm(n)'; end   % same shuffles for all measures

[r0,p0,rs,ps,pr,pp,ciLo,ciHi,ho,pho,r0Lo,r0Hi] = deal(nan(nM,1));
RX = nan(n,nM);
for j = 1:nM
    x = X(:,j);
    [r0(j),p0(j)] = corr(x, y);
    [rs(j),ps(j)] = corr(x, y, 'Type','Spearman');
    [pr(j),pp(j)] = partialcorr(x, y, covs);
    zf = atanh(pr(j)); se = 1/sqrt(n-3-k);
    ciLo(j) = tanh(zf-1.96*se); ciHi(j) = tanh(zf+1.96*se);
    r0Lo(j) = tanh(atanh(r0(j))-1.96/sqrt(n-3)); r0Hi(j) = tanh(atanh(r0(j))+1.96/sqrt(n-3));
    RX(:,j) = res(x);
    ho(j) = loo_r(RX(:,j), ry);
    nul = zeros(nPerm,1);
    for p = 1:nPerm, nul(p) = loo_r(RX(:,j), ry(perms_(:,p))); end
    pho(j) = (1 + sum(nul >= ho(j))) / (nPerm + 1);
end
q = nan(nM,1); q(pre) = fdr_bh(pp(pre));
% Cross-species comparison uses PLAIN Pearson r in both species (no covariates in either),
% because age and sex are not available for the chimpanzees.
% Do the two species differ? independent-samples Fisher z test on chimp r vs human r
zdiff = (atanh(chimpR(:)) - atanh(r0)) ./ sqrt(1/(chimpN-3) + 1/(n-3));
pdiff = 2*normcdf(-abs(zdiff));

fprintf('\n%-16s %22s | %6s %6s | %6s | %8s %16s %7s %6s | %14s | %s\n', 'measure','chimp r [95% CI]', ...
        'r','p','rho','partial','95% CI','p','q','held-out r (p)','chimp vs human p');
for j = 1:nM
    if isnan(q(j)), qs = '     -'; else, qs = sprintf('%6.3f', q(j)); end
    fprintf('%-16s %+6.2f [%+.2f, %+.2f] | %+6.2f %6.3f | %+6.2f | %+8.2f [%+.2f, %+.2f] %7.4f %s | %+.2f (%.4f) | %.3f\n', ...
        labels{j}, chimpR(j), chimpLo(j), chimpHi(j), r0(j), p0(j), rs(j), pr(j), ciLo(j), ciHi(j), pp(j), qs, ho(j), pho(j), pdiff(j));
end
has = ~isnan(chimpR);
[rx, px] = corr(chimpR(has)', r0(has));
fprintf('\nChimp r vs human r (both plain Pearson) across the %d shared measures: r = %+.2f (p = %.3f; measures are not independent, descriptive only)\n', ...
        sum(has), rx, px);

%% ---- colours ----
cSig = [0.17 0.37 0.54];   % significant (q < 0.05 for Fig 2 measures, p < 0.05 otherwise)
cNS  = [0.55 0.58 0.63];
cCh  = [0.78 0.44 0.16];   % chimpanzee
isSig = pp < 0.05; isSig(pre) = q(pre) < 0.05;

%% ---- Figure 1: partial-residual scatter ----
show = [1 2 3 4 5 6 7 9];
figure('Color','w','Position',[40 40 1500 760]);
for s = 1:numel(show)
    j = show(s); x = RX(:,j);
    subplot(2,4,s); hold on
    if isSig(j), col = cSig; else, col = cNS; end
    mdl = fitlm(x, ry);
    xx = linspace(min(x), max(x), 60)';
    [yp, yci] = predict(mdl, xx);
    fill([xx; flipud(xx)], [yci(:,1); flipud(yci(:,2))], col, 'FaceAlpha', 0.15, 'EdgeColor','none');
    scatter(x, ry, 18, col, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor','w');
    plot(xx, yp, 'Color', col, 'LineWidth', 2);
    if ~isnan(q(j)), st = sprintf('q = %.3f', q(j)); else, st = sprintf('p = %.3f', pp(j)); end
    if ~isnan(chimpR(j)), ct = sprintf('chimp r = %+.2f', chimpR(j)); else, ct = 'no chimp value'; end
    if isSig(j), fw = 'bold'; else, fw = 'normal'; end
    title({labels{j}, sprintf('partial r = %+.2f, %s | %s', pr(j), st, ct)}, 'FontSize', 9, 'FontWeight', fw);
    xlabel([labels{j} ' (adj. age, sex)']);
    if s == 1 || s == 5, ylabel('PMAT24 (adj. age, sex)'); end
    box off; set(gca,'TickDir','out','FontSize',9);
end
sgtitle(sprintf(['Humans, HCP 7T (n = %d), %s: PMAT24 vs dynamical measures. ' ...
    'Colour: q < 0.05 (FDR over the 5 Fig 2 measures) or p < 0.05 (others)'], n, bandNames{b}), 'FontSize', 10);

%% ---- Figure 2: forest plot, both species ----
figure('Color','w','Position',[80 80 760 520]); hold on
yy = nM:-1:1;
off = 0.15;                                        % humans slightly above, chimps below
for j = 1:nM
    plot([r0Lo(j) r0Hi(j)], [yy(j) yy(j)]+off, '-', 'Color', cSig, 'LineWidth', 2);
    h1 = plot(r0(j), yy(j)+off, 'o', 'MarkerFaceColor', cSig, 'MarkerEdgeColor','w', 'MarkerSize', 8);
    if ~isnan(chimpR(j))
        plot([chimpLo(j) chimpHi(j)], [yy(j) yy(j)]-off, '-', 'Color', cCh, 'LineWidth', 2);
        h2 = plot(chimpR(j), yy(j)-off, 'd', 'MarkerFaceColor','w', 'MarkerEdgeColor', cCh, 'LineWidth', 2, 'MarkerSize', 8);
        if pdiff(j) < 0.05, text(0.66, yy(j), '\ne', 'Color', cNS, 'FontSize', 11, 'Clipping','off'); end
    end
end
xline(0, 'k:');
set(gca, 'YTick', fliplr(yy), 'YTickLabel', fliplr(labels), 'TickDir','out', 'FontSize', 10);
xlim([-0.7 0.7]); ylim([0.4 nM+0.6]); box off
xlabel('correlation with cognition (95% CI)');
legend([h1 h2], {'Humans: r with PMAT24 (HCP 7T, n = 181), 95% CI', ...
                 'Chimpanzees: r with PCTB (model, n = 43), 95% CI'}, 'Location','southeast', 'Box','off');
title({'Same measures, both species: plain Pearson r, no covariates', '(\ne: species differ, p < 0.05)'});

%% ---- Figure 3: chimp r vs human r ----
figure('Color','w','Position',[120 120 560 500]); hold on
patch([-0.7 0 0 -0.7], [-0.45 -0.45 0 0], cSig, 'FaceAlpha', 0.06, 'EdgeColor','none');
patch([0 0.7 0.7 0],   [0 0 0.45 0.45],   cSig, 'FaceAlpha', 0.06, 'EdgeColor','none');
xline(0,'k:'); yline(0,'k:');
idx = find(has);
dy = [-0.025 -0.03 0.03 0 0 0.025 0];              % label nudges to avoid overlaps
for t = 1:numel(idx)
    j = idx(t);
    plot([chimpR(j) chimpR(j)], [r0Lo(j) r0Hi(j)], '-', 'Color', cSig, 'LineWidth', 1);     % human CI
    plot([chimpLo(j) chimpHi(j)], [r0(j) r0(j)], '-', 'Color', cCh, 'LineWidth', 1);       % chimp CI
    plot(chimpR(j), r0(j), 'o', 'MarkerFaceColor', cSig, 'MarkerEdgeColor','w', 'MarkerSize', 8);
    text(chimpR(j)+0.02, r0(j)+dy(t), labels{j}, 'FontSize', 9);
end
xlim([-0.7 0.7]); ylim([-0.45 0.45]); box off; set(gca,'TickDir','out');
xlabel('chimpanzee r (model, PCTB), 95% CI'); ylabel('human r (empirical, PMAT24), 95% CI');
title(sprintf('Chimp r vs human r, plain Pearson (%d measures): r = %+.2f', sum(has), rx));

%% ---- supplementary: all 10 tests, partial r (age, sex), for transparency ----
testNames = {'PMAT24','ProcSpeed','CardSort','ListSort','Flanker','VSPLOT','IWRD','PicSeq','ReadEng','PicVocab'};
Tk = tests(keep,:);
Rall = nan(nM, size(Tk,2));
for c = 1:size(Tk,2)
    for j = 1:nM, Rall(j,c) = partialcorr(X(:,j), Tk(:,c), covs, 'Rows','complete'); end
end
fprintf('\nSupplement: partial r (age, sex) of each measure with each test\n%-16s', '');
fprintf('%9s', testNames{:}); fprintf('\n');
for j = 1:nM, fprintf('%-16s', labels{j}); fprintf('%+9.2f', Rall(j,:)); fprintf('\n'); end

save('results_HCP7T_PMAT_vs_chimp.mat', 'labels','names','n','chimpR','r0','p0','rs','ps', ...
     'pr','pp','ciLo','ciHi','q','ho','pho','chimpLo','chimpHi','pdiff','r0Lo','r0Hi','rx','px','Rall','testNames');

%% ---- local functions ----
function r = loo_r(x, y)
    % exact leave-one-out prediction of y from x (with intercept), via the hat matrix
    A = [ones(numel(x),1) x];
    H = A * pinv(A'*A) * A';
    e = y - H*y;
    yhat = y - e ./ (1 - diag(H));
    r = corr(y, yhat);
end

function q = fdr_bh(p)
    p = p(:); m = numel(p);
    [ps_, order] = sort(p);
    qs_ = ps_ .* m ./ (1:m)';
    qs_ = flipud(cummin(flipud(qs_)));
    q = nan(m,1); q(order) = min(qs_, 1);
end
