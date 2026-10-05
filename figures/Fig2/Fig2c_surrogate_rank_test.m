% Fig2c_surrogate_rank_test.m
% Is the brain-cognition relationship in Fig 2 specific to the real chimpanzee
% connectomes? Compare it with connectome surrogates that keep each animal's
% degree sequence (randmio_und) or degree + strength (null_model_und_sign).
%
% Main test: rank (permutation-style) test
%   p_rank = (1 + #{surrogates with |r_null| >= |r_real|}) / (nNull + 1)
%   The smallest possible p is 1/(nNull+1): 0.048 with 20 surrogates, 0.0099 with 100.
% Descriptive only (not a p-value):
%   z      = (r_real - mean(r_null)) / std(r_null)
%   #sig   = number of surrogates whose own correlation with cognition has p < 0.05
%   r(real,null) = across animals, correlation between the real measure and the
%                  surrogate-average measure (how much of it degree/strength explains)
%
% The number of surrogates is detected from the files, so the script works
% unchanged once more surrogates (>= 100) are generated.
%
% Before submission: generate the surrogates with G* re-found for each surrogate
% (max metastability of the surrogate itself). The current files were simulated
% at the G* of the real connectome. Point SURR_DIR / patterns at the new files.

clear all; close all

%% ---------------- settings ----------------
SURR_DIR = fullfile(pwd,'surrogates_sweepfreq');   % Fig2_surrogate_connectomes.m (Sep 2026: same settings as the real run, G* of the real connectome)
families = { ...
  'Degree-preserving null (randmio\_und)',                 'functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_Degree_*.mat'; ...
  'Degree + strength-preserving null (null\_model\_und\_sign)', 'functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_Strength__*.mat'};
names  = {'Synchrony','Ignition','Edge turbulence','Entropy eFCD'};   % metastability left out (not related to cognition), 26 Sep 2026
getX   = @(RH) [RH.Synchro_sub2(:), RH.Ignition(:), mean(RH.EdgeMetaNode,2), ...
                RH.EntroEdgesMeta(:)];

%% ---------------- data and subject selection (as in Fig2_CHIMP_MVH.m) ----------------
load('meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat')           % Metasim_sub2
load('CHIMP_CWAS_DATA_FORDECO.mat')                           % AapCog_normalized
REAL = load('functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_REAL_0926_001.mat');

for ii = 1:59
    Metasim_sub(ii,:) = smooth(squeeze(mean(Metasim_sub2(ii,:,:),3)));
end
Meta = max(Metasim_sub,[],2);
indx = find(Meta>0);
indx = setdiff(indx,[23,27]);                         % SC outliers (justify in paper)
indx = indx(~isnan(AapCog_normalized(indx)));         % animal without cognition score
y    = AapCog_normalized(indx);
n    = numel(indx);
nM   = numel(names);

Xreal = getX(REAL.RH);
Xreal = Xreal(indx,:);

% |r| needed for p < 0.05 (two-sided) with this n, used to shade the plot
tcrit = tinv(0.975, n-2);
rcrit = tcrit / sqrt(n - 2 + tcrit^2);
fprintf('n = %d animals; |r| > %.2f is p < 0.05\n', n, rcrit);

%% ---------------- rank test per null family ----------------
results = struct();
figure('Color','w','Position',[50 50 1250 470])
for f = 1:size(families,1)

    % find surrogate files: only numbered ones (..._001.mat, ..._002.mat, ...)
    d = dir(fullfile(SURR_DIR, families{f,2}));
    keep = ~cellfun(@isempty, regexp({d.name}, '_\d{3}\.mat$', 'once'));
    d = d(keep);
    nNull = numel(d);
    if nNull == 0, error('No surrogate files found for: %s', families{f,2}); end

    Xnull = nan(n, nM, nNull);
    for s = 1:nNull
        S = load(fullfile(d(s).folder, d(s).name));
        tmp = getX(S.RH);
        Xnull(:,:,s) = tmp(indx,:);
    end

    fam = strrep(families{f,1},'\','');
    fprintf('\n===== %s: %d surrogates (min possible p = %.3f)\n', fam, nNull, 1/(nNull+1));
    fprintf('%-16s %7s %7s %20s %8s %7s %6s %10s\n', ...
        'measure','r_real','p_real','null r mean [range]','p_rank','z','#sig','r(real,null)');

    for k = 1:nM
        [rr, pr] = corr(Xreal(:,k), y);
        [rn, pn] = corr(squeeze(Xnull(:,k,:)), y);    % nNull x 1
        p_rank = (1 + sum(abs(rn) >= abs(rr))) / (nNull + 1);
        z      = (rr - mean(rn)) / std(rn);
        nsig   = sum(pn < 0.05);
        simRN  = corr(Xreal(:,k), mean(Xnull(:,k,:),3));

        fprintf('%-16s %+7.2f %7.3f %+7.2f [%+.2f,%+.2f] %8.3f %+7.1f %3d/%-3d %+8.2f\n', ...
            names{k}, rr, pr, mean(rn), min(rn), max(rn), p_rank, z, nsig, nNull, simRN);

        results(f,k).family  = fam;
        results(f,k).measure = names{k};
        results(f,k).r_real  = rr;
        results(f,k).p_real  = pr;
        results(f,k).r_null  = rn;
        results(f,k).p_rank  = p_rank;
        results(f,k).z       = z;
        results(f,k).n_sig_surrogates = nsig;
        results(f,k).r_real_vs_nullmean = simRN;

        % ---- plot: surrogate r (grey dots) vs real r (blue line) ----
        ax = subplot(size(families,1), nM, (f-1)*nM + k); hold on
        patch([-0.6 -rcrit -rcrit -0.6], [0 0 1 1], [0.95 0.95 0.96], 'EdgeColor','none');
        patch([ rcrit  0.6  0.6  rcrit], [0 0 1 1], [0.95 0.95 0.96], 'EdgeColor','none');
        xline(0, 'Color', [0.8 0.8 0.8]);
        jit = 0.5 + 0.35*(rand(nNull,1) - 0.5);
        scatter(rn, jit, 36, [0.61 0.64 0.69], 'filled', 'MarkerEdgeColor','w');
        xline(rr, 'Color', [0.17 0.37 0.54], 'LineWidth', 2.2);
        if rr < 0, ha = 'left'; dx = 0.03; else, ha = 'right'; dx = -0.03; end
        text(rr+dx, 0.92, sprintf('r = %+.2f', rr), 'HorizontalAlignment', ha, 'FontWeight','bold', 'FontSize', 9);
        text(rr+dx, 0.82, sprintf('p_{rank} = %.3f', p_rank), 'HorizontalAlignment', ha, 'FontSize', 8, 'Color', [0.42 0.45 0.5]);
        xlim([-0.6 0.6]); ylim([0 1]);
        set(ax, 'YTick', [], 'XTick', [-0.4 0 0.4], 'YColor', 'none', 'TickDir', 'out', 'Box', 'off');
        if f == 1, title(names{k}); end
        if k == 1, text(-0.75, 0.5, strtok(fam,'('), 'Rotation', 90, 'HorizontalAlignment','center', 'FontSize', 9); end
        if f == size(families,1), xlabel('r with cognition (PCTB)'); end
    end
end
sgtitle(sprintf('Real connectomes (blue) vs surrogates (grey); shaded: p < 0.05 (|r| > %.2f), n = %d', rcrit, n), 'FontSize', 10);

% ---- FDR of the rank p across measures, within each null family ----
fprintf('\nFDR (Benjamini-Hochberg) of p_rank across the %d measures:\n', nM);
for f = 1:size(families,1)
    pr_ = [results(f,:).p_rank]; qr_ = fdr_bh(pr_);
    for k = 1:nM, results(f,k).q_rank = qr_(k); end
    fprintf('  %-45s q = %s\n', strrep(families{f,1},'\',''), sprintf('%.3f ', qr_));
end
save('results_Fig2c_surrogate_rank_test.mat', 'results', 'indx', 'n', 'rcrit');

function q = fdr_bh(p)
    p = p(:); m = numel(p); [ps_, o] = sort(p);
    qs_ = ps_ .* m ./ (1:m)'; qs_ = flipud(cummin(flipud(qs_)));
    q = nan(m,1); q(o) = min(qs_, 1);
end
