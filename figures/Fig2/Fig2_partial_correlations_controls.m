% Fig2_partial_correlations_controls.m
% Controls for Figure 2: do the model-inferred functional measures relate to
% cognition beyond SC density, SC strength and the optimal coupling G*?
% Uses the same files and the same subject selection as Fig2_CHIMP_MVH.m.
% p-values are FDR-corrected (Benjamini-Hochberg) across the 5 measures in each analysis.

clear all; close all

load('meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat')           % Metasim_sub2 (59 x G x reps)
load('CHIMP_CWAS_DATA_FORDECO.mat')                           % CJ_CHIMP, AapCog_normalized
load('functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_REAL_0926_001.mat')  % RH

num_mammals = 59;
Grange = [0:0.1:10]/1000;

% ---- optimal working point G* (identical to Fig2_CHIMP_MVH.m) ----
for ii = 1:num_mammals
    Metasim_sub(ii,:) = smooth(squeeze(mean(Metasim_sub2(ii,:,:),3)));
end
[Meta, Gs2] = max(Metasim_sub,[],2);
Gstar = Grange(Gs2)';                 % G* value (linear in the index, so same result)

% ---- subject selection (identical to Fig2_CHIMP_MVH.m) ----
indx = find(Meta>0);
indx = setdiff(indx,[23,27]);         % SC outliers (to be justified in the paper)
indx = indx(~isnan(AapCog_normalized(indx)));   % drop the animal with no cognition score
fprintf('n = %d\n', numel(indx));

% ---- SC summary covariates ----
dens  = nan(num_mammals,1);
stren = nan(num_mammals,1);
for s = 1:num_mammals
    C = squeeze(CJ_CHIMP(:,:,s));
    C(isnan(C)) = 0;
    dens(s)  = mean(C(:) > 0);        % SC density: fraction of non-zero connections
    stren(s) = mean(C(:));            % SC strength: mean connection weight (zeros included)
end

% ---- functional measures used in Fig 2 ----
X = [RH.Integration(:), RH.Ignition(:), mean(RH.EdgeMetaNode,2), RH.EntroEdgesMeta(:), RH.Qs(:)];
names = {'meanFC','Ignition','EdgeTurbulence','EntropyeFCD','Modularity'};

y  = AapCog_normalized(indx);
Z1 = dens(indx);                                  % control 1: density only
Z3 = [dens(indx), stren(indx), Gstar(indx)];      % control 2: density + strength + G*

nM = numel(names);
[r0,p0,r1,p1,r3,p3] = deal(nan(nM,1));
for k = 1:nM
    x = X(indx,k);
    [r0(k),p0(k)] = corr(x, y);
    [r1(k),p1(k)] = partialcorr(x, y, Z1);   % partialcorr uses the correct df = n-2-k
    [r3(k),p3(k)] = partialcorr(x, y, Z3);
end

% ---- FDR correction (Benjamini-Hochberg) across the 5 measures, per analysis ----
q0 = fdr_bh(p0);  q1 = fdr_bh(p1);  q3 = fdr_bh(p3);

fprintf('\n%-15s | %6s %6s %6s | %6s %6s %6s | %8s %6s %6s\n', ...
    'measure','r','p','q','r|dens','p','q','r|d,s,G*','p','q');
for k = 1:nM
    fprintf('%-15s | %+6.2f %6.3f %6.3f | %+6.2f %6.3f %6.3f | %+8.2f %6.3f %6.3f\n', ...
        names{k}, r0(k),p0(k),q0(k), r1(k),p1(k),q1(k), r3(k),p3(k),q3(k));
end
fprintf('q = Benjamini-Hochberg FDR-adjusted p-value across the %d measures.\n', nM);

% How strongly each measure tracks SC density (the confound we control for)
fprintf('\nCorrelation of each measure with SC density:\n');
for k = 1:numel(names)
    fprintf('%-15s r = %+.2f\n', names{k}, corr(X(indx,k), dens(indx)));
end

% ---- local function: Benjamini-Hochberg adjusted p-values (no toolbox needed) ----
function q = fdr_bh(p)
    p = p(:);  n = numel(p);
    [ps, order] = sort(p);
    qs = ps .* n ./ (1:n)';
    qs = flipud(cummin(flipud(qs)));   % enforce monotonicity
    q = nan(n,1);
    q(order) = min(qs, 1);
end
