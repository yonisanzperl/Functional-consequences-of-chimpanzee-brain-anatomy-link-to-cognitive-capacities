% Fig5a_working_point.m
% Fig 5 (human validation, part 1): working point of the human model (HCP 7T group SC,
% Schaefer 100) run with EXACTLY the chimpanzee settings (0.02-0.03 Hz, noise 0.04, C/mean,
% heterogeneous omega), in the style and colours of Fig 2.
%   a  metastability vs G (mean +/- sd over the 20 runs); grey band = empirical metastability
%      (mean +/- sd across subjects, same band); dashed line = G* (max of the smoothed curve,
%      as for the chimpanzees)
%   b  fit to the empirical group FC vs G: SSIM (left axis) and distance (right axis);
%      triangles = optimum of each criterion
% (26 Sep 2026: the FC-correlation curve and the FC-matrix panel were removed; the FC correlation
%  is still printed in the Command Window. Previous version: Fig5a_working_point_withFCcorr_OLD.m)
% Input, whichever exists (the first is preferred):
%   results_model_Humans_HCP7T_asCHIMP_bandFC.mat  <- model_Humans_HCP7T_asChimp.m (fixed version):
%        FC fit against the empirical FC in the SAME band, with FC matrices (panel c)
%   results_model_Humans_HCP7T_asCHIMP.mat         <- run of 24 Sep: FC fit against the 0.008-0.08 Hz
%        empirical FC (labelled as such); empirical metastability taken from
%        HCP7T_measures_for_chimp_comparison.mat (0.02-0.03 Hz); panel c left empty
% Needs cmocean.m (this folder). Saves Fig5a_working_point.png (300 dpi) and .svg.

clear all; close all

%% ---- data ----
if exist('results_model_Humans_HCP7T_asCHIMP_bandFC.mat', 'file')
    R = load('results_model_Humans_HCP7T_asCHIMP_bandFC.mat');
    metaEmp = R.Metasim_sub22(:);
    fitLabel = 'empirical FC, 0.02-0.03 Hz';
    haveFC = true;
else
    R = load('results_model_Humans_HCP7T_asCHIMP.mat');
    E = load('HCP7T_measures_for_chimp_comparison.mat');          % M(:,8,1) = metastability, 0.02-0.03 Hz
    metaEmp = E.M(:,8,1);
    fitLabel = 'empirical FC, 0.008-0.08 Hz';
    haveFC = false;
    warning('Using the 24 Sep sweep: FC fit is against the 0.008-0.08 Hz empirical FC. Rerun model_Humans_HCP7T_asChimp.m for the band-matched version.');
end
G = R.Grange(:)';
metaEmp = metaEmp(~isnan(metaEmp));
mM = mean(R.mMetansub);  sM = std(R.mMetansub);
mF = mean(R.fitFCnusb);  sF = std(R.fitFCnusb);
mS = mean(R.ssimFCnusb); sS = std(R.ssimFCnusb);
mD = mean(R.distFCnusb); sD = std(R.distFCnusb);

mMsm = smooth(mM);  [~, iG] = max(mMsm);  Gstar = G(iG);          % as for the chimpanzees
Gsub = zeros(size(R.mMetansub,1),1);                             % G* of each run (spread)
for r = 1:size(R.mMetansub,1), [~, k] = max(smooth(R.mMetansub(r,:))); Gsub(r) = G(k); end
[~, iF] = max(mF); [~, iS] = max(mS); [~, iD] = min(mD);

fprintf('G* = %.4f (runs: %.4f +/- %.4f, range %.4f-%.4f)\n', Gstar, mean(Gsub), std(Gsub), min(Gsub), max(Gsub));
fprintf('metastability: model at G* %.3f +/- %.3f | empirical %.3f +/- %.3f (n = %d)\n', ...
    mM(iG), sM(iG), mean(metaEmp), std(metaEmp), numel(metaEmp));
fprintf('FC fit (%s): r max at G = %.4f (%.3f); SSIM max at G = %.4f (%.3f); distance min at G = %.4f (%.4f)\n', ...
    fitLabel, G(iF), mF(iF), G(iS), mS(iS), G(iD), mD(iD));
fprintf('at G*: r = %.3f, SSIM = %.3f, distance = %.4f\n', mF(iG), mS(iG), mD(iG));

Col = cmocean('thermal', 44);
cMeta = Col(16,:); cR = Col(4,:); cS = Col(22,:); cD = Col(34,:);
gx = G*1000;                                                     % plot G in units of 10^-3

%% ---- figure ----
figure
xSize = 13; ySize = 6.5;
set(gcf,'PaperUnits','centimeters','PaperPosition',[(21-xSize)/2 (30-ySize)/2 xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

% a: metastability
ax = axes('position', [.10 .19 .34 .70]); hold on
e1 = mean(metaEmp) + [-1 1]*std(metaEmp);
patch([gx(1) gx(end) gx(end) gx(1)], [e1(1) e1(1) e1(2) e1(2)], [.8 .8 .8], 'EdgeColor','none', 'FaceAlpha', .6);
plot([gx(1) gx(end)], mean(metaEmp)*[1 1], '-', 'Color', [.5 .5 .5], 'LineWidth', 1);
fill([gx fliplr(gx)], [mM+sM fliplr(mM-sM)], cMeta, 'EdgeColor','none', 'FaceAlpha', .3);
plot(gx, mM, '-', 'Color', cMeta, 'LineWidth', 1.8);
xline(Gstar*1000, 'k--', 'LineWidth', 1);
text(Gstar*1000, max(mM+sM)*1.08, sprintf(' G* = %.1f', Gstar*1000), 'FontSize', 7.5);
text(.40, .20, {sprintf('model at G*: %.2f', mM(iG)), sprintf('empirical: %.2f \\pm %.2f', mean(metaEmp), std(metaEmp))}, ...
    'units','normalized','fontsize',7)
xlabel('G (\times10^{-3})'); ylabel('Metastability')
xlim([gx(1) gx(end)]); ylim([0 max(mM+sM)*1.15])
set(ax,'TickDir','out','FontSize',8,'LineWidth',.8,'Box','off')
text(-.24, 1.07, 'a', 'units','normalized','fontsize',12,'fontweight','bold')
axis('square')
% b: FC fit criteria
ax = axes('position', [.58 .19 .30 .70]); hold on
yyaxis left
fill([gx fliplr(gx)], [mS+sS fliplr(mS-sS)], cS, 'EdgeColor','none', 'FaceAlpha', .25);
h2 = plot(gx, mS, '-', 'Color', cS, 'LineWidth', 1.6);
plot(gx(iS), mS(iS), 'v', 'MarkerFaceColor', cS, 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
ylabel('FC SSIM'); ylim([0 max(mS+sS)*1.15])
yyaxis right
fill([gx fliplr(gx)], [mD+sD fliplr(mD-sD)], cD, 'EdgeColor','none', 'FaceAlpha', .2);
h3 = plot(gx, mD, '-', 'Color', cD, 'LineWidth', 1.6);
plot(gx(iD), mD(iD), '^', 'MarkerFaceColor', cD, 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
ylabel('FC distance'); ylim([0 max(mD+sD)*1.05])
ax.YAxis(1).Color = [.2 .2 .2]; ax.YAxis(2).Color = [.2 .2 .2];
xline(Gstar*1000, 'k--', 'LineWidth', 1);
xlabel('G (\times10^{-3})'); xlim([gx(1) gx(end)])
legend([h2 h3], {'SSIM','distance'}, 'Location','southeast', 'Box','off', 'FontSize', 6.5)
title(['fit to ' fitLabel], 'FontWeight','normal', 'FontSize', 7.5)
set(ax,'TickDir','out','FontSize',8,'LineWidth',.8,'Box','off')
text(-.26, 1.07, 'b', 'units','normalized','fontsize',12,'fontweight','bold')
axis('square')

print('-dpng','-r300','Fig5a_working_point.png'); print('-dsvg','Fig5a_working_point.svg');
