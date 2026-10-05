% model_Humans_HCP7T_asChimp.m
% Human working point (group SC, HCP 7T) with EXACTLY the chimpanzee model settings:
% 0.02-0.03 Hz, noise 0.04, C/mean normalisation, heterogeneous omega, Grange 0-0.02.
% Changes vs the version run on 24 Sep (kept as model_Humans_HCP7T_asChimp_ORIGINAL.m):
%   1) empirical target (FC and metastability) now in the SAME band as the model (0.02-0.03 Hz),
%      from HCP7T_empirical_FC_chimpband.m  (was 0.008-0.08 Hz from results_HCP_7T_empiricalmeasures.mat)
%   2) G* = maximum of the smoothed mean metastability curve (as for the chimps), no hard-coded values
%      (was vline(1.3) and FCsim_meanG(12,:,:))
%   3) the simulated FC at every G (mean over runs, FCsim_Gmean) and the empirical metastability are saved, so the figure
%      (Fig5a_working_point.m) can draw the FC matrices without rerunning
%   4) output: results_model_Humans_HCP7T_asCHIMP_bandFC.mat (the previous results file is not overwritten)
% The simulation itself is unchanged.
clear all
close all


load('SC_schaefer100_17Networks_32fold_groupconnectome_2mm_symm.mat')
E = load('HCP7T_empirical_FC_chimpband.mat');   % <- HCP7T_empirical_FC_chimpband.m (0.02-0.03 Hz)
FCemp = E.FCemp; Metasim_sub22 = E.Metaemp(:)';

N=100;



TR=1;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.02;                    % lowpass frequency of filter (Hz)
fhi = 0.03;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

dt=0.1*TR/2; % ms
sig=0.04;
dsig = sqrt(dt)*sig;
a=-0.02;
f_diff= 0.025*(rand(1,N)*0.1+1);
omega = repmat(2*pi*f_diff(1:N)',1,2); omega(:,1) = -omega(:,1);

%% G optimo
Grange=[0:0.1:20]/1000;

Tmax=800;

cut=10;;

Isubdiag = find(tril(ones(N),-1));
IsubdiagT = find(tril(ones(Tmax-2*cut+1),-1));

FCemp_global=squeeze(nanmean(FCemp,1));

C=SC;

for nsub=1:20
    nG=1;
    nsub
    for G=Grange
        C=C/mean(C(:));
        wC = G*C;
        
        %  C = randmio_und(C,1000); % rand preserving degree
        %C= null_model_und_sign(C); % rand preserving degree, strength and distribution

        for trials=1:10

            sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
            xs=zeros(Tmax,N);
            z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
            nn=0;
            % discard first 3000 time steps
            for t=0:dt:1000
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
            end
            % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
            for t=0:dt:((Tmax-1)*TR)
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
                if abs(mod(t,TR))<0.01
                    nn=nn+1;
                    xs(nn,:)=z(:,1)';
                end
            end

            %%%%
            BOLD=xs';
            clear signal_filt;
            for seed=1:N
                ts33(seed,:)=detrend(BOLD(seed,:)-nanmean(BOLD(seed,:)));
                signal_filt(seed,:) =filtfilt(bfilt,afilt,ts33(seed,:));

            end
            signal_filt=signal_filt(:,cut:end-cut);
            ts2=signal_filt;

            for seed=1:N
                Xanalytic = hilbert(ts2(seed,:)-mean(ts2(seed,:)));
                Amplitude(seed,:)=abs(Xanalytic);
                PhasesS(seed,:) = angle(Xanalytic);
            end
            FCsim2(trials,:,:)=corrcoef(ts2');
            gKOMsim=nansum(complex(cos(PhasesS),sin(PhasesS)))/N;
            enstrophy1sim=abs(gKOMsim);
            Metasim_sub22_sim(trials)=nanstd(enstrophy1sim(:)); %metastability
        end
        FCsim_m=squeeze(mean(FCsim2));
        mMeta(nG)=squeeze(mean(Metasim_sub22_sim));



        fitFC(nG)=corr(FCsim_m(Isubdiag),FCemp_global(Isubdiag))
        distFC(nG)=nanmean(nanmean((FCsim_m-FCemp_global).^2))
        ssimFC(nG)=ssim(FCsim_m,FCemp_global)
        FCsim_m_G(nG,:,:)=FCsim_m;
        
        nG=nG+1;

    end
    FCsim_G(nsub,:,:,:)=FCsim_m_G;
    fitFCnusb(nsub,:)=fitFC;
    distFCnusb(nsub,:)=distFC;
    ssimFCnusb(nsub,:)=ssimFC;
    mMetansub(nsub,:)=mMeta;

end

%% ---- working point: G* = max of the smoothed metastability curve (as for the chimps) ----
mMetaSm = smooth(mean(mMetansub));
[~, iG] = max(mMetaSm);  Gstar = Grange(iG);
[~, iFit]  = max(mean(fitFCnusb));
[~, iDist] = min(mean(distFCnusb));
[~, iSsim] = max(mean(ssimFCnusb));
mMetaemp = nanmean(Metasim_sub22);  sMetaemp = nanstd(Metasim_sub22);
fprintf('G* (max metastability) = %.4f, metastability model %.3f vs empirical %.3f +/- %.3f\n', ...
    Gstar, mean(mMetansub(:,iG)), mMetaemp, sMetaemp);
fprintf('best FC fit: r at G = %.4f, distance at G = %.4f, SSIM at G = %.4f\n', ...
    Grange(iFit), Grange(iDist), Grange(iSsim));

figure
subplot(2,2,1); plot_mean_std(Grange, mean(fitFCnusb), std(fitFCnusb)); xline(Gstar,'k--'); title('FC correlation')
subplot(2,2,2); plot_mean_std(Grange, mean(distFCnusb), std(distFCnusb),'r'); xline(Gstar,'k--'); title('FC distance')
subplot(2,2,3); plot_mean_std(Grange, mean(mMetansub), std(mMetansub)); hold on
plot_mean_std(Grange, repmat(mMetaemp,1,numel(Grange)), repmat(sMetaemp,1,numel(Grange)), 'g');
xline(Gstar,'k--'); title('Metastability (green: empirical, same band)')
subplot(2,2,4); plot_mean_std(Grange, mean(ssimFCnusb), std(ssimFCnusb)); xline(Gstar,'k--'); title('FC SSIM')

% FCs comparison at G*: upper triangle empirical, lower triangle model
FCsim_meanG = squeeze(mean(FCsim_G));
FCsim_meanGopt = squeeze(FCsim_meanG(iG,:,:));
FCemp_mean = FCemp_global;
Cplot = triu(FCemp_mean,1) + tril(FCsim_meanGopt,-1);
figure; imagesc(Cplot); axis square; colorbar;
title(sprintf('upper: empirical, lower: model at G* = %.4f', Gstar))

FCsim_Gmean = FCsim_meanG;          % mean over the 20 runs (G x N x N), small enough to save
save results_model_Humans_HCP7T_asCHIMP_bandFC.mat fitFCnusb distFCnusb ssimFCnusb mMetansub Grange ...
    FCsim_Gmean FCemp_global Metasim_sub22 Gstar iG
