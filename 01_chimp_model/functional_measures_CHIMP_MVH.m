clear all
close all

addpath('/Users/cbc/Documents/sanz/Mami/')
addpath('/Users/cbc/Documents/sanz/utils')

load('meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat')

load('CHIMP_CWAS_DATA_FORDECO.mat')

% metastability computation

N=114;
num_mammals=59;
indexN=1:N;

for s=1:num_mammals

    Metasim_sub(s,:)=smooth(squeeze(mean(Metasim_sub2(s,:,:),3)));

end


Grange=[0:0.1:20]/1000;


Meta = nan(num_mammals,1);
Gs = nan(num_mammals,1);


[Meta Gs2] = max(Metasim_sub,[],2);
Gs(:,1)= Grange(Gs2);



TR=0.72;  % Repetition Time (seconds)
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
%f_diff=0.025*ones(1,N);

f_diff= 0.025*(rand(1,N)*0.1+1);


omega = repmat(2*pi*f_diff(1:N)',1,2); omega(:,1) = -omega(:,1);

%% G optimo
Grange=[0:0.1:10]/1000;

Tmax=1200;
N=114;

cut=10;;

Isubdiag = find(tril(ones(N),-1));
IsubdiagT = find(tril(ones(Tmax-2*cut+1),-1));

for repe=1:1
    for nsub=1:num_mammals

        nsub
        C=squeeze(CJ_CHIMP(:,:,nsub));
        if length(find(isnan(C)))==(114*114)
            nsub
        end
        if length(find(isnan(C)))<(114*114)
            C(isnan(C))=0;

            C=C/mean(C(:));
          %  C = randmio_und(C,1000); % rand preserving degree
            %C= null_model_und_sign(C); % rand preserving degree, strength and distribution

            G=Grange(Gs2(nsub));



            for trials=1:200
                trials
                wC = G*C;
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
                    tise =signal_filt(seed,:);
                    ev1=tise>std(tise);   %% >std(tise)
                    ev2=[0 ev1(1:end-1)];
                    events(seed,:)=(ev1-ev2)>0;
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
                Metasim_sub22(trials)=nanstd(enstrophy1sim(:)); %metastability
                Synchro_sub22(trials)=nanmean(enstrophy1sim(:)); % sincro
                    
                ts22=zscore(ts2');
                for t=1:size(ts22,1)
                    tmp1=ts22(t,:)'*ts22(t,:);
                    EdgesL(:,t)=tmp1(Isubdiag)';
                    GBCEdge(:,t)=mean(tmp1);
                end

                FCD=(EdgesL'*EdgesL)./(vecnorm(EdgesL)'*vecnorm(EdgesL));
                EdgesMeta2(trials)=std(EdgesL(:)); %Edge metastability
                EdgeMetaNode2(trials,:)=std(GBCEdge,[],2);

                EntroEdgesMeta2(trials)=0.5*(log(2*pi*var(FCD(IsubdiagT))))+0.5; % entropy mesta?

                for seed=1:N
                    ig2(seed)=corr2(Amplitude(seed,1:end-2),enstrophy1sim(3:end));
                end

                Ignition2(trials)=mean(ig2);
                IgStd2(trials)=std(ig2);
                [~,Cmean2(trials)] = LZc_measure(ts2);
                %     [ LZcom22(trials), ~, ~] = calc_lz_complexity(events(:), 'exhaustive', true);

                MaxLag =3;
                GCN= granger(ts2,MaxLag);
                m_granger(trials)=mean(GCN(:));
                s_granger(trials)=std(GCN(:));
                stdevokedintegS=ignition(events,114,Tmax);
                m_ignition(trials)=mean(stdevokedintegS);
                s_ignition(trials)=std(stdevokedintegS);

            end

            RH.Metasim_sub2(nsub)=mean(Metasim_sub22);
            RH.Synchro_sub2(nsub)=mean(Synchro_sub22);

            RH.Ignition(nsub)=mean(Ignition2);
            RH.IgStd(nsub)=std(IgStd2);
            RH.EntroEdgesMeta(nsub)=mean(EntroEdgesMeta2);
            RH.EdgesMeta(nsub)=mean(EdgesMeta2);
            RH.EdgeMetaNode(nsub,:)=mean(EdgeMetaNode2);
            % measuers over the FC
            FCsim= squeeze(nanmean(FCsim2));
            RH.FClh(nsub)=mean(mean(FCsim(1:39,1:39)));
            RH.FCrh(nsub)=mean(mean(FCsim(40:end,40:end)));

            RH.FCsim(nsub,:,:) = FCsim;
            RH.FCSCcorr(nsub) = corr(FCsim(:),C(:));

            RH.FCintehe(nsub)=(mean(mean(FCsim(1:39,40:end)))+mean(mean(FCsim(40:end,1:39))))*0.5;

            [Ci Q] = modularity_und(FCsim);
            RH.Qs(nsub)=Q;
            RH.Integration(nsub)=nanmean(nanmean(FCsim));
            RH.stdGBC(nsub)=std(nanmean(FCsim));
            RH.Cmean(nsub)=mean(Cmean2);
            %    RH.LZcom(nsub)=mean(LZcom22);

            RH.m_granger(nsub)=mean(m_granger);
            RH.s_granger(nsub)=mean(s_granger);

            RH.m_ignition(nsub)=mean(m_ignition);
            RH.s_ignition(nsub)=mean(s_ignition);
        end
    end
    RH.sigma=sig;
    save(sprintf('functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_REAL_0926_%03d.mat',repe),'RH');


end