clear all
load('CHIMP_CWAS_DATA_FORDECO.mat')

N=114;
NSUB=59;
indexN=1:N;





% Parameters of the data
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
f_diff= 0.025*(rand(1,N)*0.1+1);
omega = repmat(2*pi*f_diff(1:N)',1,2); omega(:,1) = -omega(:,1);

cut=10;;


%% G optimo
Grange=[0:0.1:20]/1000;

Tmax=1200;

%% Group
for nsub=1:NSUB
    nsub
    C=squeeze(CJ_CHIMP(:,:,nsub));
    if length(find(isnan(C)))==(114*114)
        nsub
    end
    if length(find(isnan(C)))<(114*114)
        C(isnan(C))=0;
        %C=C/max(max(C))*0.2;
        C=C/nanmean(C(:));
        nnG=1;
        for G=Grange
            G
            for nrep=1:20
                for trials=1:2
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

                    end
                    signal_filt=signal_filt(:,cut:end-cut);


                    for seed=1:N
                        Xanalytic = hilbert(demean(signal_filt(seed,:)));
                        PhasesS(seed,:) = angle(Xanalytic);
                    end
                    FCsim2(trials,:,:)=corrcoef(signal_filt');
                    gKOMsim=nansum(complex(cos(PhasesS),sin(PhasesS)))/N;
                    enstrophy1sim=abs(gKOMsim);
                    Metasim_sub22(trials)=nanstd(enstrophy1sim(:)); %metastability
                end
                Metasim_sub2(nsub,nnG,nrep)=mean(Metasim_sub22);
                stdMetasim_sub2(nsub,nnG,nrep)=std(Metasim_sub22);
        %        FCsim(nsub,nnG,:,:)=squeeze(mean(FCsim2,1));
  %              FCsim_vec=FCsim(nsub,nnG,:,:);
  %              SCFCcorr(nsub,nnG)=corr2(C(:),FCsim_vec(:));
            end
            nnG=nnG+1;
        end
    end
end

%FCemp=squeeze(mean(FC));

 save meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat Metasim_sub2 stdMetasim_sub2


