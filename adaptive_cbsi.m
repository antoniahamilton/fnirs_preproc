
%% this code does adaptive CBSI

clear all
close all

%% Create syntethic Data
freq=1; % Sampling Frequency
duration=650*freq; % Signal length 600 s
tt = freq:freq:(duration*freq);  %% generate time course

onsetsA = [ 121   297   356   428   528];   %% some fixed values for simplicity
onsetsB = [  50   179   250   496   568];
durA = [25 25 25 25 25];
durB = [25 25 25 25 25];

stim=zeros(duration,2);  %% generate box car
for i=1:length(onsetsA)
    stim(onsetsA(i) : onsetsA(i)+durA(i),1)=1;
    stim(onsetsB(i) : onsetsB(i)+durB(i),2)=1;
end

% convolve with HRF
hh = spm_hrf(freq,[6,16,1,1,6,0,32]);  %% spm12 hrf with default params

cfA = conv(stim(:,1),hh);
cfB = conv(stim(:,2),hh);
cfA=cfA(1:duration);
cfB=cfB(1:duration);

cfA=cfA./max(cfA); % 1 a.u.
cfB=cfB./max(cfB); % 1 a.u.

%% fix the ground-truth beta values
gtbeta = [0.1,0.2];   %% Stim B is double stim A.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nppt = 25;   %% number of simulated participants
nrep = 20;  %% number of replications of each expt

%% mean noise level and inter-signal SD of noise
snrlevels = [0.2, 0;   %% constant low noise
    0.4, 0;        %% constant high noise
    0.2, 0.05;
    0.4, 0.05;
    0.2, 0.2;      %% variable low noise
    0.4, 0.2];      %% variable high noise

for i=1:size(snrlevels,1)
    snrnames{i} = ['m=',num2str(snrlevels(i,1)),' s=',num2str(snrlevels(i,2))];
end

allsnrlevels = kron(snrlevels,ones(nrep,1));

deoxysc = -1/3;  %% scaling on deoxy data relative to oxy

for kk=1:length(allsnrlevels)   %% for each simulated experiment
    
    snr = repmat(allsnrlevels(kk,1),nppt,2) + ...  %% mean noise
     randn(nppt,2)*allsnrlevels(kk,2);  %% std noise

    Mnoise = 0.1;  %% weighting on Mayer-wave
        
    for ppt = 1:nppt
        
        % Create Mayer Wave - [0.08-0.15] Hz
        x=Mnoise*randn(duration,1);
        [D,C]=butter(4,[0.08 0.15]./(freq/2));
        Mayer=filtfilt(D,C,x);
        
        %Gaussian Noise on HBO
        oNoise= snr(ppt,1)*randn(duration,1);        
        % Synthetic HBO signal ( Task component + Noise + Mayer wave )
        HBO_raw= gtbeta(1)*cfA+ gtbeta(2)*cfB +Mayer+oNoise; % Stim B= 2 a.u.
        
        %Gaussian Noise on HBR
        dNoise=snr(ppt,2)*randn(duration,1);  %% new noise so it is not the same
        % Synthetic HBR signal ( Task component + Noise)
        HBR_raw=(  gtbeta(1)*cfA+ gtbeta(2)*cfB+dNoise)*(deoxysc)  ; % Stim B= 2 a.u.
        
        %Band pass filter: 0.008 0.2
        [B,A] = butter(2,[0.008 0.2]./(freq/2));
        HBO_f(:,ppt)=filtfilt(B,A,HBO_raw);
        HBR_f(:,ppt)=filtfilt(B,A,HBR_raw);
        
        %% the synthetic data has been generated
        %--------------------------------
        %% now apply different processing options
                        
        %% calculate alpha in CBSI so we can store it.
        sd_oxy = std(HBO_f(:,ppt),0,1);
        sd_deoxy = std(HBR_f(:,ppt),0,1);
        
        alpha(kk,ppt) = sd_oxy/sd_deoxy;
        
    end
    
    %% find the average alpha value
    adapt_alpha = mean(alpha(kk,:));
    %adapt_alpha = 3;
    
    %% generate the design matrix
    desX = [cfA, cfB, ones(length(cfA),1)];   %% design matrix
        
    %% now analyse each ppt
    for ppt = 1:nppt
        
        %% HDiff measure
        hdiff = HBO_f(:,ppt) - HBR_f(:,ppt);
        
        %% trad CBSI
        [HB_cbsi]=Cbsi(HBO_f(:,ppt),HBR_f(:,ppt));
    
        %% adaptive CBSI
        Hada = HBO_f(:,ppt) - adapt_alpha * HBR_f(:,ppt);
        
        %% now fit standard GLM for both
        
        %% fit to filtered data
        [boxy] = regress(HBO_f(:,ppt),desX);
        [bdxy] = regress(HBR_f(:,ppt),desX);
        
        %%%-------------------------------------
        %% fit to cbsi data
        [bcbsi] = regress(HB_cbsi,desX);
        
        %% now do H-diff measure
        bhdiff = regress(hdiff,desX);
        
        %% now do H-adaptive measure
        bhada = regress(Hada,desX);
      
        %% save all betas
        b(1,kk).name = 'oxy-filt';
        b(1,kk).vals(:,ppt) = boxy;
        b(2,kk).name = 'dxy-filt';
        b(2,kk).vals(:,ppt) = bdxy;
        b(3,kk).name = 'CBSI-Hb';
        b(3,kk).vals(:,ppt) = bcbsi;
        b(4,kk).name = 'Hdiff';
        b(4,kk).vals(:,ppt) = bhdiff;
        b(5,kk).name = 'Hadapt';
        b(5,kk).vals(:,ppt) = bhada;

        drawnow
        
    end
    kk
end

save adapt_cbsi_results b nrep deoxysc snrnames alpha

%% calculate the t-values across all simulations
for kk=1:length(b)
    for i=1:5
        
        [h,prob,ci,stats] = ttest(b(i,kk).vals(1,:)',b(i,kk).vals(2,:)');  %% ttest between StimA and StimB for this simulation
        tval(i,kk) = stats.tstat;
        meanbetaA(i,kk) = mean(b(i,kk).vals(1,:));
        meanbetaB(i,kk) = mean(b(i,kk).vals(2,:));
        tnames{i} = b(i,kk).name;
    end
end

cols = {[0.8,0,0],[0,0,0.8],[0.8,0,0.8],[0,0.8,0],[0,0.5,1]}

%% plot the t values over all the simulations
figure(4), clf
for i=1:6
    ind = (1:nrep)+(i-1)*nrep;  %% get row indices for this simulation
    dat = abs(tval(:,ind))';
    
    subplot(3,2,i)
    md = mean(dat);
    sd = std(dat);
    errorbar(3:5,md(3:5),sd(3:5),'k.')
    hold on
    for j=3:5
        h=bar(j,md(j));
        set(h(1),'FaceColor',cols{j},'EdgeColor',cols{j})
    end
    errorbar(3:5,md(3:5),sd(3:5),'k.')

    set(gca,'XTick',3:5,'XTickLabel',tnames(3:5),'FontSize',14,'XLim',[2,6])
  %  ylabel('absolute T value')
    title(['Noise level: ',snrnames{i}],'FontSize',14)
    set(gca,'YLim',[0,40])
    
       tbl{i,1} = snrnames{i};
    for j=3:5
        tbl{i,j-1} = [num2str(md(j),'%0.3g'),' (',num2str(sd(j),'%0.3g'),')'];
    end
    
    %% do t-tests
    pp = [3,4; 3,5; 4,5];
    for j=1:length(pp)
        [h,sig(i,j),tstat] = ttest(dat(:,pp(j,1))-dat(:,pp(j,2)))
        tbl{i,j+4} = sig(i,j);
    end
    
    
end

xlswrite('cbsi_adapt_results.xlsx',tbl)


%% plot the adaptive alpha value
salpha = mean(alpha,2);  %% average over ppt
clear dat

figure(5), clf
for i=1:6
    ind = (1:nrep)+(i-1)*nrep;
    dat(:,i) = salpha(ind,:);
end
plot(dat','b.')
hold on
plot(mean(dat),'k-')
set(gca,'XLim',[0,7],'XTick',1:6)
xlabel('noise levels')
ylabel('value of alpha')

    
    


