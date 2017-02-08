
%% this code should generate all figures for the CBSI paper

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

for kk=1:length(allsnrlevels)
    
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
        
        %% the synthetic data has been generated
        %--------------------------------
        %% now apply different processing options
        
        % Signal preprocessing
        %Band pass filter: 0.008 0.2
        [B,A] = butter(2,[0.008 0.2]./(freq/2));
        HBO_f=filtfilt(B,A,HBO_raw);
        HBR_f=filtfilt(B,A,HBR_raw);
        
        % CBSI
        [HB_cbsi]=Cbsi(HBO_f,HBR_f);
        
        %% calculate alpha in CBSI so we can store it.
        sd_oxy = std(HBO_f,0,1);
        sd_deoxy = std(HBR_f,0,1);
        
        alpha(kk,ppt) = sd_oxy/sd_deoxy;
        
        %% generate the design matrix
        desX = [cfA, cfB, ones(length(cfA),1)];   %% design matrix
        
        %% HDiff measure
        hdiff = HBO_f - 3* HBR_f;
        
        %% illustrate data generation process
        if(ppt==1 && rem(kk,nrep)==0)
            
            ylimit=sum(allsnrlevels(kk,:))*3;
            
            figure(1), clf
            subplot(4,2,1)
            patch(tt,stim(:,1),[0,0.8,0])
            hold on
            patch(tt,stim(:,2),[0.8,0,0.8])
            title('Stimulus boxcars')
            set(gca,'XLim',[0,max(tt)])
            
            subplot(4,2,2)
            plot(tt,cfA,'g-','Color',[0,0.8,0])
            hold on
            plot(tt,cfB,'-','Color',[0.8,0,0.8])
            title('Convolved with HRF')
            set(gca,'XLim',[0,max(tt)])
            
            subplot(4,2,3)
            plot(tt,oNoise,'-','Color',[0.7,0,0])
            title('noise on oxy-Hb')
            set(gca,'XLim',[0,max(tt)],'YLim',[-ylimit,ylimit])
            
            subplot(4,2,4)
            plot(tt,dNoise,'-','Color',[0,0,0.7])
            title('noise on dxy-Hb')
            set(gca,'XLim',[0,max(tt)],'YLim',[-ylimit,ylimit])
            
            subplot(4,2,5)
            plot(tt,Mayer,'k-')
            title('Mayer wave noise')
            set(gca,'XLim',[0,max(tt)],'YLim',[-ylimit,ylimit])
            
            subplot(4,2,6), cla
            plot(HBO_raw,'r-')
            hold on
            set(gca,'XLim',[0,max(tt)],'YLim',[-ylimit,ylimit])
            plot(tt,desX(:,1)*gtbeta(1),'k')
            plot(tt,desX(:,2)*gtbeta(2),'-.k')
            title('oxy-Hb')
            
            subplot(4,2,7), cla
            plot(HBR_raw,'b-')
             hold on
            set(gca,'XLim',[0,max(tt)],'YLim',[-ylimit,ylimit])
            plot(tt,desX(:,1)*gtbeta(1),'k')
            plot(tt,desX(:,2)*gtbeta(2),'-.k')
            title('dxy-Hb')
            
            %% illustrate data preprocessing stream
           % ttl = {'Raw signals','filtered signals','CBSI-Hb','H-diff','H-total'};
            figure(2), clf
            for i=1:9
                subplot(3,3,i)
                plot(tt,desX(:,1)*gtbeta(1),'w')
                hold on
                plot(tt,desX(:,2)*gtbeta(2),'-.w')
                switch i
                      case 1
                        plot(tt,HBO_raw,'r-','Color',[0.8,0,0])
                        ttl = 'raw HbO_2'
                    case 4
                        plot(tt,HBR_raw,'c-','Color',[0,0,0.8])
                        ttl = 'raw HHb'
                   case 2
                        plot(tt,HBO_f,'r-','Color',[0.8,0,0])
                        ttl = 'filtered HbO_2'
                    case 5
                        plot(tt,HBR_f,'c-','Color',[0,0,0.8])
                        ttl = 'filtered HHb'
                    case 3
                        plot(tt,HB_cbsi,'g-','Color',[0.8,0,0.8])
                        ttl = 'H CBSI'
                    case 6
                        plot(tt,hdiff,'k-','Color',[0,0.8,0])
                        ttl = 'Hdiff'
                    case 9
                        plot(tt,HBO_f + HBR_f, 'c-','Color',[1,0.5,0])
                        ttl = 'Htotal'
                    otherwise
                        axis off
                end
               plot(tt,desX(:,1)*gtbeta(1),'k','LineWidth',2)
                hold on
                plot(tt,desX(:,2)*gtbeta(2),'-.k','LineWidth',2)
           % xlabel('time (s)')
                set(gca,'XLim',[0,max(tt)/2],'YLim',[-ylimit,ylimit])                
                title(ttl)
            end
            
            set(1,'Position',[50,50,800,900])
            set(2,'Position',[50,50,900,700],'PaperOrientation','Landscape')
          
            pause(0.1)
            
            pfn = ['sample_generate',num2str(kk),'.tiff'];
            print(1,'-dtiff',pfn)
                 
            pfn = ['sample_process',num2str(kk),'.tiff'];
            print(2,'-dtiff',pfn)
            
        end
        
        %% now get param estimates for all
        %% now fit standard GLM for both
        
        %% fit to filtered data
        [boxy] = regress(HBO_f,desX);
        [bdxy] = regress(HBR_f,desX);
        
        %%%-------------------------------------
        %% fit to cbsi data
        [bcbsi] = regress(HB_cbsi,desX);
        
        %% now do H-diff measure
        bhdiff = regress(hdiff,desX);
        
        %% now do H-total
        htot = HBO_f + HBR_f;
        bhtot = regress(htot,desX);
        
        %% save all betas
        b(1,kk).name = 'oxy-filtered';
        b(1,kk).vals(:,ppt) = boxy;
        b(2,kk).name = 'dxy-filtered';
        b(2,kk).vals(:,ppt) = bdxy;
        b(3,kk).name = 'CBSI-Hb';
        b(3,kk).vals(:,ppt) = bcbsi;
        b(4,kk).name = 'Hdiff';
        b(4,kk).vals(:,ppt) = bhdiff;
        b(5,kk).name = 'Htotal';
        b(5,kk).vals(:,ppt) = bhtot;

        drawnow
        
    end
    kk
end

save all_cbsi_results b nrep deoxysc snrnames alpha

