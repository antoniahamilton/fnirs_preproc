
clear all, 
%%%     plot CBSI results nicely

load all_cbsi_results

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

tnames = {'Hb0_2','HHb','H-CBSI','H-diff','H-total'}

cols = {[0.8,0,0],[0,0,0.8],[0.8,0,0.8],[0,0.8,0],[1,0.5,0]}

%% plot the t values over all the simulations
figure(3), clf
for i=1:6
    ind = (1:nrep)+(i-1)*nrep;  %% get row indices for this simulation
    dat = abs(tval(:,ind))';
    
    subplot(3,2,i)
    md = mean(dat);
    sd = std(dat);
    errorbar(1:5,md,sd,'k.')
    hold on
    for j=1:5
        h=bar(j,md(j));
        set(h(1),'FaceColor',cols{j},'EdgeColor',cols{j})
    end
    
    set(gca,'XTick',1:5,'XTickLabel',tnames)
    set(gca,'YLim',[0,30])
   % ylabel('absolute T value')
    title(['Noise level: ',snrnames{i}])
    set(gca,'FontSize',14)
    
    tbl{i,1} = snrnames{i};
    for j=1:5
        tbl{i,j+1} = [num2str(md(j),'%0.3g'),' (',num2str(sd(j),'%0.3g'),')'];
    end
    
    %% do t-tests
    pp = [1,2; 1,3; 1,4; 1,5; 2,3; 2,4; 2,5; 3,4; 3,5; 4,5]
    for j=1:length(pp)
        [h,sig(i,j),tstat] = ttest(dat(:,pp(j,1))-dat(:,pp(j,2)))
        tbl{i,j+7} = sig(i,j);
    end
    
end

%% do a table of results

xlswrite('cbsi_basic_results.xlsx',tbl)


