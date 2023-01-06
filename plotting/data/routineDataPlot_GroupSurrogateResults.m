function routineDataPlot_GroupSurrogateResults(R)
close all
[bigStore,subjects,typeLabels,featBank,surFeatBank,burstWaveFormStore,burstSurWaveFormStore] = loadABCGeneric(R,{'bigStore','subjects','typeLabels',...
    'featBank','surFeatBank','burstWaveFormStore','burstSurWaveFormStore'});
typeLabels = {'ISI','Move prep.','Move exec.','Move imag.'};

%% Now do the plotting
cmap = linspecer(5);
figure(2643)
% waveform plots
for task = 1:numel(typeLabels)
    subplot(2,2,task)
    %% Retrieve burst waveforms
    D = burstWaveFormStore{task};
    D(D==0) = nan;
    realN = sum(~isnan(squeeze(D(2,end,:))));
    X = nanmean(D,3);
    XS =  nanstd(D,[],3)/sqrt(realN);
    
    DS = burstSurWaveFormStore{task};
    DS(DS==0) = nan;
    XSur = nanmean(DS,3);
    XSSur = nanstd(DS,[],3)/sqrt(realN);
    
    [hl(1),hp(1)] =  boundedline(X(1,:),X(3,:),XS(3,:));
    hp(1).FaceAlpha = 0; hp(1).FaceColor = cmap(task,:);
    hl(1).Color = cmap(task,:); hl(1).LineWidth = 1.5;
    
    hold on
    [hl(2),hp(2)] =  boundedline(XSur(1,:),XSur(3,:),XSSur(3,:));
    hp(2).FaceAlpha = 0; hp(2).FaceColor = cmap(task,:);
    hl(2).Color = cmap(task,:); hl(2).LineWidth = 1.5; hl(2).LineStyle = '--';
    
    plot(X(1,:),X(2,:),'Color',cmap(task,:))
    hold on
    plot(XSur(1,:),XSur(2,:),'Color',cmap(task,:),'LineStyle','--')
    
    plot([550 650],[-1 -1],'k','LineWidth',4)
    xlabel('Time from onset (ms)'); ylabel('Amplitude');
    xlim([-150 750]); ylim([-1.5 1.5])
    
    permSigStat(X(1,:),squeeze(D(3,:,:)),squeeze(DS(3,:,:)),500,1.25,0,cmap(task,:),0.1,5)
    
    a = gca;
    a.PlotBoxAspectRatioMode = 'auto';
end
set(gcf,'Position',[ 488.0000  381.8000  870.6000  380.2000]);

R.data.feat_xscale{1} = R.frqz;
for feat = 1:3
    for task = 1:numel(typeLabels)
        %% This retrieves the features
        inds = 1:numel(R.data.feat_xscale{feat});
        % Plot the features
        X = nanmean([featBank{task,feat}],2); %
        XS = nanstd([featBank{task,feat}],[],2)./sqrt(size(featBank{task,feat},2));
        nX = sum(sum(isnan([featBank{task,feat}]))==0); % check for n
        % Surrogate features
        Y = nanmean([surFeatBank{task,feat}],2); %
        YS = nanstd([surFeatBank{task,feat}],[],2)./sqrt(size(surFeatBank{task,feat},2));
        nY = sum(sum(isnan([surFeatBank{task,feat}]))==0); % check for n
        
        X = X(inds);
        Y = Y(inds);
        
        tmp = 100.*([featBank{task,feat}]-[surFeatBank{task,feat}])./max([featBank{task,feat}]);
        Z = nanmean(tmp,2);
        ZS = nanstd(tmp,[],2)./sqrt(size(tmp,2));
        
        figure(1)
        subplot(1,3,feat)
        %         subplot(1,4,feat)
        plot(R.data.feat_xscale{feat}(inds),X,'color',cmap(task,:),'LineWidth',1.5)
        hold on
        plot(R.data.feat_xscale{feat}(inds),Y,'color',cmap(task,:),'LineWidth',1.5,'LineStyle','--')

        if feat == 1
            axis([6 40 0 3])
            xlabel('frequency (Hz)')
            ylabel('pdf')
            title('ECoG Spectra')
        elseif feat == 2
            axis([1 8 0 0.9])
            xlabel('burst amp. (norm)')
            ylabel('pdf')
            title('Burst amplitudes')
        elseif feat ==3
            axis([30 600 0 0.01])
            xlabel('burst dur.(ms)')
            ylabel('pdf')
            title('Burst durations')
        elseif feat ==4
            axis([0 2.5e3 0 1e-3])
            xlabel('burst int.(ms)')
            ylabel('pdf')
            title('Burst intervals')
        end
        figure(2644)
        
        subplot(2,3,feat)
        [h,p] = boundedline(R.data.feat_xscale{feat}(inds),Z,ZS);
        p.FaceColor = cmap(task,:);
        p.FaceAlpha = 0.2;
        h.Color =  cmap(task,:);
        h.LineWidth = 2;
        
        hold on
        if feat == 1
            axis([6 40 -50 40])
            xlabel('Frequency (Hz)')
            ylabel('% diff from surr.')
            title('ECoG Spectra')
        elseif feat == 2
            axis([1 8 -70 40])
            xlabel('Burst amp. (norm)')
            ylabel('% diff from surr.')
            title('Burst amplitudes')
        elseif feat ==3
            axis([30 600  -70 40])
            xlabel('Burst dur.(ms)')
            ylabel('% diff from surr.')
            title('Burst durations')
        elseif feat ==4
            axis([0 2.5e3  -70 40])
            xlabel('Burst int.(ms)')
            ylabel('% diff from surr.')
            title('Burst intervals')
        end
        box off; axis square; grid on
        
    end
    
end

%% Look at Nonlinearities also
DC = table2array(bigStore);
LABEL = DC(:,2);

featid = [14 15 16 17];
for id = 1:3
    % plot SNR
    subplot(2,3,id+3)
    X = DC(:,featid(id));
    [X,adjLABEL] = removeNormOutliers(X,LABEL,1.96);
    
    [p1 t stats] = kruskalwallis(X,adjLABEL,'off');
    [B,V,E] = scatterBar(X,adjLABEL,cmap)
    
    for i = 1:numel(V)
        [p0(i,id),h,stat] = signrank(X(adjLABEL==i),1);
        
        for j = 1:numel(V)
            [pt(i,j),h] = ranksum((X(adjLABEL==j)),(X(adjLABEL==i)));
            
            if pt(i,j)<0.001
                text(i,max(X)+(0.05*max(X))+(0.02*j*max(X)),'***','Color',cmap(j,:))
            elseif pt(i,j)<0.01
                text(i,max(X)+(0.05*max(X))+(0.02*j*max(X)),'**','Color',cmap(j,:))
            elseif pt(i,j)<0.05
                text(i,max(X)+(0.05*max(X))+(0.02*j*max(X)),'*','Color',cmap(j,:))
            end
        end
    end
    a = gca;
    a.XTickLabel = typeLabels;
    a.XTickLabelRotation = -45;
    ylabel('Fit to lin surrogate (R2)');
    box off; axis square; grid on
    title(sprintf('ANOVA P = %.4f',p1),'FontSize',6)
    ylim([0 1.5]);
end
a = [];
set(gcf,'Position',[ 488.0000  263.4000  920.2000  498.6000])

