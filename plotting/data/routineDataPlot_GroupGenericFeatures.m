function routineDataPlot_GroupGenericFeatures(R)
close all
[bigStore,subjects,typeLabels,featBank] = loadABCGeneric(R,{'bigStore','subjects','typeLabels','featBank'});
typeLabels = {'ISI','Move. prep.','Move. exec.','Move imag.'};

%% Now do the plotting
cmap = linspecer(5);
R.data.feat_xscale{1} = R.frqz;
for feat = 1:size(featBank,2)
    for task = 1:numel(typeLabels)
        % Select by SNR
        N = [featBank{task,feat}]; N = sum(~isnan(N(1,:)));
        % Plot the features
        X = nanmean([featBank{task,feat}],2); %
        XS = nanstd([featBank{task,feat}],[],2)./sqrt( N);
        figure(1)
        if feat<4
            figure(1)
            subplot(2,3,feat)
        else
            figure(223); %SI
            subplot(2,2,1);
        end
        [B(task),P(task)] = boundedline((R.data.feat_xscale{feat}),X,XS);
        P(task).FaceAlpha = 0.2;
        P(task).FaceColor  = cmap(task,:);
        B(task).Color = cmap(task,:);
        B(task).LineWidth = 1.5;
        hold on
        a = gca;
        
        if feat == 1
            a.YAxis.Scale = 'linear';
            a.XAxis.Scale = 'linear';
        else
            %             a.YAxis.Scale = 'log'
            %             a.XAxis.Scale = 'log'
        end
        if feat == 1
            axis([6 40 0 3])
            xlabel('Frequency (Hz)')
            ylabel('Power (a.u.)')
            title('ECoG Spectra')
        elseif feat == 2
            axis([1 8 0 0.9])
            xlabel('Burst amp. (Z)')
            ylabel('Probability density')
            title('Burst amplitudes')
        elseif feat ==3
            axis([30 600 0 0.01])
            xlabel('Burst dur.(ms)')
            ylabel('Probability density')
            title('Burst durations')
        elseif feat ==4
            axis([0 2.5e3 0 6e-4])
            xlabel('Burst int.(ms)')
            ylabel('Probability density')
            title('Burst intervals')
        end
        box off; axis square; grid on
    end
end
legend(B,typeLabels);
set(gcf,'Position',[650.6000  191.4000  712.8000  579.2000])


figure(2)
DC = table2array(bigStore);
LABEL = DC(:,2);
featname = {'Narrow-band SNR (dB)','Mean dur. (ms)','Mean amp (Z)','Wide-band SNR (dB)','Peak freq (hz)','Mean Int (ms)'};
featid = [6 8 10 5 7 12];
ylimz = [5 30;
    100 400;
    2 4;
    10 30;
    0 40;
    300 900];

for id = 1:6
    if id<4
        figure(1)
        subplot(2,3,id+3)
    else
        figure(223);
        subplot(2,2,1+id-3);
    end
    
    % plot SNR
    X = DC(:,featid(id));
    [X,adjLABEL] = removeNormOutliers(X,LABEL,1.96);
    [p1 t stats] = anova1(log10(X),adjLABEL,'off');
    [B,V,E] = scatterBar(X,adjLABEL,cmap)
    
    for i = 1:numel(V)
        for j = 1:numel(V)
            [h,pt(i,j),ci(:,i,j)] = ttest2(log10(X(adjLABEL==j)),log10(X(adjLABEL==i)));
            
            if pt(i,j)<=0.001
                text(i,median(X)+(0.1*median(X))+(0.05*j*max(X)),'***','Color',cmap(j,:))
            elseif pt(i,j)<=0.01
                text(i,median(X)+(0.1*median(X))+(0.05*j*max(X)),'**','Color',cmap(j,:))
            elseif pt(i,j)<=0.05
                text(i,median(X)+(0.1*median(X))+(0.05*j*max(X)),'*','Color',cmap(j,:))
            end
        end
    end
    a = gca;
    a.XTickLabel = typeLabels;
    a.XTickLabelRotation = -45;
    ylabel(featname{id});
    box off; axis square; grid on
    title(sprintf('ANOVA P = %.4f',p1),'FontSize',6)
    ylim([ylimz(id,:)])
    xlim([0.5 numel(unique(adjLABEL))+.5])
end
figure(1)
set(gcf,'Position',[425.0000  221.0000  988.8000  541.0000])
figure(223)
set(gcf,'Position',[623.4000  112.2000  654.4000  633.6000])
