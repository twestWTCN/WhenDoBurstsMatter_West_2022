function routineModelAnalysisPlot_GroupLevelFitStats(R)
close all
dttag = {'type A','type B','diff'};
cmap = linspecer(5);
[FitAnalysisData,mergeLabels] = loadABCGeneric(R,{'FitAnalysisDataGroupLevel','mergeLabels'});

%% First is the plot comparing the R2 between fit types
figure(1)
featnames = {'All','Spectra','Burst dur.','Burst amp'};
for dtype = 1:3
    subplot(2,3,dtype)
    if dtype<3
        for typ = 1:4
            % throw out failed simulations
            sucIndx = ~any(isinf(FitAnalysisData.r2{typ,dtype}));
            
            r2sum{dtype}(:,typ) = nanmedian(FitAnalysisData.r2{typ,dtype}(:,sucIndx),2);
            r2std{dtype}(:,typ) = nanstd(FitAnalysisData.r2{typ,dtype}(:,sucIndx),[],2);
            r2sem{dtype}(:,typ) = nanstd(FitAnalysisData.r2{typ,dtype}(:,sucIndx),[],2)./sqrt(sum(sucIndx));
        end
    else
        for typ = 1:4
            % throw out failed simulations
            sucIndx = ~any(isinf(FitAnalysisData.r2{typ,1}) | isinf(FitAnalysisData.r2{typ,2}));
            
            r2sum{dtype}(:,typ) = nanmedian(FitAnalysisData.r2{typ,1}(:,sucIndx),2) - nanmedian(FitAnalysisData.r2{typ,2}(:,sucIndx),2); % difference between dtype 1 and 2
            r2std{dtype}(:,typ) = nanstd(FitAnalysisData.r2{typ,1}(:,sucIndx)-FitAnalysisData.r2{typ,2}(:,sucIndx),[],2);
            r2sem{dtype}(:,typ) = nanstd(FitAnalysisData.r2{typ,1}(:,sucIndx)-FitAnalysisData.r2{typ,2}(:,sucIndx),[],2)./sqrt(sum(sucIndx));
        end
        
    end
    
    [B,E] = makeGroupBarPlotError(r2sum{dtype}(1:4,:),r2std{dtype}(1:4,:),cmap);
    a= gca;
    a.XTickLabel = featnames;
    a.XTickLabelRotation = 45;
    ylabel('R2')
    if dtype<3
        ylim([0.4 1]);
    else
        ylim([-0.5 0.05]);
    end
    xlim([.5 4.5])
    title(dttag{dtype})
    if dtype == 2
        L = legend(B,mergeLabels);
        L.EdgeColor = 'none';
        L.Box = 'off';
        L.Color = 'none';
        L.NumColumns = 4;
    end
end

set(gcf,'Position',[488.0000  412.2000  927.4000  349.8000])

%% Plot of the errors in estimation
ampError = []; ampErrorStd = [];
durError = []; durErrorStd = [];
for dtype = 1:2
    for typ = 1:4
        X = FitAnalysisData.ampMeanSTD{typ,dtype}(1,:)-FitAnalysisData.ampStatEmp{typ,dtype}(1,:);
        ampError(typ,dtype) = X;
        ampErrorStd(typ,dtype) = stdDiff(FitAnalysisData.ampMeanSTD{typ,dtype}(2,:),FitAnalysisData.ampStatEmp{typ,dtype}(2,:))./sqrt(256);
        X = FitAnalysisData.durMeanSTD{typ,dtype}(1,:)-FitAnalysisData.durStatEmp{typ,dtype}(1,:);
        durError(typ,dtype) = X;
        durErrorStd(typ,dtype) = stdDiff(FitAnalysisData.durMeanSTD{typ,dtype}(2,:),FitAnalysisData.durStatEmp{typ,dtype}(2,:))./sqrt(256);
    end
end
subplot(2,3,4)
[B,E] = makeGroupBarPlotError(ampError',ampErrorStd',cmap);
a = gca;
a.XTickLabel = dttag(1:2);
a.XTickLabelRotation = 45;
ylabel('Error in amplitude (Z)')
subplot(2,3,5)
[B,E] = makeGroupBarPlotError(durError',durErrorStd',cmap);
a = gca;
a.XTickLabel = dttag(1:3);
a.XTickLabelRotation = 45;
ylabel('Error in duration (ms)')

%% PLOT OF MLE
% Plot Scatter of R2
for dtype = 1:3
    for typ = 1:4
        
        if dtype<3
            MLE{typ,dtype} = FitAnalysisData.MLE{typ,dtype};
            MLE{typ,dtype}(MLE{typ,dtype}>10) = nan;
            MLE{typ,dtype}(MLE{typ,dtype}<-10) = nan;
            
            fMLE_Mean(typ,dtype) = nanmedian(MLE{typ,dtype});
            fMLE_Std(typ,dtype) = nanstd(MLE{typ,dtype})./sqrt(size(MLE{typ,dtype},1));
        else
            fMLE_Mean(typ,dtype) = nanmedian(MLE{typ,2}) - nanmedian(MLE{typ,1});
            fMLE_Std(typ,dtype) = stdDiff(nanstd(MLE{typ,1}),nanstd(MLE{typ,2}))./sqrt(size(MLE{typ,1},1));           
        end
    end
end
subplot(2,3,6)

[B,E] = makeGroupBarPlotError(fMLE_Mean',fMLE_Std',cmap);
a = gca;
a.XTickLabel = dttag(1:3);
a.XTickLabelRotation = 45;
ylabel('max. Lyapunov exponent')
set(gcf,'Position',[ 652   182   938   580])
L.Position = [0.3224    0.0277    0.4147    0.0345];

%% FINAL PLOT CONFIGURATIONS
figure(200)
for dtype = 1:2
    subplot(2,2,2+(dtype-1))
    X = repmat(1:4,4,1)+linspace(-.2,.2,4)'; Y  = r2sum{dtype}(1:4,:)'; YE = r2std{dtype}(1:4,:)'; XC = repmat(1:4,1,4);
    errorbar(X(:),Y(:),YE(:),'LineStyle','none','Color','k','CapSize',0,'LineWidth',1)
    hold on
    scatter(X(:),Y(:),50,cmap(XC(:),:),'filled')
    a= gca;
    a.XTick = 1:4;
    a.XTickLabel = featnames;
    a.XTickLabelRotation = 45;
    ylabel('R2')
    ylim([0.65 1]);xlim([0.5 4.5])
    title(dttag{dtype})
    box off
end

subplot(2,2,4)
[B,E] = makeGroupBarPlotError(r2sum{3}(1:4,:),r2sem{3}(1:4,:),cmap);
a= gca;
a.XTickLabel = featnames;
a.XTickLabelRotation = 45;
ylabel('R2')
ylim([-0.20 0.02]);
xlim([.5 4.5])
title(dttag{3})
L = legend(B,mergeLabels);
L.EdgeColor = 'none';
L.Box = 'off';
L.Color = 'none';
L.NumColumns = 4;
box off
