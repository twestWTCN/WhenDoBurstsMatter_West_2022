function routineDataPlot_ExampleTimeSeries(R)
close all
[exampDatStore] = loadABCGeneric(R,{'exampDatStore'});
mergeLabels = {'ISI','PreMove','PostMove','Imagery'};
datSel = [1 2 3 6];
fsamp = 1000;
%% Now do the plotting
cmap = linspecer(5);
figure(1)
for task = 1:numel(mergeLabels)
    X = exampDatStore{datSel(task)};
    tvec = linspace(0,numel(X)/fsamp,numel(X));
    
    plot(tvec,X-(task*6),'Color',cmap(task,:),'LineWidth',1.5)
    xlim([12.2 15.2])
    hold on
    if task == 1
    plot([14.8 15],[-3 -3],'k','LineWidth',4)
    end
end

a = gca;
a.PlotBoxAspectRatioMode = 'auto';
box off; grid off; axis off

