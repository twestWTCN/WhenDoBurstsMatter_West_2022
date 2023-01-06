function routineModelAnalysisPlot_GroupLevelFeatures(R)

% Overide axis spread
R.plot.feat(1).axlim = [4 38 0 3.2];
R.plot.feat(2).axlim = [0 6 0 1];
R.plot.feat(3).axlim = [0 600 0 0.01];

R.plot.feat(1).axtit{2} = 'Power (a.u.)';
R.plot.feat(2).axtit{2} = 'Probability density';
R.plot.feat(3).axtit{2} = 'Probability density';

closeMessageBoxes; close all
cmap = linspecer(5);
typLabels = {'ISI','Move prep.','Move exec.','Move imag.'};
dataLab = {'ISI_FT','PreMove_FT','PostMove_FT','handImg'};
xlabs = {'freq (Hz)','burst amp. (ms)','burst dur. (ms)'};
rng(9523)
for typ = 1:4
    if typ<4
        load([R.path.localDatPath '\Generic\zt\zt_' dataLab{typ} '_Data.mat'])
    else
        load([R.path.localDatPath '\Generic\jc\jc_' dataLab{typ} '_Data.mat'])
    end
    dataStore = normalize(dataStore(:));
    dataStore = dataStore(1:end)';
    dataStore = bandpass(dataStore,[4 48],1000);
    dat{1} = dataStore;
    tvec{1} = linspace(0,size(dataStore,2)/1000,size(dataStore,2));
    
    % Load Simulated Data SPEC
    R.out.dag = ['fitGroupLevelV1_' typLabels{typ}  '_spec'];
    [permMod,Rmod,mmod] = loadABCGeneric(R,{'modeProbs','R','modelspec'});
    feat_sim{2} = computeFeatureAverage(Rmod,permMod.feat_rep);
    pMAP = permMod.MAP;
    [r2,~,~,~,xsims_gl,~,Rout] = computeSimData_160620(Rmod,mmod,[],pMAP,64,0);
    dat{2} = xsims_gl{1}(4,:);
    dat{2} = bandpass(dat{2}',[4 48],1000)';
    dat{2} = normalize(dat{2}(1:end-1));
    tvec{2} = Rout.IntP.tvec_obs;
    
    % Load Simulated Data ALL
    R.out.dag = ['fitGroupLevelV1_' typLabels{typ}  '_all'];
    [permMod,Rmod,mmod] = loadABCGeneric(R,{'modeProbs','R','modelspec'});
    feat_sim{3} = computeFeatureAverage(Rmod,permMod.feat_rep);
    pMAP = permMod.MAP;
    [r2,~,~,~,xsims_gl,~,Rout] = computeSimData_160620(Rmod,mmod,[],pMAP,64,0);
    dat{3} = xsims_gl{1}(4,:);
    dat{3} = bandpass(dat{3}',[4 48],1000)';
    dat{3} = normalize(dat{3}(1:end-1));
    tvec{3} = Rout.IntP.tvec_obs;
    
    % Reload in empirical data
    feat_sim{1} = Rmod.data.feat_emp;
    
    axlimSave = [
        16.6937   17.6937  -17.8917    5.5106
        2.8378    3.8378  -18.4233    3.4541
        15.7387   16.7387  -18.8090    3.2249
        19.8025   20.8025  -19.7211    2.8720
        ];
    
    figure(15)
    for i = 1:3
        %Plot time series
        if i == 1 % modify coloration
            modC = 0.5; ofc = 0;
        elseif i== 2
            modC = 0.7; ofc = -6;
        elseif i ==3
            modC = 1.0; ofc = -12;
        end
        subplot(4,4,typ)
        if i == 1
            plot(tvec{i},dat{i}+ofc,'Color','k','LineStyle','-','LineWidth',1); hold on
        else
            plot(tvec{i},dat{i}+ofc,'Color',cmap(typ,:).*modC,'LineWidth',1); hold on
        end
        axis(axlimSave(typ,:))
        if i == 3
            plot([10.2 10.4],[-16 -16],'Color','k','LineWidth',2)
        end
        box off; axis off; grid off
        
        %Plot data features
        for feat = 1:3
            subplot(4,4,typ+(4*feat))
            if i == 1 %empirical
                plot(Rmod.data.feat_xscale{feat},squeeze(feat_sim{i}{feat}),'Color','k','LineStyle','--','LineWidth',1); hold on
            else
                [hl, hp] = boundedline(Rmod.data.feat_xscale{feat},feat_sim{i}{feat}(:,1),feat_sim{i}{feat}(:,2),'cmap',cmap(typ,:).*modC,'alpha','transparency',0.2);
                hl.LineWidth = 1;
            end
            axis(R.plot.feat(feat).axlim);
            xlabel(xlabs{feat});         ylabel(R.plot.feat(feat).axtit{2}); box off
        end
    end
end
set(gcf,'Position',[ 41.0000   69.0000  680.8000  637.6000])
