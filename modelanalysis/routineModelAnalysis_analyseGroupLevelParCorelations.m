function routineModelAnalysis_analyseGroupLevelParCorelations(R,fresh)

% Load in FitAnalysis
FitAnalysisData = loadABCGeneric(R,{'FitAnalysisDataGroupLevel','mergeLabels'});

dtype = 2; % all only
parMod = linspace(-2,2,25);
NRep = 8;

if fresh
    %% Now analyse the fits
    figure(100)
    clf
    for typ = 1:4
        parBase = spm_unvec(FitAnalysisData.parMAPSave{typ,dtype},FitAnalysisData.parStruc.exampP);
        pMu=  FitAnalysisData.parStruc.pMu;
        pMuMap = [spm_vec(pMu.int{1}.G); spm_vec(pMu.int{1}.C)];
        
        MMod = FitAnalysisData.simStruc.mmod{typ,dtype};
        RMod = FitAnalysisData.simStruc.Rmod{typ,dtype};
        
        RMod.obs.gainmeth = {'unitvar','leadfieldOneSignal'};
        corrStore = []; lyapStore = []; statStore =[];
        
        %% get base power
        peakSizeBase = []; wflag = [];
        for rep = 1:NRep % average over realizations
            [~,~,feat_sim,~,~,wflag(:,rep)] = computeSimData_160620(RMod,MMod,[],parBase,48,0);
            if wflag(1,rep) ~= 1
                % peak
                fsamp = 1/RMod.IntP.dt;
                peakSizeBase(rep) = findBetaPeak(RMod.data.feat_xscale{1},squeeze(feat_sim{1}(1,1,1,1,:)),[14 30]);
            end
        end
        peakSizeBase = mean(peakSizeBase);
        
        for par = 1:numel(pMuMap)
            parfor alphai = 1:size(parMod,2)
                % Setup parameters
                parSet = parBase;
                pMAP_flat = spm_vec(parSet);
                pMAP_flat(pMuMap(par)) = parMod(alphai);
                parSet = spm_unvec(pMAP_flat,parSet);
                
                statPar = []; wflag= []; modelError = []; errorVec = []; r2linSave = []; lyap = [];
                % Run through each parameter
                for rep = 1:NRep
                    [modelError(rep),~,feat_sim,xsims,xsims_gl,wflag(:,rep),~,errorVec(:,rep),J] = computeSimData_160620(RMod,MMod,[],parSet,48,0);
                    if wflag(1,rep) ~= 1
                        % peak
                        fsamp = 1/RMod.IntP.dt;
                        [peakSize,betapeak] = findBetaPeak(RMod.data.feat_xscale{1},squeeze(feat_sim{1}(1,1,1,1,:)),[14 30]);
                        
                        statPar(rep,1,:) = betapeak; % peak freq;
                        statPar(rep,2,:) = 100*((peakSize-peakSizeBase)./(peakSizeBase)); % change in beta power
                        
                        dataX = xsims_gl{1}(4,:)';
                        dataX = bandpass(dataX,RMod.obs.trans.bursts.frq,fsamp);
                        dataX = (dataX-mean(dataX))./std(dataX);
                        minbs = (2/(RMod.obs.trans.bursts.frq(2)))*fsamp;
                        
                        [~,~,~,amp] = burstAmpHist(dataX',fsamp,R.data.feat_xscale{2},minbs);
                        statPar(rep,3,:) = mean(amp);
                        [~,~,~,dur] = burstDurHist(dataX',fsamp,R.data.feat_xscale{3},minbs);
                        statPar(rep,4,:) = mean(dur);
                        [~,~,~,int] = burstIntHist(dataX',fsamp,R.data.feat_xscale{4},minbs);
                        statPar(rep,5,:) = mean(int);
                        
                        % Lyapun
                        A = sort(real(eig(J{1})));
                        lyap(rep) = A(end);
                    else
                        statPar(rep,1:4,:) = nan(1,4)
                    end
                end
                
                if ~all(isnan(wflag))
                    lyapStore(:,alphai) = nanmedian(lyap);
                    statStore(:,alphai) = nanmedian(statPar);
                    r2Store(:,alphai) =  nanmedian(errorVec,2);
                else
                    lyapStore(:,alphai) = nan(1);
                    statStore(:,alphai) = nan(4,1);
                    r2Store(:,alphai) =  nan(5,1);
                end
            end % end of loop around modulations
            
            Y = statStore(3,:);
            X = parMod;
            p = []; index = [];
            [RCf(1),p(1),~,index{1}] = restrictedCorrelation(X,Y,8); %
            
            Y = statStore(3,:);
            X = parMod;
            [RCf(2),p(2),~,index{2}] = restrictedCorrelation(X,Y,8);
            
            [~,indSel] = max(abs(RCf)); %choose the range that yields the best correlation
            r2StoreSave{typ} = r2Store;
            statStoreSave{typ} = statStore;
            statIndSelSave{typ} = index{indSel};
            lyapStoreSave{typ} = lyapStore;
            
            % Put together the correlations with stats
            % this loops through: (1) peak freq; (2) % power change; (3)
            % burst amplitude; (4) burst duration; (5) burst in
            for i = 1:5
                X = statStore(i,:);
                [C,p] = corr(parMod(index{indSel})',X(index{indSel})','type','Pearson');
                corrStore(:,i,par) = [C(1,2) p(1,2) numel(index{indSel}) range(X(index{indSel}))];
            end
            
            % this looks at model fit (R2)
            [C,p] = corr(X(index{indSel})',r2Store(2,index{indSel})','type','Pearson');
            corrStore(:,6,par) = [C(1,1) p(1,1) numel(index{indSel}) range(X(index{indSel}))];
            
        end
        corrStoreSave(:,:,:,typ) = corrStore;
        % corrStoreSave has dimension 4 x 10 x 18 x 4
        % dim1 = [C p N range]
        % dim2 = [peakhz %power amp dur int R2 R2LinSpec R2Lamp R2Ldur
        %                                                       R2Lint
        % dim3 = [parlist]
        % dim4 = type
        
        R.out.dag = 'model_ParCorrAnalysis';
        saveABCGeneric(R,{'corrStoreSave','statStoreSave','statIndSelSave','r2StoreSave','lyapStoreSave'},{corrStoreSave,statStoreSave,statIndSelSave,r2StoreSave,lyapStoreSave});
    end
end

