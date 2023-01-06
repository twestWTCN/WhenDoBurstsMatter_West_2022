function routineModelAnalysis_analyseGroupLevelExampleParCorrs(R,fresh)
cmap = linspecer(5);
%Typelist
%% Sublist
FitAnalysisData = loadABCGeneric(R,{'FitAnalysisDataGroupLevel','mergeLabels'});

dtype = 2; % all only
parMod = linspace(-2,2,25);

parSelName = {'int MMC G7','int MMC G10','int MMC G4'}; % SP->SP; DP->DP; II->II
NRep = 3; % number of reps per alphai

if fresh
    for parN = 1:3
        statPar = nan(4,NRep,size(parMod,2)); lyap = nan(NRep,size(parMod,2),4); lyapS = nan(4,NRep,size(parMod,2),4); conjL = nan(8,18); r2linSave = nan(4,NRep,size(parMod,2),4);
        XEq = {}; lambdaReal = []; lambdaImag = [];
        %% Now analyse the fits
        for typ = 1:4
            parBase = spm_unvec(FitAnalysisData.parMAPSave{typ,dtype},FitAnalysisData.parStruc.exampP);
            
            MMod = FitAnalysisData.simStruc.mmod{typ,dtype};
            RMod = FitAnalysisData.simStruc.Rmod{typ,dtype};
            RMod.obs.gainmeth = {'unitvar','leadfieldOneSignal'}; % get rid of other criteria
            MFit = FitAnalysisData.mfitSave{typ,dtype};
            RMod.plot.flag = 0;
            RMod = setSimTime(RMod,32);
            [pInd,pMu,pSig] = parOptInds_110817(RMod,MFit.MAP,MMod.m); % in structure form
            
            [pInd,pMu,pSig,pIndMap,pMuMap] = parOptInds_110817(RMod,parBase,MMod.m); % in structure form
            pfldnm = getParFieldNames(pMu,MMod); %convertNameListToReadable(R,getParFieldNames(pMu,MMod));
            ParIdx = pMuMap(find(ismember(pfldnm, parSelName{parN})));
            for alphai = 1:numel(parMod)
                SP = []; feat_sim = [];
                pMAP = parBase;
                % Setup parameters
                pMAP_flat = spm_vec(pMAP);
                pMAP_flat(ParIdx) = pMAP_flat(ParIdx)+parMod(alphai);
                pMAP = spm_unvec(pMAP_flat,pMAP);
                [~,~,~,xsims] = computeSimData_160620(RMod,MMod,[],pMAP,32,0);
                MModtmp = MMod;
                MModtmp.x{1} = mean(xsims{1}(:,end-(2/RMod.IntP.dt):end),2)';
                
                RModtmp = RMod;
                RModtmp.IntP.Utype = 'constant';
                RModtmp.obs.gainmeth = {}; % get rid of other criteria
                
                [~,~,~,xsims,xsims_gl,~,RJ,~,J] = computeSimData_160620(RModtmp,MModtmp,[],pMAP,0,0);
                X = xsims{1}(4,end-(2/RMod.IntP.dt):end);
                XD = diff(X);
                tmp = unique(round(X(abs(XD)<0.1),1));
                XEq{alphai,typ} = double(tmp); 
                lambdaReal(:,alphai,typ) = sort(real(eig(J{1})));
                lambdaImag(:,alphai,typ) = sort(imag(eig(J{1})));
                eigSave(:,alphai,typ) = eig(J{1});
                XDet{alphai,typ} = xsims_gl{1};
                
                for rep = 1:NRep
                    % Run through each parameter going back to noise
                    [modelError(alphai),~,feat_sim{rep},xsims,xsims_gl{alphai},~,~,~,J] = computeSimData_160620(RMod,MMod,[],pMAP,32,0);
                    A = sort(real(eig(J{1})));
                    lyap(rep,alphai,typ) = max(A);
                    lyapS(:,rep,alphai,typ) = A;
                    XStoc{alphai,typ} = xsims_gl{alphai}{1};
                   
                    % compute lin surr
                    n_surr = 25; %5;
                    r2ln = bddlSurrogate(xsims_gl{alphai}{1}(RMod.siminds,:),n_surr,RMod,feat_sim{rep});
                    r2linSave(:,rep,alphai,typ) = mean(r2ln,2);
                    fsamp = 1/RMod.IntP.dt;
                    [md id] = max(squeeze(feat_sim{rep}{1}(1,1,1,1,RMod.data.feat_xscale{1}>=14 & RMod.data.feat_xscale{1} <=30)));
                    SP(1) = 11+ RMod.data.feat_xscale{1}(id); % peak freq;
                    dataX = xsims_gl{alphai}{1}(1,:)';
                    dataX = bandpass(dataX,RMod.obs.trans.bursts.frq,fsamp);
                    dataX = (dataX-mean(dataX))./std(dataX);
                        minbs = (2/(RMod.obs.trans.bursts.frq(2)))*fsamp;
                    
                    [~,~,~,amp] = burstAmpHist(dataX',fsamp,RMod.data.feat_xscale{2},minbs);
                    SP(2) = mean(amp);
                    [~,~,~,dur] = burstDurHist(dataX',fsamp,RMod.data.feat_xscale{3},minbs);
                    SP(3) = mean(dur);
                    [~,~,~,int] = burstIntHist(dataX',fsamp,RMod.data.feat_xscale{4},minbs);
                    SP(4) = mean(int);
                    statPar(:,rep,alphai,typ) = SP;
                    disp([alphai rep])
                end
                feat_sim_avRep{alphai,typ} = averageFeatCells(feat_sim);
                feat_sim_Rep{alphai,typ} = feat_sim;
            end
            
            corrExampSave.XEq{parN} = XEq;
            corrExampSave.lambdaRealXEq{parN} = lambdaReal;
            corrExampSave.lambdaImagXEq{parN} = lambdaImag;
            corrExampSave.eigSave{parN} = eigSave;
            corrExampSave.parName = parSelName;
            corrExampSave.lyapPar{parN} = lyap;
            corrExampSave.lyapSPar{parN} = lyapS;
            corrExampSave.r2linParSave{parN} = r2linSave;
            corrExampSave.statParPar{parN} = statPar;
            corrExampSave.feat_sim_RepPar{parN} = feat_sim_Rep;
            corrExampSave.feat_sim_avRepPar{parN} = feat_sim_avRep;
            corrExampSave.parMod = parMod;
            simTracesStoc{parN} = XStoc;
            simTracesDet{parN} = XDet;
            
            R.out.dag = 'model_ParCorrAnalysis';
            saveABCGeneric(R,{'corrExampSave'},{corrExampSave})
            saveABCGeneric(R,{'simTracesDet','simTracesStoc'},{simTracesDet,simTracesStoc})
        end
    end
end
