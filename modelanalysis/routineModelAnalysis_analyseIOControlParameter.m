function routineModelAnalysis_analyseIOControlParameter(R,fresh)
[FitAnalysisData,mergeLabels] = loadABCGeneric(R,{'FitAnalysisDataGroupLevel','mergeLabels'});

dtype = 2; % all only
parMod = linspace(-2,2,25);
parSelName = {'int MMC G1','int MMC G5','int MMC G7','int MMC G9','int MMC G11'}; 
% MP self, MP->II, SP self, II->DP, SP->DP

NRep = 6; % number of reps per alphai

if fresh
    %% Now analyse the fits
    for parN = 1:5
        for typ = 1:4
            parBase = spm_unvec(FitAnalysisData.parMAPSave{typ,dtype},FitAnalysisData.parStruc.exampP);
            
            MMod = FitAnalysisData.simStruc.mmod{typ,dtype};
            RMod = FitAnalysisData.simStruc.Rmod{typ,dtype};
            RMod.obs.gainmeth = {'unitvar'  'obsnoise','leadfieldOneSignal'};
            MFit = FitAnalysisData.mfitSave{typ,dtype};
            RMod.plot.flag = 0;
            RMod = setSimTime(RMod,32);
            [pInd,pMu,pSig] = parOptInds_110817(RMod,MFit.MAP,MMod.m); % in structure form
            
            [pInd,pMu,pSig,pIndMap,pMuMap] = parOptInds_110817(RMod,parBase,MMod.m); % in structure form
            pfldnm = getParFieldNames(pMu,MMod); %convertNameListToReadable(R,getParFieldNames(pMu,MMod));
            ParIdx = pMuMap(find(ismember(pfldnm, parSelName{parN})));
            
            parfor alphai = 1:numel(parMod)
                for rep = 1:NRep
                    pMAP = parBase;
                    % Setup parameters
                    pMAP_flat = spm_vec(pMAP);
                    pMAP_flat(ParIdx) = pMAP_flat(ParIdx)+parMod(alphai);
                    pMAP = spm_unvec(pMAP_flat,pMAP);
                    [MI(rep,alphai),corrC(rep,alphai),MI_ds(rep,alphai),corrC_ds(rep,alphai),ampDur(:,rep,alphai),EI(rep,alphai),dEI(rep,alphai)] = computeSyncDesyncMI(RMod,MMod,pMAP);
                end
                disp([alphai typ parN])
            end
            PerturbResult.MI(:,:,typ,parN) = MI;
            PerturbResult.corrC(:,:,typ,parN) = corrC;
            PerturbResult.MI_ds(:,:,typ,parN) = MI_ds;
            PerturbResult.corrC_ds(:,:,typ,parN) = corrC_ds;
            PerturbResult.ampDur(:,:,:,typ,parN) = ampDur;
            PerturbResult.EI(:,:,typ,parN) = EI;
            PerturbResult.dEI(:,:,typ,parN) = dEI;
            PerturbResult.parnames = parSelName;
            PerturbResult.par = parMod;
        end
    end
    R.out.dag = 'model_ParCorrAnalysis';
    saveABCGeneric(R,{'PerturbResult'},{PerturbResult});
else
    PerturbResult = loadABCGeneric(R,{'PerturbResult'});
end