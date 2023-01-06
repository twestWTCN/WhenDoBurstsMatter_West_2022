function routineModelAnalysis_analyseGroupLevel(R,fresh)
dttag = {'spec','all','diff'};
[SAMod,mergeLabels] = loadABCGeneric(R,{'SAMod','mergeLabels'});

%% Get list of data
DL = 1;
for typ = 1:4
    for dtype = 1:2
        Rtmp =  R;
        Rtmp.out.dag = ['fitGroupLevelV1_' mergeLabels{typ}  '_' dttag{dtype}];
        
        try
            loadABCGeneric(Rtmp,{'modelfit'});
            daglist{DL} = Rtmp.out.dag;
            dagCode(:,DL) = [typ dtype];
            frqstat(:,DL) =  SAMod.peakFrq(typ);
            
            [ampMeanSTDEmp(1),ampMeanSTDEmp(2)] = pdfExpVal(SAMod.featStore(typ).feat_emp{2}, SAMod.featStore(typ).feat_xscale{2}');
            ampStatEmp(:,DL) = ampMeanSTDEmp;
            [durMeanSTDEmp(1),durMeanSTDEmp(2)] = pdfExpVal(SAMod.featStore(typ).feat_emp{3}, SAMod.featStore(typ).feat_xscale{3}');
            durStatEmp(:,DL) = durMeanSTDEmp;
            
            DL = DL + 1;
        catch
            disp('Not found!')
        end
    end
end


%% Get Model Draws (Bayesian confidence intervals)
if fresh
    R.analysis.dagtype = 'arbitrary';
    R.analysis.modEvi.N = 512;
    R.analysis.comptype = 'switch'; % this will use the feature weighting in R not Rmod
    R.objfx.featweight = [1 1 1 0];
    R.SimAn.RealzRep = 3; % 3 realizations per process
    modID = modelCompMaster_160620(R,1:numel(daglist),[],daglist);
end

%% Retrieve LE Spectrum
for DL = 1:numel(daglist)
    Rtmp = R;
    Rtmp.out.dag = daglist{DL};
    [permMod] = loadABCGeneric(Rtmp,{'modeProbs'});
    MLE(:,DL) = permMod.lyap;
end

%% Now analyse the fits
for DL = 1:numel(daglist)
    R.out.dag = daglist{DL};
    typ =  dagCode(1,DL);
    dtype =  dagCode(2,DL);
    [permMod,Rmod,mmod,mfit] = loadABCGeneric(R,{'modeProbs','R','modelspec','modelfit'});
    Rmod.path = R.path;
    
    % Save the fits
    inds = ~isnan(permMod.r2rep);
    r2{typ,dtype} = squeeze(permMod.errorVec_rep(1,:,inds));
    parBestSave{typ,dtype} = spm_vec(permMod.bestP); % this is the best fitting draw
    parMAPSave{typ,dtype} = spm_vec(permMod.MAP);    % this is the MAP estimate
    [parPermodAverage{typ,dtype} parPermodVec{typ,dtype}] = saveParSets(permMod.par_rep,inds); % this is average over draws
    
    % Retrieve burst stats
    ampMeanSTD{typ,dtype} = reconstructBurstStats(permMod.feat_rep,R.data.feat_xscale,2);
    durMeanSTD{typ,dtype} = reconstructBurstStats(permMod.feat_rep,R.data.feat_xscale,3);
    
    ampStatEmpSave{typ,dtype} = ampStatEmp(:,DL);
    durStatEmpSave{typ,dtype} = durStatEmp(:,DL);
    frqStatSave(typ,dtype) = frqstat(DL);
    MLESave{typ,dtype} = MLE(:,DL);
    
    RModSave{typ,dtype} = Rmod;
    mmodSave{typ,dtype} = mmod;
    mfitSave{typ,dtype} = mfit;
    
end

% get common parameters (due to reduced free parameters compared to ISI)
btype = 2; % baseline state
[pInd,pMu,pSig] = parOptInds_110817(RModSave{btype,2}, mfitSave{btype,2}.MAP,mmodSave{btype,2}.m); % in structure form
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.C'};
[pIndNew,pMuNew,pSigNew] = parOptInds_110817(RModSave{1,2}, mfitSave{btype,2}.MAP,mmodSave{btype,2}.m); % in structure form
[~,IA] = intersect(spm_vec(pMu),spm_vec(pMuNew));

% Get parameter stuctures
[pInd,pMu,pSig] = parOptInds_110817(Rmod,permMod.bestP,mmod.m); % in structure form

% Form descriptives
pMuMap = spm_vec(pMu);
pSigMap = spm_vec(pSig);
pfldnm = getParFieldNames(pMu,mmod);

% Put together
FitAnalysisData.r2 = r2;
FitAnalysisData.parBestSave = parBestSave;
FitAnalysisData.parMAPSave = parMAPSave;
FitAnalysisData.parPermodAverage = parPermodAverage;
FitAnalysisData.parPermodVec = parPermodVec;

FitAnalysisData.frq = frqstat;
FitAnalysisData.ampStatEmp = ampStatEmpSave;
FitAnalysisData.durStatEmp = durStatEmpSave;
% FitAnalysisData.r2lnStatEmp = r2lnStatEmp;
FitAnalysisData.ampMeanSTD = ampMeanSTD;
FitAnalysisData.durMeanSTD = durMeanSTD;

FitAnalysisData.MLE = MLESave;

FitAnalysisData.daglist = daglist;
FitAnalysisData.dagCode = dagCode;
FitAnalysisData.parStruc.exampP = permMod.bestP;
FitAnalysisData.parStruc.exampm = mmod.m;
FitAnalysisData.parStruc.pMuMap = pMuMap;
FitAnalysisData.parStruc.pSigMap = pSigMap;
FitAnalysisData.parStruc.pInd = pInd;
FitAnalysisData.parStruc.pMu = pMu;
FitAnalysisData.parStruc.pSig = pSig;
FitAnalysisData.parStruc.pfldnm = pfldnm;

FitAnalysisData.simStruc.mmod = mmodSave;
FitAnalysisData.simStruc.Rmod = RModSave;

FitAnalysisData.mfitSave = mfitSave;
R.out.dag = 'dataprep_beta';
saveABCGeneric(R,{'FitAnalysisDataGroupLevel'},{FitAnalysisData});


