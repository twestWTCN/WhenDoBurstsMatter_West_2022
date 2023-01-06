function KMGenericPrepareData(R,fresh)
close all
load('ns_1k_1_300_filt','nsfilt') % amplifier roll-off function (in amplitude)

R.data.feat_xscale{1} = R.frqz; % fill in
% Add KM Toolbox
addpath(genpath(fullfile(R.path.KM_root,'ctmr')));
addpath(genpath(fullfile(R.path.KM_root,'toolbox')));

%Typelist
typeLabels = {'ISI_FT','PreMove_FT','PostMove_FT','ISI_FT_T2','PreMove_FT_T2','PostMove_FT_T2','handImg','handRImg','tongueImg','tongueRImg'};
typeCodes = 1:numel(typeLabels);
mergeLabels = {'ISI','Move prep.','Move exec.','Move imag.'};
mergeCodes = [1 2 3 1 2 3 4 1 4 1]; % maps data to the four motor states
%% Sublist
subjects= ls([R.path.localDatPath '\Generic']);
subjects = cellstr(subjects);
subjects = subjects(find(strncmp(subjects,'bp',2)):end);
subCodes = 1:numel(subjects);
bandnames = {'beta','alpha'};
if fresh
    for BAND = 1
        bigStore = []; featBank = [];
        for sub = 1:size(subjects,1)
            subCodeList = [];
            for typ = 1:numel(typeLabels)
                file = [R.path.localDatPath '\Generic\' subjects{sub} '\' subjects{sub} '_' typeLabels{typ} '_Data.mat'];
                if isfile(file)
                    load(file,'dataStore','chDist','acceptflag')
                    sub
                    subCodeList = [subCodeList mergeCodes(typ)]; % this prevents duplicate states being calculated
                else
                    dataStore = [];
                end
                
                if ~isempty(dataStore) & acceptflag == 1
                    fsamp = 1000;
                    %Start computing data features
                    tp = computeAmpSpec(dataStore,nsfilt,fsamp); %
                    wbPower = 10*log10(mean(tp(4:48))./mean(tp([40:55 65:125]))); %wide band power
                    
                    [Pxx,hz] =pwelch(dataStore,fsamp/2,0,fsamp,fsamp);
                    Pxx = mean(Pxx,2);
                    
                    % find SNR and peak frequency
                    if BAND == 1
                        [betaSNR] = computeBetaSNR_AbsPeak(hz,Pxx,[14 30]);
                        [peakSize,betapeak] = findBetaPeak(hz,Pxx,[11 30]);
                        frq  = [betapeak-3 betapeak+3];
                    end
                    
                    % Burst features
                    R.data.winmark = getWinMark(dataStore);
                    dataStore = normalize(dataStore(:));
                    [durMeanSTD,ampMeanSTD,intMeanSTD,feat] = getBurstStats(R,dataStore(:),frq,fsamp,0);
                    
                    for ft = 1:4
                        if ft == 1
                            featBank{mergeCodes(typ),ft}(:,sub) = squeeze(feat{ft}(1,1,1,1,:));
                        else
                            featBank{mergeCodes(typ),ft}(:,sub) = squeeze(feat{ft});
                        end
                    end
                    
                    %% Compute Average Burst Waveforms
                    [burstSelIndsout,XF,XEnv,XPhi,epsAmp,segLout,segAmpout] = simpleBurstDefineMMC(dataStore',fsamp,frq,3,75);
                    [tvec,XXF,XXE] = getTimeLockedBurstWaveform(burstSelIndsout{1},XF,XEnv,[-0.150 0.600],fsamp);
                    burstWaveFormStore{mergeCodes(typ)}(:,:,sub) = [tvec; XXF; XXE];
                    
                    %% Compute Linear Surrogates
                    n_surr = 50; %5;
                    R.obs.trans.bursts.frq = frq;
                    [r2ln,bddls,surrFeat,surrBurstWave] = bddlSurrogate(dataStore,n_surr,R,feat);
                    r2ln = median(r2ln,2);
                    r2lnSave(:,typ,sub) = median(r2ln,2);
                    for ft = 1:4
                        surFeatBank{mergeCodes(typ),ft}(:,sub) =surrFeat{ft};
                    end
                    burstSurWaveFormStore{mergeCodes(typ)}(:,:,sub) = surrBurstWave; % store surrogate burst waveforms
                    % example data
                    exampDatStore{typ} = dataStore;
                    
                    % Save to bigStore
                    bigStore = [bigStore; typeCodes(typ) mergeCodes(typ) subCodes(sub) max(chDist) wbPower...
                        betaSNR betapeak durMeanSTD ampMeanSTD intMeanSTD r2ln'];
                    
                else % if not file
                    disp('Nothing loaded')
                    try
                        featBank{mergeCodes(typ),ft}(:,sub)
                    catch
                        for ft = 1:4
                            featBank{mergeCodes(typ),ft}(:,sub) = nan(size(R.data.feat_xscale{ft}));
                            surFeatBank{mergeCodes(typ),ft}(:,sub) = nan(size(R.data.feat_xscale{ft}));
                        end
                    end
                end % load file
            end % loop types
            
            %% Save Example Data
            if sub == 1
                R.out.dag  = ['dataprep_' bandnames{BAND}];
                saveABCGeneric(R,{'exampDatStore'},{exampDatStore})
            end
        end % loop subjects
        subjectList = subjects(bigStore(:,3));
        bigStore = array2table(bigStore,'VariableNames',{'type','merge','sub','chDist','wideSNR',...
            'betaSNR','betapeak','durMean','durSTD','ampMean','ampSTD','intMean','intSTD','r2lnSpec','r2lnAmp','r2lnDur','r2lnInt'});
        
        R.out.dag  = ['dataprep_' bandnames{BAND}];
        saveABCGeneric(R,{'bigStore','subjects','typeLabels','featBank',...
            'r2lnSave','mergeLabels','mergeCodes','surFeatBank','burstWaveFormStore','burstSurWaveFormStore'},...
            {bigStore,subjectList,typeLabels,featBank,...
            r2lnSave,mergeLabels,mergeCodes,surFeatBank,burstWaveFormStore,burstSurWaveFormStore})
        
    end
end
