function motorBasicPrepareData(R,fresh)
close all
experimentName = 'motor_basic';
R.path.KM_exp = fullfile(R.path.KM_root,experimentName);

%% set defaults
srate=1000; % sample rate
warning('off','signal:psd:PSDisObsolete'); %annoying

%%
subjects={
    'bp',...  % mot good; im good
    'fp',... %       % mot bad; im bad
    'hh',...%        % mot good; im good
    'jc',...
    'jm',...%
    'rh',...
    'rr'...
    };

typeCodes = {'ISI_FT_T2','PreMove_FT_T2','PostMove_FT_T2','TransWindow_FT_T2'};
%%
for k=1:length(subjects); % 13 15 is HARD
    %% load data and electrodes
    disp(['Subject ' subjects{k}])
    load([R.path.KM_exp '/data/' subjects{k} '_mot_t_h'],'data','stim');
    load([R.path.KM_exp '/locs/' subjects{k} '_electrodes'],'electrodes');
    
    % clear previous
    for epc = 1:4; delete([R.path.localDatPath '/Generic/' subjects{k} '/' subjects{k} '_' typeCodes{epc} '_Data.mat']); end
    
    %% Common Preprocessing
    [data,fsamp] = commonPreProc(data);
    load('ns_1k_1_300_filt','nsfilt') % amplifier roll-off function (in amplitude)
    
    blocksize = 1000 + (floor(srate/2)+1);
    clear a*
    
    %% INITIAL EPOCHING (YOU REDEFINE LATER)
    % find rest conditions specific to behavior that preceded them
    stim(find(and(...
        [zeros(blocksize,1); stim(1:(length(stim)-blocksize)) ]==12, ...
        stim==0 ...
        )))=120; %hand rest block condition
    
    stim(find(and(...
        [zeros(blocksize,1); stim(1:(length(stim)-blocksize)) ]==11, ...
        stim==0 ...
        )))=110; %tongue rest block condition
    
    %% catalogue trials
    trialnr=0*stim; tr_sc=0; %initialize
    trtemp=1; trialnr(1)=trtemp; %start variable
    for n=2:length(stim)
        if stim(n)~=stim(n-1),
            trtemp=trtemp+1;
            tr_sc=[tr_sc stim(n)];
        end
        trialnr(n)=trtemp;
    end
    clear n trtemp
    trials=unique(trialnr);
    
    %% calculate spectra for all trials
    num_chans=length(data(1,:));
    mean_PSD=zeros(200,num_chans); % to calculate the mean spectra over the entire period
    all_PSD=[]; dataSave = [];
    for cur_trial=1:max(trials),
        % isolate relevant data
        curr_data=data(find((trialnr == cur_trial)),:);
        % keep only .5 to block end of data, for each block
        if length(curr_data(:,1))>=blocksize
            curr_data=curr_data((floor(srate/2)+1):blocksize,:);
            dataSave(:,:,cur_trial) = curr_data;
        else
            dataSave(:,:,cur_trial) = nan(blocksize-floor(srate/2),size(data,2));
        end
    end % session
    
    clear cur* *temp* p
    
    %% find motor cortex
    [bp_data,chloc,indM1,dist] = findClosest2ROI(data,electrodes);
    if ~isempty(indM1)
        clear wbPower epochPxx betaSNR betapeak
        for mind = 1:numel(indM1)
            for epc =1:4
                if epc == 1
                    code = 11; ls = 'r'; % tongue movement
                elseif epc == 2
                    code = 110; ls = 'k'; % tongue rest
                elseif epc == 3
                    code = 12; ls = 'b'; % hand movement
                elseif epc == 4
                    code = 120; ls = 'g';% hand rest
                end
                
                
                %% average spectra for different conditions
                % Run through epochs
                epochDat = squeeze(dataSave(:,indM1(mind),tr_sc==code));
                epochDat(:,any(isnan(epochDat))) = [];
                [tp,hz,bgrnd] = computeAmpSpec(epochDat,nsfilt,fsamp);
                wbPower(mind,epc) = mean(tp);
                epochPxx(:,mind,epc) = tp;
                
                % find SNR and peak frequency
                [betaSNR(mind,epc),betapeak(mind,epc)] = computeBetaSNR_AbsPeak(hz,epochPxx(:,mind,epc),[14 30]);
                
                %% Plot
                subplot(1,2,1)
                plot(hz,tp,ls); xlim([4 48]);
                hold on
                title(betaSNR(mind,epc))
            end % loop through epochs
            if any(betaSNR(mind,[2 4])>0)
                disp((betaSNR(mind,[2 4])))
            else
                disp('not accepted:')
                disp((betaSNR(mind,[2 4])))
            end
            clf
        end % loop through channels
        legend({'tonguemove','tonguerest','Handmove','Handrest'})
        % find beta desync from rest to movement
        lr = [1 3]; powDsync = []; dsyncFrq = [];
        for seg = 1:2 % either tongue or hand
            [powDsync(:,seg),dsyncFrq(:,seg)] = findDesync(hz,epochPxx,[lr(seg)+1 lr(seg)]); % [ premove move]
        end
        
        % choose contact in M1 with maximum desync
        [chsel(1),maxds(1),snr(1,:),acceptflagCat(1)] = selectDataCriteria(squeeze(betaSNR(:,2))',squeeze(powDsync(:,1))','powmod',5);
        if maxds(1,:)>-10; acceptflagCat(1) = 0; end
        [chsel(2),maxds(2),snr(2,:),acceptflagCat(2)] = selectDataCriteria(squeeze(betaSNR(:,4))',squeeze(powDsync(:,2))','powmod',5);
        if maxds(2)>-10; acceptflagCat(2) = 0; end
        
        
        %% Setup trial definitions
        pmint = 1; %0.75; % period
        pofs = 0.3; % offset from movement onset (ensures no cross over)
        
        % Find rest trials
        isiTrial = SplitVec(find(stim == 120),'Consecutive');
        isiInds = [];
        for seg = 1:numel(isiTrial)
            isiOnset = isiTrial{seg}(1);
            isiInds(seg,:) = [isiOnset+fix(pofs*fsamp) isiOnset+fix(pmint*fsamp)+fix(pofs*fsamp)];
        end
        
        % Redefine the trials using beta
        % Define Premovement Onset from Beta
        moveTrial = SplitVec(find(stim == 12),'Consecutive');
        s = 0; moveOnset = [];
        for seg = 1:numel(moveTrial)
            segInds = moveTrial{seg};
            if any(segInds<1)
                continue
            end
            try
                % load subject specific reaction time
                load([R.path.localDatPath '/Generic/' subjects{k} '/' subjects{k} '_PreMove_FT_Data'],'stidnSave')
                mOnsTmp = median(stidnSave);
            catch
                % load group mean
                load([R.path.localDatPath '/Generic/AcceptSaveFingerTap'],'medianResponseTime')
                mOnsTmp = medianResponseTime;
            end
            if ~isempty(mOnsTmp)&  (segInds(1) > (pmint+pofs)*fsamp)
                s = s+1;
                moveOnset(s) = segInds(1)+fix(mOnsTmp*fsamp);
            end
        end
        
        %
        digInds_premm = moveOnset';
        digInds_premm = [digInds_premm-fix(pmint*fsamp)-fix(pofs*fsamp) digInds_premm-fix(pofs*fsamp)];
        
        digInds_pstmm = moveOnset';
        digInds_pstmm = [digInds_pstmm digInds_pstmm+(pmint*fsamp)];
        
        digInds_transWin = moveOnset';
        digInds_transWin = [digInds_transWin-fix(2*fsamp) digInds_transWin+fix(1.5*fsamp)];
        
        % ISI Pre Post Trans
        lr = [2 2 2 2]; % all are from the hand
        for epc = 1:4
            if epc == 1
                ind = isiInds; ls = 'r'; % you only need this once or else repeat across digits
            elseif epc == 2
                ind = digInds_premm; ls = 'k';
            elseif epc == 3
                ind = digInds_pstmm; ls = 'b';
            elseif epc == 4
                ind = digInds_transWin; ls = '--';
            end
            
            
            if ~acceptflagCat(lr(epc))
                [~,ip] = max(snr(lr(epc),:));
                dataStore = squeeze(dataSave(:,indM1(ip),tr_sc==code));
                dataStore(:,any(isnan(dataStore))) = [];
                
                [Pxx,hz] =pwelch(dataStore,pmint*fsamp,0,fsamp,fsamp);
                Pxx = mean(Pxx,2);
                plot(hz,Pxx); xlim([8 48]);
                hold on
                title('Not accepted')
                dataStore = []; chDist = []; chSelOverall(epc) = nan;
            else
                dataStore = unpack(data(:,indM1(chsel(lr(epc)))),ind);
                dataStore(:,any(isnan(dataStore))) = [];
                chDist = dist(chsel(lr(epc)));
                chSelOverall(epc) = indM1(chsel(lr(epc)));
                [Pxx,hz] =pwelch(dataStore,pmint*fsamp,0,fsamp,fsamp);
                Pxx = mean(Pxx,2);
                plot(hz,Pxx); xlim([8 48]);
                hold on
                title('Accepted')
            end
            if fresh
                acceptflag = acceptflagCat(lr(epc));
                AcceptSave(k,epc) = acceptflag;
            else
                load([R.path.localDatPath '/Generic/AcceptSaveMotorBasic'],'AcceptSave')
                acceptflag = AcceptSave(k,epc);
            end
            pause(0.2)
            indStruc = [];
            indStruc.isiInds = isiInds;
            indStruc.digInds_premm = digInds_premm;
            indStruc.digInds_pstmm = digInds_pstmm;
            indStruc.digInds_transWin = digInds_transWin;
            
            mkdir([R.path.localDatPath '/Generic/' subjects{k} '/'])
            save([R.path.localDatPath '/Generic/' subjects{k} '/' subjects{k} '_' typeCodes{epc} '_Data.mat'],'dataStore','chDist','acceptflag','chSelOverall','indStruc')
            save([R.path.localDatPath '/Generic/epochdetail/' subjects{k} '_' typeCodes{epc} '_EpochDetail.mat'],'chDist','acceptflag','chSelOverall','indStruc')
        end
    end % loop epochs
    clf
end % check for sensorss
if fresh
    save([R.path.localDatPath '/Generic/AcceptSaveMotorBasic'],'AcceptSave')
end
