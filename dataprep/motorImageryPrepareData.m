function motorImageryPrepareData(R,fresh)
close all
experimentName = 'imagery_basic';
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
typeCodes = {'tongueImg','tongueRImg','handImg','handRImg'};
%%
for k=1:length(subjects)
    
    %% load data and electrodes
    disp(['Subject ' subjects{k}])
    load([R.path.KM_exp '/data/' subjects{k} '_im_t_h'],'data','stim');
    load([R.path.KM_exp '/locs/' subjects{k} '_electrodes'],'electrodes');
    
    %% Common Preprocessing
    [data,fsamp] = commonPreProc(data);
    load('ns_1k_1_300_filt','nsfilt') % amplifier roll-off function (in amplitude)
    
    %% reshape stimulus variable
    blocksize = 1000 + (floor(srate/2)+1);
    clear a*
    
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
            pmint = 1;
            [Pxx,hz] =pwelch(epochDat(:),floor(srate*pmint),0,srate,srate);
            Pxx = mean(Pxx,2);
            epochPxx(:,mind,epc) = Pxx;
            
            % find SNR and peak frequency
            [betaSNR(mind,epc),betapeak(mind,epc)] = computeBetaSNR_AbsPeak(hz,epochPxx(:,mind,epc),[14 30]);
            
            %% Plot
            subplot(1,2,1)
            plot(hz,Pxx,ls); xlim([4 48]);
            hold on
            title(betaSNR(mind,epc))
        end % loop through epochs
        clf
    end % loop through channels
    legend({'tongueImg','tongueRImg','handImg','handRImg'})
    % find beta desync from rest to movement
    lr = [1 3]; powDsync = []; dsyncFrq = [];
    for seg = 1:2 % either tongue or hand
        [powDsync(:,seg),dsyncFrq(:,seg)] = findDesync(hz,epochPxx,[lr(seg)+1 lr(seg)]);
    end
    
    % choose contact in M1 with maximum desync
    [chsel(1),maxds(1),snr(1,:),acceptflagCat(1)] = selectDataCriteria(squeeze(betaSNR(:,2))',squeeze(powDsync(:,1))','powmod',5);
    if maxds(1,:)>-10; acceptflagCat(1) = 0; end
    [chsel(2),maxds(2),snr(2,:),acceptflagCat(2)] = selectDataCriteria(squeeze(betaSNR(:,4))',squeeze(powDsync(:,2))','powmod',5);
    if maxds(2)>-10; acceptflagCat(2) = 0; end
    
    lr = [1 1 2 2];
    for epc = 4:-1:1
        if epc == 1
            code = 11; % tongue imag
        elseif epc == 2
            code = 110; % tongue rest
        elseif epc == 3
            code = 12;  % hand imag
        elseif epc == 4
            code = 120; % hand rest
        end
        if ~acceptflagCat(lr(epc))
            [~,ip] = max(snr(lr(epc),:));
            dataStore = squeeze(dataSave(:,indM1(ip),tr_sc==code));
            dataStore(:,any(isnan(dataStore))) = [];
            [Pxx,hz] =pwelch(dataStore,pmint*fsamp,0,fsamp,fsamp);
            Pxx = mean(Pxx,2);
            plot(hz,Pxx); xlim([4 48]);
            hold on
            title('Not accepted')
            dataStore = []; chDist = [];
        else
            dataStore = squeeze(dataSave(:,indM1(chsel(lr(epc))),tr_sc==code));
            dataStore(:,any(isnan(dataStore))) = [];
            chDist = dist(chsel(lr(epc)));
            [Pxx,hz] =pwelch(dataStore,pmint*fsamp,0,fsamp,fsamp);
            Pxx = mean(Pxx,2);
            plot(hz,Pxx); xlim([4 48]);
            hold on
            title('Accepted')
        end
        
        if fresh
            acceptflag = acceptflagCat(lr(epc));
            AcceptSave(k,epc) = acceptflag; 
        else
            load([R.path.localDatPath '/Generic/AcceptSaveMotorImagery'],'AcceptSave')
            acceptflag = AcceptSave(k,epc);
        end
        
        mkdir([R.path.localDatPath '/Generic/' subjects{k} '/'])
        save([R.path.localDatPath '/Generic/' subjects{k} '/' subjects{k} '_' typeCodes{epc} '_Data.mat'],'dataStore','chDist','acceptflag')
        
    end % loop epochs
end
if fresh
    save([R.path.localDatPath '/Generic/AcceptSaveMotorImagery'],'AcceptSave')
end