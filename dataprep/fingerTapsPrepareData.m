function fingerTapsPrepareData(R,fresh)
close all

% Get experiment path names
experimentName = 'fingerflex';
R.path.KM_exp = fullfile(R.path.KM_root,experimentName);

% Add Kai Miller's Toolbox
addpath(genpath(fullfile(R.path.KM_root,'ctmr')));
addpath(genpath(fullfile(R.path.KM_root,'toolbox')));

%% Sublist
subjects=[
    'bp';...
    'cc';...
    'zt';...
    'jp';...
% %     'ht';... % bad quality data
% %     'mv';... % missing move data
% %     'wc';... % bad quality data
    'wm';...
    'jc';...
    ];
typeCodes = {'ISI_FT','PreMove_FT','PostMove_FT','TransWindow_FT'};
for q=1:size(subjects,1)
    subject=subjects(q,:);
    %% Identify Epochs
    rng(7342); % needed for consistent choice of the ISI samples
    [st_pts,end_pts,isi_pts,flex,stidnSave]= gen_inv_pts_TW(R,subject);  % this function gets the events
    
    % Key: pts=[movementstart movementinversion movementlabel]
    % Movement label: 0 ISI; 1-5 Thumb to little finger
    % Load in Experimental Data
    load([R.path.KM_exp '/data/' subject '/' subject '_fingerflex'],'data','locs','flex')
        
    % Clear any existing data
    for epc = 1:3; delete([R.path.localDatPath '/Generic/' subject '/' subject '_' typeCodes{epc} '_Data.mat']); end
    %% Common Preprocessing
    [data,fsamp] = commonPreProc(data);
    [~,~,indM1,dist] = findClosest2ROI(data,locs); % list of indices close to M1
    if ~isempty(indM1)
        % Loop through digits (different locations for each)
        clear snr acceptflag chsel isiIndsSave digInds_premmSave digInds_execSave digInds_pstmmSave digInds_transWinSave
        for dig = 1:5
            % Premovement period
            pmint = 0.75; % period
            pofs = 0.25; % offset from movement onset (ensures no cross over)
            
            % ISI Inds
            isiInds = isi_pts(:,1);
            isiInds = [isiInds isiInds+(pmint*fsamp)];
            isiIndsSave{dig} = isiInds;
            
            % get movement onset
            digInds_premm = st_pts(st_pts(:,2) == dig,1);
            digInds_premm = [digInds_premm-fix(pmint*fsamp)-fix(pofs*fsamp) digInds_premm-fix(pofs*fsamp)];
            digInds_premmSave{dig} = digInds_premm;
            
            % get movement execution
            digInds_exec = st_pts(st_pts(:,2) == dig,1);
            digInds_exec = [digInds_exec+fix(pofs*fsamp) digInds_exec+(pmint*fsamp)+fix(pofs*fsamp)];
            digInds_execSave{dig} = digInds_exec;
            
            
            % get movement offset
            digInds_pstmm = end_pts(end_pts(:,2) == dig,1);
            digInds_pstmm = [digInds_pstmm digInds_pstmm+(pmint*fsamp)];
            digInds_pstmmSave{dig} = digInds_pstmm;
            
            %Get large window across transition
            digInds_transWin = st_pts(st_pts(:,2) == dig,1);
            digInds_transWin = [digInds_transWin-fix(2*fsamp) digInds_transWin+fix(1.5*fsamp)];
            digInds_transWinSave{dig} = digInds_transWin;
            
            clear dsyncFrq betaSNR betapeak wbPower epochPxx
            for epc = 1:3 % loop between epochs
                if epc == 1
                    ind = isiInds; ls = 'r';
                elseif epc == 2
                    ind = digInds_premm; ls = 'k';
                elseif epc == 3
                    ind = digInds_exec; ls = 'b';
                end
                for mind = 1:numel(indM1) % loop through channels
                    % Plot
                    epochDat = unpack(data(:,indM1(mind)),ind);
                    epochDat(:,any(isnan(epochDat))) = [];
                    load('ns_1k_1_300_filt','nsfilt') % amplifier roll-off function (in amplitude)
                    [tp,hz,bgrnd] = computeAmpSpec(epochDat,nsfilt,fsamp);
                    wbPower(mind,epc) = mean(tp);
                    disp([mind,epc,dig])
                    [Pxx,hz] =pwelch(epochDat,pmint*fsamp,0,fsamp,fsamp);
                    Pxx = mean(Pxx,2);
                    epochPxx(:,mind,epc) = Pxx;
                    
                    % find SNR and peak frequency
                    [betaSNR(mind,epc),betapeak(mind,epc)] = computeBetaSNR_AbsPeak(hz,epochPxx(:,mind,epc),[14 30]);
                    
                    %% Plot
                    subplot(1,2,1)
                    plot(hz,Pxx,ls); xlim([4 48]);
                    hold on
                    subplot(1,2,2)
                    flexDat = unpack(flex(:,dig),ind);
                    plot(linspace(0,1,size(flexDat,1)),mean(flexDat,2),ls)
                    hold on
                    clf
                end % loop through channels
            end% Loop through epochs
            
            % This is selecting the single channel per digit
            [powDsync,dsyncFrq] = findDesync(hz,epochPxx,[2 3]); % find max desync between pre and move
            [chsel(dig),maxds(dig),snr(:,dig),acceptflagDig(dig)] = selectDataCriteria(betaSNR(:,[1 2])',powDsync,'powmod',5); % 
            
            set{1} = 1:size(locs,1);
            set{2} = indM1';
            set{3} = indM1(betaSNR(:,1)>18)';
            set{4} =chsel(dig);
            roi = [37 -25 62; -37 -25 62];
            % plotSensorLocations(locs,roi,set,powDsync,betaSNR) % Uncomment to retrieve sensor locations
        end % loop through digits

        % Identify the digit that exhibits the maximum desync
        [dsyncSel,digSelMax] = min(maxds);
        digSel = find(maxds<-5);
        
        %% Now re-epoch the data using the selected channel
        for epc = 1:4
            ind = [];
            for dig = digSel
                if epc == 1
                    ind{1} = [isiIndsSave{dig}]; ls = 'r'; % you only need this once or else repeat across digits
                elseif epc == 2
                    ind{dig} = [digInds_premmSave{dig}]; ls = 'k';
                elseif epc == 3
                    ind{dig} = [digInds_execSave{dig}]; ls = 'b';
                elseif epc == 4
                    ind{dig} = [digInds_transWinSave{dig}]; ls = '--';
                end
            end
            chSelOverall = []; 
            if ~acceptflagDig(digSelMax)
                % This plots the bad channels
                [~,ip] = max(betaSNR(:,epc));
                if epc<4
                    dataStore = unpack(data(:,indM1(ip)),ind);
                else
                    dataStore = unpack(data(:,indM1(ip)),ind);
                end
                [Pxx,hz] =pwelch(dataStore,pmint*fsamp,0,fsamp,fsamp);
                Pxx = mean(Pxx,2);
                PxxS(:,1,epc) = Pxx;
                plot(hz,Pxx,ls); xlim([4 48]);
                hold on
                title('Not accepted')
                dataStore = []; chDist = []; chSelOverall(epc) = nan;
                
            else
                % This plots the bad channels
                dataStore = []; flexStore = [];
                for dig = digSel
                    if epc== 1
                        dataStore = unpack(data(:,indM1(chsel(digSelMax))),ind{1});
                    elseif epc <4
                        dataStore = [dataStore unpack(data(:,indM1(chsel(dig))),ind{dig})];
                    else
                        dataStore = [dataStore unpack(data(:,indM1(chsel(digSel))),ind{dig})];
                        flexStore = [flexStore unpack(flex(:,digSel),ind{dig})];
                    end
                end
                
                if epc == 4
                    dataStore(:,:,1) = dataStore;
                    dataStore(:,:,2) = flexStore;
                end
                
                chSelOverall{epc} = indM1(chsel(digSel));
                chDist = dist(chsel(digSel));
                [Pxx,hz] =pwelch(squeeze(dataStore(:,:,1)),pmint*fsamp,0,fsamp,fsamp);
                Pxx = mean(Pxx,2);
                plot(hz,Pxx,ls); xlim([4 48]);
                PxxS(:,1,epc) = Pxx;
                hold on
                title('Accepted')
            end
            
            if fresh
                acceptflag = acceptflagDig(digSelMax);
                AcceptSave(q,epc) = acceptflag;
            else
                load([R.path.localDatPath '/Generic/AcceptSaveFingerTap'],'AcceptSave')
                acceptflag = AcceptSave(q,epc);
            end
            pause(0.2)
            
            mkdir([R.path.localDatPath '/Generic/' subject '/'])
            save([R.path.localDatPath '/Generic/' subject '/' subject '_' typeCodes{epc} '_Data'],'dataStore','chDist','acceptflag','chSelOverall','stidnSave')
            stidnSaveGroup{q} = stidnSave; % save so you can get group stat
        end
        clf;
    end % if no channels
end % loop through subjects

X = [stidnSaveGroup{sum(AcceptSave==1,2)==4}]; % select only the accepted
medianResponseTime = median(X);
if fresh
    save([R.path.localDatPath '/Generic/AcceptSaveFingerTap'],'AcceptSave','medianResponseTime')
end

