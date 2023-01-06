function investigateNLTransBehaviour(R)
subjects= ls([R.path.localDatPath '\Generic']);
subjects = cellstr(subjects);
subjects = subjects(5:end);
fsamp = 1000;
SL = 0;
for sub = 1:size(subjects,1)
    subCodeList = [];
    file1 = [R.path.localDatPath '\Generic\' subjects{sub} '\' subjects{sub} '_TransWindow_FT_Data.mat'];
    file2 = [R.path.localDatPath '\Generic\' subjects{sub} '\' subjects{sub} '_TransWindow_FT_T2_Data.mat'];
    if isfile(file1)
        load(file1,'dataStore','chDist','acceptflag')
    elseif isfile(file2)
        load(file2,'dataStore','chDist','acceptflag')
    else
        dataStore = [];
    end


    %         load dataTmp
    if ~isempty(dataStore) & acceptflag == 1 %& sum(subCodeList==mergeCodes(typ))<2
        SL = SL+1;
        tvec = linspace(-2,(size(dataStore,1)/fsamp)-2,size(dataStore,1));
        try
            acc = squeeze(dataStore(:,:,2));
            acc = normalize(acc);
        catch
            acc = nan;
        end
        ecog = squeeze(dataStore(:,:,1));

        %Start computing data features
        [Pxx,hz] =pwelch(ecog,fsamp/2,0,fsamp,fsamp);
        Pxx = mean(Pxx,2);

        % find SNR and peak frequency
        [betaSN,betapeak] = computeBetaSNR_AbsPeak(hz,Pxx,[14 30]); % do BETA (14-30)
        frq  = [betapeak-3 betapeak+3]; %[21 30];

        % Now compute over sliding window
        % Burst features
        ind = 1:size(ecog,1);
        winsize = 0.7.*fsamp;
        slideInds = slideWindow(ind,winsize,fix(winsize*0.9));
        swtvec = tvec(fix(median(slideInds)));
        for swin = 1:size(slideInds,2)-1
            Rtmp = R;
            X = ecog(slideInds(:,swin),:);
            Rtmp.data.winmark = getWinMark(X);
            XF = X.*hamming(size(X,1));
            [bp(:,swin),fz{swin}] = pwelch(XF(:),fsamp/2,0,fsamp,fsamp);
            X = normalize(X);
            [durMeanSTD(:,swin),ampMeanSTD(:,swin),intMeanSTD(:,swin),feat{swin}] = getBurstStats(Rtmp,X(:),frq,fsamp,0);
            n_surr =25;
            Rtmp.obs.trans.bursts.frq = frq;
            [r2ln(:,:,swin),bddl] = bddlSurrogate(X,n_surr,Rtmp,feat{swin});
            disp(swin)
        end
        acc_save(:,SL) = median(acc,2);
        bp_save(:,:,SL) = bp;
        durMeanSTD_save(:,:,SL) = durMeanSTD;
        ampMeanSTD_save(:,:,SL) = ampMeanSTD;
        intMeanSTD_save(:,:,SL) = intMeanSTD;
        r2ln_save(:,:,:,SL) = r2ln;
        fzSave = fz;
        
        figure
        subplot(3,1,1)
        imagesc(swtvec,fz{1},bp)
        ylim([4 30])
        a = gca;
        a.YDir = 'normal';
        xlabel('Time to movement onset (s)')
        ylabel('Frequency')
        yyaxis right
        plot(tvec,median(acc,2),'k','LineWidth',2)
        ylabel('Movement')
        axis normal
        xlim([-1.5 1.25])

        subplot(3,1,3)
        X = durMeanSTD(1,:);
        plot(swtvec(1:end-1),smooth(X,5)','LineWidth',1.5); hold on
        X = ampMeanSTD(1,:);
        yyaxis right
        plot(swtvec(1:end-1),smooth(X,5)','LineWidth',1.5); hold on


        axis normal

        subplot(3,1,2)
        X = squeeze(median(r2ln,2));
        plot(swtvec(1:end-1),smooth(X(3,:),4)','LineWidth',1.5); hold on
        plot(swtvec(1:end-1),smooth(X(2,:),4)','LineWidth',1.5)
        xlabel('Time to movement onset (s)')
        ylabel('Nonlinearity')
        yyaxis right
        plot(tvec,mean(acc,2))
        axis normal
        xlim([-1.5 1.25])
        ylabel('Movement')
        
        figure
        x = R.data.feat_xscale{3};
        for sw = 1: size(slideInds,2)-1
            if rem(sw,3) == 0
                Y = feat{sw}{3};
                Y = Y.*10;
                %                    Y = Y/10;
                plot(Y+swtvec(sw),x);
                hold on
                ylim([0 400])
            end
        end
        
        
        
        
        
        
        %         r2ln
        %         bp
    end
end
% durMeanSTD_save(:,:,SL) = durMeanSTD;
% ampMeanSTD_save(:,:,SL) = ampMeanSTD;
% intMeanSTD_save(:,:,SL) = intMeanSTD;
% r2ln_save(:,:,:,SL) = r2ln;

R.out.dag = 'NLTrans';
saveABCGeneric(R,{'durMeanSTD_save','ampMeanSTD_save','intMeanSTD_save','r2ln_save','bp_save','tvec','swtvec','fzSave','acc_save'},{durMeanSTD_save,ampMeanSTD_save,intMeanSTD_save,r2ln_save,bp_save,tvec,swtvec,fzSave,acc_save})
