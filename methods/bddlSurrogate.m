function [r2ln,bddls,surrFeat,surrBurstWave] = bddlSurrogate(dataStore,n_surr,R,featEmp)
frq =     R.obs.trans.bursts.frq;
fsamp = R.IntP.dt^-1;

[surr,params] = surrogate(dataStore(:)', n_surr, 'IAAFT2', 0, fsamp);

% % Plot example
% % tvec = linspace(0,size(surr,2)/fsamp,size(surr,2));
% % plot(tvec,dataStore(:))
% % hold on
% % plot(tvec,surr(1,:)-8)
% % xlim([10 15])

r2ln = [];
R.data.winmark ={[]};
XXF = nan(n_surr,751);
XXE = nan(n_surr,751);

for i_surr = 1:n_surr
    [durMeanSTDSur,ampMeanSTDSur,intMeanSTDSur,featSur,~,~,wflag] = getBurstStats(R,surr(i_surr,:),frq,fsamp,[]);
    if wflag == 0
        for ft = 1:4
            if ft == 1
                XEmp = squeeze(featEmp{ft}(1,1,1,1,:));
                XSur = squeeze(featSur{ft}(1,1,1,1,:));
            else
                XEmp = squeeze(featEmp{ft});
                XSur = squeeze(featSur{ft});
            end
            r2ln(ft,i_surr) = fxPooledR2(XEmp,XSur);
            bddls(ft,i_surr) = sum((XEmp - XSur).^2)/mean(XSur)^2;
            surrFeatTmp{ft}(i_surr,:) = XSur;
        end
        
    else
        for ft = 1:4
            r2ln(ft,i_surr) = nan;
            bddls(ft,i_surr) = nan;
            surrFeatTmp{ft}(i_surr,:) = nan(size(R.data.feat_xscale{ft}));
        end
    end
    
    % get burst waveforms
    [burstSelIndsout,XF,XEnv] = simpleBurstDefineMMC(surr(i_surr,:),fsamp,frq,3,75);
    [tvec,XXF(i_surr,:),XXE(i_surr,:)] = getTimeLockedBurstWaveform(burstSelIndsout{1},XF,XEnv,[-0.150 0.600],fsamp);
   
end

% average waveforms
surrBurstWave = [tvec; nanmean(XXF); nanmean(XXE)];
for ft = 1:4
    surrFeat{ft} = nanmean(  surrFeatTmp{ft});
end
