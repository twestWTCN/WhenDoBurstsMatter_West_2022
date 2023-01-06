function [chsel,desync,snr,acceptflag] = selectDataCriteria(betaSNR,powDsync,type,critSnr,critDsync)
% BetaSNR must be in (dB) 10 log10 SNR. critSnr is the index of the channel on
% which there must be a minimal SNR - usually the movement block
switch type
    case 'powmod' % threshold on SNR then select for max modulation
        critfilt = (betaSNR>critSnr); %Find all samples above SNR
        if size(critfilt,1)<2
            powDsync(~critfilt) = nan; % if only one signal;
        else
            powDsync(any(~critfilt)) = nan;
        end
        [maxds,d] = min((powDsync),[],2);
        d(isnan(maxds)) = nan;
        desync = maxds; % this is the criterion
        if  ~isnan(desync)
            acceptflag = 1;
        else
            acceptflag = 0;
        end
    case 'modpow' % threshold on MOD then select for max SNR
        critfilt = (powDsync<-critDsync); 
        betaSNR(:,~critfilt) = nan;
        [maxsnr,d] = max(betaSNR(critSnr,:));
        d(isnan(maxsnr)) = nan;
        desync = maxsnr; % this is the criterion
        if desync>0
            acceptflag = 1;
        else
            acceptflag = 0;
        end
    case 'maxpow'
         [maxsnr,d] = max(betaSNR);
        desync = maxsnr; % this is the criterion
        if desync>0
            acceptflag = 1;
        else
            acceptflag = 0;
        end    
    case 'maxdsync'
         [maxdsync,d] = min(critDsync);
        desync = maxdsync; % this is the criterion
        if desync>0
            acceptflag = 1;
        else
            acceptflag = 0;
        end          
end

chsel = d;
if ~isnan(d)
    snr = betaSNR(:,d)';
else
    snr = nan(1,size(betaSNR,1));
end
