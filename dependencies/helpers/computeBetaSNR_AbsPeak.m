function [betaSNR,betapeak,peakSize] = computeBetaSNR_AbsPeak(hz,Pxx,frqz)
betals = hz>=frqz(1) & hz<=frqz(2);
[~,bp] = max(Pxx(betals));

[peakSize,betapeak] = findBetaPeak(hz,Pxx,[12 30]); % beta band is slightly lower to allow cross over into alpha when peak freq is quite low
betaSNR = 10.*log10(mean(Pxx(hz>=betapeak-3 & hz<=betapeak+3))./mean(Pxx(hz>48 & hz<98)));


