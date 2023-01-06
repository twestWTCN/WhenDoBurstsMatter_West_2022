function [tp,Hz,bgrnd] = computeAmpSpec(epochDat,nsfilt,fsample)
tp=zeros(300,1);
wsize = size(epochDat,1);
h_env=hann(wsize); %hann window for edge effects
ns = size(epochDat,2);
for q=1:ns
    tm=abs(fft(h_env.*epochDat(:,q))); % amplitude spectra
    tp=tp+((tm(1:300)./nsfilt(1:300)').^2); % correct for amplifier roll-off
end

Hz = 1:(fsample/2);
Hz = Hz(1:300);
bgrnd = ((nsfilt(1:300)').^2)/ns;
tp = tp/ns;

