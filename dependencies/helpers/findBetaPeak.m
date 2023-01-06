function [peakSize,betapeak] = findBetaPeak(hz,Pxx,frqz)
betals = hz>=frqz(1) & hz<=frqz(2);
bhz = hz(betals);
Pxx = Pxx(betals);

[pks,loc] = findpeaks(Pxx);

% get maximum
loc = loc(pks==max(pks));
pks = pks(pks==max(pks));

betapeak = bhz(loc);
peakSize = pks;

% plot(bhz,Pxx);
% hold on
% scatter(bhz(loc),pks,'filled','Marker','^')
% clf
if isempty(loc)
    betapeak = 22;
    peakSize = Pxx(bhz==22);
end
