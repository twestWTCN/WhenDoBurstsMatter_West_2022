function [powDsync,dsyncFrq] = findDesync(hz,epochPxx,dsInd)
for i = 1:size(epochPxx,2)
    bandhz = hz>=14 & hz<=30;
    hzlist = hz(bandhz);
    pds =100.*(epochPxx(bandhz,i,dsInd(2))-epochPxx(bandhz,i,dsInd(1)))./epochPxx(bandhz,i,dsInd(1));
    [~,ds] = min(pds); % find the frequency with minimum desync
    dsyncFrq(i) = hzlist(ds);
    powDsync(i) = 100.*(mean(epochPxx(bandhz,i,dsInd(2)))-mean(epochPxx(bandhz,i,dsInd(1))))./mean(epochPxx(bandhz,i,dsInd(1)));
end