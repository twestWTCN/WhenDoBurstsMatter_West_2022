function pQ = getModelPriorsNH_emp(m)
for i = 1:m.m
    switch m.dipfit.model(i).source
        case 'MMC'
            pQ(i).G  = [12 12 12 12 12 6 12 12 6 3 6 12 12 6];   % intrinsic connections (mV.s)
            pQ(i).T  = [25 20 15 15];                            % synaptic time constants [mp sp ii dp] (ms)
            pQ(i).C  = [20 20 20 20]; % noise in firing rate (should be on same scale as F s^-1))
            pQ(i).Mn  = [67 64 131 44]; % maximal firing rate
            pQ(i).Sn  = [0.1 0.2 0.4 0.1]; % FI slope
            pQ(i).Bn  = [18 5 20 10];   % baseline firing rate
    end
    % Convert to seconds
    pQ(i).T = pQ(i).T./1000;
    
end
%--------------------------------------------------------------------------
% G(:,1)  mp -> mp (-ve self)  4
% G(:,2)  mp -> sp (+ve rec )  4
% G(:,3)  ii -> mp (-ve rec )  4
% G(:,4)  ii -> ii (-ve self)  4
% G(:,5)  mp -> ii (+ve rec )  4
% G(:,6)  dp -> ii (+ve rec )  2
% G(:,7)  sp -> sp (-ve self)  4
% G(:,8)  sp -> mp (+ve rec )  4
% G(:,9)  ii -> dp (-ve rec )  2
% G(:,10) dp -> dp (-ve self)  1
% G(:,11) sp -> dp (+ve rec)  2
% G(:,12) ii -> sp (-ve rec)  4
% G(:,13) sp -> ii (+ve rec)  4
% G(:,14) dp -> sp (+ve rec)  2