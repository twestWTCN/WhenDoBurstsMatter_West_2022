function [burstSelIndsout,XF,XEnv,XPhi,epsAmp,segLout,segAmpout] = simpleBurstDefineMMC(X,fsamp,frqz,minper,eps,NSamp)
if nargin<5
    eps = 75;
end
% Filter Data
[dum,bpk,apk] = ft_preproc_bandpassfilter(X, fsamp,frqz,4,'but','twopass');
XF =filtfilt(bpk,apk,X')';

% Analytic Signal
XEnv = []; XPhi = [];
for ch = 1:size(X,1)
    XEnv(ch,:) = abs(hilbert(XF(ch,:)));
    XPhi(ch,:) = angle(hilbert(XF(ch,:)));
    % Set Threshold
    epsAmp(ch) = prctile(XEnv(ch,:),eps);
end

for ch = 1:size(X,1)
    % Define Bursts
    ThreshX = double(XEnv(ch,:) > epsAmp(ch));
    minS = (minper/frqz(1))*fsamp; % set minimum to 3 periods of slowest
    betaBurstInds = SplitVec(find(ThreshX),'consecutive'); % Split up data based upon the target threshold
    
    if eps == 0
        warning('Rewriting the burst definitions as eps = 0')
        for epoch = 1:NSamp(1)
            sind = randi(size(X,2)-NSamp(2),1);
            betaBurstInds{epoch} = sind:sind+NSamp(2)-1;
        end
    end
    
    segL = cellfun('length',betaBurstInds); % Find burst lengths
    burstSelInds = segL>minS; % Select bursts with above min length
    burstSelIndsout{ch} = betaBurstInds(burstSelInds);
    %
    
    segL = cellfun('length',burstSelIndsout{ch}); % Find burst lengths
    segLout{ch} = log10(segL);
    segAmp = cellfun(@(x) mean(XEnv(1:size(X,1),x),2),burstSelIndsout{ch},'UniformOutput',0);
    segAmpout{ch} = [segAmp{:}];
end

plotopexamp = 0;
if plotopexamp
    bind = 9
    clf
    % bind = 5;
    tvec = linspace(0,size(XF,2)/fsamp,size(XF,2));
    plot(tvec,X(4,:))
    
    hold on
    plot(tvec,XF(4,:)-5e-6)
    hold on
    plot(tvec,XEnv(4,:)-5e-6)
    plot([0 500],[epsAmp(4) epsAmp(4)]-5e-6)
    plot(tvec([burstSelIndsout{bind}(1) burstSelIndsout{bind}(1)]),[5e-6 -20e-6])
    plot(tvec([burstSelIndsout{bind}(end) burstSelIndsout{bind}(end)]),[5e-6 -20e-6])
    
    plot(tvec,1e-6.*XPhi(1,:)-10e-6)
    plot(tvec,1e-6.*XPhi(4,:)-10e-6)
    
    xlim(tvec([burstSelIndsout{bind}(1)-(0.35*fsamp) burstSelIndsout{bind}(end)+(0.35*fsamp)]))
    ylim([-1.5e-5 0.3e-5])
    % pause
    
    ind = burstSelIndsout{bind}(1):5:burstSelIndsout{bind}(end);
    polarhistogram(XPhi(1,ind)-XPhi(4,ind),24)
end