function [nb,Apdf,normpdf,amp,burstinds] = burstAmpHistHMM(dataX,fsamp,bins,minbs,winmark,dataC)


dataX = decimate(dataX,5);
XH = (abs(hilbert(dataX)));

if nargin>6
dataC = decimate(dataC,5);
end
fsamp = fsamp/5;
minbs = fix(minbs/5);
%% Resample WinMark
winMarkInds = ceil(find(winmark)/5); %these are the transitions
winMarkRS = zeros(size(dataX)); %these are the transitions
winMarkRS(winMarkInds) = 1;
winmark = winMarkRS;
clear winMarkRS winMarkInds

%% HMM
options = struct();
options.initrep = 5;
options.K = 2;
options.standardise = 0;
options.verbose = 0;
options.Fs = fsamp;
options.useMEX = 1;
options.zeromean = 0;
options.dropstates = 1;
options.order = 0;
% options.DirichletDiag = 1e10; % set very large as we don't have much data

% HMM inference - we only store the Gamma time-course of posterior probabilities
T = length(XH);
[hmm, Gamma_emb,~,vpath] = hmmmar(XH',T,options);

% find ind of high power state
[~,selInd] = max([max(XH(vpath==1))    max(XH(vpath==2))]);
    
burstinds = SplitVec(find(vpath==selInd),'consecutive');

if nargin>4
    if ~isempty(winmark)
        burstinds = edgeCorrect(burstinds,winmark);
    end
else
    winmark = [];
end

segL = 1000*(cellfun('length',burstinds)/fsamp);

burstinds(segL<minbs) = [];
amp = (cellfun(@(a) max(XH(a)),burstinds));


if numel(burstinds)>2
    [Apdf,nb] =  ksdensity(amp,bins);
else
    Apdf = nan;
    nb = nan;
end

if nargout>2
    if numel(amp)>2
        [normpdf] = fitdist(amp,'normal');
    else
        normpdf = nan;
        warning('Not enough data to fit a distribution!')
    end
end

if nargin>6
    plotBurstConstruction(fsamp,dataC,dataX,XH,burstinds)
    figure
    H = histogram(amp,bins,'Normalization','pdf'); hold on
    H.FaceColor = [131 131 131]/256;
    plot(nb,Apdf,'LineWidth',1.5,'Color',H.FaceColor.*0.4)
    
    a = gca;
        a.YTickLabel = {};
   a.Color = 'none';
   box off; axis square
    xlabel('Amplitude (Z)')
    ylabel('pdf')
end