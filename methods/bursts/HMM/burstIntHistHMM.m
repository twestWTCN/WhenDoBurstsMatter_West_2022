function [nb,Apdf,normpdf,bint] = burstIntHistHMM(dataX,fsamp,bins,minbs,winmark,burstinds)

XH = abs(hilbert(dataX));
dataX = decimate(dataX,5);
fsamp = fsamp/5;

%% Resample WinMark
winMarkInds = ceil(find(winmark)/5); %these are the transitions
winMarkRS = zeros(size(dataX)); %these are the transitions
winMarkRS(winMarkInds) = 1;
winmark = winMarkRS;
clear winMarkRS winMarkInds


% %% HMM
% options = struct();
% options.initrep = 10; % Set to be large as this is a short simulation
% options.K = 2;
% options.standardise = 0;
% options.verbose = 1;
% options.Fs = fsamp;
% options.useMEX = 1;
% options.zeromean = 0;
% options.dropstates = 1;
% options.order = 0;
% options.DirichletDiag = 1e10; % set very large as we don't have much data
% 
% % HMM inference - we only store the Gamma time-course of posterior probabilities
% T = length(XH);
% [hmm, Gamma_emb,~,vpath] = hmmmar(XH',T,options);
% 
% % find ind of high power state
% [~,selInd] = max([mean(XH(vpath==1))    mean(XH(vpath==2))]);
%     
% burstinds = SplitVec(find(vpath==selInd),'consecutive');

if nargin>4
    if ~isempty(winmark)
        burstinds = edgeCorrect(burstinds,winmark);
    end
end
segL = 1000*(cellfun('length',burstinds)/fsamp);

burstinds(segL<minbs) = [];
bint = [];
for i = 1:numel(burstinds)-1
    bint(i) = burstinds{i+1}(1)-burstinds{i}(end) - 2;
end

if numel(burstinds)>2
    [Apdf,nb] =  ksdensity(bint,bins);
else
    warning('Not enough data to fit a distribution!')
    Apdf = nan;
    nb = nan;
end

if nargout>2
    if numel(segL)>2
        [normpdf] = fitdist(segL','normal');
    else
        normpdf = nan;
        warning('Not enough data to fit a distribution!')
    end
end