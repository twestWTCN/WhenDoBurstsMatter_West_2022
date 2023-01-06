function [nb,kpdf,lnpdf,segL] = burstDurHistHMM(dataX,fsamp,bins,minbs,winmark,burstinds,dataC)

dataX = decimate(dataX,5);
XH = (abs(hilbert(dataX)));

if nargin>6
dataC = decimate(dataC,5);
end
fsamp = fsamp/5;

%% Resample WinMark
winMarkInds = ceil(find(winmark)/5); %these are the transitions
winMarkRS = zeros(size(dataX)); %these are the transitions
winMarkRS(winMarkInds) = 1;
winmark = winMarkRS;
clear winMarkRS winMarkInds


% % %% HMM
% % options = struct();
% % options.initrep = 5; % Set to be large as this is a short simulation
% % options.K = 2;
% % options.standardise = 0;
% % options.verbose = 1;
% % options.Fs = fsamp;
% % options.useMEX = 1;
% % options.zeromean = 0;
% % options.dropstates = 1;
% % options.order = 0;
% % % options.DirichletDiag = 1e10; % set very large as we don't have much data
% % 
% % % HMM inference - we only store the Gamma time-course of posterior probabilities
% % T = length(XH);
% % [hmm, Gamma_emb,~,vpath] = hmmmar(XH',T,options);
% % 
% % % find ind of high power state
% % [~,selInd] = max([mean(XH(vpath==1))    mean(XH(vpath==2))]);
% %     
% % burstinds = SplitVec(find(vpath==selInd),'consecutive');

if nargin>4
    if ~isempty(winmark)
        burstinds = edgeCorrect(burstinds,winmark);
    end
end
segL = 1000*(cellfun('length',burstinds)/fsamp);
segL(segL<minbs) = [];
if numel(segL)>2
    [kpdf,nb] =  ksdensity(segL,bins);
else
    kpdf = nan;
    nb = nan;
end

if nargout>2
    if numel(segL)>2
        [raypdf] = fitdist(segL','rayleigh');
        %         [lnpdf] = fitdist(segL','Lognormal');
        [lnpdf] = fitdist(segL','inversegaussian');
    else
        raypdf = nan;
        lnpdf = nan;
        warning('Not enough data to fit a distribution!')
    end
end

if nargin>7
    plotBurstConstruction(fsamp,dataC,dataX,XH,burstinds)
    figure
    H = histogram(segL,bins,'Normalization','pdf'); hold on
    H.FaceColor = [65 138 179]/256;
    plot(nb,kpdf,'LineWidth',1.5,'Color',H.FaceColor.*0.4)
    
    a = gca;
    a.YTickLabel = {};
    a.Color = 'none';
    box off; axis square
    xlabel('Duration (ms)')
    ylabel('pdf')
    xlim([0 800])
end