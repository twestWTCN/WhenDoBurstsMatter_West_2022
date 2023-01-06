function [MI,corrC,MI_ds,corrC_ds,SP,EI,dEI] = computeSyncDesyncMI(RMod,MMod,parSet)
uc1 = innovate_timeseries(RMod,MMod,parSet);
[modelError,~,feat_sim{1},xsims1,xsims_gl] = computeSimData_160620(RMod,MMod,uc1,parSet,0,0);

%% Compute EI
pQ = getModelPriorsNH_emp_nonGains(MMod);
D = recallIntDelayTableWS_emp(RMod,parSet,MMod);
D = D{1};

Q = parSet.int{1};
G = pQ(1).G.*exp(Q.G);
%     Mn = pQ(1).Mn.*exp(Q.Mn);
%     Sn = pQ(1).Sn.*exp(Q.Sn);
%     Bn = pQ(1).Bn.*exp(Q.Bn);
x = xsims1{1}';
 
% % G(:,1)  mp -> mp (-ve self)  4
% % G(:,2)  mp -> sp (+ve rec )  4
% % G(:,3)  ii -> mp (-ve rec )  4
% % G(:,4)  ii -> ii (-ve self)  4
% % G(:,5)  mp -> ii (+ve rec )  4
% % G(:,6)  dp -> ii (+ve rec )  2
% % G(:,7)  sp -> sp (-ve self)  4
% % G(:,8)  sp -> mp (+ve rec )  4
% % G(:,9)  ii -> dp (-ve rec )  2
% % G(:,10) dp -> dp (-ve self)  1
% % G(:,11) sp -> dp (+ve rec)  2
% % G(:,12) ii -> sp (-ve rec)  4
% % G(:,13) sp -> ii (+ve rec)  4
% % G(:,14) dp -> sp (+ve rec)  2
% 
% ei_tmp = mean(abs(xsims1{1}),2);
% EI = sum(ei_tmp([1 2 4]))./ei_tmp(3);

% Get base sync and desync
dataX = bandpass(xsims1{1}(4,:),RMod.obs.trans.bursts.frq,1/RMod.IntP.dt);
dataOrg = dataX;

dataXN = (dataX-mean(dataX))./std(dataX);
XH = abs(hilbert(dataXN));
BInd = find(XH(1:end-1)>prctile(XH,75));
baseSyncs = zeros(1,size(uc1{1}{1},1));
baseSyncs(BInd) = ones(size(BInd));

BInd = find(XH(1:end-1)<prctile(XH,75));
baseDeSyncs = zeros(1,size(uc1{1}{1},1));
baseDeSyncs(BInd) = ones(size(BInd));

BInd = BInd(BInd>max(D)); % ensure all are longer than delay


% Within burst
E_wib = G(11).*abs(x(BInd-D(11),2));
I_wib = G(10).*abs(x(BInd-D(10),4)) + abs(G(9).*x(BInd-D(9),3));
EI_wib = mean(E_wib)./mean(I_wib);

% Out of burst
BInd_oob = setdiff(100:size(x,1),BInd);
E_oob = G(11).*abs(x(BInd_oob-D(11),2));
I_oob = G(10).*abs(x(BInd_oob-D(10),4)) + abs(G(9).*x(BInd_oob-D(9),3));
EI_oob = mean(E_oob)./mean(I_oob);

dEI = EI_wib-EI_oob;
EI = mean(I_wib)./mean(E_oob);

%% This does the synchronizations
% Get some random intervals
gateInt = 0.2 + 0.15.*randn(1,12000);
breakInt = 0.8 + 0.4.*randn(1,12000);
uvec = 0; k = 1; dt = RMod.IntP.dt;
sclx = 5;
while size(uvec,2)<=size(uc1{1}{1},1)
    uvec = [uvec ones(1,fix(breakInt(k)./dt)) repmat(sclx,1,fix(gateInt(k)./dt))];
    k = k + 1;
end
uvec  = uvec(1,1:size(uc1{1}{1},1));

% Simulate model
LH = 1;
L = 1;
uc2 = uc1;

uc2{1}{1}(:,LH(L)) = uc1{1}{1}(:,LH(L)).*uvec';
% Setup parameters
[modelError,~,feat_sim{2},xsims2,xsims_gl] = computeSimData_160620(RMod,MMod,uc2,parSet,0,0);

dataX = bandpass(xsims2{1}(4,:),RMod.obs.trans.bursts.frq,1/RMod.IntP.dt);
dataXN = (dataX-mean(dataX))./std(dataX);
XH = abs(hilbert(dataXN));
BInd = find(XH(1:end-1)>prctile(XH,75));
BurstDig = zeros(size(uvec));
BurstDig(BInd) = ones(size(BInd));

fsamp = 1/RMod.IntP.dt;
minbs = (2/25)*fsamp;

[~,~,~,amp] = burstAmpHist(dataX,fsamp,RMod.data.feat_xscale{2},minbs);
SP(1) = mean(amp);
[~,~,~,dur] = burstDurHist(dataX,fsamp,RMod.data.feat_xscale{3},minbs);
SP(2) = mean(dur);


MI = nan; %% YOU NEED A TOOLBOX for this mutualinfo(uvec,BurstDig);
corrC = corr(uvec',BurstDig','type','Spearman');
% %                 figure(201)
% %                 p(1) = plot(RMod.IntP.tvec,uvec);
% %                 hold on
% %                 p(2) = plot(RMod.IntP.tvec,dataOrg(1:end-1));
% %                 p(3) = plot(RMod.IntP.tvec,dataX(1:end-1));
% %                 p(4) = plot(RMod.IntP.tvec,BurstDig);
% % 
% % 
% %                 a(L) = gca;
% %                 a(L).PlotBoxAspectRatioMode = 'auto';
% %                 xlim([15.5 18.5]); ylim([-38 38])
% %                 set(gcf,'Position',[305.8000  335.4000  940.8000  344.8000])
% %                 legend('uvec','org','pert','burst')
% % end

%% This does the antisynchronizations
% Get some random intervals
gateInt = 0.2 + 0.15.*randn(1,12000);
breakInt = 0.8 + 0.4.*randn(1,12000);
uvec = 0; k = 1; dt = RMod.IntP.dt;
sclx = 0.01;
while size(uvec,2)<=size(uc1{1}{1},1)
    uvec = [uvec ones(1,fix(breakInt(k)./dt)) repmat(sclx,1,fix(gateInt(k)./dt))];
    k = k + 1;
end
uvec  = uvec(1,1:size(uc1{1}{1},1));

% Simulate model
LH = 1; %[2 1 4]; % i.e., middle only (main thalamic input)
L = 1;
% for L = 1:3
uc2 = uc1;

uc2{1}{1}(:,LH(L)) = uc1{1}{1}(:,LH(L)).*uvec';
% Setup parameters
[modelError,~,feat_sim{2},xsims2,xsims_gl] = computeSimData_160620(RMod,MMod,uc2,parSet,0,0);

dataX = bandpass(xsims2{1}(4,:),RMod.obs.trans.bursts.frq,1/RMod.IntP.dt);
dataXN = (dataX-mean(dataX))./std(dataX);
XH = abs(hilbert(dataXN));
BInd = find(XH(1:end-1)<prctile(XH,25));
BurstDig = ones(size(uvec));
BurstDig(BInd) = zeros(size(BInd));

MI_ds = nan; %% YOU NEED A TOOLBOX for this mutualinfo(uvec,BurstDig);
corrC_ds = corr(uvec',BurstDig','type','Spearman');
% %                                 figure(201)
% %                                 p(1) = plot(RMod.IntP.tvec,uvec);
% %                                 hold on
% %                                 p(2) = plot(RMod.IntP.tvec,dataOrg(1:end-1));
% %                                 p(3) = plot(RMod.IntP.tvec,dataX(1:end-1));
% %                                 p(4) = plot(RMod.IntP.tvec,BurstDig);
% %
% %
% %                                 a(L) = gca;
% %                                 a(L).PlotBoxAspectRatioMode = 'auto';
% %                                 xlim([15.5 18.5]); ylim([-38 38])
% %                                 set(gcf,'Position',[305.8000  335.4000  940.8000  344.8000])
% % end