function plotNLTransBehaviour(R)
close all
% Load in data
R.out.dag = 'NLTrans';
[durMeanSTD_save,ampMeanSTD_save,intMeanSTD_save,r2ln_save,bp_save,tvec,swtvec,fzSave,acc_save] = loadABCGeneric(R,{'durMeanSTD_save','ampMeanSTD_save','intMeanSTD_save','r2ln_save','bp_save','tvec','swtvec','fzSave','acc_save'})


AmpSmooth = []; DurSmooth = []; r2LnAmpSmooth = []; r2LnDurSmooth = []; bpSubNorm = []; ACCSubNorm = [];
for sub = 1:size(r2ln_save,4)
    r2Sub = squeeze(r2ln_save(:,:,:,sub));
    ampSub = squeeze(ampMeanSTD_save(1,:,sub));
    durSub = squeeze(durMeanSTD_save(1,:,sub));
    bpSub = squeeze(bp_save(:,:,sub));
    accSub =  squeeze(acc_save(:,sub));
    % baseline
    bpSubNorm(:,:,sub) = baselineNormSmooth(bpSub,1:8,0);
    
    ACCSubNorm(:,sub) = baselineNormSmooth(accSub',1:8,5);
    
    [AmpSmooth(:,sub),AmpSmoothSur(:,sub)] = baselineNormSmooth(ampSub,1:4,5);
    
    [DurSmooth(:,sub),DurSmoothSur(:,sub)] =  baselineNormSmooth(durSub,1:4,5);
    [r2LnAmpSmooth(:,sub),r2LnAmpSmoothSur(:,sub)] = baselineNormSmooth(squeeze(median(r2Sub(2,:,:),2))',1:4,5);
    [r2LnDurSmooth(:,sub),r2LnDurSmoothSur(:,sub)] =  baselineNormSmooth(squeeze(median(r2Sub(3,:,:),2))',1:4,5);
end

figure(22)
subplot(3,1,1)
imagesc(swtvec,fzSave{1},squeeze(median(bpSubNorm,3)))
ylim([4 48]); caxis([-4 2]);
c = colorbar;
c.Location = 'eastoutside';
a = gca;
a.YDir = 'normal';
xlabel('Time to movement onset (s)')
ylabel('Frequency')
yyaxis right
plot(tvec,nanmedian(acc_save,2),'k','LineWidth',2)
ylabel('Movement')
axis normal
xlim([-1.5 1])


[clusters, p_values, t_sums, permutation_distribution ] = permutest(AmpSmooth,AmpSmoothSur,false,0.05,1000,2);

subplot(3,1,2)
plot(swtvec(2:end),median(AmpSmooth,2))
hold on
ylabel('Peak burst amp. (Z)')
yyaxis right
[clusters, p_values, t_sums, permutation_distribution ] = permutest(DurSmooth,DurSmoothSur,false,0.05,1000,2);

plot(swtvec(2:end),median(DurSmooth,2))

 ylabel('Peak burst dur. (Z)')
axis normal
xlim([-1.5 1])
xlabel('Time to movement onset (s)')


subplot(3,1,3)
[clusters, p_values, t_sums, permutation_distribution ] = permutest(r2LnAmpSmooth,r2LnAmpSmoothSur,false,0.05,1000,2);

% boundedline(swtvec(2:end),median(r2LnAmpSmooth,2),std(r2LnAmpSmooth,[],2)/sqrt(size(r2LnAmpSmooth,2)*25))
plot(swtvec(2:end),median(r2LnAmpSmooth,2))
% hold on
% plot(swtvec(2:end),median(r2LnAmpSmoothSur,2))

ylabel('R2Ln Amp(Z)')
hold on
plot(swtvec(2:end),median(r2LnDurSmooth,2))
% hold on
% plot(swtvec(2:end),median(r2LnDurSmoothSur,2))

% boundedline(swtvec(2:end),median(r2LnDurSmooth,2),std(r2LnDurSmooth,[],2)/sqrt(size(r2LnDurSmooth,2)*25))
ylabel('R2Ln Dur(Z)')
axis normal
xlim([-1.5 1])
xlabel('Time to movement onset (s)')

function [XN,XR] = baselineNormSmooth(X,ind,smoothFlag)
XB = X(:,ind);
XN = (X-mean(XB,2))./std(XB,[],2);
XNB = XN(ind);


if smoothFlag >0
    XNB = smooth(XNB,smoothFlag);
    XN = smooth(XN,smoothFlag);
end
if nargout==2
    XR = mean(XNB)+ std(XNB).*randn(size(XN)); % surrugate match
end




