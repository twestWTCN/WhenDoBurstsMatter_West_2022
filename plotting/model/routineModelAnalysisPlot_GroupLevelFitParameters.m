function routineModelAnalysisPlot_GroupLevelFitParameters(R)
close all
cmap = linspecer(5);

[FitAnalysis,mergeLabels] = loadABCGeneric(R,{'FitAnalysisDataGroupLevel','mergeLabels'});
figure(1)
%% Now Plot Parameter differences from ISI
pMu = FitAnalysis.parStruc.pMu;
pSig = FitAnalysis.parStruc.pSig;
MMod = FitAnalysis.simStruc.mmod{1,1};
RMod = FitAnalysis.simStruc.Rmod{1,1};
pMuMap = [spm_vec(pMu.int{1}.G); spm_vec(pMu.int{1}.C)];
pSigMap = [spm_vec(pSig.int{1}.G); spm_vec(pSig.int{1}.C)];
pfldnm = getParFieldNames(pMu,MMod);

% Reduce parameter set
[~,IA] = intersect(spm_vec(pMu),pMuMap);
pfldnmSel = pfldnm(IA);
pfldnmSel = convertNameListToReadable(RMod,pfldnmSel); % makes parameters more readable
XP = FitAnalysis.parMAPSave(:,2);

X = []; XS = []; 
for typ = 1:4
    %use the raw trials
    X(:,typ) = XP{typ}(pMuMap,:);
    XS(:,typ) = sqrt(XP{typ}(pSigMap,:));
end

[h,p1] = paramZTest(X(:,2)',XS(:,2)',X(:,1)',XS(:,1)');
[h,p2] = paramZTest(X(:,2)',XS(:,2)',X(:,3)',XS(:,3)');
[h,p3] = paramZTest(X(:,2)',XS(:,2)',X(:,4)',XS(:,4)');

subplot(2,1,1)
[B,E] = makeGroupBarPlotError(X,XS,cmap);
ylabel('Parameter scaling');
alpha = 0.05; %/18;
text(B(1).XEndPoints(p1<alpha),X(p1<alpha).*0.75 + sign(X(p1<alpha)),sprintfc('*',p1(p1<alpha)),'HorizontalAlignment','center','FontSize',18,'Color',cmap(1,:).*0.95)
text(B(3).XEndPoints(p2<alpha),X(p2<alpha).*0.75 + sign(X(p2<alpha)),sprintfc('*',p2(p2<alpha)),'HorizontalAlignment','center','FontSize',18,'Color',cmap(3,:).*0.95)
text(B(4).XEndPoints(p3<alpha),X(p3<alpha).*0.75 + sign(X(p3<alpha)),sprintfc('*',p3(p3<alpha)),'HorizontalAlignment','center','FontSize',18,'Color',cmap(4,:).*0.95)
a = gca;
a.XTick = 1:1:size(X,1);
a.XTickLabel = pfldnmSel(1:1:size(X,1));
a.XTickLabelRotation = -45;
a.XAxisLocation = 'Top';
axis([0.5 size(X,1)+.5 -2 3])
axis normal; box off

L = legend(B,mergeLabels);
L.NumColumns = 4;
L.Box = 'off';
L.Position = [0.1423 0.5185 0.7768 0.0476];
set(gcf,'Position',[681.0000   52.2000  559.0000  730.4000])


%% PLOT PARAMETER CORRELATIONS
% Load the correlations
R.out.dag = 'model_ParCorrAnalysis';
corrStoreSave = loadABCGeneric(R,{'corrStoreSave'});
% corrStoreSave has dimension 4 x 10 x 18 x 4
% dim1 = [C p N range]
% dim2 = [peakhz power amp dur int R2 R2LinSpec R2Lamp R2Ldur R2Lint
% dim3 = [parlist]
% dim4 = type
featlist = {'Peak freq','% power','burst amp','burst dur'};

%% Get power as baseline first

minN = 8;

parSelInds = 1:18; % excludes input gains

subplot(2,1,2)
% Plot elements correlating with burst duration and not power/freq
symblis = {'p','o','s','d'}
for feat= 1:4
    %     subplot(3,1,feat)
    N = squeeze(corrStoreSave(3,feat,parSelInds,:));
    COR = squeeze(corrStoreSave(1,feat,parSelInds,:)).*(N>=minN);
    COR(COR==0) = nan;
    [h,critP] = fdr_bh(corrStoreSave(2,feat,:,:),0.01,'pdep','yes')
    CORSig = COR.*(squeeze(corrStoreSave(2,feat,parSelInds,:)<critP));
    CORSig(CORSig==0) = nan;

    % Plot positives
    X = CORSig>0; X = double(X);
    X(X==0) = nan;
    
    ofs = linspace(-.25,.25,4);
    for typ = 1:4
        S(feat) = scatter((1:size(X,1))+ofs(typ),X(:,typ)*feat,40,cmap(typ,:),'filled','Marker',symblis{feat})
        hold on
    end
    
    % Plot negatives
    X = CORSig<0; X = double(X);
    X(X==0) = nan;
    
    for typ = 1:4
        scatter((1:size(X,1))+ofs(typ),X(:,typ)*-feat,50,cmap(typ,:),'filled','Marker',symblis{feat})
        hold on
    end
end
a = gca;
a.XTick = 1:1:size(X,1);
a.XTickLabel = pfldnmSel(1:1:size(X,1));
a.XTickLabelRotation = 45;
a.XAxisLocation = 'Bottom';
a.YTick = [-4 -3 -2 -1  1 2 3 4];
a.YTickLabel = {'Dur','Amp','Pow','Frq','Frq','Pow','Amp','Dur'}; % {'R-' 'R+'};
ylabel('Parameter correlation')
axis([0.5 size(X,1)+.5 -5 5])
axis normal; box off

L2 = legend(S,featlist);
L2.NumColumns = 4;
L2.Box = 'off';
L2.Position = [0.2796 0.4594 0.4708 0.0285];
set(gcf,'Position',[815    19   904   806])

% Support functions
function [h,p,z_vector] = paramZTest(X1,S1,X2,S2)
% this must take mean X and standard deviation S
z_vector= (X1-X2)./ (sqrt(S1.^2 + S2.^2));
p = 2 * normcdf(-abs(z_vector),0,1);
h = p<0.05;
