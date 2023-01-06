function routineDataPlot_GroupHDS_SVM(R)
close all
[bigStore,subjects,typeLabels,featBank] = loadABCGeneric(R,{'bigStore','subjects','typeLabels','featBank'});
typeLabels = {'ISI','PreMove','PostMove','Imagery'};

%% Now do the plotting
cmap = linspecer(5);
R.data.feat_xscale{1} = R.frqz;
% Select by SNR
DC = table2array(bigStore(:,8:13));
label = table2array(bigStore(:,2));
label(isnan(DC(:,1))) = [];
DC(isnan(DC(:,1)),:) = [];
subjects(isnan(DC(:,1))) = [];
DCw = prewhiten(DC);

GX = lda(DCw,label,2);
subplot(2,2,3)
[genError,S,ROC,AUCCV] = computePlotSVM(GX,label,1:4,cmap,1);
axis([-4 4 -4 4]); ylabel('LDA 1'); ylabel('LDA 2');
title('SVM classification from burst features')
legend(S,typeLabels)

subplot(2,2,4)
plotROCCurve(ROC,cmap)

DC = table2array(bigStore(:,5:7));
DCw = prewhiten(DC);
label = table2array(bigStore(:,2));

GX = lda(DCw,label,2);
subplot(2,2,1)
[genError,S,ROC,AUCCV] = computePlotSVM(GX,label,1:4,cmap,1);
axis([-4 4 -4 4]); ylabel('LDA 1'); xlabel('LDA 2')
title('SVM classification from burst features')
subplot(2,2,2)
plotROCCurve(ROC,cmap)

