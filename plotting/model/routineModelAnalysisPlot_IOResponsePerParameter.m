function routineModelAnalysisPlot_IOResponsePerParameter(R)
cmap = linspecer(5);

close all
R.out.dag = 'dataprep_beta';
[mergeLabels] = loadABCGeneric(R,{'mergeLabels'});
R.out.dag = 'model_ParCorrAnalysis';
PerturbResult = loadABCGeneric(R,{'PerturbResult'});

MI = PerturbResult.corrC;
MI_ds = PerturbResult.corrC_ds;

EI = PerturbResult.EI;
dEI = PerturbResult.dEI;
ampDur = PerturbResult.ampDur;
parMod = PerturbResult.par;
parnames = PerturbResult.parnames;

parSelName = {'int MMC G1','int MMC G5','int MMC G7','int MMC G9','int MMC G11'};
% MP self, MP->II, SP self, II->DP, SP->DP

parSel = [3 1];
% replicate across dimension
for parN = 1:2
    for typ = 1:4
        figure(1)
        subplot(3,4,(parN-1)*4 + typ)
        yyaxis left
        MISel = squeeze(MI(:,:,typ,parSel(parN)));
        Y = mean(MISel); YS = std(MISel)./sqrt(size(MISel,1));
        [hl, hp] = boundedline(parMod,Y,YS);
        hl.Color = cmap(typ,:);
        hl.LineWidth = 1.25;
        hp.FaceColor = cmap(typ,:);
        hp.FaceAlpha = 0.3;
        if parN<3
            axis([-3 3 -0.05 0.4])
        else
            axis([-3 3 -0.05 0.5])
        end
        if typ == 1
            ylabel('IO fidelity')
        end
        if parN == 1
            title(mergeLabels{typ})
        end
        if parN == 3
            xlabel('gain')
        end
        gca; a.YColor = 'k'; 
        
        yyaxis right
        ampDurSel = squeeze(ampDur(2,:,:,typ,parN));
        Y = mean(ampDurSel); YS = std(ampDurSel)./sqrt(size(ampDurSel,1));
        E = errorbar(parMod(1:3:25),Y(1:3:25),YS(1:3:25));
        E.LineStyle = 'none';
        E.Marker = 'o';
        E.CapSize = 0;
        E.MarkerFaceColor = cmap(typ,:).*0.9;
        E.MarkerSize = 2.5;
        E.Color = cmap(typ,:).*0.9;
        axis([-3 3 120 200])
        if typ == 4
            ylabel('Burst dur')
        end
        a = gca;
        a.YTick = [150 200];
        box off
        gca; a.YColor = 'k';
        
        figure(2)
        subplot(3,4,(parN-1)*4 + typ)
        MISel = squeeze(MI(:,:,typ,parSel(parN)));
        Y = mean(MISel); YS = std(MISel)./sqrt(size(MISel,1));
        [hl, hp] = boundedline(parMod,Y,YS);
        hl.Color = cmap(typ,:);
        hl.LineWidth = 1.25;
        hp.FaceColor = cmap(typ,:);
        hp.FaceAlpha = 0.3;
        if parN<3
            axis([-2 2 -0.05 0.4])
        else
            axis([-2 2 -0.05 0.5])
        end
        if typ == 1
            ylabel('IO fidelity')
        end
        if parN == 1
            title(mergeLabels{typ})
        end
        if parN == 3
            xlabel('gain')
        end        
        gca; a.YColor = 'k';
        
        
        yyaxis right
        EISel = squeeze(EI(:,:,typ,parSel(parN)));
        Y = mean(EISel); YS = std(EISel)./sqrt(size(EISel,1));
        E = errorbar(parMod(1:3:25),Y(1:3:25),YS(1:3:25));
        E.LineStyle = 'none';
        E.Marker = 'o';
        E.CapSize = 0;
        E.MarkerFaceColor = cmap(typ,:).*0.9;
        E.MarkerSize = 2.5;
        E.Color = cmap(typ,:).*0.9;
        if typ == 4
            ylabel('E:I ratio')
        end
        a = gca;
        box off;        gca; a.YColor = 'k';
        ylim([0 4])
        xlim([-2 2])
    end
end
figure(1)
set(gcf,'Position',[ 527.4000  261.0000  872.8000  501.0000])
figure(2)
set(gcf,'Position',[ 527.4000  261.0000  872.8000  501.0000])