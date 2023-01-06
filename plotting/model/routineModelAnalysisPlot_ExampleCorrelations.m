function routineModelAnalysisPlot_ExampleCorrelations(R)
close all
R.out.dag = 'model_ParCorrAnalysis';
corrExampSave = loadABCGeneric(R,{'corrExampSave'});
[simTracesDet,simTracesStoc] = loadABCGeneric(R,{'simTracesDet','simTracesStoc'});
cmap = linspecer(5);
for typ = 1:4
    for parN = 1:3
        %         plot traces
        figure(1+parN)
        ofset = 8;
        i = 0;
        for alphai = [5 13 20]
            i = i +1;
            tn = size(simTracesStoc{parN}{alphai,typ},2);
            tvec = linspace(0,tn*R.IntP.dt,tn);
            subplot(2,4,2+((typ-1)*2))
            X = simTracesStoc{parN}{alphai,typ}(4,:)';
            X = normalize(X);
            plot(tvec,X-(ofset*i)'); xlim([3 5])
            hold on
            
            subplot(2,4,1+((typ-1)*2))
            plot(R.frqz,corrExampSave.feat_sim_avRepPar{parN}{alphai,typ}{1}(1,:));
            hold on
        end
        
        figure(parN*10)
        parMod = linspace(-2,2,25);
        %(:,rep,alphai,typ) - [frq amp dur int]
        statParPar = corrExampSave.statParPar;
        lyapPar = corrExampSave.lyapPar;
        
        % Plot 1
        subplot(3,4,typ)
        XEq = corrExampSave.XEq{parN}(:,typ);
        REig =  corrExampSave.lambdaRealXEq{parN}(:,:,typ);
        IEig =  corrExampSave.lambdaImagXEq{parN}(:,:,typ);
        Eig = corrExampSave.eigSave{parN}(:,:,typ);
        
        eigfrq = angle(Eig);
        for i = 1:numel(XEq)
            scatter(repmat(parMod(i),1,size(XEq{i},2)),XEq{i},30,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',cmap(typ,:))
            hold on
        end
        box off;
        ylabel('X Equilibrium'); axis([-3 3 0 35])
        
        figure(50+parN)
        subplot(2,4,typ)
        plot(parMod,REig','Color',cmap(typ,:),'LineWidth',1.5)
        ylabel('Real part of LE')
        axis([-3 3 -0.35 0])
        box off;
        subplot(2,4,typ+4)
        plot(parMod,IEig','Color',cmap(typ,:),'LineWidth',1.5)
        ylabel('Imag part of LE')
        axis([-3 3 -0.1 0.1])
        box off;
        
        figure(parN*10)
        % Plot 3
        subplot(3,4,typ+4)
        for feat = 2:3
            if feat == 2
                yyaxis right
            else
                yyaxis left
            end
            Y = squeeze(statParPar{parN}(feat,:,:,typ));
            X = repmat(parMod,size(Y,1),1);
            
            YM = mean(Y,1);
            YS = std(Y,1)./sqrt(size(Y,1));
            S = scatterErrorPlot(parMod,YM,YS, cmap(typ,:));
            if feat == 3
                S.Marker = 'o';
                ylabel('Burst Duration'); axis([-3 3 120 220])
            else
                S.Marker = 's';
                S.MarkerFaceColor = cmap(typ,:)*0.5;
                ylabel('Burst Amplitude'); axis([-3 3 2.2 2.8])
            end
            hold on
        end
        % Plot 3
        r2linParSave = corrExampSave.r2linParSave;
        subplot(3,4,typ+8)
        Y = squeeze(r2linParSave{parN}(3,:,:,typ));
        X = repmat(parMod,size(Y,1),1);
        YM = mean(Y,1);
        YS = std(Y,1)./sqrt(size(Y,1));
        scatterErrorPlot(parMod,YM,YS, cmap(typ,:));
        hold on
        ylabel('IAFFT2 Fit'); axis([-3 3 0.4 1])
        
        figure(1+parN)
        set(gcf,'Position',[ 552         453        1184         525]);
        figure(parN*10)
        set(gcf,'Position',[ 624          54        1108         798]);
        figure(50+parN)
        set(gcf,'Position',[  629         371        1117         442])
    end
end
function S = scatterErrorPlot(X,YM,YS,cmap)

S = scatter(X,YM,30);
S.MarkerFaceAlpha = 0.75;
S.MarkerFaceColor = cmap;
S.MarkerEdgeColor = 'none';
hold on;
if ~isempty(YS)
    E = errorbar(X,YM,YS);
    E.LineStyle = 'none';
    E.Marker = 'o';
    E.Color = cmap;
    E.MarkerFaceColor = 'none';
    E.MarkerEdgeColor = 'none';
    E.CapSize = 0;
    E.LineWidth = 1.5;
end
box off