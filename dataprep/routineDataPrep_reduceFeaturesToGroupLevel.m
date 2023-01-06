function routineDataPrep_reduceFeaturesToGroupLevel(R)
cmap = linspecer(5);
close all
[SMod,mergeLabels] = loadABCGeneric(R,{'SMod','mergeLabels'});
K = [0 0 0 0];
for sub = 1:numel(SMod)
    for typ = 1:4 % premove first
        try
            SMod(sub).featStore(typ);
        catch
            continue
        end

        % Unpack data
        if isempty(SMod(sub).featStore(typ).feat_emp)
            continue
        else
            [betaSNR,betapeak] = computeBetaSNR_AbsPeak(SMod(sub).featStore(typ).feat_xscale{1},squeeze(SMod(sub).featStore(typ).feat_emp{1}),[11 30])
            K(typ) = K(typ) + 1;
        end
        feat_emp_mean{typ}{K(typ)} = SMod(sub).featStore(typ).feat_emp;
        feat_xscale = SMod(sub).featStore(typ).feat_xscale;
    end
end

for typ = 1:4 % premove first

    X = vertcat(feat_emp_mean{typ}{:});
    XFeat = [];
    for ft = 1:4
        XM = [];
        for K = 1:numel(X(:,ft))
            if ft == 1
                XM(1,1,1,1,:,K) = X{K,ft};
            else
                XM(:,K) = X{K,ft};
            end
        end

        XFeat{ft} = mean(XM,ndims(XM));
        XSFeat{ft} = std(XM,[],ndims(XM))./sqrt(K);

        if ft == 1
            [betaSNR,betapeak] = computeBetaSNR_AbsPeak(feat_xscale{ft},squeeze(XFeat{ft}(1,1,1,1,:)),[11 30]);
        end
    end
    SAMod.featStore(typ).feat_emp = XFeat;
    SAMod.featStore(typ).feat_xscale = feat_xscale;
    SAMod.peakFrq(typ) = betapeak;

    for ft = 1:4
        subplot(2,2,ft)
        [L(ft,typ),B(ft,typ)] = boundedline(SAMod.featStore(typ).feat_xscale{ft},squeeze(SAMod.featStore(typ).feat_emp{ft}),squeeze(XSFeat{ft} ));
        L(ft,typ).Color = cmap(typ,:);
        L(ft,typ).LineWidth = 2;
        B(ft,typ).FaceAlpha = 0.2;
        B(ft,typ).FaceColor = cmap(typ,:);
        hold on
    end
end
saveABCGeneric(R,{'SAMod'},{SAMod});