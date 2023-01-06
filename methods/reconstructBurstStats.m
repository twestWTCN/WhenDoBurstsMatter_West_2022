function stat = reconstructBurstStats(feat,support,featInd)
stat = nan(2,numel(feat));
for i = 1:numel(feat)
    if ~isnan(feat{i}{1})
        [ampMeanSTD(1),ampMeanSTD(2)] = pdfExpVal(feat{i}{featInd},support{featInd}');
        stat(:,i) = ampMeanSTD;
    else
        stat(:,i) = nan;
    end
end

stat = nanmean(stat,2);