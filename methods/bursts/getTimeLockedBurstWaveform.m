function  [tvec,XXF,XXE,XXRaw,bstRaw,bstUC] = getTimeLockedBurstWaveform(burstSelIndsout,XF,XEnv,twin,fsamp,XRaw,UC)
newBsInds = fix(1+(fsamp*twin(1))):fix(1+(fsamp*twin(2)));
tvec = linspace(-150,750,numel(newBsInds));
bstXF = []; bstXE = [];
% get average envelopes
for bst = 1:numel(burstSelIndsout)
    bsInds = burstSelIndsout{bst};
    % realign to last zero crossing
    %     init = bsInds(find(diff(sign(XF(bsInds)))>1,1,'first'));
    % or the peak!
    [~,maxInd] = max(XF(bsInds));
    init = bsInds(maxInd);
    
    newBsInds = fix(init+(fsamp*twin(1))):fix(init+(fsamp*twin(2)));
    tvec = linspace(twin(1)*fsamp,twin(2)*fsamp,numel(newBsInds));
    
    if (newBsInds(end)<numel(XF)) && (~any(newBsInds<=0)) && (newBsInds(end)<numel(XEnv))
        bstXF(bst,:) = XF(newBsInds);
        bstXE(bst,:) = XEnv(newBsInds);
        if nargout>4
            bstRaw(bst,:,:) = XRaw(:,newBsInds)';
            bstUC(bst,:,:) = UC(newBsInds,:);
        end
    else
        bstXF(bst,:) = nan(size(tvec));
        bstXE(bst,:) = nan(size(tvec));
        if nargout>4
            bstRaw(bst,:,:) = nan(size(XRaw,1),size(tvec,2))';
            bstUC(bst,:,:) = nan(1,size(tvec,2),size(XRaw,1));
        end
    end
end

if ~isempty(bstXF)
    XXF = nanmean(bstXF);
    XXE = nanmean(bstXE);
    if nargout>4
        XXRaw = nanmean(bstRaw);
    end
else
    XXF = nan(size(tvec));
    XXE = nan(size(tvec));
    if nargout>4
        XXRaw = nan(size(tvec));
    end
end