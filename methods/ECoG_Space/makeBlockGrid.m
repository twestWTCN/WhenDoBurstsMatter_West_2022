function [blockStack,Xedges,Yedges,blockStackSmth,XedgesSm,YedgesSm] =  makeBlockGrid(X,Y,binX,binY,data)
% plotting
[N,Xedges,Yedges,binXInd,binYInd] = histcounts2(X,Y,binX,binY);

blockStack = zeros(size(data,1),numel(Xedges),numel(Yedges));
blockStackZero = zeros(size(data,1),numel(Xedges),numel(Yedges));
for i = 1:numel(Xedges)
    for j = 1:numel(Yedges)
        inds = (binXInd == i) & (binYInd == j);
        if any(inds)
            blockStack(:,i,j) = mean(data(:,inds),2);
             blockStackZero(:,i,j) = mean(data(:,inds),2);
        else
            blockStack(:,i,j) = nan;
        end
    end
end

% Interpolate
[Xgrid,Ygrid] = meshgrid(Xedges,Yedges);

smfac = 2;
XedgesSm = linspace(Xedges(1),Xedges(end),fix(numel(Xedges)*smfac));
YedgesSm = linspace(Yedges(1),Yedges(end),fix(numel(Xedges)*smfac));
[XgridSm,YgridSm] = meshgrid(XedgesSm,YedgesSm);
for df = 1:size(data,1)
    XS= interp2(Xgrid,Ygrid,squeeze(blockStack(df,:,:))',XgridSm,YgridSm)';
%     XS= interp2(Xgrid,Ygrid,squeeze(blockStackZero(df,:,:))',XgridSm,YgridSm)';
%     XS(XS==0) = nan;
    blockStackSmth(df,:,:) = XS;
end