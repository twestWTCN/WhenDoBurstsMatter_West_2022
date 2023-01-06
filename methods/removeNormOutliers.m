function [z,newlab] = removeNormOutliers(X,LABEL,zthresh)
labList = unique(LABEL);
z = []; newlab=[];
for lb = 1:numel(labList)
    labId = labList(lb);
    XSel = X(LABEL==labId);
    XSel = XSel(abs(XSel)<(prctile(XSel,95)));
    XSel = XSel(abs(XSel)>(prctile(XSel,5)));
    newlab = [newlab; repmat(labId,size(XSel))];
    z = [z; XSel];
end
