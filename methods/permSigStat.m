function permSigStat(X,base,fy,permnum,sigh,ofset,cmap,alpha,minclust)
A = (base);
B = (fy);

B(:,isnan(A(1,:))) = [];
A(:,isnan(A(1,:))) = [];

% rem
[clusters, p_values, t_sums, permutation_distribution ] = permutest(A,B,0,alpha,permnum,1,5);

if numel(clusters)>0
    for cl = 1:numel(clusters)
        if (numel(clusters{cl})>minclust) && p_values(cl)<0.1
            sigX = X(clusters{cl});
            a = sigh + (ofset);
            sigY = repmat(a,1,numel(sigX));
            plot(sigX,sigY,'Color',cmap,'LineWidth',2)
        end
    end
end