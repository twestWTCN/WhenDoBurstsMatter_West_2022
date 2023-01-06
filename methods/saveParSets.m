function [par_mean,par_sel] = saveParSets(par_rep,inds)
if islogical(inds)
    inds = find(inds);
end
i = 0;
for p = inds
    i = i+1;
    par_sel(:,i) = spm_vec(par_rep(p));
end

par_mean = median(par_sel,2);