function feat_sim_mean = averageFeatCellsBUV2(feat_sim)

for feat = 1:size(feat_sim{1},2)
    feat_sim_mean_tmp = []; L = 0;
    for N = 1:size(feat_sim,2)
        if ~any(isnan(feat_sim{N}{1}))
            if ndims(feat_sim{N}{feat})>2
                L = L+1;
                feat_sim_mean_tmp(:,:,:,:,:,L) = feat_sim{N}{feat};
                %             feat_sim_mean_tmp(:,N) = squeeze(feat_sim{N}{feat}(1,1,1,2,:));
            else
                L = L+1;
                feat_sim_mean_tmp(:,L) = feat_sim{N}{feat};
            end
        end
    end
    feat_sim_mean{feat} = squeeze(mean(feat_sim_mean_tmp,ndims(feat_sim_mean_tmp)));
end


if numel(feat_sim_mean)<2
    feat_sim_mean{1} = nan;
end