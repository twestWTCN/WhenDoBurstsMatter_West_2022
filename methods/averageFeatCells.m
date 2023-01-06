function feat_sim_mean = averageFeatCells(feat_sim)

for feat = 1:size(feat_sim{1},2)
    feat_sim_mean_tmp = [];
    for N = 1:size(feat_sim,2)
        if ndims(feat_sim{N}{feat})>2
            feat_sim_mean_tmp(:,:,:,:,:,N) = feat_sim{N}{feat};
%             feat_sim_mean_tmp(:,N) = squeeze(feat_sim{N}{feat}(1,1,1,2,:));
        else
            feat_sim_mean_tmp(:,N) = feat_sim{N}{feat};
        end
    end
    feat_sim_mean{feat} = squeeze(mean(feat_sim_mean_tmp,ndims(feat_sim_mean_tmp)));
end


