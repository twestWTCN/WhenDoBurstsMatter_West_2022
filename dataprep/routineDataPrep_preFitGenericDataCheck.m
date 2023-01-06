function routineDataPrep_preFitGenericDataCheck(R)
close all
closeMessageBoxes;
cmap = linspecer(5);
% % If looping through burst thresholds
% % epsList = [50 60 70 75 80 85 90 95];
% % for eps = 1:numel(epsList)
% %     R.out.dag  = ['dataprep_beta' num2str(eps)];
% Otherwise just do once
for eps = 1
    R.out.dag  = 'dataprep_beta';
    [bigStore,subjectList,mergeLabels,featBank] = loadABCGeneric(R,{'bigStore','subjects','mergeLabels','featBank'});
    bigSub = unique(subjectList);
    
    typKeepStore = [];
    R.plot.holdop = 1; % hold plots
    S = 0;
    for sub = 1:size(bigSub,1)
        close all;
        peakFrq = []; featStore = []; typSubKeep = [];
        for typ = 1:numel(mergeLabels)
            % Unpack data
            typSub = find(((bigStore.merge == typ) + strcmp(bigSub(sub),subjectList))==2); % this finds which column of bigstore/feat you need.
            if numel(typSub)>1
                [~,ind] = max(bigStore.betaSNR(typSub));
                typSub = typSub(ind);
            end
            typSubKeep = [typSubKeep; typSub];
            
            if isempty(typSub) && typ== 2
                warning('Cant load data')
                continue
            end
            if isempty(typSub) && typ== 1
                NoISIFlag = 1;
            else
                NoISIFlag = 0;
            end
            
            if ~isempty(typSub)
                peakFrq(typ) = bigStore.betapeak(typSub);
                subID = bigStore.sub(typSub);
                for feat = 1:numel(R.data.datatype)
                    if feat == 1
                        R.data.feat_emp{feat}(1,1,1,1,:) = featBank{typ,feat}(:,subID);
                    else
                        R.data.feat_emp{feat} = featBank{typ,feat}(:,subID);
                    end
                end
                R.data.feat_xscale{1} = R.frqz;
                R.plot.featcolor = cmap(typ,:);
                R.plot.outFeatFx({R.data.feat_emp},{},R.data.feat_xscale,R,1,[])
                
                featStore(typ).feat_xscale = R.data.feat_xscale;
                featStore(typ).feat_emp = R.data.feat_emp;
            end
        end
        figure(11); shg
        S = S+1;
        SMod(S).sub = bigSub{sub};
        SMod(S).ID = subID;
        SMod(S).NoISIFlag = NoISIFlag;
        SMod(S).peakFrq = peakFrq;
        SMod(S).featStore = featStore;
        typKeepStore = [typKeepStore; typSubKeep];
        drawnow; pause(0.2)
    end
    bigStoreEdit = bigStore(typKeepStore,:);
    saveABCGeneric(R,{'SMod','bigStoreEdit'},{SMod,bigStoreEdit})
end