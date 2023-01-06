function routineModelComp_compareModelsFit(R,WML)
if nargin<2
    spreadSession(R);
else
    spreadSession(R,WML);
end
closeMessageBoxes;
dttag = {'spec','all'};
close all
%Typelist
%% Sublist
[SAMod,mergeLabels] = loadABCGeneric(R,{'SAMod','mergeLabels'});

%% Get list of data
typlist = [2 1 3 4];

modFit = 0
if modFit
    typ =2
    for mod = 1:2
        % Unpack data
        peakFrq = SAMod.peakFrq(typ);
        R.data.feat_xscale = SAMod.featStore(typ).feat_xscale;
        R.data.feat_emp = SAMod.featStore(typ).feat_emp;
        X =R.data.feat_xscale{1}; Y = squeeze(R.data.feat_emp{1})';
        plotop = 1;

        R.data.feat_emp{1}(1,1,1,1,:) = reduceSpectraToGaussians(X,Y,3,[1 34],0.2,plotop);
        close all
        % Remove previously loaded fits
        try; R = rmfield(R,'Mfit'); end

        %Model Spec
        R.nmsim_name = {'MMC'}; %modules (fx) to use. These must match those listed in the fx_compile function
        R.chdat_name = {'DP'}; % observed channels (redundant)
        R.chsim_name = {'MP'  'SP'  'II'  'DP'};
        R.siminds = 4;

        switch mod
            case 1
                % Model specification
                modelspec = eval(['@MS_burstABC_MMC_NHModel' num2str(1)]);
                R.IntP.intFx = @ABC_fx_compile_NH;
                % Model Inversion
                R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.Sn','.int{src}.Bn','.int{src}.C','.int{src}.alpha'}; %,'.int{src}.Mn'

            case 2
                % Model specification
                modelspec = eval(['@MS_burstABC_MMC_JRModel' num2str(1)]);
                R.IntP.intFx = @ABC_fx_compile_JR;
                % Model Inversion
                R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S',                '.int{src}.C','.int{src}.alpha'};
        end

        R.obs.gainmeth = {'unitvar','symmetrical','flattail','leadfieldOneSignal'}; %,'leadfieldOneSignal'};
        R.IntP.Utype = 'coloured';

        for dtype = 1:2
            [R,p,m] = modelspec(R); % load model

            %             setup priors
            if typ == 2 % i.e., if PreMove use default priors
                [R,p,m] = modelspec(R); % load model
                if dtype == 2
                    Rtmp = R;
                    Rtmp.out.dag =  ['fitGroupLevelV1_' mergeLabels{typ}  '_' dttag{1} '_modtype' num2str(mod)];
                    [R,m,p] = loadABCData_160620(Rtmp); % empirical Bayes
                end
                R.SimAn.jitter = 1;
            else% use an empirical update
%                 datatmp = R.data;
%                 Rtmp = R;
%                 Rtmp.out.dag = ['fitGroupLevelV1_' mergeLabels{2} '_' dttag{dtype}];
%                 [R,m,p] = loadABCData_160620(Rtmp); % empirical Bayes
% 
%                 % Reshape Opt List so that priors match
%                 [pInd,pMu,pSig] = parOptInds_110817(R,p,m.m); % in structure form
%                 R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.C'}; %,'.int{src}.alpha'}; %,'.int{src}.Mn',,'.int{src}.Bn',
%                 [pIndNew,pMuNew,pSigNew] = parOptInds_110817(R,p,m.m); % in structure form
%                 [~,IA] = intersect(spm_vec(pMu),spm_vec(pMuNew));
%                 R.Mfit.Sigma = R.Mfit.Sigma(IA,IA);
%                 R.Mfit.Mu = R.Mfit.Mu(IA,:);
%                 R.Mfit.xf = R.Mfit.xf(IA,:);
%                 R.Mfit.ks = R.Mfit.ks(IA,:);
%                 R.Mfit.Rho = R.Mfit.Rho(IA,IA);%.*(~eye(numel(IA),numel(IA)).*0.5)+eye(numel(IA),numel(IA));
%                 R.data = datatmp; % As you reload 'R' - make sure data is replaced with task specficic features
%                 clear Rtmp datatmp;
%                 disp('Loading PreMove priors (spec fit only)')
            end

            % Weight Features Differently
            if dtype == 1
                R.objfx.featweight = [1 0 0 0];
            elseif dtype == 2
                R.objfx.featweight = [1 1 1 0];
            end

            % Data transforms
            R.Bcond = 1;
            R.obs.trans.gausSm = 0;
            R.SimAn.minRankLambda = 2.5;
            R = setSimTime(R,48);
            R.SimAn.scoreweight = [1 1e-4];
            R.SimAn.rep = 512;

            % set repetitions (for stochastics)
            if dtype == 1
                R.SimAn.RealzRep = 1; R.SimAn.convIt.dEps = 1e-4; % spectra
            else
                R.SimAn.RealzRep = 4; R.SimAn.convIt.dEps = 1e-5; %7.5e-5; % bursts
            end

            frq = peakFrq;
            R.obs.trans.bursts.frq = [frq-3 frq+3]; % this is baked into KMFingerTapsJoin - dont change!
            R.plot.outFeatFx({R.data.feat_emp},{},R.data.feat_xscale,R,1,[])
            R.out.dag = ['fitGroupLevelV1_' mergeLabels{typ}  '_' dttag{dtype} '_modtype' num2str(mod)];

            try
                %                 SimAn_ABC_201120_MRealz(R,p,m);
                SimAn_ABC_201120(R,p,m);
                %   R.IntP.intFx
                %   R.SimAn.pOptList
                % a= 1;
            catch
                warning('Fitting failed')
                try
                    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ERLIST'])
                    ABCErrorFX(R)
                catch
                    ABCErrorFX(R)
                end
            end

        end
        %         disp([sub typ dtype])
    end
end

R.analysis.dagtype = 'arbitrary';
daglist = {
    ['fitGroupLevelV1_' mergeLabels{2}  '_' dttag{1} '_modtype1'],...
    ['fitGroupLevelV1_' mergeLabels{2}  '_' dttag{2} '_modtype1'],...
    ['fitGroupLevelV1_' mergeLabels{2}  '_' dttag{1} '_modtype2'],...
    ['fitGroupLevelV1_' mergeLabels{2}  '_' dttag{2} '_modtype2']
    };


modnames = {'WC Spec','WC All','JR Spec','JR All'}
R.objfx.featweight = [1 1 1 0];
R.analysis.comptype = 'switch';
% modID = modelCompMaster_160620(R,1:4,[],daglist)

R.plot.featsel = [1 1 1 0]; % these are the feats to be printed
R.modcomp.modN = [1:4];
R.modcompplot.NPDsel = 1:4
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
cmap = cmap(end:-1:1,:);
R.data.datatype = R.data.datatype(1:3);
plotModComp_310520BURSTSPECIFIC(R,cmap,daglist,modnames)

set(gcf,'Position',[913    41   751   721])




% % %% Reset priors
% % for mod = [1 3]
% %     % Wilson Cowan
% %     % retrieve WC prior
% %     Rtmp = R;
% %     Rtmp.out.dag = daglist{mod};
% %     % Load Config
% %     RMod = loadABCData_160620(Rtmp);
% %     orgPrior = RMod.Mfit.prior;
% %
% %     % Update Prior
% %     Rtmp.out.dag = daglist{mod+1};
% %     RMod = loadABCData_160620(Rtmp);
% %     RMod.Mfit.prior = orgPrior; % Remember this isnt the real prior used in fitting, as it was an empirical update
% %     [KL,DKL] = KLDiv(RMod,RMod.Mfit,[],[],1);
% %     RMod.Mfit.DKL = DKL;
% %     saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' Rtmp.out.dag '\R_' R.out.tag '_' Rtmp.out.dag '.mat'],RMod)
% % end
% % clear Rtmp Rmod KL DKL