function R = projectABCConfigurations(R)
%
%% DATA SPECIFICATION
R.data.datatype{1} = 'CSD'; %%'NPD'
R.data.datatype{2} = 'ENVPDF';
R.data.datatype{3} = 'DURPDF';
R.data.datatype{4} = 'INTPDF';

R.frqz = [8:.2:38];
% % R.frqzfull = [1:.2:200]; % legacy (only needed for SPM Structured Innovations)
R.data.feat_xscale{2} = 0:0.1:8;
R.data.feat_xscale{3} = 10:25:1400;
R.data.feat_xscale{4} = 10:25:2500;
R.nmsim_name = {'MMC','THAL'}; %modules (fx) to use. These must match those listed in the fx_compile function
R.chdat_name = {'SP'}; % observed channels (redundant)
% R.datinds = 1:4; % Specify this when you deal with the data - ensure its not the wrong order!
R.chsim_name = {'MP','SP','II','DP','THAL'}; % simulated channel names (must be same as list of models)
R.datinds = 1; % Maps from data to Sim
R.siminds = 1; % Maps from sim output to data i.e., should be the chsimname that maps to chdat_name
R.condnames = {'LIN'};
% Spectral characteristics
R.obs.csd.df = 0.5;
R.obs.csd.reps = 32; %96;

%% INTEGRATION
% Main dynamics function
R.IntP.intFx = @ABC_fx_compile_150121;
R.IntP.dt = .001;
R.IntP.Utype = 'white_covar_perPop'; %OR 'white_covar' OR 'DCM_Str_Innov'
R.IntP.bufferExt = ceil(0.050*(1/R.IntP.dt)); % buffer for extrinsic delays
R.IntP.bufferInt = ceil(0.01*(1/R.IntP.dt)); % buffer for  intrinsic delays

N = R.obs.csd.reps; % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.obs.SimOrd = 8;
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

dfact = fsamp/(2*2^(R.obs.SimOrd));
fprintf('The target simulation df is %.2f Hz',R.obs.csd.df);
fprintf('The actual simulation df is %.2f Hz',dfact);

%% OBSERVATION
% observation function (you reset lots of these within the fitting scripts
R.obs.obsFx = @observe_data;
% % R.obs.gainmeth = {'unitvar','obsnoise'}; % NOW set within the mainfitting functions
R.obs.glist =0; %linspace(-5,5,12);  % gain sweep optimization range [min max listn] (log scaling)
R.obs.brn =2; % 2; % burn in time
LF = [1]*10; % Fit visually and for normalised data
R.obs.LF = LF;
R.obs.Cnoise = [0.2]; % Sensor Noise SNR prior i.e. 1/x signal to noise ratio

% Data Features
% fx to construct data features
R.obs.transFx = @constructGenCrossMatrix;
% These are options for transformation (NPD)
R.obs.logscale = 0;
R.obs.trans.zerobase = 1;
R.obs.trans.norm = 1;
R.obs.trans.normcat = 0;
R.obs.trans.logscale = 0;
R.obs.trans.logdetrend = 0;
R.obs.trans.gauss3 = 0;
R.obs.trans.gausSm = 0; % This is off but is switched on to 1 Hz at data processing stage
R.obs.trans.interptype = 'linear';
R.obs.trans.npdscalar = 1; % optional NPD rescaling
%% OBJECTIVE FUNCTION
R.objfx.compFx = @compareData_031121;
R.objfx.errorFx = @fxPooledR2;
R.objfx.feattype = 'magnitude'; %%'ForRev'; %
R.objfx.specspec = 'npd'; %'npd'; %%'auto'; % which part of spectra to fit
R.objfx.featweight = [1 1 1 1];
%% OPTIMISATION
% % R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.DInt','.int{src}.S','.int{src}.C','.obs.Cnoise'}; 
R.SimAn.pOptBound = [-12 12];
R.SimAn.pOptRange = R.SimAn.pOptBound(1):.1:R.SimAn.pOptBound(2);
R.SimAn.searchMax = 200;
R.SimAn.convIt.dEps = 7e-4; %3e-2;
R.SimAn.convIt.eqN = 5;
R.analysis.modEvi.N  = 500;
R.SimAn.scoreweight = [1 0]; %1e-8];
R.SimAn.rep = 512; % Repeats per temperature
% R.SimAn.saveout = 'xobs1';
R.SimAn.jitter = 1; % Global precision
R.SimAn.minRankLambda = 2;
%% PLOTTING
R.plot.outFeatFx = @genplotter_200420; 
R.plot.save = 'False';
R.plot.distchangeFunc = @plotDistChange_KS;

R.plot.feat(1).axlim = [0 45 0 5];
R.plot.feat(1).axtit = {'Freq (Hz)','Norm Power'};
R.plot.feat(2).axlim = [0 6 0 1.6];
R.plot.feat(2).axtit = {'Envelope amplitude (Norm.)','p.d.f.'};
R.plot.feat(3).axlim = [0 1000 0 0.015];
R.plot.feat(3).axtit = {'Burst duration (ms)','p.d.f.'};
R.plot.feat(4).axlim = [0 1000 0 0.015];
R.plot.feat(4).axtit = {'Burst interval (ms)','p.d.f.'};






