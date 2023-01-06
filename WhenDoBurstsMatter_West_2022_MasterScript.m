 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%      Master Script for Simulations/Figures in West et al. (2022)      %%
%          "When do bursts matter in the primary motor cortex?            %  
%          Investigating changes in the intermittencies of beta           %  
%               rhythms associated with movement states."                 %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Timothy West; MRC Brain Networks Dynamics Unit,Nuffield Department of   %
% Clinical Neurosciences, University of Oxford;                           %
% Wellcome Centre for Human Neuroscience, University College London.      %
% 2020-2022                                                               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Setup Paths
clear; close all
restoredefaultpath
repopath = cd;
cd(repopath)
setGraphicsDefaults; % sets the graphics defaults

% Set up project paths
routname = 'mainproject';
R = setupPaths(repopath,routname);

% Setup configurations for ABC fits
R = projectABCConfigurations(R);
R.out.tag  = 'burstcortex';
R.out.dag = 'dataprep_beta';

%% Data Preparation
% These functions do some basic preprocessing, epoching, and channel
% selection
fingerTapsPrepareData(R,1); % Fingertap prepare data
motorBasicPrepareData(R,1); % Motor basic prepare data
motorImageryPrepareData(R,1); % Motor imagery prepare data
KMGenericPrepareData(R,1); % Merges altogether

%% Data Analysis
routineDataPlot_ExampleTimeSeries(R)
routineDataPlot_GroupGenericFeatures(R)
routineDataPlot_GroupHDS_SVM(R)
routineDataPlot_GroupSurrogateResults(R)

%% Modelling
% Preprocessing
routineDataPrep_preFitGenericDataCheck(R)
% Reduce data features to group level
routineDataPrep_reduceFeaturesToGroupLevel(R)
% Model fitting
routineModelFit_fitGroupLevel(R,[])
% Model Analysis
routineModelAnalysis_analyseGroupLevel(R,1)
routinModelAnalysis_analyseGroupLevelParCorelations(R,1)
routineModelAnalysis_analyseGroupLevelExampleParCorrs(R,1); % this does the example correlations
routineModelAnalysis_analyseIOControlParameter(R,1)

%% Analysis of Models
% Plot outcomes of model analysis
routineModelAnalysisPlot_GroupLevelFeatures(R); % This plots SI Fit Features
routineModelAnalysisPlot_GroupLevelFitStats(R)
routineModelAnalysisPlot_GroupLevelFitParameters(R)
routineModelAnalysisPlot_ExampleCorrelations(R)
routineModelAnalysisPlot_IOResponsePerParameter(R)
