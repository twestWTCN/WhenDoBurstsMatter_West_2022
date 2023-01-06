function R = setupPaths(repopath,routname)
switch getenv('computername')
    case {'DESKTOP-0HO6J14','OPAL'}
        gitpath = 'Z:\TimWest\GITHUB';
        spmpath = 'Z:\TimWest\GITHUB\spm12';
        R.path.KM_root = 'Z:\TimWest\DATA\Data\K_Miller_Data';
% Creat new case for your computer
% %     case #PC name#
% %         gitpath = #path for github dependencies#
% %         spmpath = #path to SPM#
% %         R.path.KM_root = #path to Miller data#
% %         R.path.localDatPath = #write path for local data storage#;
end
R.path.root = [repopath];
R.path.rootn = R.path.root; 
R.path.projectn = routname;
R.path.projpath =  [R.path.root '\Projects\' R.path.projectn];
R.path.datapath =  [R.path.projpath filesep 'data'];
R.path.localDatPath  = R.path.datapath;
pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(spmpath, pathCell));
if ~onPath; addpath(spmpath); spm eeg; close all; end


addpath(genpath([R.path.rootn '\Methods']))
addpath(genpath([R.path.rootn '\dependencies']))

%% You need these dependencies
addpath(genpath([gitpath '\ABCNeuralModellingToolbox\ABC_dependencies'])); % https://github.com/twestWTCN/ABCNeuralModellingToolbox
addpath(genpath([gitpath '\ABCNeuralModellingToolbox\sim_machinery']))
addpath(genpath([gitpath '\Violinplot-Matlab'])); % https://github.com/bastibe/Violinplot-Matlab

addpath(genpath(repopath));
addpath(genpath(R.path.projpath))


