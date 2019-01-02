% Startup script for an fmi.LAB project
% slightly modified by luschana -- 2017-12-27

% -------------------------------------------------------------------------
try
    if isdeployed
        disp('luschana: startup.m skipped, something is deployed...?!?')
        return
    else
        disp('luschana: called startup.m')
        clear all;
        % start AVL fmi.LAB environment -- does not too much: sets some version
        artelabrc;
        clear ACI_channel_list;
    end
catch exc
        warning(message('ülülüü:', exc.identifier, exc.message));
end

%%open required libraries
%addpath('E:\AVL\buildEnv\lib');
%open_MyLibs

%%%%%%%%%%%%%%%%%%%%%%%%
%%  Model spec. defs %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Define file names
modelfile  = 'CombEng';
% Define model step size
Ts=1e-4; % [s]

checkfiles = {''}; %names script(s) to be run with "modelname" as single argument

% Define global variables
global MDL

% load parameter set defined in ini file
disp(['Reading parameter file from: "' modelfile '.ini" ...']);
[err,matfile]=get_matfilename_from_inifile([modelfile '.ini']);
if ~err
  load(matfile);
else
    disp(['Failed to load parameter file defined in"' modelfile '.ini" ...']);
end

% Set options for ArteLabModelInstaller, see help create_setup_file -- but only for rta!?!
%ArteLabModelInstallerOptions = {};
%ArteLabModelInstallerOptions = { 'Setup:EptName='; 'Setup:ParameterFile='; 'Setup:ParameterSource='}; % no parameters, no EPT
%ArteLabModelInstallerOptions = { 'Setup:EptName='; 'Setup:ParameterSource=PumaParameterPool'}; % via RMP block without EPT 
ArteLabModelInstallerOptions = { 'Setup:EptName='; 'Setup:ParameterSource=ParameterFile'};     % via MAT file without EPT

% settings for RTZ package
%ArteLabModelStartupArguments = ['-f3 -b3 -w..\cfg -i' modelfile '.ini -l2'];
ArteLabModelStartupArguments = '';
%ArteLabModelMemorySettings = {'PoolMin=200000';'PoolMax=-1';'VsegSize=80000000';'ObjDirSize=200'};
ArteLabModelMemorySettings = {''};


% Open model file
disp('Opening model file...');
eval(modelfile);