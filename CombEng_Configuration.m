Ts=1e-4; % [s]

% Add more parameter checkers/updaters if required
checkfiles = {strcat(gcs, '_ParamCheck')};

% Set options for ArteLabModelInstaller, see help create_setup_file
% Note: parameterization via RMP block (ParameterSource=PumaParameterPool) is only supported for PUMA Open 1.5.1 or higher
%ArteLabModelInstallerOptions = {};
%ArteLabModelInstallerOptions = { 'Setup:EptName=';                'Setup:ParameterFile='; 'Setup:ParameterSource='};                  % no parameters, no EPT
%ArteLabModelInstallerOptions = {['Setup:EptName=' checkfiles{1}];                         'Setup:ParameterSource=PumaParameterPool'}; % via RMP block with EPT
%ArteLabModelInstallerOptions = {['Setup:EptName=' checkfiles{1}];                         'Setup:ParameterSource=ParameterFile'};     % via MAT file with EPT (default)
%ArteLabModelInstallerOptions = { 'Setup:EptName=';                                        'Setup:ParameterSource=PumaParameterPool'}; % via RMP block without EPT 
ArteLabModelInstallerOptions = { 'Setup:EptName=';                                        'Setup:ParameterSource=ParameterFile'};     % via MAT file without EPT

% settings for RTZ package
%ArteLabModelStartupArguments = ['-f3 -b3 -w..\cfg -i' modelfile '.ini -l2'];
ArteLabModelStartupArguments = '';
%ArteLabModelMemorySettings = {'PoolMin=200000';'PoolMax=-1';'VsegSize=80000000';'ObjDirSize=200'};
ArteLabModelMemorySettings = {''};
 
