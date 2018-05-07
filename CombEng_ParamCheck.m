function res=CombEng_ParamCheck(varargin)
% CombEng_PARAMCHECK Checks CombEng parameters for consistency and
%   updates previous versions to the current version using the following rules:
%    - Creates a default parameter set when MAT file does not exist.
%    - Creates a CombEng parameter set in a MAT file which does not contain CombEng parameters.
%    - Updates CombEng parameters in a MAT file, when the CombEng parameters are from a previous version.
%    - Checks CombEng parameters for correctness.
%    - Calculates derived parameters which are typically not visible to the user.
%   
%   CombEng_PARAMCHECK(matfilename) or [MDL] = CombEng_PARAMCHECK(MDL) 
%      operates on the parameter set MDL.CombEng. Nothing will be done in case the field CombEng 
%      does not exist in the structure MDL. In case the MAT-file does not exist a default
%      parameter file containing MDL.CombEng will be created.
%
%   CombEng_PARAMCHECK(matfilename, OPT) or [MDL] = CombEng_PARAMCHECK(MDL, OPT) 
%      operates on the parameter set MDL.CombEng which the following options defined by OPT:
%
%      OPT='-a' adds the structure field 'CombEng' in case it does not exist in the structure MDL.
%               In case the matfile does not exist a default parameter MAT-file containing 
%               MDL.CombEng will be created.
%
%      OPT='-d' deletes the structure filed 'CombEng' in case it exists in the structure MDL.

vers_ParChk2 = '3';
sub_ParChk2  = '0';

disp([mfilename ': parameter checker/updater for model "CombEng v' vers_ParChk2 '.' sub_ParChk2 '".']);

%--- Determine save options ---------------------------------------------------
s         = version;
isMatlab6 = s(1)<'7';

%--- Check parameter options --------------------------------------------------
bOptAdd=0;
bOptDel=0;

%--- Create empty MDL structure -----------------------------------------------
MDL = MPE_Folder ('1.0', 'MDL', '', ''); 

%--- Load argument file -------------------------------------------------------
err  = 0;
done = 0;

%--- Check command line arguments ---------------------------------------------
if nargin>0

  %-- Process options ---------------------------------------------------------
  for ii=2:nargin
    if strcmp(char(varargin{ii}),'-a')
      bOptAdd=1;
    else
      if strcmp(char(varargin{ii}),'-d')
        bOptDel=1;
      end
    end
  end

  if bOptAdd
    bOptDel=0;
  end 

  %-- Process MDL structure ---------------------------------------------------
  if ischar(varargin{1})
    %-- The argument is a file name -------------------------------------------
    fname = char(varargin{1});
 
    %--- Check if filename is either a valid mat-file or not a file at all ----
    [filepath,filename,fileext] = fileparts(fname);

    if length(fileext)==0
      fileext='.mat';
    elseif ~strcmp(fileext,'.mat')
      error(['Invalid file extension: "' fileext '".']);
    end

    if length(filename)==0
      error(['Invalid file name: "' fileext '".']);
    end

    fullfilename=fullfile(filepath,[filename fileext]);

    if exist(fullfilename,'file')
      load(fullfilename,'MDL');
    else
      %--- Create default parameter set ---------------------------------------
      MDL.CombEng = CombEng_Param_V3_0;
      done=1;
    end
  
  else
    %-- The argument is the structure MDL -------------------------------------
    MDL = varargin{1}{1};
  end
  
  %--- Add default parameter set ----------------------------------------------
  if ~isfield(MDL,'CombEng')
    if bOptAdd
      % add component only if option is set
      MDL.CombEng = CombEng_Param_V3_0;
    end
    done=1;
  else
    if bOptDel
      MDL = rmfield(MDL,'CombEng');
      done=1;  
    end
  end
  
else
  error('### Missing argument such as file name or structure variable.')
end

%---- Parameter update --------------------------------------------------------
if done==0

  %--- Check if version property exists----------------------------------------
  if ~isfield(MDL.CombEng,'Type')
    MDL.CombEng
    error('### Invalid data structure, component "MDL.CombEng.Type" does not exist.');  
  end
    
  if ~isfield(MDL.CombEng,'Version')
    MDL.CombEng
    error('### Invalid data structure, component "MDL.CombEng.Version" does not exist.');  
  end
  
  if (length(MDL.CombEng.Version)<3)
    MDL.CombEng.Version
    error('### Invalid version string in component "MDL.CombEng.Version".');  
  end

  switch compare_version_strings(MDL.CombEng.Version,[vers_ParChk2 '.' sub_ParChk2])  % check if equal

    %--- Update previous parameter sets to current version --------------------
    case -1

      % ...Insert here your parameter updater coding, e.g. from v1.0 to v1.1
      %    Be careful with the use of model version numbers to allow the
      %    model to be updated later on ! Use only the string variables
      %    vers_ParChk2 and sub_ParChk2 or the script name CombEng_Param_V...

      error(['### Parameter file is an older version than the current version. ' ... 
             'Don''t know how to update from version v' MDL.CombEng.Version ' to version v' vers_ParChk2 '.' sub_ParChk2 '.']);  

    %--- Update meta data in current version ----------------------------------
    case 0

      MDL_act = MDL.CombEng;
      MDL_def = CombEng_Param_V3_0;

      [err,MDL.CombEng] = MPE_parameter_update(MDL_def, MDL_act);

      if err
        error(['### Error during parameter property update of v' vers_ParChk2 '.' sub_ParChk2 '.']);
      end

    %--- Downdate newer parameter sets to current version ---------------------
    case 1

      % ...Insert here your parameter downdater coding, e.g. from v1.1 to v1.0
      %    Be careful with the use of model version numbers to allow the
      %    model to be updated later on ! Use only the string variables
      %    vers_ParChk2 and sub_ParChk2 or the script name CombEng_Param_V...

      error(['### Parameter file is a newer version than the current version. ' ...
             'Don''t know how to downdate from version v' MDL.CombEng.Version ' to version v' vers_ParChk2 '.' sub_ParChk2 '.']);  

    otherwise
      error(['### Error during parameter update from version v' MDL.CombEng.Version ' to version v' vers_ParChk2 '.' sub_ParChk2 '.']);
  end

  %--- Check parameters -------------------------------------------------------
  % ...Insert here your parameter checking alorithm
  % [err,MDL.CombEng]=CombEng_checker(MDL.CombEng); 
  % if err
  %   error('### Error during parameter checking.');
  % end


  %--- Preprocess parameters --------------------------------------------------
  % ...Insert here your parameter preprocessing for deried parameters
  % [err,MDL.CombEng]=CombEng_preproc(MDL.CombEng);    
  % if err
  %   error('### Error during parameter preprocessing.');
  % end
end

%--- Save MDL to mat-file or return MDL structure -----------------------------
if ischar(varargin{1})
  MDL.MetaVersion = MPE_MetaVersion;
  if isMatlab6
    save(fullfilename, 'MDL');
  else
    save(fullfilename, 'MDL', '-V6');
  end
  res=0;
else
  res={MDL};
end


%------------------------------------------------------------------------------
% Compares main and subversion
% sv1,sv2: strings of type '1.0.00.000', i.e. '<mainver>.<subver1>.<subver2>'
% r:     1 if v1  > v2
%        0 if v1 == v2
%       -1 if v1  < v2
%------------------------------------------------------------------------------
function r = compare_version_strings(sv1, sv2)

r=0;
v1 = get_version(sv1);
v2 = get_version(sv2);

if v1(1)== v2(1),
  if v1(2) == v2(2),
    r = 0;
  elseif v1(2) > v2(2),
    r = 1;
  else
    r = -1;
  end
elseif v1(1) > v2(1)
  r = 1;
else
  r = -1;
end

function v = get_version(s)

v = zeros(1,4);
dots = [0 find(s=='.') size(s,2)+1];
for ii=1:size(dots,2)-1,
  v(ii) = str2num(s(dots(ii)+1:dots(ii+1)-1));
end
