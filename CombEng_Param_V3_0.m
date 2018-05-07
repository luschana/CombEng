function CombEng=CombEng_Param_V3_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CombEng = MPE_Folder ('3.0.0.0000', 'CombEng', '', ''); % version, comment, help, options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CombEng.Gain1 = MPE_Scalar('1.0', 'Gain1', '', '', 'Gain1', '-', [0 100], 4.0, 'DOUBLE', 2, 1);

