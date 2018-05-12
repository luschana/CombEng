dataStruct=importdata(datFile);
cp_meas=dataStruct.data(:,1);
t=dataStruct.data(:,2)/1000;

fitParams=polyfit(t,cp_meas,3);

%t_Sh=linspace(min(t),max(t),25);
t_Sh=linspace(0.2,2,50);
cp_Sh=polyval(fitParams,t_Sh);
figure
plot(t,cp_meas,'o')
hold on
plot(t_Sh,cp_Sh)

% this were the "old values for octane -- completly different to
% fitted/measured
%oldPoly=[351.455, 279.288];
%cp_Old=polyval(oldPoly,t_Sh);
%plot(t_Sh,cp_Old, '--');

legend('cp_meas','cp_Sh')

formatSpec='ShDataEntry(%.2f, %.2f, %f, %f, %f, %f, %1.1f, %1.1f, %1.1f, %1.1f)';
%ShDataEntry(298.0, 6000.0, 351.455, 279.288, 0.0, 0.0, 0.0, 0.0, 0.0, -208700.0)
fprintf('//%s\n',datFile)
fprintf('const ShDataEntry ShData_Fuel[2] = {\n')
fprintf(formatSpec, min(t)*1000, max(t)*1000, fitParams(4), fitParams(3), fitParams(2), fitParams(1), 0, 0, 0, 0)
fprintf(',\n')
fprintf(formatSpec, max(t)*1000, 6000, cp_meas(length(cp_meas)), 0, 0, 0, 0, 0, 0, 0)
fprintf('};\n')

%eof