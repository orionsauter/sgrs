function AnalyzeStudy(study_name, xlab, varargin)
params = readtable(strcat(study_name,'/params.txt'));
if size(params,2) > 2
    param_values = params.(strcat('param_values_',string(varargin{1})));
else
    param_values = params.param_values;
end
for i=params.idx'
    matfile = dir(strcat(study_name,'/',study_name,string(i),'_*/TMrel_*.mat'));
    data(i) = load(strcat(matfile.folder,'/',matfile.name));
end

figure;
plot(param_values, arrayfun(@(x) max(abs(x.xTM(:,1))), data), '.-');
xlabel(xlab);
ylabel('max x disp. [m]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",false);
close;

figure;
plot(param_values, arrayfun(@(x) max(abs(x.xTM(:,2))), data), '.-');
xlabel(xlab);
ylabel('max y disp. [m]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) max(abs(x.xTM(:,3))), data), '.-');
xlabel(xlab);
ylabel('max z disp. [m]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) max(abs(x.bTM(:,1))), data), '.-');
xlabel(xlab);
ylabel('max psi disp. [rad]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) max(abs(x.bTM(:,2))), data), '.-');
xlabel(xlab);
ylabel('max theta disp. [rad]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) max(abs(x.bTM(:,3))), data), '.-');
xlabel(xlab);
ylabel('max phi disp. [rad]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) x.vTM(end,1), data), '.-');
xlabel(xlab);
ylabel('final x vel. [m/s]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) x.vTM(end,2), data), '.-');
hold on;
if max(arrayfun(@(x) x.vTM(end,2), data))/4.5e-6 > 0.1
    yline(4.5e-6, '--r');     % Max residual velocity from "ResVel_2DOF.m"
    legend('', 'limit');
end
xlabel(xlab);
ylabel('final y vel. [m/s]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) x.vTM(end,3), data), '.-');
xlabel(xlab);
ylabel('final z vel. [m/s]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) x.wTM(end,1), data), '.-');
xlabel(xlab);
ylabel('final psi vel. [rad/s]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) x.wTM(end,2), data), '.-');
xlabel(xlab);
ylabel('final theta vel. [rad/s]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;

figure;
plot(param_values, arrayfun(@(x) x.wTM(end,3), data), '.-');
xlabel(xlab);
ylabel('final phi vel. [rad/s]');
exportgraphics(gcf,strcat(study_name,'/Results.pdf'),"Append",true);
close;
end
