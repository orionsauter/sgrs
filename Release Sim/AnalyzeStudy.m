function AnalyzeStudy(study_name, xlab)
params = readtable(strcat(study_name,'/params.txt'));
% figure(1);
% hold on;
y_vels = zeros([size(params,1),1]);
for i=params.idx'
    matfile = dir(strcat(study_name,'/',study_name,string(i),'_*/TMrel_*.mat'));
    load(strcat(matfile.folder,'/',matfile.name), 't','vTM');
    y_vels(i) = max(abs(vTM(:,2)));
    % plot(t, vTM(:,2), "-");
end
% legend(string(params.param_values'));
% xlabel('time [s]');
% ylabel('TM y vel. [m/s]');
plot(params.param_values, y_vels, ".-");
xlabel(xlab);
ylabel('max y vel. [m/s]');
saveas(figure(1),strcat(study_name,'/y_vel.png'));
end
