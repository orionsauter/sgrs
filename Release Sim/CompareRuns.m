function [meandiffs] = CompareRuns(path1, path2)
t1 = load(path1,'t').t;
t2 = load(path2,'t').t;
meandiffs = struct('xTM', [], 'bTM', [], 'vTM', [], 'wTM', []);
vars = fieldnames(meandiffs);
for vari = 1:length(vars)
    var = vars{vari};
    var1 = load(path1,char(var)).(var);
    var2 = load(path2,char(var)).(var);
    tmax = min([t1(end),t2(end)]);
    var2i = [interp1(t2,var2(:,1),t1(t1<=tmax)),...
             interp1(t2,var2(:,2),t1(t1<=tmax)),...
             interp1(t2,var2(:,3),t1(t1<=tmax))];
    var1 = var1(t1<=tmax);
    % se = ((var2i-var1)./var1).^2;
    se = ((var2i-var1)).^2;
    meandiffs.(var) = sqrt(mean(se(isfinite(se))));
end
meandiffs
end
