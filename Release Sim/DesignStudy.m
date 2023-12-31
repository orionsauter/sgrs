function DesignStudy(study_name, param_name, param_values, xlab)
config = readstruct('config.xml');
mkdir(study_name);
for i=1:length(param_values)
    this_run = strcat(study_name, string(i));
    config.run_name = strcat("'",this_run,"'");
    config.(param_name) = mat2str(param_values(i,:));
    config.plots = 0;
    config_file = strcat(study_name, '/', this_run, '.xml');
    writestruct(config, config_file, 'StructNodeName', 'config');
end
if size(param_values,2) > 1
    dp = find(abs(diff(param_values,1,1))>0,1);
    analyze_cmd = strcat('sbatch --dependency=afterok:$SLURM_ARRAY_JOB_ID RunAnalyze.sh "', ...
        study_name, '" "', xlab,'" ',dp);
else
    analyze_cmd = strcat('sbatch --dependency=afterok:$SLURM_ARRAY_JOB_ID RunAnalyze.sh "', ...
        study_name, '" "', xlab,'"');
end
slurm = [...
    '#!/bin/bash --login',...
    strcat('#SBATCH --job-name=', study_name),...
    '#SBATCH --mail-type=NONE',...
    '#SBATCH --account=conklin',...
    '#SBATCH --qos=conklin',...
    '#SBATCH --cpus-per-task=2',...
    '#SBATCH --mem=16gb',...
    '#SBATCH --time=0-01:00:00',...
    strcat('#SBATCH --output=', study_name, '_%j.log'),...
    strcat('#SBATCH --array=1-', string(length(param_values))),...
    'pwd; hostname; date',...
    'cd /blue/conklin/orionsauter/SGRS/',...
    'module load matlab',...
    strcat('/usr/bin/time -v matlab -nodisplay -r "TMrel_HG(''', ...
        study_name, '/', study_name, '$SLURM_ARRAY_TASK_ID.xml''); exit;"'), ...
    'if [ $SLURM_ARRAY_TASK_ID -eq 1 ]', 'then', ...
    analyze_cmd, 'fi'];
writelines(slurm,strcat('Run_', study_name, '.sh'));
idx = (1:length(param_values))';
writetable(table(idx,param_values), strcat(study_name, '/params.txt'));
end
