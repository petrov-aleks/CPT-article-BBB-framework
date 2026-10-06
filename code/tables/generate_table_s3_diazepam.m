%% Generate the Diazepam column of Table S3 from the latest optimization

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(repoRoot,'code','utilities'));

baseFile = fullfile(repoRoot,'data','model_parameters','table_s3_diazepam_base_parameters.csv');
tableS1File = fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
fitFile = fullfile(repoRoot,'results','model_parameters','diazepam_optimized_permeability.csv');
settingsFile = fullfile(repoRoot,'data','model_parameters','diazepam_simulation_settings.csv');
outputFile = fullfile(repoRoot,'results','tables','table_s3_diazepam.csv');

tableS3Diazepam = build_table_s3_diazepam(baseFile,tableS1File,fitFile,settingsFile);
writetable(tableS3Diazepam,outputFile);
disp(tableS3Diazepam)
fprintf('Saved %s\n',outputFile);
