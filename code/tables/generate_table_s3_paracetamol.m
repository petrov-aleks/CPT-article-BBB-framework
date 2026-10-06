%% Generate the Paracetamol column of Table S3 from the latest optimization

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(repoRoot,'code','utilities'));

baseFile = fullfile(repoRoot,'data','model_parameters','table_s3_paracetamol_base_parameters.csv');
tableS1File = fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
fitFile = fullfile(repoRoot,'results','model_parameters','paracetamol_optimized_permeability.csv');
settingsFile = fullfile(repoRoot,'data','model_parameters','paracetamol_simulation_settings.csv');
outputFile = fullfile(repoRoot,'results','tables','table_s3_paracetamol.csv');

tableS3Paracetamol = build_table_s3_paracetamol(baseFile,tableS1File,fitFile,settingsFile);
writetable(tableS3Paracetamol,outputFile);
disp(tableS3Paracetamol)
fprintf('Saved %s\n',outputFile);
