%% Generate the Mannitol column of Table S3 from the latest optimization
clear; clc;
scriptDir=fileparts(mfilename('fullpath')); repoRoot=fileparts(fileparts(scriptDir));
addpath(fullfile(repoRoot,'code','utilities'));
baseFile=fullfile(repoRoot,'data','model_parameters','table_s3_mannitol_base_parameters.csv');
tableS1File=fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
fitFile=fullfile(repoRoot,'results','model_parameters','mannitol_optimized_permeability.csv');
settingsFile=fullfile(repoRoot,'data','model_parameters','mannitol_simulation_settings.csv');
outputFile=fullfile(repoRoot,'results','tables','table_s3_mannitol.csv');
tableS3Mannitol=build_table_s3_mannitol(baseFile,tableS1File,fitFile,settingsFile);
writetable(tableS3Mannitol,outputFile); disp(tableS3Mannitol); fprintf('Saved %s\n',outputFile);
