%% Generate the Ibuprofen column of Table S3 from the latest optimization
clear; clc;
scriptDir=fileparts(mfilename('fullpath')); repoRoot=fileparts(fileparts(scriptDir));
addpath(fullfile(repoRoot,'code','utilities'));
baseFile=fullfile(repoRoot,'data','model_parameters','table_s3_ibuprofen_base_parameters.csv');
tableS1File=fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
fitFile=fullfile(repoRoot,'results','model_parameters','ibuprofen_optimized_permeability.csv');
settingsFile=fullfile(repoRoot,'data','model_parameters','ibuprofen_simulation_settings.csv');
outputFile=fullfile(repoRoot,'results','tables','table_s3_ibuprofen.csv');
tableS3Ibuprofen=build_table_s3_ibuprofen(baseFile,tableS1File,fitFile,settingsFile);
writetable(tableS3Ibuprofen,outputFile); disp(tableS3Ibuprofen); fprintf('Saved %s\n',outputFile);
