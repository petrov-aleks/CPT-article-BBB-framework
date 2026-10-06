%% Generate the Indomethacin column of Table S3 from the latest optimization
clear; clc;
scriptDir=fileparts(mfilename('fullpath')); repoRoot=fileparts(fileparts(scriptDir));
addpath(fullfile(repoRoot,'code','utilities'));
baseFile=fullfile(repoRoot,'data','model_parameters','table_s3_indomethacin_base_parameters.csv');
tableS1File=fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
fitFile=fullfile(repoRoot,'results','model_parameters','indomethacin_optimized_permeability.csv');
settingsFile=fullfile(repoRoot,'data','model_parameters','indomethacin_simulation_settings.csv');
outputFile=fullfile(repoRoot,'results','tables','table_s3_indomethacin.csv');
tableS3Indomethacin=build_table_s3_indomethacin(baseFile,tableS1File,fitFile,settingsFile);
writetable(tableS3Indomethacin,outputFile); disp(tableS3Indomethacin); fprintf('Saved %s\n',outputFile);
