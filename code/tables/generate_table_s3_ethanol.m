%% Generate the Ethanol column of Table S3

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
inputFile = fullfile(repoRoot, 'data', 'model_parameters', ...
    'table_s3_ethanol_base_parameters.csv');
tableS1File = fullfile(repoRoot, 'data', 'model_parameters', ...
    'table_s1_base_parameters.csv');
settingsFile = fullfile(repoRoot, 'data', 'model_parameters', ...
    'ethanol_simulation_settings.csv');
fitResultFile = fullfile(repoRoot, 'results', 'model_parameters', ...
    'ethanol_optimized_permeability.csv');
outputDir = fullfile(repoRoot, 'results', 'tables');
outputFile = fullfile(outputDir, 'table_s3_ethanol.csv');

requiredFiles = {inputFile,tableS1File,settingsFile,fitResultFile};
for i = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{i})
        error(['Required input was not found: %s\nRun the ethanol ' ...
            'Figure 4 simulation workflow before generating Table S3.'], ...
            requiredFiles{i});
    end
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

addpath(fullfile(repoRoot, 'code', 'utilities'));
tableS3Ethanol = build_table_s3_ethanol( ...
    inputFile, tableS1File, fitResultFile, settingsFile);
writetable(tableS3Ethanol, outputFile);

disp(tableS3Ethanol)
fprintf('Saved %s\n', outputFile);
