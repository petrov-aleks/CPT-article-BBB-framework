%% Generate Table S4: clinical studies and modeled dosing regimens

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
inputFile = fullfile(repoRoot, 'data', 'clinical', ...
    'table_s4_clinical_studies.csv');
outputDir = fullfile(repoRoot, 'results', 'tables');
outputFile = fullfile(outputDir, 'table_s4.csv');

if ~isfile(inputFile)
    error('Required Table S4 input was not found: %s', inputFile);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

opts = detectImportOptions(inputFile);
opts = setvartype(opts, ...
    {'Drug','DosingRegimen','AgeYears','Indication','SampleCollection', ...
     'Reference','ModelDoseUnit'}, 'string');
tableS4 = readtable(inputFile, opts);

requiredVariables = {'Drug','DosingRegimen','NSubjects','AgeYears', ...
    'MalePercent','Indication','SampleCollection','Reference', ...
    'ModelDoseValue','ModelDoseUnit','ModelAdministrationDurationMin'};
if ~all(ismember(requiredVariables, tableS4.Properties.VariableNames))
    error('Table S4 input does not contain all required variables.');
end

writetable(tableS4, outputFile);

disp(tableS4)
fprintf('Saved %s\n', outputFile);
