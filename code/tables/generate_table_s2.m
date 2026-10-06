%% Generate Table S2: drug-specific systemic PK parameters

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
inputDir = fullfile(repoRoot, 'data', 'model_parameters');
parameterFile = fullfile(inputDir, 'table_s2_plasma_pk_parameters.csv');
sourceFile = fullfile(inputDir, 'table_s2_compound_sources.csv');
outputDir = fullfile(repoRoot, 'results', 'tables');

requiredFiles = {parameterFile, sourceFile};
for i = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{i})
        error('Required Table S2 input was not found: %s', requiredFiles{i});
    end
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

opts = detectImportOptions(parameterFile);
opts = setvartype(opts, opts.VariableNames, 'string');
parameters = readtable(parameterFile, opts);
sources = readtable(sourceFile, 'TextType', 'string');

compounds = ["Ethanol","Diazepam","Paracetamol", ...
    "Ibuprofen","Indomethacin","Mannitol"];
tableS2 = parameters(:, {'Category','Parameter','Description','Unit'});

for compound = compounds
    valueColumn = compound + "_value";
    cvColumn = compound + "_CV_percent";
    values = parameters.(valueColumn);
    cvs = parameters.(cvColumn);
    displayed = values + " (" + cvs + ")";
    tableS2.(compound) = displayed;
end

% Table S2 reports compact rounded values, while its source CSV retains the
% unrounded transformed mannitol parameters used by the simulations.
tableS2.Mannitol(parameters.Parameter == "Vpone") = "13 (91)";
tableS2.Mannitol(parameters.Parameter == "Qcpone") = "36.1 (75)";
tableS2.Mannitol(parameters.Parameter == "CL") = "7.3 (74)";

writetable(tableS2, fullfile(outputDir, 'table_s2.csv'));
writetable(sources, fullfile(outputDir, 'table_s2_sources.csv'));

disp(tableS2)
fprintf('Saved %s\n', fullfile(outputDir, 'table_s2.csv'));
fprintf('Saved %s\n', fullfile(outputDir, 'table_s2_sources.csv'));
