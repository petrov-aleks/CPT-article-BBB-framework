%% Generate the permeability values reported in Table 1
% Reads the source permeability data and compound descriptors, calculates
% the two in silico estimates, and converts all estimates to PS_BBB (L/h).

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
inputDir = fullfile(repoRoot, 'data', 'permeability');
outputDir = fullfile(repoRoot, 'results', 'tables');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

descriptors = readtable(fullfile(inputDir, 'compound_descriptors.csv'), ...
    'TextType', 'string');
experimental = readtable(fullfile(inputDir, 'experimental_permeability.csv'), ...
    'TextType', 'string');
modelFile = fullfile(outputDir, 'table_s5.csv');
if ~isfile(modelFile)
    error(['Missing %s. Run code/analysis/fit_tsinman_pampa_model.m and ' ...
        'code/tables/generate_table_s5.m first.'], modelFile);
end
models = readtable(modelFile, 'TextType', 'string');

% Use the exact same machine-readable equations reported in Table S5.
smithModel = selectModel(models, "Smith et al. 2024");
tsinmanModel = selectModel(models, "Tsinman-based model");
smithPS0 = evaluateModel(smithModel, descriptors.logP, descriptors.MW_g_per_mol);
tsinmanP0 = evaluateModel(tsinmanModel, descriptors.logP, descriptors.MW_g_per_mol);

n = height(descriptors);
predicted = table( ...
    repelem(descriptors.Compound, 2), ...
    repmat(["in silico perfusion"; "in silico PAMPA"], n, 1), ...
    reshape([smithPS0, NaN(n,1)].', [], 1), ...
    reshape([NaN(n,1), tsinmanP0].', [], 1), ...
    repmat([smithModel.Reference; tsinmanModel.Reference], n, 1), ...
    'VariableNames', experimental.Properties.VariableNames);

allSources = [experimental; predicted];

% Keep the compound and source order used in the manuscript.
compoundOrder = ["Ibuprofen", "Indomethacin", "Diazepam", ...
    "Ethanol", "Paracetamol", "Mannitol"];
sourceOrder = ["in situ perfusion", "ex vivo UC", "in vitro MDCK", ...
    "in vitro Caco-2", "in vitro PAMPA", ...
    "in silico perfusion", "in silico PAMPA"];
compoundRank = zeros(height(allSources), 1);
sourceRank = zeros(height(allSources), 1);
for i = 1:height(allSources)
    compoundRank(i) = find(compoundOrder == allSources.Compound(i), 1);
    sourceRank(i) = find(sourceOrder == allSources.Source(i), 1);
end

[~, rowOrder] = sortrows([compoundRank, sourceRank], [1 2]);
allSources = allSources(rowOrder, :);

% Display units used in Table 1.
PS0_1e4_mL_per_g_per_s = allSources.PS0_mL_per_g_per_s .* 1e4;
P0_1e6_cm_per_s = allSources.P0_cm_per_s .* 1e6;

% Translate normalized rat brain perfusion data and permeability data to
% the human BBB permeability-surface area product used by the PBPK model.
referenceBrainWeight_g = 1450;
humanBBBSurfaceArea_m2 = 15;
PS0_to_PSBBB = 1e-4 * 1e-3 * 3600 * referenceBrainWeight_g;
P0_to_PSBBB = 1e-6 * 1e-2 * 3600 * humanBBBSurfaceArea_m2 * 1000;

% Preserve the two-stage rounding used in the original analysis: PS0 and
% P0 are rounded before their conversion to PS_BBB.
PS0_1e4_mL_per_g_per_s = roundForOriginalAnalysis( ...
    PS0_1e4_mL_per_g_per_s);
P0_1e6_cm_per_s = roundForOriginalAnalysis(P0_1e6_cm_per_s);

figure3PSBBB = PS0_1e4_mL_per_g_per_s .* PS0_to_PSBBB;
useP0 = isnan(figure3PSBBB);
figure3PSBBB(useP0) = P0_1e6_cm_per_s(useP0) .* P0_to_PSBBB;
figure3PSBBB = roundForOriginalAnalysis(figure3PSBBB);

table1 = table(allSources.Compound, allSources.Source, ...
    PS0_1e4_mL_per_g_per_s, P0_1e6_cm_per_s, figure3PSBBB, ...
    allSources.Reference, ...
    'VariableNames', {'Compound','Source','PS0_1e4_mL_per_g_per_s', ...
    'P0_1e6_cm_per_s','PS_BBB_L_per_h','Reference'});
figure3Data = table(table1.Compound, table1.Source, ...
    table1.PS_BBB_L_per_h, table1.Reference, ...
    'VariableNames', {'Compound','Source','PS_BBB_L_per_h','Reference'});

writetable(table1, fullfile(outputDir, 'table_1.csv'));
writetable(figure3Data, fullfile(outputDir, 'figure_3_permeability.csv'));

disp(table1)
fprintf('Saved %s\n', fullfile(outputDir, 'table_1.csv'));
fprintf('Saved %s\n', fullfile(outputDir, 'figure_3_permeability.csv'));

function row = selectModel(models, modelName)
row = models(models.Model == modelName, :);
assert(height(row) == 1, 'Expected exactly one Table S5 row for %s.', modelName);
required = ["Intercept","LogP_Coefficient","Log10_MW_Coefficient", ...
    "MW_over_100_Coefficient"];
assert(all(ismember(required, string(row.Properties.VariableNames))), ...
    'Table S5 lacks machine-readable coefficients. Regenerate it first.');
end

function prediction = evaluateModel(model, logP, molecularWeight)
exponent = model.Intercept + model.LogP_Coefficient .* logP + ...
    model.Log10_MW_Coefficient .* log10(molecularWeight) + ...
    model.MW_over_100_Coefficient .* (molecularWeight ./ 100);
prediction = 10.^exponent;
end

function values = roundForOriginalAnalysis(values)
for i = 1:numel(values)
    x = values(i);
    if ~isnan(x) && x ~= 0
        if abs(x) >= 10
            values(i) = round(x);
        else
            decimals = -floor(log10(abs(x))) + 1;
            values(i) = round(x, decimals);
        end
    end
end
end
