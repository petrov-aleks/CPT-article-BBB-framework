%% Calculate the GMR and GMFE values used in Figure 5

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
clinicalDataDir = fullfile(repoRoot, 'data', 'clinical', 'figure_4');
simulationDataDir = fullfile(repoRoot, 'results', 'simulations', ...
    'main_figures', 'figure_4');
outputDir = fullfile(repoRoot, 'results', 'derived', ...
    'main_figures', 'figure_5');
outputFile = fullfile(outputDir, 'fold_error_summary.csv');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

compounds = {'Ethanol', 'Paracetamol', 'Mannitol', ...
    'Ibuprofen', 'Diazepam', 'Indomethacin'};
tolerance = 1e-6;
summaryRows = {};

for c = 1:numel(compounds)
    compound = compounds{c};
    clinicalFile = fullfile(clinicalDataDir, ...
        sprintf('%s_clinical_data.csv', compound));
    simulationFile = fullfile(simulationDataDir, ...
        sprintf('%s_simulations.csv', compound));

    if ~isfile(clinicalFile)
        error('Required clinical data file was not found: %s', clinicalFile);
    end
    if ~isfile(simulationFile)
        error(['Required generated simulation file was not found: %s\n' ...
            'Run the Figure 4 simulation workflows first.'], simulationFile);
    end

    clinical = readtable(clinicalFile);
    opts = detectImportOptions(simulationFile);
    opts = setvartype(opts, 3, 'char');
    simulations = readtable(simulationFile, opts);

    % Figure 5 evaluates CNS observations; plasma observations are excluded.
    clinical = clinical(~strcmp(clinical.site, 'plasma'), :);
    clinical.site(strcmp(clinical.site, 'brain mass')) = {'brain ECF'};
    simulations = simulations(~strcmp(simulations.site, 'plasma'), :);

    % Harmonize the two cell-based assay labels used by the simulations.
    cellBased = strcmp(simulations.source, 'in vitro MDCK/UC') | ...
        strcmp(simulations.source, 'in vitro Caco-2');
    simulations.source(cellBased) = {'in vitro cell-based'};

    permeabilityValues = unique(simulations.PS);
    sites = unique(simulations.site);

    for i = 1:numel(permeabilityValues)
        ps = permeabilityValues(i);
        simulationsAtPS = simulations( ...
            abs(simulations.PS - ps) < tolerance, :);

        for j = 1:numel(sites)
            site = sites{j};
            simulatedSite = simulationsAtPS( ...
                strcmp(simulationsAtPS.site, site), :);
            observedSite = clinical(strcmp(clinical.site, site), :);

            if isempty(simulatedSite) || isempty(observedSite) || ...
                    height(simulatedSite) < 2
                continue
            end

            predicted = interp1(simulatedSite.time, simulatedSite.conc, ...
                observedSite.time, 'linear', 'extrap');
            observed = observedSite.conc;

            valid = predicted > 0 & observed > 0;
            logRatios = log(predicted(valid) ./ observed(valid));
            if isempty(logRatios)
                continue
            end

            GMR = exp(mean(logRatios));
            GMFE = exp(mean(abs(logRatios)));
            summaryRows(end+1,:) = {compound, GMR, GMFE, ...
                simulatedSite.source{1}, ps}; %#ok<SAGROW>
        end
    end
end

foldErrorSummary = cell2table(summaryRows, 'VariableNames', ...
    {'Compound', 'GMR', 'GMFE', 'Source', 'PS'});
writetable(foldErrorSummary, outputFile);

disp(foldErrorSummary)
fprintf('Saved %s\n', outputFile);
