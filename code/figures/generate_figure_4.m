%% Generate Figure 4: concentration-time profiles

clear; clc; close all;

% -------------------------------------------------------------------------
% GLOBAL FIGURE STYLE PARAMETERS (CPT:PSP journal compliant)
% -------------------------------------------------------------------------

FIG_WIDTH  = 17.8;   % cm
FIG_HEIGHT = 18;     % cm

AXIS_LINE_WIDTH = 1.2;

LINE_WIDTH_MAIN = 1.5;
LINE_WIDTH_SIM  = 1.5;

MARKER_SIZE = 5;
DATA_LINE_WIDTH = 1.2;

X_TICK_FONT = 7;
Y_TICK_FONT = 7;

AXIS_LABEL_FONT = 7;
TITLE_FONT = 7;
PLASMA_LEGEND_FONT = 4;
MAIN_LEGEND_FONT = 7;

% -------------------------------------------------------------------------
% Legend location per compound (plasma panels)
% -------------------------------------------------------------------------

legendLocationMap = containers.Map( ...
    {'Ethanol','Paracetamol','Mannitol','Ibuprofen','Diazepam','Indomethacin'}, ...
    {'northeast','southwest','southeast','southeast','southeast','southeast'} );

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
clinicalDataDir = fullfile(repoRoot, 'data', 'clinical', 'figure_4');
simulationDataDir = fullfile(repoRoot, 'results', 'simulations', ...
    'main_figures', 'figure_4');
outputDir = fullfile(repoRoot, 'results', 'figures', 'main');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------

figure
set(gcf,'Units','centimeters','Position',[2 2 FIG_WIDTH FIG_HEIGHT])

tl = tiledlayout(6,3,'TileSpacing','compact','Padding','compact');

% -------------------------------------------------------------------------
% Style configuration
% -------------------------------------------------------------------------

colorMap = containers.Map( ...
    {'optimized','perfusion-limited','in situ perfusion','in vitro cell-based','in vitro PAMPA', ...
     'in silico perfusion','in silico PAMPA'}, ...
    {[0 0 0], [0 0 0], [0.9255 0.0039 0.0039], [0 0.4235 0.9451], [0 0.4745 0], ...
     [0.9255 0.0039 0.0039], [0 0.4745 0]} );

lineStyleMap = containers.Map( ...
    {'optimized','perfusion-limited','in situ perfusion','in vitro cell-based','in vitro PAMPA', ...
     'in silico perfusion','in silico PAMPA'}, ...
    {'-','--','-','-','-','--','--'} );

legendOrder = { ...
    'optimized','perfusion-limited','in situ perfusion', ...
    'in vitro cell-based','in vitro PAMPA','in silico perfusion','in silico PAMPA'};

expMarkers = {'o','^','s'};

% -------------------------------------------------------------------------
% Compound configuration
% -------------------------------------------------------------------------

compounds = {
'Ethanol'      'g/L'     9    1e-2   1e1    3   'brain mass'  0:2:8     [1e-2 1e-1 1e0 1e1]
'Diazepam'     'ng/mL'   25   1e-2   1e3    1   'scsf'        0:4:24    [1e-2 1e0 1e2 1e4 1e6]
'Paracetamol'  'mg/L'    13   1e-2   1e2    3   'scsf'        0:2:12    [1e-2 1e0 1e2]
'Ibuprofen'    'ng/mL'   25   1e-2   1e4    1   'scsf'        0:4:24    [1e-2 1e0 1e2 1e4 1e6]
'Indomethacin' 'ng/mL'   25   1e-2   1e4    1   'scsf'        0:4:24    [1e-2 1e0 1e2 1e4 1e6]
'Mannitol'     'g/L'     7    1e-4   1e1    1   'scsf'        0:2:6     [1e-3 1e-1 1e1]
};

sim_sites = {'plasma','brain ECF','sCSF'};
titles    = {'Plasma','Brain ECF','CSF'};

tile_id = 1;

% -------------------------------------------------------------------------
% Main loops
% -------------------------------------------------------------------------

for c = 1:size(compounds,1)

    compound      = compounds{c,1};
    unit          = compounds{c,2};
    xmax          = compounds{c,3};
    ymin          = compounds{c,4};
    ymax          = compounds{c,5};
    nExp          = compounds{c,6};
    brainDataSite = compounds{c,7};
    xticks_vals   = compounds{c,8};
    yticks_vals   = compounds{c,9};

    dataFile = fullfile(clinicalDataDir, ...
        sprintf('%s_clinical_data.csv', compound));
    simFile = fullfile(simulationDataDir, ...
        sprintf('%s_simulations.csv', compound));

    if ~isfile(dataFile)
        error('Required clinical data file was not found: %s', dataFile);
    end
    if ~isfile(simFile)
        error(['Required generated simulation file was not found for %s: %s\n' ...
            'Run the corresponding code/simulations/generate_figure_4_* script first.'], ...
            compound, simFile);
    end

    %% Read data

    Tdata = readtable(dataFile);

    opts = detectImportOptions(simFile);
    opts = setvartype(opts,3,"char");
    Tsim = readtable(simFile,opts);

    if ismember('source',Tsim.Properties.VariableNames)
        Tsim.source(strcmp(Tsim.source,'in vitro MDCK/UC')) = {'in vitro cell-based'};
    end

    %% Loop sites

    for s = 1:3

        ax = nexttile(tile_id);
        hold(ax,'on')

        sim_site = sim_sites{s};

        %% Experimental data

        plotExp = false;

        if strcmp(sim_site,'plasma')
            plotExp  = true;
            dataSite = 'plasma';
        end

        if strcmp(sim_site,'brain ECF') && strcmp(brainDataSite,'brain mass')
            plotExp  = true;
            dataSite = 'brain mass';
        end

        if strcmp(sim_site,'sCSF') && strcmp(brainDataSite,'scsf')
            plotExp  = true;
            dataSite = 'scsf';
        end

        if plotExp

            Tbl = Tdata(strcmpi(Tdata.site,dataSite),:);
            expSources = unique(Tbl.source);

            for k = 1:min(length(expSources),nExp)

                R = Tbl(strcmp(Tbl.source,expSources{k}),:);

                if ismember('sd_down',Tbl.Properties.VariableNames)

                    errorbar(ax,...
                        R.time,...
                        R.conc,...
                        R.sd_down,...
                        R.sd_up,...
                        'k',...
                        'LineStyle','none',...
                        'LineWidth',DATA_LINE_WIDTH,...
                        'HandleVisibility','off');

                end

                plot(ax,...
                    R.time,...
                    R.conc,...
                    'Marker',expMarkers{k},...
                    'MarkerSize',MARKER_SIZE,...
                    'MarkerFaceColor','none',...
                    'MarkerEdgeColor','k',...
                    'LineStyle','none',...
                    'LineWidth', DATA_LINE_WIDTH, ...
                    'DisplayName',expSources{k});

            end

        end

        %% Simulations

        Sim = Tsim(strcmpi(Tsim.site,sim_site),:);

        if strcmp(sim_site,'plasma')

            plot(ax,...
                Sim.time,...
                Sim.conc,...
                'k',...
                'LineWidth',LINE_WIDTH_MAIN,...
                'HandleVisibility','off');

        else

            for l = 1:length(legendOrder)

                src = legendOrder{l};
                SrcTbl = Sim(strcmp(Sim.source,src),:);

                if isempty(SrcTbl)
                    continue
                end

                col = colorMap(src);
                ls  = lineStyleMap(src);

                PS_values = unique(SrcTbl.PS);

                for p = 1:length(PS_values)

                    R = SrcTbl(SrcTbl.PS==PS_values(p),:);

                    plot(ax,...
                        R.time,...
                        R.conc,...
                        'Color',col,...
                        'LineStyle',ls,...
                        'LineWidth',LINE_WIDTH_SIM,...
                        'HandleVisibility','off');

                end

            end

        end

        %% Legend (custom per compound)

        if strcmp(sim_site,'plasma')

            lg_loc = legendLocationMap(compound);

            lg = legend(ax,...
                'Location', lg_loc, ...
                'FontSize',PLASMA_LEGEND_FONT,...
                'Box','on');

            if strcmp(compound,'Ethanol')
                lg.Position(1) = lg.Position(1)*1.05;
                lg.Position(2) = lg.Position(2)*0.94;
            elseif strcmp(compound,'Paracetamol')
                lg.Position(1) = lg.Position(1)*0.97;
                lg.Position(2) = lg.Position(2)*0.94 ;
            end

        end

        %% Axis formatting

        set(ax,'YScale','log','LineWidth',AXIS_LINE_WIDTH)

        ax.XAxis.FontSize = X_TICK_FONT;
        ax.YAxis.FontSize = Y_TICK_FONT;
        ax.YMinorTick = 'off';

        xlim([0 xmax])
        ylim([ymin ymax])

        xticks(ax,xticks_vals)
        yticks(ax,yticks_vals)

        if strcmp(compound,'Mannitol')
            xlabel(ax,'Time (h)','FontSize',AXIS_LABEL_FONT)
        end
        if strcmp(sim_site,'plasma')
            ylabel(ax,{compound; ['Concentration (' unit ')']},...
                'FontSize',AXIS_LABEL_FONT)
        end

        if strcmp(compound,'Ethanol')
            title(ax,titles{s},...
                'FontSize',TITLE_FONT,'FontWeight','normal')
        end

        box(ax,'on')

        tile_id = tile_id + 1;

    end

end

% -------------------------------------------------------------------------
% Global legend
% -------------------------------------------------------------------------

legendHandles = gobjects(1,length(legendOrder));

for l = 1:length(legendOrder)

    src = legendOrder{l};

    legendHandles(l) = plot(nan,nan,...
        'Color',colorMap(src),...
        'LineStyle',lineStyleMap(src),...
        'LineWidth',LINE_WIDTH_SIM);

end

lgd = legend(legendHandles,legendOrder,...
    'Orientation','horizontal',...
    'FontSize',MAIN_LEGEND_FONT,...
    'Box','on');

lgd.ItemTokenSize = [20 5];   % compact legend lines
lgd.NumColumns = 4;
lgd.Layout.Tile = 'north';

% -------------------------------------------------------------------------
% Export PDF (vector)
% -------------------------------------------------------------------------

set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 FIG_WIDTH FIG_HEIGHT])
set(gcf,'PaperSize',[FIG_WIDTH FIG_HEIGHT])

pdfFile = fullfile(outputDir, 'Figure_4.pdf');
pngFile = fullfile(outputDir, 'Figure_4.png');
print(gcf, pdfFile, '-dpdf', '-vector')
exportgraphics(gcf, pngFile, 'Resolution', 600)

fprintf('Saved %s\n', pdfFile);
fprintf('Saved %s\n', pngFile);
